import csv
import json
import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import logging
import tempfile
from typing import Dict, List, Tuple
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report
import joblib

# Define the base paths for the reference_crispr directory
try:
    BASE_PATH = Path(__file__).resolve().parent.parent.parent
except NameError:
    BASE_PATH = Path.cwd()

REFERENCE_CRISPR_PATH = BASE_PATH / 'reference_crispr'
MAP_JSON = REFERENCE_CRISPR_PATH / 'map.json'
SPACERS_FOLDER = REFERENCE_CRISPR_PATH / 'spacers_csv'


class CRRFinder:
    def __init__(self, genome_fasta: Path, output_folder: str, identity_threshold: float = 95.0,
                 coverage_threshold: float = 95.0, match_threshold: float = 95.0):
        self.genome_fasta = Path(genome_fasta)
        self.map_json = MAP_JSON
        self.spacers_folder = SPACERS_FOLDER
        self.output_folder = Path(output_folder) / 'CRR_finder' / self.genome_fasta.stem
        self.output_folder.mkdir(parents=True, exist_ok=True)
        self.identity_threshold = identity_threshold
        self.coverage_threshold = coverage_threshold
        self.match_threshold = match_threshold
        self.spacers = {}
        self.groups = self.load_groups()
        self.blast_output = tempfile.NamedTemporaryFile(delete=False, suffix='_blast_results.txt').name
        self.spacer_fasta = tempfile.NamedTemporaryFile(delete=False, suffix='_spacers.fasta').name
        self.total_spacers_found = 0
        self.setup_logging()

    def setup_logging(self) -> None:
        """Setup logging configuration."""
        logging.basicConfig(filename=self.output_folder / 'genome_analyzer.log',
                            filemode='a',
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            level=logging.INFO)
        logging.info("GenomeAnalyzer initialized.")

    def load_spacers(self, csv_file: Path) -> Dict[str, str]:
        """Load spacer sequences from a CSV file."""
        spacers = {}
        try:
            with open(csv_file, mode='r') as infile:
                reader = csv.reader(infile)
                next(reader)  # skip header
                for rows in reader:
                    spacers[rows[0]] = rows[1]
        except Exception as e:
            logging.error(f"Error loading spacers from {csv_file}: {e}")
            raise
        return spacers

    def load_groups(self) -> Dict[str, Dict]:
        """Load group definitions from a JSON file."""
        try:
            with open(self.map_json, 'r') as f:
                groups = json.load(f)
        except Exception as e:
            logging.error(f"Error loading groups from {self.map_json}: {e}")
            raise
        return groups

    def create_spacer_fasta(self) -> None:
        """Create a FASTA file from the spacer sequences."""
        try:
            with open(self.spacer_fasta, 'w') as f:
                for spacer_id, sequence in self.spacers.items():
                    f.write(f">{spacer_id}\n{sequence}\n")
        except Exception as e:
            logging.error(f"Error creating spacer FASTA file: {e}")
            raise

    def run_blast(self) -> None:
        """Run BLAST to find spacer sequences in the genome."""
        try:
            makeblastdb_cmd = [
                'makeblastdb',
                '-in', self.spacer_fasta,
                '-dbtype', 'nucl'
            ]
            subprocess.run(makeblastdb_cmd, check=True)

            blastn_cmd = [
                'blastn',
                '-query', self.genome_fasta,
                '-db', self.spacer_fasta,
                '-out', self.blast_output,
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
            ]
            subprocess.run(blastn_cmd, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running BLAST: {e}")
            raise

    def parse_blast_results(self) -> Tuple[Dict[str, int], List[List[str]]]:
        """Parse BLAST results and determine the presence of each spacer."""
        presence_matrix = {spacer_id: 0 for spacer_id in self.spacers.keys()}
        blast_results = []

        try:
            with open(self.blast_output, 'r') as infile:
                reader = csv.reader(infile, delimiter='\t')
                for row in reader:
                    blast_results.append(row)
                    spacer_id = row[1]
                    pident = float(row[2])
                    length = int(row[3])
                    if spacer_id in self.spacers:
                        spacer_length = len(self.spacers[spacer_id])
                        coverage = (length / spacer_length) * 100
                        if pident >= self.identity_threshold and coverage >= self.coverage_threshold:
                            presence_matrix[spacer_id] = 1
                            self.total_spacers_found += 1
        except Exception as e:
            logging.error(f"Error parsing BLAST results: {e}")
            raise

        return presence_matrix, blast_results

    def save_presence_matrix(self, presence_matrix: Dict[str, int], output_csv: Path) -> None:
        """Save the presence/absence matrix to a CSV file."""
        try:
            df = pd.DataFrame(list(presence_matrix.items()), columns=['Spacer_ID', 'Presence'])
            df.to_csv(output_csv, index=False)
        except Exception as e:
            logging.error(f"Error saving presence matrix to {output_csv}: {e}")
            raise

    def save_blast_results(self, blast_results: List[List[str]], output_csv: Path) -> None:
        """Save the raw BLAST results to a CSV file."""
        try:
            blast_df = pd.DataFrame(blast_results, columns=[
                'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
            ])
            blast_df.to_csv(output_csv, index=False)
        except Exception as e:
            logging.error(f"Error saving BLAST results to {output_csv}: {e}")
            raise

    def identify_groups(self, presence_matrix: Dict[str, int], crr_type: str) -> Dict[str, Dict]:
        """Identify the presence of groups based on the presence matrix."""
        identified_groups = {}
        if crr_type not in self.groups:
            logging.warning(f"{crr_type} not found in the map JSON.")
            return identified_groups

        for group, subgroups in self.groups[crr_type].items():
            identified_groups[group] = {}
            for subgroup, spacers in subgroups.items():
                total_spacers = len(spacers)
                present_spacers = sum(presence_matrix.get(spacer, 0) for spacer in spacers)
                present_percentage = (present_spacers / total_spacers) * 100

                if self.total_spacers_found > 0:
                    total_spacers_found_percentage = (present_spacers / self.total_spacers_found) * 100
                else:
                    total_spacers_found_percentage = 0

                identified_groups[group][subgroup] = {
                    'present_percentage': present_percentage,
                    'total_spacers': total_spacers,
                    'present_spacers': present_spacers,
                    'total_spacers_found_in_genome_percentage': total_spacers_found_percentage,
                    'missing_spacers': [spacer for spacer in spacers if presence_matrix.get(spacer, 0) == 0]
                }

        # Sort groups by confidence score
        for group in identified_groups:
            identified_groups[group] = dict(sorted(
                identified_groups[group].items(),
                key=lambda item: (item[1]['total_spacers_found_in_genome_percentage'], item[1]['present_percentage']),
                reverse=True
            ))

        return identified_groups

    def save_identified_groups(self, identified_groups: Dict[str, Dict], output_json: Path) -> None:
        """Save the identified groups to a JSON file."""
        try:
            with open(output_json, 'w') as f:
                json.dump(identified_groups, f, indent=4)
        except Exception as e:
            logging.error(f"Error saving identified groups to {output_json}: {e}")
            raise

    def analyze_genome(self) -> None:
        """Perform the full analysis workflow."""
        all_presence_matrices = []
        all_labels = []

        spacer_files = list(self.spacers_folder.glob('*.csv'))
        for spacer_file in tqdm(spacer_files, desc="Processing spacer files"):
            crr_type = spacer_file.stem.split('_')[-1]
            logging.info(f"Processing {crr_type} spacers from {spacer_file}")
            self.spacers = self.load_spacers(spacer_file)

            output_csv = self.output_folder / f"{crr_type}_incidence_matrix.csv"
            blast_output_csv = self.output_folder / f"{crr_type}_incidence_matrix_blast_results.csv"

            self.create_spacer_fasta()

            print(f"Running BLAST for {crr_type}...")
            self.run_blast()

            print(f"Parsing BLAST results for {crr_type}...")
            presence_matrix, blast_results = self.parse_blast_results()

            print(f"Saving presence matrix for {crr_type}...")
            self.save_presence_matrix(presence_matrix, output_csv)

            print(f"Saving BLAST results for {crr_type}...")
            self.save_blast_results(blast_results, blast_output_csv)

            print(f"Identifying groups for {crr_type}...")
            identified_groups = self.identify_groups(presence_matrix, crr_type)

            print(f"Saving identified groups for {crr_type}...")
            output_json = self.output_folder / f"{crr_type}_groups.json"
            self.save_identified_groups(identified_groups, output_json)

            print(f"Analysis complete for {crr_type}.")
            logging.info(f"Analysis complete for {crr_type}.")

            # Collect presence matrices and labels for ML training
            all_presence_matrices.append(list(presence_matrix.values()))
            # Collect the group labels for each genome (assuming each genome has a single group)
            for group, subgroups in identified_groups.items():
                for subgroup, details in subgroups.items():
                    all_labels.append(f"{group}_{subgroup}")

            # Reset total spacers found for the next group
            self.total_spacers_found = 0

            # Clean up temporary files
            self.cleanup_temp_files()

            # Train the ML model
        self.train_ml_model(all_presence_matrices, all_labels)

    def train_ml_model(self, presence_matrices: List[List[int]], labels: List[str]) -> None:
        """Train a machine learning model to predict group/subgroup based on presence matrices."""
        # Ensure all presence matrices have the same length
        max_length = max(len(matrix) for matrix in presence_matrices)
        for matrix in presence_matrices:
            matrix.extend([0] * (max_length - len(matrix)))

        X_train, X_test, y_train, y_test = train_test_split(presence_matrices, labels, test_size=0.2,
                                                            random_state=42)

        model = RandomForestClassifier(n_estimators=100, random_state=42)
        model.fit(X_train, y_train)

        y_pred = model.predict(X_test)
        accuracy = accuracy_score(y_test, y_pred)
        report = classification_report(y_test, y_pred)

        print(f"Model accuracy: {accuracy}")
        print("Classification report:")
        print(report)

        # Save the model to disk
        model_path = self.output_folder / 'crr_ml_model.pkl'
        joblib.dump(model, model_path)
        print(f"Model saved to {model_path}")

    def cleanup_temp_files(self):
        """Clean up temporary files."""
        try:
            Path(self.blast_output).unlink()
            Path(self.spacer_fasta).unlink()
        except Exception as e:
            logging.error(f"Error cleaning up temporary files: {e}")
            raise



if __name__ == '__main__':
    genomes_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/crispr_test/genomes/exact_genomes'
    output_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/crispr_test/output'

    try:
        for genome_file in Path(genomes_folder).glob('*.fasta'):
            genome_fasta = genome_file
            crr_finder = CRRFinder(genome_fasta, output_folder)
            crr_finder.analyze_genome()
    except Exception as e:
        logging.error(f"An error occurred during the analysis: {e}")
        raise

