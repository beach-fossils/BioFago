import subprocess
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import os
import json
import csv
import logging
import tempfile
from typing import Dict, List, Tuple

# Define the base paths for the reference_crispr directory
BASE_PATH = Path(__file__).resolve().parent.parent.parent
REFERENCE_CRISPR_PATH = BASE_PATH / 'reference_crispr'
MAP_JSON = REFERENCE_CRISPR_PATH / 'map.json'
SPACERS_FOLDER = REFERENCE_CRISPR_PATH / 'spacers_csv'


class CRRFinder:
    def __init__(self, genome_fasta: Path, output_folder: str, identity_threshold: float = 95.0,
                 coverage_threshold: float = 90.0, match_threshold: float = 95.0):
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
                        logging.info(
                            f"Spacer {spacer_id}: pident={pident}, length={length}, spacer_length={spacer_length}, coverage={coverage}")
                        if pident >= self.identity_threshold and coverage >= self.coverage_threshold:
                            presence_matrix[spacer_id] = 1
                            self.total_spacers_found += 1
                        else:
                            logging.info(
                                f"Spacer {spacer_id} did not meet thresholds: pident={pident}, coverage={coverage}")
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

        spacers_in_genome = [spacer for spacer, presence in presence_matrix.items() if presence == 1]

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

                spacers_in_subgroup = [spacer for spacer in spacers if presence_matrix.get(spacer, 0) == 1]
                spacers_in_genome_not_in_subgroup = [spacer for spacer in spacers_in_genome if
                                                     spacer not in spacers_in_subgroup]

                identified_groups[group][subgroup] = {
                    'present_percentage': present_percentage,
                    'total_spacers': total_spacers,
                    'present_spacers': present_spacers,
                    'total_spacers_found_in_genome_percentage': total_spacers_found_percentage,
                    'missing_spacers': [spacer for spacer in spacers if presence_matrix.get(spacer, 0) == 0],
                    'spacers_in_genome_not_in_subgroup': spacers_in_genome_not_in_subgroup
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

            # Reset total spacers found for the next group
            self.total_spacers_found = 0

            # Clean up temporary files
            self.cleanup_temp_files()

    def cleanup_temp_files(self) -> None:
        """Clean up temporary files created during the analysis."""
        try:
            if Path(self.spacer_fasta).exists():
                Path(self.spacer_fasta).unlink()
            if Path(self.blast_output).exists():
                Path(self.blast_output).unlink()
        except Exception as e:
            logging.error(f"Error cleaning up temporary files: {e}")
            raise


def get_best_type(group_data):
    best_type = None
    best_value = None
    best_score = -1

    for key, value in group_data.items():
        # Calculate composite score as the average of present_percentage and total_spacers_found_in_genome_percentage
        score = (value["present_percentage"] + value["total_spacers_found_in_genome_percentage"]) / 2
        if score > best_score:
            best_score = score
            best_type = key
            best_value = value
            best_value["composite_score"] = score

    return best_type, best_value


def get_best_group(data):
    best_group = None
    best_type = None
    best_value = None
    best_score = -1

    for group, group_data in data.items():
        type_key, type_value = get_best_type(group_data)
        if type_value and type_value["composite_score"] > best_score:
            best_score = type_value["composite_score"]
            best_group = group
            best_type = type_key
            best_value = type_value

    return best_group, best_type, best_value


def process_directory(root_dir):
    results = []

    for subdir in os.listdir(root_dir):
        subdir_path = os.path.join(root_dir, subdir)
        if os.path.isdir(subdir_path):
            crr1_json_path = os.path.join(subdir_path, 'CRR1_groups.json') # here manually change...
            if os.path.exists(crr1_json_path):
                try:
                    with open(crr1_json_path, 'r') as file:
                        data = json.load(file)
                        best_group, best_type, best_value = get_best_group(data)
                        if best_value:
                            result = {
                                "Genome": subdir,
                                "Best Group": best_group,
                                "Best Type": best_type,
                                "present_percentage": best_value["present_percentage"],
                                "total_spacers": best_value["total_spacers"],
                                "present_spacers": best_value["present_spacers"],
                                "total_spacers_found_in_genome_percentage": best_value[
                                    "total_spacers_found_in_genome_percentage"],
                                "missing_spacers": ", ".join(best_value["missing_spacers"]),
                                "spacers_in_genome_not_in_subgroup": ", ".join(
                                    best_value.get("spacers_in_genome_not_in_subgroup", [])),
                                "composite_score": best_value["composite_score"]
                            }
                            results.append(result)
                except Exception as e:
                    print(f"Error processing file {crr1_json_path}: {e}")
                    continue

    # Write results to a CSV file
    output_csv = os.path.join(root_dir, 'CRR1_best_results.csv') # here manually change...
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ["Genome", "Best Group", "Best Type", "present_percentage", "total_spacers",
                      "present_spacers", "total_spacers_found_in_genome_percentage",
                      "missing_spacers", "spacers_in_genome_not_in_subgroup", "composite_score"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


if __name__ == '__main__':
    genomes_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/crispr_test/genomes/genomesEa'
    output_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/crispr_test/output'

    CRRFinder_folder = output_folder + '/CRR_finder'
    # #
    try:
        for genome_file in Path(genomes_folder).glob('*.fasta'):
            genome_fasta = genome_file
            crr_finder = CRRFinder(genome_fasta, output_folder)
            crr_finder.analyze_genome()
    except Exception as e:
        logging.error(f"An error occurred during the analysis: {e}")
        raise
    process_directory(CRRFinder_folder)
