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
from typing import Dict, List, Tuple, Any


BASE_PATH = Path(__file__).resolve().parent.parent.parent
REFERENCE_CRISPR_PATH = BASE_PATH / 'reference_crispr'
MAP_JSON = REFERENCE_CRISPR_PATH / 'map.json'
SPACERS_FOLDER = REFERENCE_CRISPR_PATH / 'spacers_csv'


class CRRFinder:
    IDENTITY_THRESHOLD = 95.0
    COVERAGE_THRESHOLD = 90.0
    MATCH_THRESHOLD = 95.0

    def __init__(self, genome_fasta: Path, genomes_folder: Path):
        self.genome_fasta = Path(genome_fasta)
        self.genomes_folder = Path(genomes_folder)
        self.output_folder = self.genomes_folder / 'CRR_finder' / self.genome_fasta.stem
        self.output_folder.mkdir(parents=True, exist_ok=True)
        self.spacers: Dict[str, str] = {}
        self.groups = self._load_groups()
        self.blast_output = tempfile.NamedTemporaryFile(delete=False, suffix='_blast_results.txt').name
        self.spacer_fasta = tempfile.NamedTemporaryFile(delete=False, suffix='_spacers.fasta').name
        self.total_spacers_found = 0
        self._setup_logging()


    def _setup_logging(self) -> None:
        """Setup logging configuration."""
        logging.basicConfig(
            filename=self.output_folder / 'genome_analyzer.log',
            filemode='a',
            format='%(asctime)s - %(levelname)s - %(message)s',
            level=logging.INFO
        )
        logging.info("CRRFinder initialized.")

    def _load_spacers(self, csv_file: Path) -> Dict[str, str]:
        """Load spacer sequences from a CSV file."""
        spacers = {}
        try:
            with open(csv_file, mode='r') as infile:
                reader = csv.reader(infile)
                next(reader)  # skip header
                spacers = {rows[0]: rows[1] for rows in reader}
        except FileNotFoundError:
            logging.error(f"Spacer file not found: {csv_file}")
        except csv.Error as e:
            logging.error(f"Error reading CSV file {csv_file}: {e}")
        return spacers

    def _load_groups(self) -> Dict[str, Any]:
        """Load group definitions from a JSON file."""
        try:
            with open(MAP_JSON, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            logging.error(f"Group definition file not found: {MAP_JSON}")
        except json.JSONDecodeError as e:
            logging.error(f"Error parsing JSON from {MAP_JSON}: {e}")
        return {}

    def _create_spacer_fasta(self) -> None:
        """Create a FASTA file from the spacer sequences."""
        try:
            with open(self.spacer_fasta, 'w') as f:
                for spacer_id, sequence in self.spacers.items():
                    f.write(f">{spacer_id}\n{sequence}\n")
        except IOError as e:
            logging.error(f"Error creating spacer FASTA file: {e}")

    def _run_blast(self) -> None:
        """Run BLAST to find spacer sequences in the genome."""
        try:
            subprocess.run([
                'makeblastdb',
                '-in', self.spacer_fasta,
                '-dbtype', 'nucl'
            ], check=True)

            subprocess.run([
                'blastn',
                '-query', self.genome_fasta,
                '-db', self.spacer_fasta,
                '-out', self.blast_output,
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
            ], check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running BLAST: {e}")
            raise

    def _parse_blast_results(self) -> Tuple[Dict[str, int], List[List[str]]]:
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
                        if pident >= self.IDENTITY_THRESHOLD and coverage >= self.COVERAGE_THRESHOLD:
                            presence_matrix[spacer_id] = 1
                            self.total_spacers_found += 1
                        else:
                            logging.info(
                                f"Spacer {spacer_id} did not meet thresholds: pident={pident}, coverage={coverage}")
        except FileNotFoundError:
            logging.error(f"BLAST output file not found: {self.blast_output}")
        except csv.Error as e:
            logging.error(f"Error parsing BLAST results: {e}")

        return presence_matrix, blast_results

    def _save_presence_matrix(self, presence_matrix: Dict[str, int], output_csv: Path) -> None:
        """Save the presence/absence matrix to a CSV file."""
        try:
            df = pd.DataFrame(list(presence_matrix.items()), columns=['Spacer_ID', 'Presence'])
            df.to_csv(output_csv, index=False)
        except IOError as e:
            logging.error(f"Error saving presence matrix to {output_csv}: {e}")

    def _save_blast_results(self, blast_results: List[List[str]], output_csv: Path) -> None:
        """Save the raw BLAST results to a CSV file."""
        try:
            blast_df = pd.DataFrame(blast_results, columns=[
                'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
            ])
            blast_df.to_csv(output_csv, index=False)
        except IOError as e:
            logging.error(f"Error saving BLAST results to {output_csv}: {e}")

    def _identify_groups(self, presence_matrix: Dict[str, int], crr_type: str) -> Dict[str, Dict]:
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

                total_spacers_found_percentage = (
                                                             present_spacers / self.total_spacers_found) * 100 if self.total_spacers_found > 0 else 0

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

    def _save_identified_groups(self, identified_groups: Dict[str, Dict], output_json: Path) -> None:
        """Save the identified groups to a JSON file."""
        try:
            with open(output_json, 'w') as f:
                json.dump(identified_groups, f, indent=4)
        except IOError as e:
            logging.error(f"Error saving identified groups to {output_json}: {e}")

    def analyze_genome(self) -> None:
        """Perform the full analysis workflow."""
        spacer_files = list(SPACERS_FOLDER.glob('*.csv'))
        for spacer_file in tqdm(spacer_files, desc="Processing spacer files"):
            crr_type = spacer_file.stem.split('_')[-1]
            logging.info(f"Processing {crr_type} spacers from {spacer_file}")
            self.spacers = self._load_spacers(spacer_file)

            output_csv = self.output_folder / f"{crr_type}_incidence_matrix.csv"
            blast_output_csv = self.output_folder / f"{crr_type}_incidence_matrix_blast_results.csv"

            self._create_spacer_fasta()

            print(f"Running BLAST for {crr_type}...")
            self._run_blast()

            print(f"Parsing BLAST results for {crr_type}...")
            presence_matrix, blast_results = self._parse_blast_results()

            print(f"Saving presence matrix for {crr_type}...")
            self._save_presence_matrix(presence_matrix, output_csv)

            print(f"Saving BLAST results for {crr_type}...")
            self._save_blast_results(blast_results, blast_output_csv)

            print(f"Identifying groups for {crr_type}...")
            identified_groups = self._identify_groups(presence_matrix, crr_type)

            print(f"Saving identified groups for {crr_type}...")
            output_json = self.output_folder / f"{crr_type}_groups.json"
            self._save_identified_groups(identified_groups, output_json)

            print(f"Analysis complete for {crr_type}.")
            logging.info(f"Analysis complete for {crr_type}.")

            # Reset total spacers found for the next group
            self.total_spacers_found = 0

            # Clean up temporary files
            self._cleanup_temp_files()


    def _cleanup_temp_files(self) -> None:
        """Clean up temporary files created during the analysis."""
        try:
            for file in [self.spacer_fasta, self.blast_output]:
                if Path(file).exists():
                    Path(file).unlink()
        except OSError as e:
            logging.error(f"Error cleaning up temporary files: {e}")

    def process_directory(self):
        results = []
        csv_folder = SPACERS_FOLDER

        for csv_file in csv_folder.glob('*.csv'):
            crr_type = csv_file.stem.split('_')[-1]
            crr_json_path = self.output_folder / f'{crr_type}_groups.json'
            if crr_json_path.exists():
                try:
                    with open(crr_json_path, 'r') as file:
                        data = json.load(file)
                        best_group, best_type, best_value = self._get_best_group(data)
                        if best_value:
                            result = {
                                "Genome": self.genome_fasta.stem,
                                "CRR Type": crr_type,
                                "Best Group": best_group,
                                "Best Subgroup": best_type,
                                "present_percentage": round(best_value["present_percentage"], 2),
                                "total_spacers": best_value["total_spacers"],
                                "present_spacers": best_value["present_spacers"],
                                "total_spacers_found_in_genome_percentage": round(
                                    best_value["total_spacers_found_in_genome_percentage"], 2),
                                "missing_spacers": ", ".join(best_value["missing_spacers"]),
                                "spacers_in_genome_not_in_subgroup": ", ".join(
                                    best_value.get("spacers_in_genome_not_in_subgroup", [])),
                                "composite_score": round(best_value["composite_score"], 2)
                            }
                            results.append(result)
                except Exception as e:
                    logging.error(f"Error processing file {crr_json_path}: {e}")
                    continue

        # Write results to a CSV file for each genome
        output_csv = self.output_folder / f'{self.genome_fasta.stem}_CRR_best_results.csv'
        fieldnames = ["Genome", "CRR Type", "Best Group", "Best Subgroup", "present_percentage", "total_spacers",
                      "present_spacers", "total_spacers_found_in_genome_percentage",
                      "missing_spacers", "spacers_in_genome_not_in_subgroup", "composite_score"]
        try:
            with open(output_csv, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(results)
        except IOError as e:
            logging.error(f"Error writing results to CSV: {e}")

    def _get_best_type(self, group_data: Dict[str, Any]) -> Tuple[str, Dict[str, Any]]:
        best_type = max(
            group_data.items(),
            key=lambda x: (x[1]["present_percentage"] + x[1]["total_spacers_found_in_genome_percentage"]) / 2
        )
        best_type[1]["composite_score"] = (best_type[1]["present_percentage"] + best_type[1][
            "total_spacers_found_in_genome_percentage"]) / 2
        return best_type

    def _get_best_group(self, data: Dict[str, Any]) -> Tuple[str, str, Dict[str, Any]]:
        best_group = max(
            ((group, self._get_best_type(group_data)) for group, group_data in data.items()),
            key=lambda x: x[1][1]["composite_score"]
        )
        return best_group[0], best_group[1][0], best_group[1][1]

    def get_crr_summary(self) -> str:
        """Summarize the CRR best results for the genome."""
        best_results_file = self.output_folder / f'{self.genome_fasta.stem}_CRR_best_results.csv'
        if best_results_file.exists():
            try:
                df = pd.read_csv(best_results_file)
                crr_info = "CRR: "
                for _, row in df.iterrows():
                    crr_info += f"{row['CRR Type']} - Group: {row['Best Group']}, Subgroup: {row['Best Subgroup']}, Score: {row['composite_score']:.2f}; "
                return crr_info.rstrip('; ')  # Remove trailing semicolon and space
            except Exception as e:
                logging.error(f"Error reading CRR results from {best_results_file}: {e}")
                return "CRR: Error in processing results"
        return "CRR: No results"






def main():
    genome_dir = Path('/Users/josediogomoura/Documents/BioFago/BioFago/data/crispr_test/test_with_all/missing_genomes/GCA_000367565.2_ASM36756v2_genomic')
    genome_file = genome_dir / f"{genome_dir.name}.fasta"

    try:
        crr_finder = CRRFinder(genome_file, genome_dir.parent)
        crr_finder.analyze_genome()
        crr_finder.process_directory()
        crr_info = crr_finder.get_crr_summary()
        print(crr_info)
    except Exception as e:
        logging.error(f"An error occurred during the analysis of {genome_file}: {e}")


if __name__ == '__main__':
    main()
