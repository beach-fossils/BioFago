import fcntl
import shutil
import subprocess
import threading
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import os
import json
import logging
import tempfile
from collections import Counter
import multiprocessing
import csv
from io import StringIO
import re
from typing import Dict, Tuple, Optional, Any, List




BASE_PATH = Path(__file__).resolve().parent.parent.parent
REFERENCE_CRISPR_PATH = BASE_PATH / 'reference_crispr'
MAP_JSON = REFERENCE_CRISPR_PATH / 'updated_crispr.json'
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

        # Create unique temporary directory for this instance
        self.temp_dir = Path(tempfile.mkdtemp(prefix=f'crr_finder_{self.genome_fasta.stem}_'))
        self.blast_output = self.temp_dir / f'{self.genome_fasta.stem}_blast_results.txt'
        self.spacer_fasta = self.temp_dir / f'{self.genome_fasta.stem}_spacers.fasta'

        self.total_spacers_found = 0
        self._setup_logging()
        self._lock = threading.Lock()

    def _setup_logging(self) -> None:
        self.logger = logging.getLogger(f'CRRFinder_{self.genome_fasta.stem}')
        self.logger.setLevel(logging.INFO)

        # Add file handler
        log_file = self.output_folder / f'{self.genome_fasta.stem}_crr_finder.log'
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        self.logger.addHandler(file_handler)

        self.logger.info("CRRFinder initialized.")

    def _create_spacer_fasta(self) -> None:
        """Create a FASTA file from the spacer sequences with proper error handling"""
        try:
            # Ensure the parent directory exists
            self.spacer_fasta.parent.mkdir(parents=True, exist_ok=True)

            with open(self.spacer_fasta, 'w') as f:
                for spacer_id, sequence in self.spacers.items():
                    f.write(f">{spacer_id}\n{sequence}\n")

            self.logger.info(f"Created spacer FASTA file: {self.spacer_fasta}")
        except Exception as e:
            self.logger.error(f"Error creating spacer FASTA file: {e}")
            raise

    def _run_blast(self) -> None:
        """Run BLAST with proper error handling and database creation"""
        try:
            # Create BLAST database
            makeblastdb_cmd = [
                'makeblastdb',
                '-in', str(self.spacer_fasta),
                '-dbtype', 'nucl',
                '-logfile', str(self.output_folder / 'makeblastdb.log')
            ]

            subprocess.run(makeblastdb_cmd, check=True, capture_output=True)

            # Run BLAST search
            blastn_cmd = [
                'blastn',
                '-query', str(self.genome_fasta),
                '-db', str(self.spacer_fasta),
                '-out', str(self.blast_output),
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
            ]

            subprocess.run(blastn_cmd, check=True, capture_output=True)

            self.logger.info(f"BLAST search completed successfully")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Error running BLAST: {e}")
            raise
        except Exception as e:
            self.logger.error(f"Unexpected error in BLAST processing: {e}")
            raise

    def analyze_genome(self) -> None:
        """Perform the full analysis workflow with proper file handling"""
        try:
            spacer_files = list(SPACERS_FOLDER.glob('*.csv'))
            for spacer_file in tqdm(spacer_files, desc="Processing spacer files"):
                try:
                    crr_type = spacer_file.stem.split('_')[-1]
                    self.logger.info(f"Processing {crr_type} spacers from {spacer_file}")

                    with self._lock:
                        self.spacers = self._load_spacers(spacer_file)
                        if not self.spacers:
                            self.logger.warning(f"No spacers loaded from {spacer_file}")
                            continue

                    # Create output paths
                    output_csv = self.output_folder / f"{crr_type}_incidence_matrix.csv"
                    blast_output_csv = self.output_folder / f"{crr_type}_incidence_matrix_blast_results.csv"
                    output_json = self.output_folder / f"{crr_type}_groups.json"

                    # Ensure output directory exists
                    output_csv.parent.mkdir(parents=True, exist_ok=True)

                    self._create_spacer_fasta()
                    self._run_blast()

                    presence_matrix, blast_results = self._parse_blast_results()

                    if presence_matrix and blast_results:
                        self._save_presence_matrix(presence_matrix, output_csv)
                        self._save_blast_results(blast_results, blast_output_csv)

                        identified_groups = self._identify_groups(presence_matrix, crr_type)
                        self._save_identified_groups(identified_groups, output_json)

                    self.logger.info(f"Analysis complete for {crr_type}")

                except Exception as e:
                    self.logger.error(f"Error processing {crr_type}: {str(e)}")
                    continue

                finally:
                    self.total_spacers_found = 0

        except Exception as e:
            self.logger.error(f"Error in analyze_genome: {str(e)}")
            raise
        finally:
            self._cleanup_temp_files()

    def _cleanup_temp_files(self) -> None:
        """Clean up temporary files and directories"""
        try:
            if self.temp_dir.exists():
                shutil.rmtree(self.temp_dir)
                self.logger.info(f"Cleaned up temporary directory: {self.temp_dir}")
        except Exception as e:
            self.logger.error(f"Error cleaning up temporary files: {e}")

    def _save_presence_matrix(self, presence_matrix: Dict[str, int], output_csv: Path) -> None:
        """Save the presence matrix to CSV with error handling"""
        try:
            df = pd.DataFrame(list(presence_matrix.items()), columns=['Spacer_ID', 'Presence'])
            df.to_csv(output_csv, index=False)
            self.logger.info(f"Saved presence matrix to {output_csv}")
        except Exception as e:
            self.logger.error(f"Error saving presence matrix to {output_csv}: {e}")
            raise

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
                data = json.load(f)

            # Recursively strip whitespace from spacer IDs
            def strip_spacers(obj):
                if isinstance(obj, dict):
                    return {k: strip_spacers(v) for k, v in obj.items()}
                elif isinstance(obj, list):
                    return [strip_spacers(item) for item in obj]
                elif isinstance(obj, str):
                    return obj.strip()
                else:
                    return obj

            return strip_spacers(data)
        except FileNotFoundError:
            logging.error(f"Group definition file not found: {MAP_JSON}")
        except json.JSONDecodeError as e:
            logging.error(f"Error parsing JSON from {MAP_JSON}: {e}")
        return {}

    def _parse_blast_results(self) -> Tuple[Dict[str, int], List[List[str]]]:
        presence_matrix = {spacer_id: 0 for spacer_id in self.spacers.keys()}
        blast_results = []
        spacer_counter = Counter()

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
                            spacer_counter[spacer_id] += 1
                        else:
                            logging.info(
                                f"Spacer {spacer_id} did not meet thresholds: pident={pident}, coverage={coverage}")
        except FileNotFoundError:
            logging.error(f"BLAST output file not found: {self.blast_output}")
        except csv.Error as e:
            logging.error(f"Error parsing BLAST results: {e}")

        # Update presence_matrix based on spacer_counter
        for spacer_id, count in spacer_counter.items():
            presence_matrix[spacer_id] = count

        self.total_spacers_found = sum(spacer_counter.values())
        return presence_matrix, blast_results

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
        identified_groups = {}
        if crr_type not in self.groups:
            logging.warning(f"{crr_type} not found in the map JSON.")
            return identified_groups

        # Trim whitespace from keys in presence_matrix
        presence_matrix = {k.strip(): v for k, v in presence_matrix.items()}

        spacers_in_genome = set(spacer for spacer, count in presence_matrix.items() if count > 0)

        def process_group(path, group_data):
            if isinstance(group_data, dict):
                for key, value in group_data.items():
                    process_group(path + [key], value)
            elif isinstance(group_data, list):
                # Reached the spacers list
                subgroup_name = ' > '.join(path)
                spacers = set(spacer.strip() for spacer in group_data)  # Strip whitespace from spacers
                total_spacers = len(spacers)
                present_spacers = len(spacers.intersection(spacers_in_genome))
                present_percentage = (present_spacers / total_spacers) * 100 if total_spacers > 0 else 0
                total_spacers_found_percentage = min(
                    (present_spacers / len(spacers_in_genome)) * 100 if spacers_in_genome else 0, 100)

                missing_spacers = list(spacers - spacers_in_genome)
                spacers_in_genome_not_in_subgroup = list(spacers_in_genome - spacers)

                composite_score = (present_percentage + total_spacers_found_percentage) / 2

                identified_groups[subgroup_name] = {
                    'present_percentage': present_percentage,
                    'total_spacers': total_spacers,
                    'present_spacers': present_spacers,
                    'total_spacers_found_in_genome_percentage': total_spacers_found_percentage,
                    'missing_spacers': missing_spacers,
                    'spacers_in_genome_not_in_subgroup': spacers_in_genome_not_in_subgroup,
                    'composite_score': composite_score
                }
            else:
                logging.warning(f"Unexpected data type in group structure at {'.'.join(path)}: {type(group_data)}")

        process_group([], self.groups[crr_type])

        # Sort the identified groups by composite score
        identified_groups = dict(sorted(
            identified_groups.items(),
            key=lambda item: item[1]['composite_score'],
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

    def _get_best_group(self, data: Dict[str, Any]) -> Tuple[Optional[str], Optional[str], Optional[Dict[str, Any]]]:
        """
        Get the best group from analyzed data.

        Args:
            data: Dictionary containing group analysis data

        Returns:
            Tuple of (best_group, best_subgroup, best_value) where each can be None if no data is available
        """
        if not data:
            return '', '', {}  # Return empty strings and dict instead of None

        try:
            best_group_name, best_value = max(
                data.items(),
                key=lambda item: item[1]['composite_score']
            )

            # Split the subgroup path to extract group and subgroup names
            group_parts = best_group_name.split(' > ')
            best_group = group_parts[0] if group_parts else ''
            best_subgroup = ' > '.join(group_parts[1:]) if len(group_parts) > 1 else ''

            return best_group, best_subgroup, best_value
        except Exception as e:
            self.logger.error(f"Error finding best group: {str(e)}")
            return '', '', {}  # Return empty strings and dict in case of error

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


def traverse_genomes(root_dir: Path) -> List[Tuple[Path, str, str, str]]:
    """
    Traverse the directory structure and yield genome files with their clade information.
    """
    for clade in root_dir.iterdir():
        if clade.is_dir():
            for subclade in clade.iterdir():
                if subclade.is_dir():
                    for subsubclade in subclade.iterdir():
                        if subsubclade.is_dir():
                            for genome_file in subsubclade.glob('*.fn*'):  # Match both .fna and .fasta
                                yield genome_file, clade.name, subclade.name, subsubclade.name
                        else:
                            for genome_file in subclade.glob('*.fn*'):  # Match both .fna and .fasta
                                yield genome_file, clade.name, subclade.name, ''
                else:
                    for genome_file in clade.glob('*.fn*'):  # Match both .fna and .fasta
                        yield genome_file, clade.name, '', ''
def collect_results(root_dir: Path) -> List[Tuple[str, str, str, str, str]]:
    results = []
    for clade in root_dir.iterdir():
        if clade.is_dir():
            for subclade in clade.iterdir():
                if subclade.is_dir():
                    for subsubclade in subclade.iterdir():
                        if subsubclade.is_dir():
                            results.extend(process_directory(subsubclade, clade.name, subclade.name, subsubclade.name))
                        else:
                            results.extend(process_directory(subclade, clade.name, subclade.name, ''))
                else:
                    results.extend(process_directory(clade, clade.name, '', ''))
    return results


def process_directory(directory: Path, clade: str, subclade: str, subsubclade: str) -> List[
    Tuple[str, str, str, str, str]]:
    """
    Process a directory to find CRR results for genome files.

    Args:
        directory: Path object pointing to the directory containing genome files
        clade: Clade information
        subclade: Subclade information
        subsubclade: Sub-subclade information

    Returns:
        List of tuples containing (genome_name, clade, subclade, subsubclade, crr_info)
    """
    results = []
    try:
        # Ensure directory is a Path object
        directory = Path(directory)

        # Find all genome files
        for genome_file in directory.glob('*.fn*'):
            # Construct paths using joinpath or / operator
            crr_finder_dir = directory.joinpath('CRR_finder', genome_file.stem)
            results_file = crr_finder_dir.joinpath(f'{genome_file.stem}_CRR_best_results.csv')

            if results_file.exists():
                try:
                    with results_file.open('r') as f:
                        reader = csv.DictReader(f)
                        crr_info = "CRR: "
                        for row in reader:
                            crr_info += (
                                f"{row.get('CRR Type', 'Unknown')} - "
                                f"Group: {row.get('Best Group', 'Unknown')}, "
                                f"Subgroup: {row.get('Best Subgroup', 'Unknown')}, "
                                f"Score: {row.get('composite_score', '0.0')}; "
                            )
                        crr_info = crr_info.rstrip('; ')

                    results.append((
                        genome_file.stem,
                        str(clade),
                        str(subclade),
                        str(subsubclade),
                        crr_info
                    ))
                except Exception as e:
                    logging.error(f"Error processing {results_file}: {e}")
                    results.append((
                        genome_file.stem,
                        str(clade),
                        str(subclade),
                        str(subsubclade),
                        f"Error: {str(e)}"
                    ))
            else:
                results.append((
                    genome_file.stem,
                    str(clade),
                    str(subclade),
                    str(subsubclade),
                    "CRR: No results"
                ))

    except Exception as e:
        logging.error(f"Error processing directory {directory}: {e}")

    return results



# def process_genome(genome_file, genome_dir):
#     try:
#         logging.info(f"Processing genome: {genome_file.name}")
#
#         # Create a CRRFinder instance for each genome
#         crr_finder = CRRFinder(genome_file, genome_dir)
#
#         # Analyze the genome
#         crr_finder.analyze_genome()
#
#         # Process the results
#         crr_finder.process_directory()
#
#         # Get the CRR summary
#         crr_info = crr_finder.get_crr_summary()
#
#         # Extract genome accession/name (remove .fasta extension)
#         genome_accession = genome_file.stem
#
#         logging.info(f"Finished processing {genome_file.name}")
#         return genome_accession, crr_info
#     except Exception as e:
#         logging.error(f"An error occurred during the analysis of {genome_file}: {e}")
#         return genome_file.stem, f"Error: {str(e)}"

def index_genome_files(root_dir: Path) -> List[Tuple[Path, str, str, str]]:
    """
    Index all genome files and their clade information.
    """
    indexed_files = []
    for clade in root_dir.iterdir():
        if clade.is_dir():
            print(f"Checking clade: {clade.name}")
            for item in clade.iterdir():
                if item.is_dir():
                    print(f"  Checking subclade: {item.name}")
                    for subitem in item.iterdir():
                        if subitem.is_dir():
                            print(f"    Checking subsubclade: {subitem.name}")
                            for genome_file in subitem.glob('*.[ff][na]*'):  # Match .fna and .fasta
                                print(f"      Found file: {genome_file.name}")
                                indexed_files.append((genome_file, clade.name, item.name, subitem.name))
                        else:
                            if subitem.suffix.lower() in ['.fna', '.fasta']:
                                print(f"    Found file: {subitem.name}")
                                indexed_files.append((subitem, clade.name, item.name, ''))
                else:
                    if item.suffix.lower() in ['.fna', '.fasta']:
                        print(f"  Found file: {item.name}")
                        indexed_files.append((item, clade.name, '', ''))

    print(f"Total files indexed: {len(indexed_files)}")
    return indexed_files


def process_genome(genome_file: Path) -> Tuple[str, str]:
    try:
        logging.info(f"Starting to process genome: {genome_file.name}")

        # Create a CRRFinder instance for each genome
        crr_finder = CRRFinder(genome_file, genome_file.parent)
        logging.info(f"CRRFinder instance created for {genome_file.name}")

        # Analyze the genome
        crr_finder.analyze_genome()
        logging.info(f"Genome analysis completed for {genome_file.name}")

        # Process the results
        crr_finder.process_directory()
        logging.info(f"Directory processing completed for {genome_file.name}")

        # Get the CRR summary
        crr_info = crr_finder.get_crr_summary()
        logging.info(f"CRR summary obtained for {genome_file.name}")

        # Extract genome accession/name (remove file extension)
        genome_accession = genome_file.stem

        logging.info(f"Finished processing {genome_file.name}")
        return genome_accession, crr_info
    except Exception as e:
        logging.error(f"An error occurred during the analysis of {genome_file}: {e}", exc_info=True)
        return genome_file.stem, f"Error: {str(e)}"

def analyze_genomes(csv_file_path):
    warnings = []
    exceptions = []

    with open(csv_file_path, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header row

        for row in reader:
            if len(row) < 2:
                continue  # Skip empty rows

            genome = row[0]
            crr_summary = row[1]

            # Extract CRR1 and CRR2 information using regex
            crr1_match = re.search(r'CRR1 - Group: (.*?), Subgroup: .*?, Score: (\d+\.\d+)', crr_summary)
            crr2_match = re.search(r'CRR2 - Group: (.*?), Subgroup: .*?, Score: (\d+\.\d+)', crr_summary)

            if crr1_match and crr2_match:
                crr1_group = crr1_match.group(1)
                crr2_group = crr2_match.group(1)
                crr1_score = float(crr1_match.group(2))
                crr2_score = float(crr2_match.group(2))

                # Check for warnings (not 100% on both CRR1 and CRR2)
                if crr1_score < 100 or crr2_score < 100:
                    warnings.append((genome, f"CRR1: {crr1_score}%, CRR2: {crr2_score}%"))

                # Check for exceptions (different groups)
                if crr1_group != crr2_group:
                    exceptions.append((genome, f"CRR1: {crr1_group}, CRR2: {crr2_group}"))

    return warnings, exceptions



def main():
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Specify the directory containing the FASTA files
    fasta_dir = Path('/Users/josediogomoura/Documents/BioFago/BioFago/test-data/genomes')

    # Get all FASTA files in the directory
    fasta_files = list(fasta_dir.glob('*.[ff][na]*'))  # This will match .fna and .fasta files
    logging.info(f"Found {len(fasta_files)} FASTA files.")

    # Set up multiprocessing
    num_processes = multiprocessing.cpu_count()  # Use all available CPU cores
    pool = multiprocessing.Pool(processes=num_processes)

    # Process genomes in parallel
    logging.info("Processing genomes...")
    results = pool.map(process_genome, fasta_files)

    # Close the pool and wait for all processes to finish
    pool.close()
    pool.join()

    # Write results to CSV
    output_csv = fasta_dir / 'genome_crr_summary.csv'
    try:
        with open(output_csv, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(['Genome', 'CRR Summary'])  # Write header
            csv_writer.writerows(results)
        logging.info(f"Results written to {output_csv}")
    except Exception as e:
        logging.error(f"Error writing to CSV file: {e}")

    # Print results
    for genome, crr_info in results:
        print(f"{genome}: {crr_info}")

    # Usage

    # csv_data = '/Volumes/Crucial_X9/BioFago/data/ApproachFlankGenes/genomes/ErwiniaAmyl/ncbi_dataset/fastas/genome_crr_summary.csv'
    # warnings, exceptions = analyze_genomes(csv_data)
    #
    # print("Warnings (not 100% on both CRR1 and CRR2):")
    # for genome, warning in warnings:
    #     print(f"{genome}: {warning}")
    #
    # print("\nExceptions (different groups for CRR1 and CRR2):")
    # for genome, exception in exceptions:
    #     print(f"{genome}: {exception}")


if __name__ == '__main__':
    main()
