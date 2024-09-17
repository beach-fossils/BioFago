import subprocess
from pathlib import Path
from typing import List, Dict
from Bio import SeqIO
import logging
import pandas as pd

# Set up logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Dynamically get the reference_plasmids folder
current_script_path = Path(__file__).resolve()
project_root = current_script_path.parent.parent.parent  # Adjust as necessary based on script location
REFERENCE_PLASMIDS = project_root / "reference_plasmids"

class PlasmidFinder:
    def __init__(self, plasmid_reference_folder: Path = REFERENCE_PLASMIDS,
                 identity_threshold: float = 90.0,
                 coverage_threshold: float = 0.8) -> None:

        if 'REFERENCE_PLASMIDS' not in globals():
            logging.error("REFERENCE_PLASMIDS is not defined. Please check the configuration.")
            raise NameError("REFERENCE_PLASMIDS is not defined")
        self.reference_plasmids = REFERENCE_PLASMIDS

        self.plasmid_reference_folder = plasmid_reference_folder
        self.identity_threshold = identity_threshold
        self.coverage_threshold = coverage_threshold
        self.plasmid_lengths: Dict[str, int] = {}
        self.results: List[Dict[str, any]] = []
        self.present_plasmids: List[str] = []


    def combine_plasmid_sequences(self) -> None:
        """Combine all plasmid sequences into a single FASTA file and store lengths."""
        try:
            self.combined_plasmid_fasta = self.db_folder / "combined_plasmids.fasta"
            with self.combined_plasmid_fasta.open('w') as outfile:
                for plasmid_file in self.plasmid_reference_folder.glob('*.fasta'):
                    with plasmid_file.open('r') as infile:
                        seq_records = list(SeqIO.parse(infile, "fasta"))
                        for record in seq_records:
                            record.id = plasmid_file.stem
                            record.description = ""
                            self.plasmid_lengths[record.id] = len(record.seq)
                            SeqIO.write(record, outfile, "fasta")
            logger.info(f"Combined plasmid sequences into {self.combined_plasmid_fasta}")
        except Exception as e:
            logger.error(f"Error combining plasmid sequences: {e}")
            raise

    def create_blast_db(self) -> None:
        """Create a single BLAST database for all combined plasmid sequences."""
        self.blast_db_name = self.db_folder / "combined_plasmids"
        try:
            make_db_cmd = [
                'makeblastdb',
                '-in', str(self.combined_plasmid_fasta),
                '-dbtype', 'nucl',
                '-out', str(self.blast_db_name)
            ]
            subprocess.run(make_db_cmd, check=True, capture_output=True, text=True)
            logger.info(f"BLAST database created at {self.blast_db_name}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error in creating BLAST database: {e.stderr}")
            raise

    def run_blast_for_genome(self, genome_file: Path) -> List[str]:
        """Run BLAST for the genome file against the combined plasmid database."""
        self.genome_file = genome_file
        self.output_folder = genome_file.parent.parent / "plasmid_finder" / self.genome_file.stem
        self.output_folder.mkdir(parents=True, exist_ok=True)
        self.db_folder = self.output_folder / "db"
        self.db_folder.mkdir(exist_ok=True)
        self.results_folder = self.output_folder / "results"
        self.results_folder.mkdir(exist_ok=True)

        self.results = []
        self.present_plasmids = []

        try:
            self.combine_plasmid_sequences()
            self.create_blast_db()

            result_file = self.results_folder / f"{self.genome_file.stem}_blast_results.csv"
            blast_cmd = [
                'blastn',
                '-query', str(self.genome_file),
                '-db', str(self.blast_db_name),
                '-outfmt', "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                '-dust', 'no',
                '-perc_identity', str(self.identity_threshold),
                '-out', str(result_file)
            ]

            logger.info(f"Running BLAST command: {' '.join(blast_cmd)}")
            subprocess.run(blast_cmd, check=True, capture_output=True, text=True)
            logger.info(f"BLAST search completed for {self.genome_file} with results saved to {result_file}")

            self._parse_results(result_file)
            self._extract_present_plasmids()
        except subprocess.CalledProcessError as e:
            logger.error(f"Error in running BLAST search: {e.stderr}")
        except Exception as e:
            logger.error(f"Unexpected error: {e}")

        logger.info(f"Returning present plasmids: {self.present_plasmids}")
        return self.present_plasmids  # Return the list of present plasmids

    def _parse_results(self, result_file: Path) -> None:
        query_length = self._get_sequence_length(self.genome_file)
        aggregated_results = {}

        def merge_intervals(intervals: List[tuple]) -> List[tuple]:
            """Merge overlapping intervals."""
            sorted_intervals = sorted(intervals, key=lambda x: x[0])
            merged = []
            for higher in sorted_intervals:
                if not merged or merged[-1][1] < higher[0]:
                    merged.append(higher)
                else:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], higher[1]))
            return merged

        try:
            df = pd.read_csv(result_file, names=[
                'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
            ])
            logger.info(f"Number of lines in BLAST results: {len(df)}")

            # Add the plasmid_len column
            df['plasmid_len'] = df['sseqid'].map(self.plasmid_lengths)

            # Write the raw results with column names to the CSV file
            df.to_csv(result_file, index=False)

            for _, row in df.iterrows():
                try:
                    plasmid_name = row['sseqid']
                    pident = row['pident']
                    alignment_length = row['length']
                    sstart, send = sorted([row['sstart'], row['send']])
                    if plasmid_name not in aggregated_results:
                        aggregated_results[plasmid_name] = {
                            'intervals': [], 'total_identity': 0, 'alignment_lengths': 0
                        }
                    aggregated_results[plasmid_name]['intervals'].append((sstart, send))
                    aggregated_results[plasmid_name]['total_identity'] += pident * alignment_length
                    aggregated_results[plasmid_name]['alignment_lengths'] += alignment_length
                except ValueError as e:
                    logger.error(f"Error parsing row: {row} -> {e}")

            for plasmid_name, data in aggregated_results.items():
                intervals = merge_intervals(data['intervals'])
                total_coverage_bp = sum(end - start + 1 for start, end in intervals)
                plasmid_length = self.plasmid_lengths.get(plasmid_name, 0)
                if plasmid_length == 0:
                    logger.warning(f"Plasmid length not found for {plasmid_name}")
                    continue
                total_coverage = total_coverage_bp / plasmid_length
                avg_identity = data['total_identity'] / data['alignment_lengths']
                present = "Yes" if avg_identity >= self.identity_threshold and total_coverage >= self.coverage_threshold else "No"
                logger.debug(
                    f"Plasmid: {plasmid_name}, Identity: {avg_identity}, Coverage: {total_coverage}, Present: {present}")
                self.results.append({
                    'plasmid': plasmid_name,
                    'identity': round(avg_identity, 2),
                    'alignment_length': round(total_coverage_bp, 0),
                    'plasmid_length': round(plasmid_length, 0),
                    'coverage': round(total_coverage, 2),
                    'present': present
                })

            # Write processed results to a new file
            processed_result_file = self.results_folder / f"{self.genome_file.stem}_blast_processed_results.csv"
            df_results = pd.DataFrame(self.results)
            df_results.to_csv(processed_result_file, index=False)

            # Extract the present plasmids
            self._extract_present_plasmids()

        except FileNotFoundError:
            logger.error(f"Result file not found: {result_file}")
        except Exception as e:
            logger.error(f"Error parsing results from {result_file}: {e}")

    def _extract_present_plasmids(self) -> None:
        """Extract the list of present plasmids and remove the 'sequence_' prefix."""
        try:
            processed_result_file = self.results_folder / f"{self.genome_file.stem}_blast_processed_results.csv"
            df_results = pd.read_csv(processed_result_file)
            present_df = df_results[df_results['present'] == 'Yes']
            self.present_plasmids = present_df['plasmid'].str.replace('sequence_', '').tolist()
            self.present_plasmids = list(set(self.present_plasmids))  # Remove duplicates
            logger.info(f"Present plasmids: {self.present_plasmids}")
        except Exception as e:
            logger.error(f"Error extracting present plasmids: {e}")

    @staticmethod
    def _get_sequence_length(fasta_file: Path) -> int:
        try:
            with fasta_file.open('r') as f:
                for record in SeqIO.parse(f, "fasta"):
                    return len(record.seq)
        except Exception as e:
            logger.error(f"Error reading sequence length from {fasta_file}: {e}")
        return 0

    def get_results(self) -> List[Dict[str, any]]:
        return self.results

    def get_present_plasmids(self) -> List[str]:
        return self.present_plasmids

    def process_all_genomes_in_folder(self, genomes_folder: Path) -> None:
        """Process all genome files in the given folder."""
        genomes_folder_path = genomes_folder
        for genome_folder in genomes_folder_path.glob('*'):
            genome_file = genome_folder / f"{genome_folder.stem}.fasta"
            if genome_file.exists():
                try:
                    self.run_blast_for_genome(genome_file)
                    self.update_species_csv(genomes_folder_path / "species_finder" / f"{genome_folder.stem}.csv", genome_folder.stem)
                    present_plasmids = self.get_present_plasmids()
                    logger.info(f"Present plasmids for {genome_folder.stem}: {present_plasmids}")
                except Exception as e:
                    logger.error(f"Error processing genome {genome_folder.stem}: {e}")

    def update_species_csv(self, csv_file: Path, genome_name: str) -> None:
        """Update the species CSV file with the plasmid results."""
        if not csv_file.exists():
            logger.error(f"CSV file not found: {csv_file}")
            return
        try:
            df = pd.read_csv(csv_file)
            plasmid_column = ', '.join(self.present_plasmids) if self.present_plasmids else 'None'
            df.loc[df['Name'] == genome_name, 'Plasmid'] = plasmid_column
            df.to_csv(csv_file, index=False)
            logger.info(f"Updated CSV file {csv_file} for genome {genome_name} with plasmids: {plasmid_column}")
        except Exception as e:
            logger.error(f"Error updating CSV file {csv_file}: {e}")


# Example usage
if __name__ == "__main__":
    plasmid_finder = PlasmidFinder()
    plasmid_finder.process_all_genomes_in_folder(Path("/path/to/your/genomes/folder"))
