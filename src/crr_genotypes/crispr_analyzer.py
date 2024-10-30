import shutil
import threading
from io import TextIOWrapper
from Bio import SeqIO
import logging
import time
from pathlib import Path
import pandas as pd
from typing import List, Dict, TextIO
from collections import defaultdict
import os
import multiprocessing

class CRISPRAnalyzer:
    _instances = {}
    _processed_genomes = set()
    _global_lock = threading.Lock()
    _logger = None

    def __new__(cls, input_dir: Path, output_dir: Path, batch_size: int = 20):
        key = str(Path(input_dir).resolve())
        with cls._global_lock:
            if key not in cls._instances:
                instance = super(CRISPRAnalyzer, cls).__new__(cls)
                instance._lock = threading.Lock()  # Instance-specific lock
                cls._instances[key] = instance
            return cls._instances[key]

    def __init__(self, input_dir: Path, output_dir: Path, batch_size: int = 20):
        with self._lock:
            if not hasattr(self, '_initialized'):
                self.input_dir = Path(input_dir).resolve()
                self.output_dir = Path(output_dir).resolve()
                self.results_dir = self.output_dir / "Results"  # Restore results_dir
                self.crr_types = ["CRR1", "CRR2", "CRR4"]
                self.batch_size = batch_size
                self._results_cache = {}
                self._initialized = True
                self.setup_directories()
                self.setup_logging()


    def setup_directories(self):
        """Setup directories with proper error handling"""
        try:
            # Create output directory if it doesn't exist
            self.output_dir.mkdir(parents=True, exist_ok=True)
            self.results_dir.mkdir(parents=True, exist_ok=True)

        except Exception as e:
            if CRISPRAnalyzer._logger:
                CRISPRAnalyzer._logger.error(f"Error setting up directories: {str(e)}")
            raise

    def setup_logging(self):
        if CRISPRAnalyzer._logger is None:
            CRISPRAnalyzer._logger = logging.getLogger('CRISPRAnalyzer')
            CRISPRAnalyzer._logger.setLevel(logging.INFO)

            # Remove any existing handlers
            CRISPRAnalyzer._logger.handlers = []

            # Add console handler
            console_handler = logging.StreamHandler()
            console_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            CRISPRAnalyzer._logger.addHandler(console_handler)

    @staticmethod
    def find_crr12(query: str) -> int:
        crr1 = "GTGTTCCCCGCGTGAGCGGGGATAAACCG"
        crr2 = "GTGTTCCCCGCGTATGCGGGGATAAACCG"
        crr1_match = sum(q == r for q, r in zip(query, crr1))
        crr2_match = sum(q == r for q, r in zip(query, crr2))
        return crr1_match if crr1_match >= crr2_match else -crr2_match

    @staticmethod
    def find_crr4(query: str) -> int:
        crr4 = "GTTCACTGCCGTACAGGCAGCTTAGAAA"
        return sum(q == r for q, r in zip(query, crr4))

    def find_crispr_spacers(self, genome_file: Path, output_prefix: str, output_dir: Path):
        """Find CRISPR spacers with proper file handling"""
        try:
            genome_key = str(genome_file.resolve())

            with self._lock:
                if genome_key in self._processed_genomes:
                    self._logger.debug(f"Skipping already processed genome: {genome_file}")
                    return
                self._processed_genomes.add(genome_key)

            self._logger.info(f"Finding CRISPR spacers in {genome_file}...")

            # Create all required directories
            output_dir.mkdir(parents=True, exist_ok=True)

            # Initialize results file with headers
            results_file = output_dir / f"{output_prefix}Results.csv"
            try:
                # Ensure parent directory exists
                results_file.parent.mkdir(parents=True, exist_ok=True)
                with open(results_file, 'w', newline='') as f:
                    f.write("Strain ID,CRR1 Spacers,CRR2 Spacers,CRR4 Spacers,Total Spacers\n")
            except Exception as e:
                self._logger.error(f"Error creating Results file {results_file}: {e}")
                raise

            file_handlers = {}
            try:
                # Open all file handlers
                file_handlers = {
                    'all': open(output_dir / f"{output_prefix}AllCRR.fasta", "w"),
                    'crr1': open(output_dir / f"{output_prefix}CRR1.fasta", "w"),
                    'crr2': open(output_dir / f"{output_prefix}CRR2.fasta", "w"),
                    'crr4': open(output_dir / f"{output_prefix}CRR4.fasta", "w"),
                    'error': open(output_dir / f"{output_prefix}Error.fasta", "w")
                }

                file_handlers['error'].write("Potential assembly errors found in the following crr_genotypes:\n\n")

                crr_counts = defaultdict(int)
                for record in SeqIO.parse(genome_file, "fasta"):
                    self._logger.info(f"Processing record: {record.id}, length: {len(record.seq)}")

                    for strand, seq in [("+", str(record.seq)), ("-", str(record.seq.reverse_complement()))]:
                        self.process_sequence(seq, record.id, strand, crr_counts, file_handlers)

                # Write results to CSV after processing
                total_spacers = sum(crr_counts.values())
                try:
                    with open(results_file, 'a', newline='') as f:
                        f.write(
                            f"{output_prefix},{crr_counts['CRR1']},{crr_counts['CRR2']},{crr_counts['CRR4']},{total_spacers}\n"
                        )
                except Exception as e:
                    self._logger.error(f"Error writing to Results file {results_file}: {e}")
                    raise

                # Cache the results
                with self._lock:
                    self._results_cache[genome_key] = {
                        'CRR1': crr_counts['CRR1'],
                        'CRR2': crr_counts['CRR2'],
                        'CRR4': crr_counts['CRR4'],
                        'total': total_spacers
                    }

            finally:
                # Close all file handlers
                for handler in file_handlers.values():
                    try:
                        handler.close()
                    except Exception as e:
                        self._logger.error(f"Error closing file handler: {str(e)}")

        except Exception as e:
            self._logger.error(f"Error in find_crispr_spacers for {genome_file}: {str(e)}")
            raise

    def process_sequence(self, seq: str, record_id: str, strand: str,
                         crr_counts: Dict[str, int],
                         file_handlers: Dict[str, TextIO]):
        seq_length = len(seq)
        for m in range(seq_length - 29):
            if m + 28 >= seq_length:
                break
            if seq[m:m + 6] == "GTGTTC" or seq[m + 20:m + 28] == "GATAAACC":
                self.process_crr12(seq, m, record_id, strand, crr_counts, file_handlers)
            elif seq[m:m + 6] == "GTTCAC" or seq[m + 10:m + 17] == "GTACGGG":
                self.process_crr4(seq, m, record_id, strand, crr_counts, file_handlers)

    def process_crr12(self, seq: str, m: int, record_id: str, strand: str,
                      crr_counts: Dict[str, int],
                      file_handlers: Dict[str, TextIO]):
        seq_length = len(seq)
        for n in range(100):
            if m + 88 + n >= seq_length:
                break
            crr12_match = self.find_crr12(seq[m:m + 29])
            if abs(crr12_match) >= 22 and abs(self.find_crr12(seq[m + 59 + n:m + 88 + n])) >= 22:
                spacer = seq[m + 29:m + 59 + n]
                crr_type = 'CRR1' if crr12_match > 0 else 'CRR2'
                if len(spacer) > 35:
                    file_handlers['error'].write(
                        f">{record_id} {crr_type}Spacer{crr_counts[crr_type]} {m + 30}:{m + 59 + n} Strand:{strand}\n{spacer}\n")
                else:
                    crr_counts[crr_type] += 1
                    for key in ['all', crr_type.lower()]:
                        file_handlers[key].write(
                            f">{record_id} {crr_type}Spacer{crr_counts[crr_type]} {m + 30}:{m + 59 + n} Strand:{strand}\n{spacer}\n")
                    logging.debug(
                        f"{crr_type} spacer found: {record_id} {crr_type}Spacer{crr_counts[crr_type]} {m + 30}:{m + 59 + n} Strand:{strand}")
                break

    def process_crr4(self, seq: str, m: int, record_id: str, strand: str,
                     crr_counts: Dict[str, int],
                     file_handlers: Dict[str, TextIO]):
        seq_length = len(seq)
        for n in range(30, 101):
            if m + n + 28 >= seq_length:
                break
            crr4_start = self.find_crr4(seq[m:m + 28])
            crr4_end = self.find_crr4(seq[m + n:m + n + 28])
            if crr4_start >= 20 and crr4_end >= 20:
                spacer = seq[m + 28:m + n]
                if 20 <= len(spacer) <= 50:
                    crr_counts['CRR4'] += 1
                    for key in ['all', 'crr4']:
                        file_handlers[key].write(
                            f">{record_id} CRR4Spacer{crr_counts['CRR4']} {m + 29}:{m + n} Strand:{strand}\n{spacer}\n")
                    logging.info(
                        f"CRR4 spacer found: {record_id} CRR4Spacer{crr_counts['CRR4']} {m + 29}:{m + n} Length:{len(spacer)} Strand:{strand}")
                else:
                    file_handlers['error'].write(
                        f">{record_id} CRR4Spacer{crr_counts['CRR4']} {m + 29}:{m + n} Strand:{strand}\n{spacer}\n")
                    logging.debug(
                        f"CRR4 spacer length out of range: {record_id} CRR4Spacer{crr_counts['CRR4']} {m + 29}:{m + n} Length:{len(spacer)} Strand:{strand}")
                break

    def process_genome_batch(self, genome_batch: List[Path]):
        for genome_file in genome_batch:
            genome_name = genome_file.stem
            logging.info(f"Processing genome: {genome_name}")

            genome_output_dir = self.results_dir / genome_name
            genome_output_dir.mkdir(parents=True, exist_ok=True)

            self.find_crispr_spacers(genome_file, genome_name, genome_output_dir)

    def run_analysis(self):
        """Process a single genome file"""
        with self._lock:
            if not any(self.input_dir.glob('*.fasta')):
                self._logger.warning(f"No FASTA files found in {self.input_dir}")
                return

            for genome_file in self.input_dir.glob('*.fasta'):
                genome_key = str(genome_file.resolve())
                if genome_key not in self._processed_genomes:
                    self._logger.info(f"Processing genome: {genome_file.name}")
                    genome_output_dir = self.results_dir / genome_file.stem
                    genome_output_dir.mkdir(parents=True, exist_ok=True)
                    self.find_crispr_spacers(genome_file, genome_file.stem, genome_output_dir)

    def get_crispr_summary(self) -> Dict[str, str]:
        """Get summary of CRISPR analysis"""
        summaries = {}
        for genome_dir in self.results_dir.iterdir():
            if genome_dir.is_dir():
                genome_name = genome_dir.name
                results_file = genome_dir / f"{genome_name}Results.csv"

                try:
                    if results_file.exists() and results_file.stat().st_size > 0:
                        df = pd.read_csv(results_file, skipinitialspace=True)
                        if not df.empty:
                            row = df.iloc[0]
                            summaries[genome_name] = (
                                f"CRR1: {row['CRR1 Spacers']}, "
                                f"CRR2: {row['CRR2 Spacers']}, "
                                f"CRR4: {row['CRR4 Spacers']}, "
                                f"Total: {row['Total Spacers']}"
                            )
                        else:
                            summaries[genome_name] = "No spacers found"
                    else:
                        summaries[genome_name] = "No results"
                except Exception as e:
                    self._logger.error(f"Error reading results for {genome_name}: {e}")
                    summaries[genome_name] = "Error reading results"

        return summaries

    def clear_cache(self):
        """Clear the cache and processed genomes set"""
        with self._lock:
            self._processed_genomes.clear()
            self._results_cache.clear()



