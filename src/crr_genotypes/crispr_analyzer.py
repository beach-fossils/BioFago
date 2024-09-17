import os
import shutil
from io import TextIOWrapper
from Bio import SeqIO
import logging
import time
import csv
from pathlib import Path
import pandas as pd
from typing import List, Tuple, Dict
from collections import defaultdict


class CRISPRAnalyzer:
    def __init__(self, input_dir: Path, output_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir) / "CRISPR_finder"
        self.results_dir = self.output_dir / "Results"
        self.crr_types = ["CRR1", "CRR2", "CRR4"]
        self.setup_logging()
        self.setup_directories()

    def setup_logging(self):
        log_file = self.output_dir / 'crispr_analysis.log'
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            filename=str(log_file),
                            filemode='w')
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)
        logging.info(f"Logging initialized. Log file: {log_file}")

    def setup_directories(self):
        logging.info("Setting up directories...")
        if self.results_dir.exists():
            shutil.rmtree(self.results_dir)
        self.results_dir.mkdir(parents=True)
        logging.info("Directories setup complete.")

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
        logging.info(f"Finding CRISPR spacers in {genome_file}...")

        file_handlers = {
            'all': open(output_dir / f"{output_prefix}AllCRR.fasta", "w"),
            'crr1': open(output_dir / f"{output_prefix}CRR1.fasta", "w"),
            'crr2': open(output_dir / f"{output_prefix}CRR2.fasta", "w"),
            'crr4': open(output_dir / f"{output_prefix}CRR4.fasta", "w"),
            'csv': open(output_dir / f"{output_prefix}Results.csv", "w"),
            'error': open(output_dir / f"{output_prefix}Error.fasta", "w")
        }

        file_handlers['csv'].write("Strain ID, CRR1 Spacers, CRR2 Spacers, CRR4 Spacers, Total Spacers\n")
        file_handlers['error'].write("Potential assembly errors found in the following crr_genotypes:\n\n")

        for record in SeqIO.parse(genome_file, "fasta"):
            logging.info(f"Processing record: {record.id}, length: {len(record.seq)}")
            crr_counts = defaultdict(int)

            for strand, seq in [("+", str(record.seq)), ("-", str(record.seq.reverse_complement()))]:
                self.process_sequence(seq, record.id, strand, crr_counts, file_handlers)

            total_spacers = sum(crr_counts.values())
            file_handlers['csv'].write(
                f"{record.id}, {crr_counts['CRR1']}, {crr_counts['CRR2']}, {crr_counts['CRR4']}, {total_spacers}\n")
            logging.info(
                f"Spacers found for {record.id}: CRR1 {crr_counts['CRR1']}, CRR2 {crr_counts['CRR2']}, CRR4 {crr_counts['CRR4']}, total {total_spacers}")

        for handler in file_handlers.values():
            handler.close()
        logging.info("CRISPR spacer finding complete.")

    def process_sequence(self, seq: str, record_id: str, strand: str, crr_counts: Dict[str, int],
                         file_handlers: Dict[str, 'TextIOWrapper']):
        for m in range(len(seq) - 29):
            if seq[m:m + 6] == "GTGTTC" or seq[m + 20:m + 28] == "GATAAACC":
                self.process_crr12(seq, m, record_id, strand, crr_counts, file_handlers)
            elif seq[m:m + 6] == "GTTCAC" or seq[m + 10:m + 17] == "GTACGGG":
                self.process_crr4(seq, m, record_id, strand, crr_counts, file_handlers)

    def process_crr12(self, seq: str, m: int, record_id: str, strand: str, crr_counts: Dict[str, int],
                      file_handlers: Dict[str, 'TextIOWrapper']):
        for n in range(100):
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

    def process_crr4(self, seq: str, m: int, record_id: str, strand: str, crr_counts: Dict[str, int],
                     file_handlers: Dict[str, 'TextIOWrapper']):
        for n in range(30, 101):
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

    def run_analysis(self):
        logging.info(f"Starting CRISPR analysis for genomes in {self.input_dir}")
        start_time = time.time()

        for genome_file in self.input_dir.glob('*.fasta'):
            genome_name = genome_file.stem
            logging.info(f"Processing genome: {genome_name}")

            genome_output_dir = self.results_dir / genome_name
            genome_output_dir.mkdir(parents=True, exist_ok=True)

            self.find_crispr_spacers(genome_file, genome_name, genome_output_dir)

        end_time = time.time()
        logging.info(f"CRISPR analysis complete. Total time: {end_time - start_time:.2f} seconds")

    def get_crispr_summary(self) -> Dict[str, str]:
        summaries = {}
        for genome_dir in self.results_dir.iterdir():
            if genome_dir.is_dir():
                genome_name = genome_dir.name
                results_file = genome_dir / f"{genome_name}Results.csv"
                if results_file.exists():
                    try:
                        df = pd.read_csv(results_file, skipinitialspace=True)
                        crr1_total = df['CRR1 Spacers'].sum()
                        crr2_total = df['CRR2 Spacers'].sum()
                        crr4_total = df['CRR4 Spacers'].sum()
                        total_spacers = df['Total Spacers'].sum()
                        summaries[
                            genome_name] = f"CRR1: {crr1_total}, CRR2: {crr2_total}, CRR4: {crr4_total}, Total: {total_spacers}"
                    except Exception as e:
                        logging.error(f"Error processing {results_file}: {e}")
                        summaries[genome_name] = "Error processing results"
                else:
                    summaries[genome_name] = "No results file found"

        return summaries


if __name__ == "__main__":
    genome_dir = '/path/to/your/genome/directory'
    crispr_output_dir = Path(genome_dir).parent / "CRISPR_finder"
    crispr_analyzer = CRISPRAnalyzer(Path(genome_dir), crispr_output_dir)
    crispr_analyzer.run_analysis()
    crispr_info = crispr_analyzer.get_crispr_summary()
    print(crispr_info)



