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
        self.genome_name = self.input_dir.name
        self.output_dir = Path(output_dir) / "CRISPR_finder" / self.genome_name
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
        for crr in self.crr_types + ["AllCRR"]:
            (self.results_dir / crr).mkdir()
        logging.info("Directories setup complete.")

    def combine_fasta_files(self, output_filename: str) -> Path:
        logging.info(f"Combining FASTA files into {output_filename}...")
        combined_fasta = self.results_dir / f"{output_filename}.fasta"
        with combined_fasta.open("w") as outfile:
            for file in self.input_dir.glob("*.fasta"):
                SeqIO.write(SeqIO.parse(file, "fasta"), outfile, "fasta")
        logging.info(f"Combined FASTA file created: {combined_fasta}")
        return combined_fasta

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

    def find_crispr_spacers(self, input_fasta: Path, output_prefix: str):
        logging.info(f"Finding CRISPR spacers in {input_fasta}...")

        file_handlers = {
            'all': open(self.results_dir / f"AllCRR/{output_prefix}AllCRR.fasta", "w"),
            'crr1': open(self.results_dir / f"CRR1/{output_prefix}CRR1.fasta", "w"),
            'crr2': open(self.results_dir / f"CRR2/{output_prefix}CRR2.fasta", "w"),
            'crr4': open(self.results_dir / f"CRR4/{output_prefix}CRR4.fasta", "w"),
            'csv': open(self.results_dir / f"{output_prefix}Results.csv", "w"),
            'error': open(self.results_dir / f"{output_prefix}Error.fasta", "w")
        }

        file_handlers['csv'].write("Strain ID, CRR1 Spacers, CRR2 Spacers, CRR4 Spacers, Total Spacers\n")
        file_handlers['error'].write("Potential assembly errors found in the following crr_genotypes:\n\n")

        check1, check2 = "GTGTTC", "GATAAACC"
        check3, check4 = "GTTCAC", "GTACGGG"

        for record in SeqIO.parse(input_fasta, "fasta"):
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

    @staticmethod
    def match_crr(query: str, reference: str) -> List[int]:
        crr_match = sum(q == r for q, r in zip(query, reference))
        crr_rev_match = sum(q == r for q, r in zip(query[::-1], reference))
        return [crr_match, crr_rev_match]

    def align_spacers(self, input_file: Path, output_prefix: str, crr_type: str):
        logging.info(f"Aligning spacers for {crr_type}...")
        seqObj = SeqIO.parse(input_file, "fasta")
        sampleList = self.group_samples(seqObj)
        consensus = self.generate_consensus(sampleList)

        self.write_alignment_files(sampleList, consensus, output_prefix, crr_type)
        logging.info(f"Spacer alignment for {crr_type} complete.")

    def group_samples(self, seqObj):
        sampleList = []
        preName = "N"
        sample = []

        for record in seqObj:
            if preName in record.id:
                sample.append(str(record.seq))
            else:
                if sample:
                    sample.insert(0, preName)
                    sampleList.append(sample)
                preName = record.id
                sample = [str(record.seq)]

        if sample:
            sample.insert(0, preName)
            sampleList.append(sample)

        return sampleList

    def generate_consensus(self, sampleList):
        consensus = []
        for cassette in sampleList:
            if not consensus:
                consensus = cassette[1:]
            else:
                for d in range(1, len(cassette)):
                    match = False
                    for e in range(len(consensus)):
                        consensusMatch = self.match_crr(cassette[d], consensus[e])
                        if consensusMatch[0] >= (len(consensus[e]) - 2) or consensusMatch[1] >= (len(consensus[e]) - 2):
                            if len(cassette[d]) > len(consensus[e]):
                                consensus[e] = cassette[d]
                            match = True
                            break
                    if not match:
                        consensus.append(cassette[d])
        consensus.insert(0, "Consensus")
        return consensus

    def write_alignment_files(self, sampleList, consensus, output_prefix, crr_type):
        aligned_csv = self.results_dir / f"{crr_type}/{output_prefix}.{crr_type}Aligned.csv"
        aligned_fasta = self.results_dir / f"{crr_type}/{output_prefix}.{crr_type}Aligned.fasta"
        unique_fasta = self.results_dir / f"{crr_type}/{output_prefix}.{crr_type}Unique.fasta"
        binary_phy = self.results_dir / f"{crr_type}/{output_prefix}.{crr_type}Binary.phy"

        with aligned_csv.open("w") as csv_file, aligned_fasta.open("w") as fasta_file, \
                unique_fasta.open("w") as unique_file, binary_phy.open("w") as binary_file:

            csv_file.write(",".join(consensus) + "\n")
            for h in range(1, len(consensus)):
                unique_file.write(f">Consensus{h}\n{consensus[h]}\n")

            fasta_file.write(f">Consensus\n{''.join(consensus[1:])}\n")
            binary_file.write(f"{len(sampleList)} {len(consensus) - 1}\n")

            for i, sample in enumerate(sampleList):
                csv_file.write(f"{sample[0]},")
                fasta_file.write(f">{sample[0]}\n")
                binary_count = 0
                for p in range(3, min(13, len(sample[0]))):
                    binary_file.write(sample[0][p])
                    binary_count += 1
                binary_file.write(" " * (13 - binary_count))

                for k in range(1, len(consensus)):
                    hsp = False
                    for j in range(1, len(sample)):
                        consensusMatch = self.match_crr(sample[j], consensus[k])
                        if consensusMatch[0] >= len(consensus[k]) - 2 or consensusMatch[1] >= len(consensus[k]) - 2:
                            hsp = True
                            csv_file.write(f"{crr_type}S{j},")
                            binary_file.write("1")
                            if len(sample[j]) != len(consensus[k]):
                                if consensusMatch[0] < consensusMatch[1]:
                                    fasta_file.write("-" * (len(consensus[k]) - len(sample[j])))
                                fasta_file.write(sample[j])
                                if consensusMatch[0] > consensusMatch[1]:
                                    fasta_file.write("-" * (len(consensus[k]) - len(sample[j])))
                            else:
                                fasta_file.write(sample[j])
                            break
                    if not hsp:
                        csv_file.write(",")
                        binary_file.write("0")
                        fasta_file.write("-" * len(consensus[k]))
                csv_file.write("\n")
                fasta_file.write("\n")
                binary_file.write("\n")

    def combine_alignments(self, output_prefix: str):
        logging.info("Combining alignments...")
        all_aligned_fasta = self.results_dir / f"AllCRR/{output_prefix}.AllCRRAligned.fasta"
        all_aligned_csv = self.results_dir / f"AllCRR/{output_prefix}.AllCRRAligned.csv"
        all_binary_phy = self.results_dir / f"AllCRR/{output_prefix}.AllCRRBinary.phy"

        self.merge_aligned_fasta_files(output_prefix, all_aligned_fasta)
        self.combine_csv_files(output_prefix, all_aligned_csv)
        self.combine_binary_files(output_prefix, all_binary_phy)

        logging.info("Alignment combination complete.")

    def merge_aligned_fasta_files(self, output_prefix: str, all_aligned_fasta: Path):
        with all_aligned_fasta.open('w') as outfile:
            for crr_type in self.crr_types:
                crr_file = self.results_dir / f"{crr_type}/{output_prefix}.{crr_type}Aligned.fasta"
                if crr_file.exists():
                    SeqIO.write(SeqIO.parse(crr_file, "fasta"), outfile, "fasta")

    def combine_csv_files(self, output_prefix: str, all_aligned_csv: Path):
        spacers = {crr_type: [] for crr_type in self.crr_types}
        for crr_type in self.crr_types:
            csv_file = self.results_dir / f"{crr_type}/{output_prefix}.{crr_type}Aligned.csv"
            if csv_file.exists():
                with csv_file.open('r') as file:
                    spacers[crr_type] = [line.strip().split(',') for line in file]

        max_length = max(len(spacer_list) for spacer_list in spacers.values())
        with all_aligned_csv.open('w') as outfile:
            for m in range(max_length):
                row = []
                for crr_type in self.crr_types:
                    if m < len(spacers[crr_type]):
                        row.extend(spacers[crr_type][m][1:] if row else spacers[crr_type][m])
                outfile.write(','.join(row) + '\n')

    def combine_binary_files(self, output_prefix: str, all_binary_phy: Path):
        total_lines = 0
        total_columns = 0
        binary_content = []
        for crr_type in self.crr_types:
            binary_file = self.results_dir / f"{crr_type}/{output_prefix}.{crr_type}Binary.phy"
            if binary_file.exists():
                with binary_file.open('r') as infile:
                    lines = infile.readlines()
                    if lines:
                        header = lines[0].strip().split()
                        total_lines += int(header[0])
                        total_columns += int(header[1])
                        binary_content.extend(lines[1:])

        with all_binary_phy.open('w') as outfile:
            outfile.write(f"{total_lines} {total_columns}\n")
            outfile.writelines(binary_content)

    def run_analysis(self):
        logging.info(f"Starting CRISPR analysis for {self.genome_name}")
        start_time = time.time()

        combined_fasta = self.combine_fasta_files(self.genome_name)
        self.find_crispr_spacers(combined_fasta, self.genome_name)

        for crr_type in self.crr_types:
            input_file = self.results_dir / f"{crr_type}/{self.genome_name}{crr_type}.fasta"
            if input_file.exists():
                self.align_spacers(input_file, self.genome_name, crr_type)
            else:
                logging.warning(f"Input file for {crr_type} not found: {input_file}")

        self.combine_alignments(self.genome_name)

        end_time = time.time()
        logging.info(f"CRISPR analysis complete. Total time: {end_time - start_time:.2f} seconds")

    def get_crispr_summary(self) -> str:
        results_file = self.results_dir / f"{self.genome_name}Results.csv"

        if not results_file.exists():
            logging.error(f"Results file not found: {results_file}")
            return "CRISPR: No results file found"

        try:
            df = pd.read_csv(results_file, skipinitialspace=True)
            logging.info(f"Successfully read CSV file: {results_file}")
            logging.info(f"CSV columns: {df.columns}")
            logging.info(f"First few rows of the CSV:\n{df.head()}")

            expected_columns = ['CRR1 Spacers', 'CRR2 Spacers', 'CRR4 Spacers', 'Total Spacers']
            missing_columns = [col for col in expected_columns if col not in df.columns]

            if missing_columns:
                logging.error(f"Missing columns in {results_file}: {missing_columns}")
                return f"CRISPR: Error - Missing columns in results file: {', '.join(missing_columns)}"

            crr1_total = df['CRR1 Spacers'].sum()
            crr2_total = df['CRR2 Spacers'].sum()
            crr4_total = df['CRR4 Spacers'].sum()
            total_spacers = df['Total Spacers'].sum()

            logging.info(
                f"Calculated totals: CRR1: {crr1_total}, CRR2: {crr2_total}, CRR4: {crr4_total}, Total: {total_spacers}")

            return f"CRISPR: CRR1: {crr1_total}, CRR2: {crr2_total}, CRR4: {crr4_total}, Total: {total_spacers}"
        except pd.errors.EmptyDataError:
            logging.error(f"The file {results_file} is empty")
            return "CRISPR: Error - Results file is empty"
        except pd.errors.ParserError as e:
            logging.error(f"Error parsing CSV file {results_file}: {e}")
            return "CRISPR: Error - Unable to parse results file"
        except Exception as e:
            logging.error(f"Unexpected error reading CRISPR results from {results_file}: {e}")
            return f"CRISPR: Error in processing results - {str(e)}"


if __name__ == "__main__":
    # genome_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/crispr_test/test2/genomes/genomes'
    # output_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/crispr_test/test2/genomes/genomes/output_michael'
    # analyzer = CRISPRAnalyzer(genome_path, output_path)
    # analyzer.run_analysis("genomes")
    genome_dir = '/Users/josediogomoura/Documents/BioFago/BioFago/data/crispr_test/test_with_all/missing_genomes/GCA_000367565.2_ASM36756v2_genomic'
    crispr_output_dir = Path(genome_dir).parent / "CRISPR_finder"
    crispr_analyzer = CRISPRAnalyzer(Path(genome_dir), crispr_output_dir)
    crispr_analyzer.run_analysis()
    crispr_info = crispr_analyzer.get_crispr_summary()
    print(crispr_info)



