import logging
import os
import shutil
import subprocess
import time
from pathlib import Path

import pandas as pd

# Configure logging
logging.basicConfig(
    filename='ani_analysis.log',
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


class FastaStatistics:
    def __init__(self, file_path):
        self.file_path = Path(file_path)
        self.contigs = self._read_fasta()

    def _read_fasta(self):
        contigs = []
        try:
            with self.file_path.open('r') as file:
                contig_name = ''
                contig_sequence = ''
                for line in file:
                    if line.startswith('>'):
                        if contig_sequence:
                            contigs.append((contig_name, contig_sequence))
                            contig_sequence = ''
                        contig_name = line.strip()
                    else:
                        contig_sequence += line.strip()
                if contig_sequence:
                    contigs.append((contig_name, contig_sequence))
        except IOError:
            logging.error(f"Error: The file {self.file_path} could not be opened.")
        return contigs

    def get_number_of_contigs(self):
        return len(self.contigs)

    def get_largest_contig_length(self):
        if not self.contigs:
            return 0, None
        largest_contig = max(self.contigs, key=lambda x: len(x[1]))
        return len(largest_contig[1]), largest_contig[0]

    def get_total_length_of_contigs(self):
        return sum(len(contig[1]) for contig in self.contigs)

    def get_average_contig_length(self):
        num_contigs = self.get_number_of_contigs()
        return self.get_total_length_of_contigs() / num_contigs if num_contigs else 0

    def calculate_gc_content(self):
        gc_count = sum(contig[1].count('G') + contig[1].count('C') for contig in self.contigs)
        total_bases = self.get_total_length_of_contigs()
        return (gc_count / total_bases * 100) if total_bases > 0 else 0

    def calculate_n50(self):
        lengths = sorted((len(contig[1]) for contig in self.contigs), reverse=True)
        total_length = sum(lengths)
        cumulative_sum = 0
        for length in lengths:
            cumulative_sum += length
            if cumulative_sum >= total_length / 2:
                return length
        return 0

    def count_non_standard_bases(self):
        standard_bases = {'A', 'C', 'G', 'T'}
        return sum(1 for _, sequence in self.contigs for base in sequence.upper() if base not in standard_bases)

    def generate_assembly_statistics(self):
        largest_contig_length, largest_contig_name = self.get_largest_contig_length()
        return {
            'Number of Contigs': self.get_number_of_contigs(),
            'Largest Contig Length': largest_contig_length,
            'Largest Contig Name': largest_contig_name,
            'Total Length of Contigs': self.get_total_length_of_contigs(),
            'Average Contig Length': self.get_average_contig_length(),
            'GC Content (%)': self.calculate_gc_content(),
            'N50 Value': self.calculate_n50(),
            'Total Non-standard Bases': self.count_non_standard_bases()
        }

    def output_to_csv(self, output_path):
        stats = self.generate_assembly_statistics()
        df = pd.DataFrame([stats])
        try:
            df.to_csv(output_path, index=False)
            logging.info(f"Assembly statistics saved to {output_path}")
        except IOError:
            logging.error(f"Error: Could not write to file {output_path}")


class OptimizedLocalANIExecutor:
    def __init__(self, single_sequence_path, genomes_directory, results_file, threshold):
        self.single_sequence_path = Path(single_sequence_path)
        self.genomes_directory = Path(genomes_directory)
        self.results_file = Path(results_file)
        self.input_dir = self.single_sequence_path.parent
        self.output_dir = self.input_dir.parent / 'output'
        self.threshold = threshold
        self.match_found = False
        self.best_match = {'species': 'Unknown', 'ani': 0.0}

    def run_pyani_docker(self):
        logging.info("Running ANI analysis using Docker...")
        if self.output_dir.exists():
            shutil.rmtree(self.output_dir)
        time.sleep(1)

        base_dir = self.input_dir.parent.absolute()
        genome_folder_name = self.single_sequence_path.parent.name
        docker_command = (
            f"docker run --rm --platform linux/amd64 -v {base_dir}:/host_dir "
            f"leightonpritchard/average_nucleotide_identity:v0.2.9 "
            f"-i /host_dir/{genome_folder_name} -o /host_dir/output -m ANIm -g --gformat eps"
        )

        try:
            process = subprocess.run(docker_command, shell=True, capture_output=True, text=True, check=True)
            logging.info("Docker ANI calculation completed successfully.")
            return True
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to run docker command: {e.stderr}")
        except Exception as e:
            logging.error(f"Exception during Docker execution: {e}")
        return False

    def process_output(self, reference_genome):
        results_path = self.output_dir / 'ANIm_percentage_identity.tab'
        if not results_path.exists():
            logging.warning(f"Expected result file not found: {results_path}")
            return False

        df = pd.read_csv(results_path, sep='\t', index_col=0)
        write_header = not self.results_file.exists()
        df.to_csv(self.results_file, sep='\t', mode='a', header=write_header)
        logging.info(f"Results appended to {self.results_file}")

        ani_value = df.loc[reference_genome, self.single_sequence_path.stem]
        species_name = self._extract_species_from_filename(reference_genome)

        if ani_value > self.best_match['ani']:
            self.best_match = {'species': species_name, 'ani': ani_value}

        if ani_value >= self.threshold:
            self.match_found = True
            logging.info(f"Match found: {species_name} with ANI value {ani_value:.4f} (>= {self.threshold})")

        return True

    def execute(self):
        if not self.input_dir.exists():
            logging.error(f"Input directory not found: {self.input_dir}")
            return False

        erwinia_amylovora = next(self.genomes_directory.glob('Erwinia_amylovora*.fasta'), None)
        if erwinia_amylovora and self._process_genome(erwinia_amylovora):
            return True

        for genome_file in self.genomes_directory.glob('*.fasta'):
            if genome_file != erwinia_amylovora and self._process_genome(genome_file):
                return True

        logging.info(f"Best match found: {self.best_match['species']} with ANI {self.best_match['ani']:.4f}")
        return True

    def _process_genome(self, genome_file):
        try:
            shutil.copy(genome_file, self.input_dir)
            if self.run_pyani_docker():
                self.process_output(genome_file.stem)
            return self.match_found
        except Exception as e:
            logging.error(f"Error during processing {genome_file.name}: {e}")
        finally:
            (self.input_dir / genome_file.name).unlink(missing_ok=True)
            if self.output_dir.exists():
                shutil.rmtree(self.output_dir)
        return False

    @staticmethod
    def _extract_species_from_filename(filename):
        parts = filename.split('_')
        if parts[0] == "Candidatus":
            return f"Candidatus {parts[1]} {parts[2]}" if len(parts) >= 3 else "Unknown"
        return f"{parts[0]} {parts[1]}" if len(parts) >= 2 else "Unknown"


class NewSpeciesTabModifier:
    def __init__(self, tab_file):
        self.tab_file = Path(tab_file)

    def modify_tab_file(self, output_file):
        try:
            with self.tab_file.open('r') as f, Path(output_file).open('w') as out:
                for line in f:
                    if line.strip() == '':
                        out.write(line)
                        continue
                    values = line.strip().split('\t')
                    if len(values) < 2:
                        logging.warning(f"Invalid line format: {line.strip()}. Skipping...")
                        continue
                    species = self._extract_species_from_filename(values[0])
                    modified_line = f"{species}\t{values[0]}\t" + "\t".join(values[1:]) + "\n"
                    out.write(modified_line)
            logging.info(f"Modified tab file saved to {output_file}")
        except IOError as e:
            logging.error(f"Error modifying tab file: {e}")

    @staticmethod
    def _extract_species_from_filename(filename):
        parts = filename.split('_')
        if parts[0] == "Candidatus":
            return f"Candidatus {parts[1]} {parts[2]}" if len(parts) >= 3 else "Unknown"
        return f"{parts[0]} {parts[1]}" if len(parts) >= 2 else "Unknown"

    def check_species_above_threshold(self, tab_file, threshold=0.95):
        best_match = {'species': 'Unknown', 'ani': 0.0}
        input_genome = None
        try:
            with Path(tab_file).open('r') as f:
                for line in f:
                    if line.strip() == '' or line.startswith('Unknown'):
                        continue
                    values = line.strip().split('\t')
                    if len(values) < 3:
                        logging.warning(f"Invalid line format: {line.strip()}. Skipping...")
                        continue

                    species = self._extract_species_from_filename(values[0])

                    if input_genome is None:
                        input_genome = values[1]
                        continue

                    try:
                        ani_value = float(values[2])
                    except ValueError:
                        logging.warning(f"Invalid ANI value: {values[2]}. Skipping...")
                        continue

                    if ani_value > best_match['ani']:
                        best_match = {'species': species, 'ani': ani_value}

            if best_match['ani'] >= threshold:
                logging.info(f"Best match: {best_match['species']} with ANI {best_match['ani']:.4f}")
                return [best_match['species']]
            else:
                logging.warning(f"No species identified above the threshold. Best match: {best_match['species']} with ANI {best_match['ani']:.4f}")
                return ['Unknown']
        except IOError as e:
            logging.error(f"Error reading tab file: {e}")
            return ['Unknown']


def main_local():
    pass


if __name__ == "__main__":
    main_local()
