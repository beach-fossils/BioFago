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
        self.output_dir = self.input_dir.parent / 'output'  # Docker's output directory
        self.threshold = threshold
        self.match_found = False
        self.best_match = {'species': 'Unknown', 'ani': 0.0}

    def run_pyani_docker(self):
        logging.info("Running ANI analysis using Docker...")

        # Clean up any existing output directory
        if self.output_dir.exists():
            try:
                shutil.rmtree(self.output_dir)
            except PermissionError:
                logging.warning(f"Permission error cleaning up output directory: {self.output_dir}. Continuing anyway.")
            except Exception as e:
                logging.warning(f"Error cleaning up output directory: {e}. Continuing anyway.")
        time.sleep(1)

        base_dir = self.input_dir.parent.absolute()
        genome_folder_name = self.single_sequence_path.parent.name

        # Get current user ID for Docker to use
        try:
            user_id = os.getuid()
            group_id = os.getgid()
        except AttributeError:
            # Windows doesn't have getuid/getgid
            user_id = 1000
            group_id = 1000
            logging.info("Running on Windows, using default user/group IDs")

        # Split the Docker command into a list of arguments
        # Add user parameter to make files owned by the current user
        docker_command = [
            "docker", "run", "--rm",
            "--platform", "linux/amd64",
            "-v", f"{base_dir}:/host_dir:rw",
            "--user", f"{user_id}:{group_id}",  # Run as current user
            "leightonpritchard/average_nucleotide_identity:v0.2.9",
            "-i", f"/host_dir/{genome_folder_name}",
            "-o", "/host_dir/output",
            "-m", "ANIm", "--workers", "4",
            "-g", "--gformat", "eps"
        ]

        logging.info(f"Executing docker command: {' '.join(docker_command)}")

        try:
            process = subprocess.run(
                docker_command,  # Pass the command as a list
                shell=False,  # Use shell=False
                capture_output=True,
                text=True,
                check=True
            )
            logging.info("Docker ANI calculation completed successfully")
            return True
        except subprocess.CalledProcessError as e:
            logging.error(f"Docker command failed: {e.stderr}")
            return False
        except Exception as e:
            logging.error(f"Exception during Docker execution: {e}")
            return False

    def process_output(self, reference_genome):
        # Look for results in Docker's output directory
        results_path = self.output_dir / 'ANIm_percentage_identity.tab'
        if not results_path.exists():
            logging.warning(f"Expected result file not found: {results_path}")
            return False

        try:
            df = pd.read_csv(results_path, sep='\t', index_col=0)
            ani_value = df.loc[reference_genome, self.single_sequence_path.stem]
            species_name = self._extract_species_from_filename(reference_genome)

            # Update best match regardless of threshold
            if ani_value > self.best_match['ani']:
                self.best_match = {'species': species_name, 'ani': ani_value}
                logging.info(f"New best match: {species_name} with ANI value {ani_value:.4f}")

            # Set match_found if we meet the threshold
            if ani_value >= self.threshold:
                self.match_found = True
                logging.info(f"Threshold met: {species_name} with ANI value {ani_value:.4f}")

            return True

        except Exception as e:
            logging.error(f"Error processing ANI results: {e}")
            return False

    def execute(self):
        if not self.input_dir.exists():
            logging.error(f"Input directory not found: {self.input_dir}")
            return False

        try:
            # First try with Erwinia amylovora if present
            erwinia_amylovora = next(self.genomes_directory.glob('Erwinia_amylovora*.fasta'), None)
            if erwinia_amylovora:
                self._process_genome(erwinia_amylovora)

            # If no match found or ANI below threshold, try other genomes
            if not self.match_found:
                for genome_file in self.genomes_directory.glob('*.fasta'):
                    if genome_file != erwinia_amylovora:
                        self._process_genome(genome_file)
                        if self.match_found:  # Stop if we find a match above threshold
                            break

            logging.info(f"Best match found: {self.best_match['species']} with ANI {self.best_match['ani']:.4f}")
            return True

        except Exception as e:
            logging.error(f"Error during execution: {e}")
            return False

    def _process_genome(self, genome_file):
        try:
            # Copy reference genome to input directory
            shutil.copy(genome_file, self.input_dir)

            if self.run_pyani_docker():
                return self.process_output(genome_file.stem)
            return False

        except Exception as e:
            logging.error(f"Error processing {genome_file.name}: {e}")
            return False
        finally:
            # Clean up with better error handling
            try:
                input_file = self.input_dir / genome_file.name
                if input_file.exists():
                    try:
                        input_file.unlink()
                    except PermissionError:
                        logging.warning(f"Permission error removing input file: {input_file}. Continuing execution.")
                    except Exception as e:
                        logging.warning(f"Error removing input file: {e}. Continuing execution.")
                
                if self.output_dir.exists():
                    try:
                        # Try to remove problematic file first if it exists
                        nucmer_tar = self.output_dir / 'nucmer_output.tar.gz'
                        if nucmer_tar.exists():
                            try:
                                # Try chmod first to ensure we can delete it
                                os.chmod(nucmer_tar, 0o666)
                                nucmer_tar.unlink()
                                logging.info(f"Successfully removed problematic file: {nucmer_tar}")
                            except Exception as e:
                                logging.warning(f"Could not remove problematic file: {nucmer_tar}, {e}")
                        
                        # Now try to remove the directory
                        shutil.rmtree(self.output_dir)
                    except PermissionError:
                        logging.warning(f"Permission error removing output directory: {self.output_dir}. Continuing execution.")
                    except Exception as e:
                        logging.warning(f"Error removing output directory: {e}. Continuing execution.")
            except Exception as e:
                logging.error(f"Error during cleanup: {e}")
                # Continue execution despite cleanup errors

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

                    parts = line.strip().split('\t')
                    if len(parts) < 2:
                        out.write(line)
                        continue

                    # Extract species names
                    species1 = parts[0]
                    species2 = parts[1]

                    # Modify if needed
                    modified_species1 = self._extract_species_from_filename(species1)
                    modified_species2 = self._extract_species_from_filename(species2)

                    # Replace parts and join
                    parts[0] = modified_species1
                    parts[1] = modified_species2
                    out.write('\t'.join(parts) + '\n')
                
            logging.info(f"Modified tab file saved to {output_file}")
            return True
        except Exception as e:
            logging.error(f"Error modifying tab file: {e}")
            return False

    @staticmethod
    def _extract_species_from_filename(filename):
        parts = filename.split('_')
        if parts[0] == "Candidatus":
            return f"Candidatus {parts[1]} {parts[2]}" if len(parts) >= 3 else "Unknown"
        return f"{parts[0]} {parts[1]}" if len(parts) >= 2 else "Unknown"

    def check_species_above_threshold(self, tab_file, threshold=0.95):
        try:
            with Path(tab_file).open('r') as f:
                for line in f:
                    if line.strip() == '' or not '\t' in line:
                        continue

                    parts = line.strip().split('\t')
                    if len(parts) < 3:
                        continue

                    # Extract ANI value
                    try:
                        ani_value = float(parts[2])
                        if ani_value >= threshold:
                            logging.info(f"Found ANI value {ani_value} above threshold {threshold} for {parts[0]} and {parts[1]}")
                            return True, parts[0], parts[1], ani_value
                    except (ValueError, IndexError):
                        continue

            logging.info(f"No ANI values above threshold {threshold} found")
            return False, None, None, 0.0
        except Exception as e:
            logging.error(f"Error checking species above threshold: {e}")
            return False, None, None, 0.0


def main_local():
    # This is a placeholder for local testing
    pass


if __name__ == "__main__":
    main_local()
