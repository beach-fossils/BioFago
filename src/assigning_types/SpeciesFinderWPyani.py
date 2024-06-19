import json
import os
import shutil
import subprocess
import time
import logging
import pandas as pd
from src.database.RemoteDockerExecutor import RemoteServerManager
import tempfile
from pathlib import Path


# Setup logging
# logging.basicConfig(filename='ani_analysis.log', level=logging.INFO,
#                     format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
#                     datefmt='%Y-%m-%d %H:%M:%S')


class FastaStatistics:
    def __init__(self, file_path):
        self.file_path = file_path
        self.contigs = []  # This will store tuples of (name, sequence)
        self.read_fasta()

    def read_fasta(self):
        """
        Reads a FASTA file and stores contigs in a list along with their names.
        """
        try:
            with open(self.file_path, 'r') as file:
                contig_name = ''
                contig_sequence = ''
                for line in file:
                    if line.startswith('>'):
                        if contig_sequence:
                            self.contigs.append((contig_name, contig_sequence))
                            contig_sequence = ''
                        contig_name = line.strip()  # Store the contig name
                    else:
                        contig_sequence += line.strip()
                if contig_sequence:  # Add the last contig if exists
                    self.contigs.append((contig_name, contig_sequence))
        except IOError:
            print(f"Error: The file {self.file_path} could not be opened.")

    def get_number_of_contigs(self):
        """
        Returns the number of contigs in the FASTA file.
        """
        return len(self.contigs)

    def get_largest_contig_length(self):
        """
        Returns the length of the largest contig and its name.
        """
        if not self.contigs:
            return (0, None)
        largest_contig = max(self.contigs, key=lambda x: len(x[1]))
        return (len(largest_contig[1]), largest_contig[0])

    def get_total_length_of_contigs(self):
        """
        Returns the total length of all contigs.
        """
        return sum(len(contig[1]) for contig in self.contigs)

    def get_average_contig_length(self):
        """
        Returns the average length of contigs.
        """
        total_length = self.get_total_length_of_contigs()
        num_contigs = self.get_number_of_contigs()
        return total_length / num_contigs if num_contigs else 0

    def calculate_gc_content(self):
        """
        Calculates the overall GC content of the contigs.
        """
        gc_count = sum(contig[1].count('G') + contig[1].count('C') for contig in self.contigs)
        total_bases = self.get_total_length_of_contigs()
        return (gc_count / total_bases * 100) if total_bases > 0 else 0

    def calculate_n50(self):
        """
        Calculates the N50 value of the contigs.
        """
        lengths = sorted((len(contig[1]) for contig in self.contigs), reverse=True)
        total_length = sum(lengths)
        cumulative_sum = 0
        for length in lengths:
            cumulative_sum += length
            if cumulative_sum >= total_length / 2:
                return length
        return 0

    def count_non_standard_bases(self):
        """
        Counts non-standard nucleotide bases (anything other than A, C, G, T) in all contigs.
        """
        total_non_standard_count = 0

        # Set of standard nucleotide bases
        standard_bases = {'A', 'C', 'G', 'T'}

        for _, sequence in self.contigs:
            for character in sequence:
                if character.upper() not in standard_bases:
                    total_non_standard_count += 1

        return total_non_standard_count

    def generate_assembly_statistics(self):
        """
        Generates a dictionary with assembly statistics, including a simplified count of non-standard bases.
        """
        largest_contig_length, largest_contig_name = self.get_largest_contig_length()
        non_standard_count = self.count_non_standard_bases()  # This now returns only an int
        stats = {
            'Number of Contigs': self.get_number_of_contigs(),
            'Largest Contig Length': largest_contig_length,
            'Largest Contig Name': largest_contig_name,
            'Total Length of Contigs': self.get_total_length_of_contigs(),
            'Average Contig Length': self.get_average_contig_length(),
            'GC Content (%)': self.calculate_gc_content(),
            'N50 Value': self.calculate_n50(),
            'Total Non-standard Bases': non_standard_count
        }
        return stats

    def output_to_csv(self, output_path):
        """
        Outputs the assembly statistics to a CSV file.
        """
        stats = self.generate_assembly_statistics()
        df = pd.DataFrame([stats])
        try:
            df.to_csv(output_path, index=False)
        except IOError:
            print(f"Error: Could not write to file {output_path}.")


class ReferenceSpecies:
    def __init__(self, folder_w_genomes, jsonl_file, sequence_target_path=None):
        self.folder_w_genomes = folder_w_genomes
        self.jsonl_file = jsonl_file
        self.sequence_target_path = sequence_target_path
        self.genomes = self.get_genomes()
        self.species = self.get_species()
        self.species_dict = self.get_species_dict()

    def get_genomes(self):
        genomes = {}
        for folder in os.listdir(self.folder_w_genomes):
            folder_path = os.path.join(self.folder_w_genomes, folder)
            if os.path.isdir(folder_path):
                for file in os.listdir(folder_path):
                    if file.endswith(".fna"):
                        genomes[folder] = os.path.join(folder_path, file)
        return genomes

    def get_species(self):
        species = {}
        with open(self.jsonl_file, 'r') as file:
            for line in file:
                data = json.loads(line)
                species[data['accession']] = data['organism']['organismName']
        return species

    def get_species_dict(self):
        return {accession: self.species.get(accession, 'Unknown') for accession, filepath in self.genomes.items()}

    def export_species_mapping(self, output_file):
        species_mapping = {}
        for accession, species in self.species.items():
            species_mapping[accession] = species
        with open(output_file, 'w') as f:
            json.dump(species_mapping, f, indent=4)


class ANIAnalyzer:
    def __init__(self, ani_output_path, ref_species, threshold=0.95):
        self.ani_output_path = ani_output_path
        self.ref_species = ref_species  # ReferenceSpecies instance
        self.threshold = threshold

    def analyze_ani_results(self):
        matches = {}
        with open(self.ani_output_path, 'r') as file:
            headers = next(file).strip().split("\t")
            for line in file:
                parts = line.strip().split("\t")
                genome = parts[0]
                ani_values = map(float, parts[1:])
                for target_genome, ani in zip(headers[1:], ani_values):
                    if ani >= self.threshold:
                        matches[genome] = matches.get(genome, []) + [(target_genome, ani)]

        logging.info(f"Found matches for {len(matches)} genomes.")
        return matches

    def display_matches(self):
        matches = self.analyze_ani_results()
        for genome, matched_genomes in matches.items():
            genome_species = self.ref_species.get_species_by_accession(genome)
            logging.info(f"{genome} ({genome_species}) matches with:")
            for target_genome, ani in matched_genomes:
                target_genome_species = self.ref_species.get_species_by_accession(target_genome)
                logging.info(f"  {target_genome} ({target_genome_species}): {ani * 100}%")


def copy_fna_files(source_dir, target_dir, sequence_target_path=None):
    """
    Copy all .fna files from subdirectories within source_dir to target_dir.
    Additionally, copy the genomes sequence file if provided.
    """
    # Ensure genomes directory exists
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Copy genomes
    copied_files = set()
    for root, dirs, files in os.walk(source_dir):
        for file in files:
            if file.endswith('.fna') and file not in copied_files:
                shutil.copyfile(os.path.join(root, file), os.path.join(target_dir, file))
                copied_files.add(file)

    # Optionally copy the genomes sequence
    if sequence_target_path:
        shutil.copyfile(sequence_target_path, os.path.join(target_dir, os.path.basename(sequence_target_path)))
        logging.info(f"Copied genomes sequence to {target_dir}")


def run_docker(local_dir, output_path):
    password = 'l8Wn057sNWsekCz'
    docker_command = (
        "docker run -d -v /home/jmoura/pyani/remote_try_5:/host_dir "
        "leightonpritchard/average_nucleotide_identity:v0.2.9 "
        "-i /host_dir/genomes -o /host_dir/output -m ANIm -g --gformat png,pdf,eps"
    )

    remote_server = RemoteServerManager(hostname='palsson.di.uminho.pt', username='jmoura', password=password)

    try:
        remote_server.copy_files(local_dir, '/home/jmoura/pyani/remote_try_5/genomes')
        container_id = remote_server.run_command(docker_command).strip()
        logging.info(f"Container ID: {container_id}")

        while not remote_server.check_container_status(container_id):
            logging.info("Waiting for the container to finish...")
            time.sleep(60)
        logging.info("Container execution completed.")

        output_dir_remote = '/home/jmoura/pyani/remote_try_5/output'
        # check if the folder is created and has files
        if remote_server.run_command(f"ls {output_dir_remote}"):
            logging.info(f"Output directory created at {output_dir_remote}")

    except Exception as e:
        logging.error(f"Error running Docker command: {str(e)}")
        return None, None

    return container_id, output_dir_remote


def transfer_and_analyze_ani_results(container_id, remote_output_dir, local_output_dir, remote_server, ref_species):
    if container_id is None:
        logging.error("Container ID is None, skipping file transfer and analysis.")
        return

    try:
        remote_server.copy_files_to_local(remote_output_dir, local_output_dir)
        logging.info(f"ANI results transferred to {local_output_dir}")
    except Exception as e:
        logging.error(f"Error transferring ANI results: {str(e)}")
        return

    # Adjusted to consider the additional 'output' folder created by the remote server
    ani_output_path = os.path.join(local_output_dir, 'output', 'ANIm_percentage_identity.tab')

    try:
        if os.path.exists(ani_output_path):
            # Now passing ref_species to ANIAnalyzer
            ani_analyzer = ANIAnalyzer(ani_output_path, ref_species)
            ani_analyzer.display_matches()
        else:
            logging.error(f"ANI output file not found at {ani_output_path}")
    except Exception as e:
        logging.error(f"Error analyzing ANI results: {str(e)}")


class LocalANIExecutor:
    def __init__(self, single_sequence_path, genomes_directory, results_file):
        self.single_sequence_path = Path(single_sequence_path)
        self.genomes_directory = Path(genomes_directory)
        self.results_file = Path(results_file)
        self.input_dir = self.single_sequence_path.parent
        self.output_dir = self.input_dir.parent / 'output'

    def run_pyani_docker(self):
        logging.info("Running ANI analysis using Docker...")
        if self.output_dir.exists():
            shutil.rmtree(self.output_dir)
        time.sleep(1)

        base_dir = self.input_dir.parent.absolute()
        genome_folder_name = self.single_sequence_path.parent.stem
        docker_command = (
            f"docker run --rm --platform linux/amd64 -v {base_dir}:/host_dir "
            f"leightonpritchard/average_nucleotide_identity:v0.2.9 "
            f"-i /host_dir/{genome_folder_name} -o /host_dir/output -m ANIm -g --gformat png,pdf,eps"
        )

        try:
            process = subprocess.run(docker_command, shell=True, capture_output=True, text=True)
            if process.returncode == 0:
                logging.info("Docker ANI calculation completed successfully.")
                return True
            else:
                logging.error(f"Failed to run docker command: {process.stderr}")
                return False
        except Exception as e:
            logging.error(f"Exception during Docker execution: {e}")
            return False

    def process_output(self):
        results_path = self.output_dir / 'ANIm_percentage_identity.tab'
        if results_path.exists():
            df = pd.read_csv(results_path, sep='\t')
            write_header = not self.results_file.exists()
            with self.results_file.open('a') as f:
                df.to_csv(f, sep='\t', index=False, header=write_header)
                f.write("\n")
            logging.info(f"Results appended to {self.results_file}")
            return True
        else:
            logging.warning(f"Expected result file not found: {results_path}")
            return False

    def execute(self):
        if not self.input_dir.exists():
            logging.error(f"Input directory not found: {self.input_dir}")
            return

        existing_paths = list(self.input_dir.glob('*'))
        for path in existing_paths:
            if path.is_file() and path.name != self.single_sequence_path.name:
                try:
                    path.unlink()
                except PermissionError:
                    logging.error(f"Permission denied when trying to delete {path}.")

        for genome_file in self.genomes_directory.glob('*.fasta'):
            try:
                existing_files = list(self.input_dir.glob('*'))
                for file in existing_files:
                    if file.name != self.single_sequence_path.name:
                        file.unlink()

                shutil.copy(genome_file, self.input_dir)
                if self.run_pyani_docker():
                    if not self.process_output():
                        continue
                else:
                    continue

                (self.input_dir / genome_file.name).unlink()
            except Exception as e:
                logging.error(f"Error during processing {genome_file.name}: {e}")
            finally:
                if self.output_dir.exists():
                    shutil.rmtree(self.output_dir)

        try:
            existing_files = list(self.input_dir.glob('*'))
            for file in existing_files:
                if file.name != self.single_sequence_path.name:
                    file.unlink()
            logging.info("Final cleanup of input directory completed.")
        except Exception as e:
            logging.error(f"Error during final cleanup of input directory: {e}")


def modify_fasta_headers(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):

                header = line.split()[0]  # Splits on whitespace and takes the first part
                outfile.write(header + '\n')  # Writes the modified header
            else:
                outfile.write(line)  # Writes the sequence lines unchanged


class SpeciesTabModifier:
    def __init__(self, tab_file):
        self.tab_file = Path(tab_file)

    def modify_tab_file(self, output_file):
        modified_lines = []
        with open(self.tab_file, 'r') as f:
            for line in f:
                if line.strip() == '':  # Skip blank lines
                    modified_lines.append(line)
                    continue
                values = line.strip().split('\t')
                if len(values) < 3:
                    logging.warning(f"Invalid line format: {line.strip()}. Skipping...")
                    continue

                accession = values[0]
                try:
                    species = self.extract_species_from_filename(accession)
                    modified_line = f"{species}\t{accession}\t{values[1]}\t{values[2]}\n"
                    modified_lines.append(modified_line)
                except Exception as e:
                    logging.error(f"Error processing line: {line.strip()}. Error: {str(e)}")
                    continue

        with open(output_file, 'w') as f:
            f.writelines(modified_lines)

    def extract_species_from_filename(self, filename):
        parts = filename.split('_')
        if len(parts) >= 3:
            return '_'.join(parts[:-2])
        return "Unknown"

    def check_species_above_threshold(self, tab_file, threshold=0.95):
        species_list = []
        with open(tab_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.strip() == '' or line.startswith(
                        'Unknown'):  # Skip blank lines and lines starting with 'Unknown'
                    continue
                values = line.strip().split('\t')
                if len(values) < 4:
                    logging.warning(f"Invalid line format: {line.strip()}. Skipping...")
                    continue
                species, accession, value1, value2 = values  # Get the species name and the values from the line
                if float(value1) == 1.0 and float(value2) == 1.0:
                    species_list.append(species)
                elif (float(value1) > threshold and float(value1) != 1.0) or (
                        float(value2) > threshold and float(value2) != 1.0):
                    species_list.append(species)

        if len(species_list) == 0:
            logging.warning("No species identified above the threshold.")
            species_list.append('Unknown')

        # Remove duplicates
        species_list = list(set(species_list))

        # Check if the list has more than one species
        if len(species_list) > 1:
            logging.warning("Multiple species identified above the threshold.")

        return species_list


def move_genomes_to_folder(source_folder, copy_to_folder):
    # Get the list of folders in the source directory
    folders = os.listdir(source_folder)

    # Iterate over each folder in the source directory
    for folder in folders:
        folder_path = os.path.join(source_folder, folder)
        # Check if it's a directory
        if os.path.isdir(folder_path):
            # Find the 'fna' file in the genome folder
            fna_files = [file for file in os.listdir(folder_path) if file.endswith('.fna')]
            if fna_files:
                for fna_file in fna_files:
                    # Copy the 'fna' file to the destination folder
                    shutil.copy(os.path.join(folder_path, fna_file), copy_to_folder)
            else:
                print(f"'fna' file not found in {folder_path}.")

    print("Copying completed.")


def main_move_genomes():
    source_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/Pyani_2/pantoea_all_genomes/data'
    dest_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/Pyani_2/pantoea_all_genomes/genomes'

    move_genomes_to_folder(source_folder, dest_folder)


def main_local():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    stats_output = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/Pyani_2/csv_output/stats.csv'
    folder_w_genomes = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/Pyani_2/pantoea_all_genomes/genomes'
    sequence_target_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/Pyani_2/genomes/PRR_75_modified.fasta'
    results_file = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/Pyani_2/csv_output/results_file.tab'

    fasta_stats(fasta_path=sequence_target_path, output_path=stats_output)

    # Modify the fasta headers with modify_fasta_headers
    # sequence_target_path_modified = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/Pyani_2/genomes/PRR_75_modified.fasta'
    # modify_fasta_headers(sequence_target_path, sequence_target_path_modified)
    #
    # for genome_file in Path(folder_w_genomes).iterdir():
    #     if genome_file.suffix == '.fna':
    #         genome_file_modified = genome_file.parent / (genome_file.stem + '_modified.fasta')
    #         modify_fasta_headers(genome_file, genome_file_modified)

    Path(results_file).touch()

    ani_executor = LocalANIExecutor(sequence_target_path, folder_w_genomes, results_file)
    ani_executor.execute()
    logging.info("Local ANI analysis completed.")

    # results_file = results_file
    species_mapping_file = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/Pyani_2/pantoea_all_genomes/genomes/map/species_mapping.json'
    output_file = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/Pyani_2/csv_output/results_file_species.tab'
    threshold = 0.92
    #
    mapping = SpeciesTabModifier(results_file, species_mapping_file)
    mapping.modify_tab_file(output_file)

    # grab the string from the check_species_above_threshold function to add to the stats_output.csv file a new column with species

    species = mapping.check_species_above_threshold(output_file, threshold)
    print(species)

    if species is not None:
        # add to the stats_output.csv file a new column with species
        df = pd.read_csv(stats_output)
        df['Species'] = species
        df.to_csv(stats_output, index=False)

        logging.info("Species added to the results file.")


def fasta_stats(fasta_path=None, output_path=None):
    if not fasta_path:
        fasta_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/stats/genomes/PRR_75_modified.fasta'
    if not output_path:
        output_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/stats/output/stats.csv'

    fasta_stats = FastaStatistics(file_path=fasta_path)
    fasta_stats.output_to_csv(output_path=output_path)


def main():
    folder_w_genomes = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/Genomes/ncbi_dataset/ncbi_dataset/data'
    jsonl_file = '/data/recent_outputs/Genomes/ncbi_dataset/ncbi_dataset/data/assembly_data_report.jsonl'
    sequence_target_path = '/data/recent_outputs/Pyani/genomes/genomes.fna'
    ref_species = ReferenceSpecies(folder_w_genomes, jsonl_file, sequence_target_path=sequence_target_path)

    copy_to_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/Pyani/genomes'

    output_dir_local = '/data/recent_outputs/Pyani/output'
    # go on level inside the output folder, because it will copy a folder, and we will have the files inside that folder

    copy_fna_files(folder_w_genomes, copy_to_folder, sequence_target_path)

    remote_server = RemoteServerManager(hostname='palsson.di.uminho.pt', username='jmoura', password='l8Wn057sNWsekCz')
    logging.info("Starting ANI analysis...")
    logging.debug(f"Folder with genomes: {folder_w_genomes}")
    logging.debug(f"Species JSONL file: {jsonl_file}")
    container_id, output_dir_remote = run_docker(copy_to_folder, output_dir_local)

    if container_id:
        logging.info(f"Container {container_id} started successfully.")
        transfer_and_analyze_ani_results(container_id, output_dir_remote, output_dir_local, remote_server, ref_species)
    else:
        logging.error("Failed to execute Docker command properly.")


if __name__ == "__main__":
    main_local()
