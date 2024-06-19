import multiprocessing
import subprocess
import os
import logging
import time
from pathlib import Path

# Set up basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class ProkkaAnnotator:
    def __init__(self, fasta_folder, base_output_folder, prokka_path, custom_db_path, locus_tag_prefix,
                 prokka_options=None, max_cpus=None):
        self.fasta_folder = fasta_folder
        self.base_output_folder = base_output_folder
        self.prokka_path = prokka_path
        self.custom_db_path = custom_db_path  # Path to your GenBank file as custom database
        self.locus_tag_prefix = locus_tag_prefix
        self.prokka_options = prokka_options if prokka_options is not None else {}
        self.max_cpus = max_cpus or multiprocessing.cpu_count()

    def _build_command(self, fasta_file, output_folder):
        # Command setup to include GenBank file via --proteins flag
        cmd = [
            self.prokka_path,
            '--outdir', output_folder,
            '--force',
            '--proteins', self.custom_db_path,  # Using GenBank file as custom database
            '--locustag', self.locus_tag_prefix,
            '--cpus', str(self.max_cpus)
        ]

        # Adding additional Prokka options if any
        for option, value in self.prokka_options.items():
            cmd.extend([option, value])

        cmd.append(fasta_file)
        return cmd

    def _annotate_genome(self, fasta_file):
        file_name = os.path.basename(os.path.splitext(fasta_file)[0])
        output_folder = os.path.join(self.base_output_folder, file_name)
        os.makedirs(output_folder, exist_ok=True)

        fasta_path = os.path.join(self.fasta_folder, fasta_file)
        cmd = self._build_command(fasta_path, output_folder)
        logging.info(f"Running Prokka with custom database: {' '.join(cmd)}")

        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
            logging.info(result.stdout)
            if result.stderr:
                logging.error(f"Prokka generated warnings/errors:\n{result.stderr}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Prokka failed on {fasta_file} with error: {e.stderr}")

    def annotate_genomes(self):
        fasta_files = [f for f in os.listdir(self.fasta_folder) if f.endswith('.fasta')]
        with multiprocessing.Pool(processes=self.max_cpus) as pool:
            pool.map(self._annotate_genome, fasta_files)


class ProdigalTraining:
    def __init__(self, reference_genome_path, output_dir, genetic_code=11):
        """
        Initialize the ProdigalTraining class.

        Parameters:
        - reference_genome_path: Path to the reference genome file in GenBank format.
        - output_dir: Directory where the training file will be saved.
        - genetic_code: Genetic code to use. Default is 11 (Bacterial, Archaeal, and Plant Plastid Code).
        """
        self.reference_genome_path = reference_genome_path
        self.output_dir = output_dir
        self.genetic_code = genetic_code
        self.training_file_path = os.path.join(output_dir, "prodigal_training_file.trn")

    def train_prodigal(self):
        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)
        # Define the output path for the training file
        self.training_file_path = os.path.join(self.output_dir, "prodigal_training_file.trn")
        # Construct the Prodigal command with -p meta for metagenomic data
        command = [
            "prodigal",
            "-i", self.reference_genome_path,
            "-t", self.training_file_path,
            "-c", "-g", str(self.genetic_code),
            "-p", "meta"  # Use meta mode
        ]

        logging.info("Running Prodigal to generate training file...")
        try:
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
            logging.info("Prodigal training completed successfully.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Prodigal training failed: {e.stderr}")
            return False

        return True

    def get_training_file_path(self):
        """
        Return the path to the generated Prodigal training file.
        """
        return self.training_file_path


def train_prodigal():
    train_prodigal_file = ProdigalTraining(
        '/Users/josediogomoura/Documents/BioFago/BioFago/data/genomes/loci_ref_sequence.gb',
        '/Users/josediogomoura/Documents/BioFago/BioFago/data/ApproachFlankGenes/prodigal',
        11)

    train_prodigal_file.train_prodigal()


def run_w_docker():
    fasta_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/rlsA_operon/extracted_sequences'
    base_output_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/rlsA_operon/prokka'
    custom_db_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/rlsA_operon/rlsA.gb'
    locus_tag_prefix = 'rls'
    delay_seconds = 10  # Add a delay of 10 seconds between runs to avoid overloading

    # Additional Prokka options
    prokka_options = {
        '--kingdom': 'Bacteria',
        '--genus': 'Erwinia',
        '--species': 'amylovora',
        '--force': 'true'  # Add force to overwrite existing output directory
    }

    # Ensure base output folder exists
    Path(base_output_folder).mkdir(parents=True, exist_ok=True)

    # Loop through each fasta file in the directory
    for fasta_file in os.listdir(fasta_folder):
        if fasta_file.endswith('.fasta'):
            fasta_path = Path(fasta_folder) / fasta_file
            output_dir = Path(base_output_folder) / fasta_file.rstrip('.fasta')
            output_dir.mkdir(parents=True, exist_ok=True)  # Ensure each output directory is created

            prokka_command = [
                'sudo', 'docker', 'run', '--rm',
                '--platform', 'linux/amd64',
                '-v', f'{fasta_path}:/data/{fasta_file}',
                '-v', f'{output_dir}:/output',
                'staphb/prokka:latest',
                'prokka', f'/data/{fasta_file}', '-o', '/output', '--locustag', locus_tag_prefix
            ]

            # Add custom database option if provided
            if custom_db_path:
                prokka_command.extend(['--proteins', custom_db_path])

            # Add additional options
            for option, value in prokka_options.items():
                prokka_command.extend([option, value])

            # Execute Prokka command
            try:
                subprocess.run(prokka_command, check=True)
                print(f"Prokka annotation completed successfully for {fasta_file}.")
                time.sleep(delay_seconds)  # Delay to avoid overloading
            except subprocess.CalledProcessError as e:
                print(f"Error: Prokka annotation failed for {fasta_file} with exit code {e.returncode}.")
                print(e.stderr)


import os
import shutil

def move_and_rename_gbk_files(source_root, target_folder):
    # Ensure the target folder exists
    os.makedirs(target_folder, exist_ok=True)

    # Walk through the source directory
    for root, dirs, files in os.walk(source_root):
        for file in files:
            if file.endswith('.gbk'):
                # Build the source file path
                source_path = os.path.join(root, file)
                # Extract the folder name and prepare the new file name
                folder_name = os.path.basename(root)
                new_file_name = f"{folder_name}.gbk"
                # Build the target file path
                target_path = os.path.join(target_folder, new_file_name)
                # Move and rename the file
                shutil.move(source_path, target_path)
                print(f"Moved and renamed {file} to {target_path}")

def run_w_dockerv2():
    fasta_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/rpsL_analysis/blast/sequences_rpsL'
    base_output_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/rpsL_analysis/blast/prokka'
    custom_db_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/rpsL_analysis/ref_rpsL/sequence.gb'
    locus_tag_prefix = 'rpsL'

    # Additional Prokka options
    prokka_options = {
        '--kingdom': 'Bacteria',
        '--genus': 'Erwinia',
        '--species': 'amylovora',
        '--force': 'true'
    }

    # Iterate over all fasta files in the directory
    for fasta_file in os.listdir(fasta_folder):
        if fasta_file.endswith(".fasta"):
            file_path = os.path.join(fasta_folder, fasta_file)
            output_folder_name = fasta_file[:-6]  # Removing ".fasta"
            output_path = os.path.join(base_output_folder, output_folder_name)
            Path(output_path).mkdir(parents=True, exist_ok=True)  # Create the output folder

            # Configure the Prokka command for each fasta file
            prokka_command = [
                'docker', 'run', '--rm',
                '--platform', 'linux/amd64',
                '-v', f'{fasta_folder}:/data',
                '-v', f'{output_path}:/output',
                'staphb/prokka:latest',
                'prokka', f'/data/{fasta_file}', '-o', '/output', '--locustag', locus_tag_prefix
            ]

            # # Add custom database option if provided
            if custom_db_path:
                prokka_command.extend(['--proteins', custom_db_path])

            # Add additional options
            for option, value in prokka_options.items():
                prokka_command.extend([option, value])

            # Execute Prokka command
            try:
                subprocess.run(prokka_command, check=True)
                print(f"Prokka annotation for {fasta_file} completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error: Prokka annotation for {fasta_file} failed with exit code {e.returncode}.")
                print(e.stderr)


def main():
    fasta_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/blast_making_sure/extracted_plus5000'
    base_output_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/blast_making_sure/old_prokka'
    prokka_path = '/Users/josediogomoura/miniconda3/envs/new_prokka_env/bin/prokka'
    custom_db_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/genomes/loci_ref_sequences/Capsule_locus_12genes_X77921.gb'  # Path to your GenBank file
    locus_tag_prefix = 'KL'

    # Additional Prokka options
    prokka_options = {
        '--kingdom': 'Bacteria',
        '--genus': 'Erwinia',
        '--species': 'tasmaniensis',
        '--force': 'true'  # Add force to overwrite existing output directory
    }

    annotator = ProkkaAnnotator(fasta_folder, base_output_folder, prokka_path, custom_db_path, locus_tag_prefix,
                                prokka_options)
    annotator.annotate_genomes()


def new_main():
    fasta_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/blast_making_sure/extracted_plus5000'
    base_output_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/blast_making_sure/prokka'
    prokka_path = '/Users/josediogomoura/miniconda3/envs/prokka_env/bin/prokka'
    custom_db_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/genomes/loci_ref_sequences/Capsule_locus_12genes_X77921.gb'
    locus_tag_prefix = 'KL'

    # Additional Prokka options
    prokka_options = {
        '--kingdom': 'Bacteria',
        '--genus': 'Klebsiella',
        '--species': 'amylovora'}

    annotator = ProkkaAnnotator(fasta_folder, base_output_folder, prokka_path, custom_db_path, locus_tag_prefix,
                                prokka_options)
    annotator.annotate_genomes()


if __name__ == '__main__':
    # main()

    #move_and_rename_gbk_files('/Users/josediogomoura/Documents/BioFago/BioFago/data/ALLAproachV2/capsule/prokka',
                             # '/Users/josediogomoura/Documents/BioFago/BioFago/data/ALLAproachV2/capsule/all_gbk_20k')
    run_w_dockerv2()
