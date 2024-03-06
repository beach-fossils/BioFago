import multiprocessing
import subprocess
import os
import logging

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
        '/Users/josediogomoura/Documents/BioFago/BioFago/data/input/loci_ref_sequences/Capsule_locus_12genes_X77921.gb',
        '/Users/josediogomoura/Documents/BioFago/BioFago/data/ApproachFlankGenes/prodigal',
        11)

    train_prodigal_file.train_prodigal()


def main():
    fasta_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/ApproachFlankGenes/lps/extracted_seq'
    base_output_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/ApproachFlankGenes/lps/prokka'
    prokka_path = '/Users/josediogomoura/miniconda3/envs/prokka_env/bin/prokka'
    custom_db_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/input/loci_ref_sequences/LPS_locus_13genes_FN434113.gb'  # Path to your GenBank file
    locus_tag_prefix = 'OL'

    # Additional Prokka options
    prokka_options = {
        '--kingdom': 'Bacteria',
        '--genus': 'Erwinia',
        '--species': 'amylovora'
    }

    annotator = ProkkaAnnotator(fasta_folder, base_output_folder, prokka_path, custom_db_path, locus_tag_prefix,
                                prokka_options)
    annotator.annotate_genomes()


if __name__ == '__main__':
    main()
