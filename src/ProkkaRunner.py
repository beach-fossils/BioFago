import subprocess
import os
import logging
import multiprocessing

# Set up basic logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class ProkkaAnnotator:
    def __init__(self, fasta_folder, base_output_folder, prokka_path, prokka_options=None):
        self.fasta_folder = fasta_folder
        self.base_output_folder = base_output_folder
        self.prokka_path = prokka_path
        self.prokka_options = prokka_options if prokka_options is not None else {}

    def _build_command(self, fasta_file, output_folder):
        # Build the Prokka command based on the options provided
        cmd = [self.prokka_path, '--outdir', output_folder, '--force']
        for option, value in self.prokka_options.items():
            cmd.extend([option, value])
        cmd.append(fasta_file)
        return cmd

    def _annotate_genome(self, fasta_file):
        # Create a unique output folder for each annotation process
        file_name = os.path.splitext(fasta_file)[0]
        output_folder = os.path.join(self.base_output_folder, file_name)
        os.makedirs(output_folder, exist_ok=True)

        fasta_path = os.path.join(self.fasta_folder, fasta_file)
        cmd = self._build_command(fasta_path, output_folder)
        logging.info(f"Running Prokka: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
            logging.info(result.stdout)
        except subprocess.CalledProcessError as e:
            logging.error(f"Prokka failed on {fasta_file} with error: {e.stderr}")

    def annotate_genomes(self):
        # Get all FASTA files in the given folder
        fasta_files = [f for f in os.listdir(self.fasta_folder) if f.endswith('.fasta')]
        # Use a multiprocessing pool to run annotations in parallel
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            pool.map(self._annotate_genome, fasta_files)


def main():
    fasta_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/extracted_sequences_erw_amy/Cellulose_locus_8genes_NZ_CAPB01000041'
    base_output_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/prokka/Cellulose_locus_8genes_NZ_CAPB01000041'
    prokka_path = '/Users/josediogomoura/miniconda3/envs/prokka_env/bin/prokka'

    # Define any additional options for Prokka here
    prokka_options = {
        '--kingdom': 'Bacteria',
        '--genus': 'Erwinia',
        '--species': 'amylovora'
    }

    # Create an instance of the ProkkaAnnotator class
    annotator = ProkkaAnnotator(fasta_folder, base_output_folder, prokka_path, prokka_options)

    # Run the annotation process in parallel
    annotator.annotate_genomes()


if __name__ == "__main__":
    main()
