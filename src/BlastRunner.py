import subprocess
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm
import logging
import os
import csv


class BlastRunner:
    def __init__(self, ref_sequence, genomes_folder, db_folder, results_folder,
                 output_format="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"):
        """
        Initializes the BlastRunner class with the reference sequence, genomes folder, database folder, results folder, and output format.
        """
        self.ref_sequence = ref_sequence
        self.genomes_folder = genomes_folder
        self.db_folder = db_folder
        self.results_folder = results_folder
        self.output_format = output_format
        self.db_name = Path(self.ref_sequence).stem
        self.blast_results = {}

        # Initialize logging in the db_folder
        log_file = Path(self.db_folder) / 'blast_runner.log'
        logging.basicConfig(filename=log_file,
                            filemode='a',
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            level=logging.INFO)
        logging.info(f"BlastRunner initialized. Log file created at {log_file}")

    def create_blast_db(self):
        """
        Creates a BLAST database from the reference sequence.
        """
        try:
            make_db_cmd = [
                'makeblastdb',
                '-in', self.ref_sequence,
                '-dbtype', 'nucl',
                '-out', Path(self.db_folder) / self.db_name
            ]
            subprocess.run(make_db_cmd, check=True)
            logging.info(f"BLAST database created for {self.ref_sequence} at {Path(self.db_folder) / self.db_name}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error in creating BLAST database: {e}")
            raise

    def run_blast(self, genome_file, result_file_name):
        """
        Runs BLAST for the given genome file against the reference sequence database.
        """
        result_file = Path(self.results_folder) / result_file_name
        try:
            blast_cmd = [
                'blastn',
                '-query', genome_file,
                '-db', Path(self.db_folder) / self.db_name,
                '-outfmt', self.output_format,
                '-out', result_file
            ]
            subprocess.run(blast_cmd, check=True)
            logging.info(f"BLAST search completed for {genome_file} with results saved to {result_file}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error in running BLAST search: {e}")
            raise

    def run_blast_on_all_genomes(self):
        """
        Runs BLAST for all genomes in the genomes folder against the reference sequence database.
        """
        try:
            self.create_blast_db()

            # Get a list of all genome files
            genome_files = list(Path(self.genomes_folder).glob('*.fasta'))
            logging.info("Starting BLAST searches for all genomes.")

            # Set up the progress bar
            with tqdm(total=len(genome_files), unit='file', desc='Running BLAST') as pbar:
                for genome_file in genome_files:
                    result_file_name = Path(genome_file).stem + "_blast_results.txt"
                    self.run_blast(genome_file, result_file_name)
                    self.blast_results[genome_file] = result_file_name

                    # Update the progress bar
                    pbar.update(1)
            logging.info("BLAST searches completed for all genomes.")
        except Exception as e:
            logging.error(f"Error in running BLAST on all genomes: {e}")
            raise

    def convert_gb_to_fasta(self, gb_file, fasta_file):
        """
        Converts a GenBank file to FASTA format.
        Args:
        - gb_file: the path to the input GenBank file
        - fasta_file: the path to the output FASTA file
        """
        try:
            with open(gb_file, "r") as input_handle, open(fasta_file, "w") as output_handle:
                sequences = SeqIO.parse(input_handle, "genbank")
                count = SeqIO.write(sequences, output_handle, "fasta")
            logging.info(f"Converted {count} sequences from {gb_file} to FASTA format in {fasta_file}")
        except Exception as e:
            logging.error(f"Error converting {gb_file} to FASTA format: {e}")
            raise

    def process_genomes(self):
        """
        Process all genome files in the genomes folder, converting GenBank files to FASTA format.
        """
        for gb_file in Path(self.genomes_folder).glob('*.gbff'):
            fasta_file = gb_file.with_suffix('.fasta')
            self.convert_gb_to_fasta(gb_file, fasta_file)
        logging.info("Genome files processed.")


class BlastResultsCompiler:
    def __init__(self, directory_path):
        self.directory_path = directory_path

    def compile_results(self, output_csv):
        with open(output_csv, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            # Write the header row with an additional file_name field
            csv_writer.writerow(['file_name', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])

            # Iterate over each file in the directory
            for filename in os.listdir(self.directory_path):
                if filename.endswith('.txt'):
                    # Extract the file name without the '_genomic_blast_results.txt' part
                    file_name = filename.replace('_genomic_blast_results.txt', '')
                    filepath = os.path.join(self.directory_path, filename)
                    with open(filepath, 'r') as file:
                        for line in file:
                            # Split the line by tab since BLAST format 6 is tab-delimited
                            row_data = line.strip().split('\t')
                            # Insert the file name at the beginning of the row data
                            csv_writer.writerow([file_name] + row_data)



def main():

    # PART 1:

    # # Constants
    # REF_SEQ = '/Users/josediogomoura/Documents/BioFago/BioFago/data/input/loci_ref_sequences/fasta/LPS_locus_13genes_FN434113.fasta'
    # GENOMES_FOLDER = '/Users/josediogomoura/Documents/BioFago/BioFago/data/genomes_erwinia_amylovora'
    # DB_FOLDER = '/Users/josediogomoura/Documents/BioFago/BioFago/data/BLAST/LPS_locus_13genes_FN434113/db_ref'
    # RESULTS_FOLDER = '/Users/josediogomoura/Documents/BioFago/BioFago/data/BLAST/LPS_locus_13genes_FN434113/results'
    #
    # # Initialize the BlastRunner
    # blast_runner = BlastRunner(REF_SEQ, GENOMES_FOLDER, DB_FOLDER, RESULTS_FOLDER)
    #
    # # Process the genomes
    # blast_runner.process_genomes()
    #
    # # Run BLAST on all genomes
    # blast_runner.run_blast_on_all_genomes()

    # PART 2:

    # directory_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/BLAST/LPS_locus_13genes_FN434113/results'
    # output_csv = '/Users/josediogomoura/Documents/BioFago/BioFago/data/BLAST/LPS_locus_13genes_FN434113/results/csv_folder/compiled_results.csv'
    #
    # # Create an instance of the class
    # compiler = BlastResultsCompiler(directory_path)
    #
    # # Compile the results into a CSV
    # compiler.compile_results(output_csv)


if __name__ == "__main__":
    main()
