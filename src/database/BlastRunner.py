import subprocess
from Bio import SeqIO
from tqdm import tqdm
import logging
import os
import csv
import shutil
from pathlib import Path


class BlastRunner:
    def __init__(self, ref_sequence, genomes_folder, db_folder, results_folder,
                 output_format="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"):
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
        try:
            self.create_blast_db()

            # Get a list of all genome files
            genome_files = list(Path(self.genomes_folder).glob('*.fasta'))
            logging.info(f"Found {len(genome_files)} genome files.")

            if not genome_files:
                logging.warning("No genome files found for BLAST searches.")
                return

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
        try:
            with open(gb_file, "r") as input_handle, open(fasta_file, "w") as output_handle:
                sequences = SeqIO.parse(input_handle, "genbank")
                count = SeqIO.write(sequences, output_handle, "fasta")
            logging.info(f"Converted {count} sequences from {gb_file} to FASTA format in {fasta_file}")
        except Exception as e:
            logging.error(f"Error converting {gb_file} to FASTA format: {e}")
            raise

    def process_genomes(self):
        for gb_file in Path(self.genomes_folder).glob('*.gbff'):
            fasta_file = gb_file.with_suffix('.fasta')
            self.convert_gb_to_fasta(gb_file, fasta_file)

        for fna_file in Path(self.genomes_folder).glob('*.fna'):
            fasta_file = fna_file.with_suffix('.fasta')
            shutil.copy(fna_file, fasta_file)

        logging.info("Genome files processed.")


class BlastResultsCompiler:
    def __init__(self, directory_path):
        self.directory_path = directory_path

    def compile_results(self, output_csv):
        with open(output_csv, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            # Write the header row with an additional file_name field
            csv_writer.writerow(
                ['file_name', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
                 'send', 'evalue', 'bitscore'])

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


class LocusIdentifier:
    def __init__(self, ref_genbank, genomes_folder, db_folder, results_folder):
        """
        Initializes the LocusIdentifier class.
        """
        self.ref_genbank = ref_genbank
        self.genomes_folder = genomes_folder
        self.db_folder = db_folder
        self.results_folder = Path(results_folder)
        self.first_gene = None
        self.last_gene = None

    def extract_flanking_genes(self):
        """
        Extracts the sequences of the first and last gene from the reference GenBank file.
        """
        with open(self.ref_genbank, "r") as input_handle:
            for record in SeqIO.parse(input_handle, "genbank"):
                genes = [feature for feature in record.features if feature.type == "gene"]
                if genes:
                    self.first_gene = genes[0]
                    self.last_gene = genes[-1]
                    # Save the first and last gene in FASTA format
                    for gene in [self.first_gene, self.last_gene]:
                        gene_name = gene.qualifiers.get("gene", [""])[0]
                        gene_seq = gene.location.extract(record).seq
                        fasta_filename = Path(self.db_folder) / f"{gene_name}.fasta"
                        with open(fasta_filename, "w") as fasta_file:
                            fasta_file.write(f">{gene_name}\n{gene_seq}\n")
                    logging.info("Flanking gene sequences extracted and saved in FASTA format.")
                else:
                    logging.error("No genes found in the GenBank file.")

    def run_blast_for_flanking_genes(self):
        """
        Runs BLAST for the first and last gene sequences against all genomes in the folder.
        """
        for gene in [self.first_gene, self.last_gene]:
            gene_name = gene.qualifiers.get("gene", [""])[0]
            fasta_path = Path(self.db_folder) / f"{gene_name}.fasta"
            # Create separate result directories for each gene
            gene_results_folder = self.results_folder / gene_name
            gene_results_folder.mkdir(exist_ok=True)

            blast_runner = BlastRunner(ref_sequence=str(fasta_path),
                                       genomes_folder=self.genomes_folder,
                                       db_folder=self.db_folder,
                                       results_folder=str(gene_results_folder))  # Pass the gene-specific results folder
            blast_runner.create_blast_db()
            blast_runner.run_blast_on_all_genomes()


class BlastResultProcessor:
    def __init__(self, first_gene_results_folder, last_gene_results_folder, output_csv):
        self.first_gene_results_folder = Path(first_gene_results_folder)
        self.last_gene_results_folder = Path(last_gene_results_folder)
        self.output_csv = output_csv

    def aggregate_positions(self):
        genome_positions = {}

        def process_file(file_path, gene_type):
            file_name = file_path.name
            with open(file_path, 'r') as file:
                for line in file:
                    parts = line.strip().split()
                    genome_id = parts[0]
                    sstart, send = sorted([int(parts[6]), int(parts[7])])  # Ensure correct order for orientation

                    if genome_id not in genome_positions:
                        genome_positions[genome_id] = {
                            'files': set([file_name]),
                            'range': [sstart, send],
                            'start_id': None,
                            'end_id': None
                        }
                    if gene_type == 'start':
                        genome_positions[genome_id]['start_id'] = genome_id
                    elif gene_type == 'end':
                        genome_positions[genome_id]['end_id'] = genome_id
                    # Update the range irrespective of gene_type
                    current_range = genome_positions[genome_id]['range']
                    genome_positions[genome_id]['range'] = [min(current_range[0], sstart), max(current_range[1], send)]
                    genome_positions[genome_id]['files'].add(file_name)

        # Process results for both genes
        for file_path in self.first_gene_results_folder.glob('*.txt'):
            process_file(file_path, 'start')
        for file_path in self.last_gene_results_folder.glob('*.txt'):
            process_file(file_path, 'end')

        return genome_positions

    def write_to_csv(self, genome_positions):
        with open(self.output_csv, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(['Filenames', 'Start_ID', 'End_ID', 'Region_Start', 'Region_End'])

            for genome_id, info in genome_positions.items():
                filenames = ", ".join(info['files'])
                start_id = info['start_id'] if info['start_id'] else "NA"
                end_id = info['end_id'] if info['end_id'] else "NA"
                csv_writer.writerow([filenames, start_id, end_id, info['range'][0], info['range'][1]])

    def process_results(self):
        genome_positions = self.aggregate_positions()
        self.write_to_csv(genome_positions)


def collect_fasta_files(source_dir, destination_dir):
    """
    Moves all .fna files from the source directory and its subdirectories to the destination directory.

    :param source_dir: The root directory to search for .fna files.
    :param destination_dir: The directory where all .fna files will be moved.
    """
    # Create the destination directory if it doesn't exist
    Path(destination_dir).mkdir(parents=True, exist_ok=True)

    # Walk through the source directory
    for subdir, dirs, files in os.walk(source_dir):
        for filename in files:
            if filename.endswith('.fna'):
                # Construct the full file path
                file_path = Path(subdir) / filename
                # Move the file
                shutil.move(str(file_path), destination_dir)
                print(f"Moved {filename} to {destination_dir}")


def main():
    # PART 1:

    # Constants
    REF_SEQ = '/Users/josediogomoura/Documents/BioFago/BioFago/data/genomes/loci_ref_sequences/fasta/Cellulose_locus_8genes_NZ_CAPB01000041.fasta'

     # collect_fasta_files(
     #    '/Users/josediogomoura/Documents/BioFago/BioFago/data/all_genomes_erwinia_complete/fasta/ncbi_dataset/data',
     #    '/Users/josediogomoura/Documents/BioFago/BioFago/data/all_genomes_erwinia_complete/fasta/ncbi_dataset/all_fasta')

    GENOMES_FOLDER = '/Users/josediogomoura/Documents/BioFago/BioFago/data/all_genomes_erwinia_complete/fasta/ncbi_dataset/all_fasta'
    #
    DB_FOLDER = '/Users/josediogomoura/Documents/BioFago/BioFago/data/BLAST/Cellulose_locus_8genes_NZ_CAPB01000041/all_genomes_complete'
    RESULTS_FOLDER = '/Users/josediogomoura/Documents/BioFago/BioFago/data/BLAST/Cellulose_locus_8genes_NZ_CAPB01000041/all_genomes_complete/results'

    # Initialize the BlastRunner
    blast_runner = BlastRunner(REF_SEQ, GENOMES_FOLDER, DB_FOLDER, RESULTS_FOLDER)

    # Process the genomes
    blast_runner.process_genomes()

    # Run BLAST on all genomes
    blast_runner.run_blast_on_all_genomes()

    # PART 2:
    #
    directory_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/BLAST/Cellulose_locus_8genes_NZ_CAPB01000041/all_genomes_complete/results'
    output_csv = '/Users/josediogomoura/Documents/BioFago/BioFago/data/BLAST/Cellulose_locus_8genes_NZ_CAPB01000041/all_genomes_complete/csv_folder/results_blast.csv_results'

    # Create an instance of the class
    compiler = BlastResultsCompiler(directory_path)

    # Compile the results into a CSV
    compiler.compile_results(output_csv)


def main2():
    # Initialize the class with paths to your reference GenBank file, genomes folder, database folder, and results folder
    # locus_identifier = LocusIdentifier(
    #     ref_genbank="/Users/josediogomoura/Documents/BioFago/BioFago/data/genomes/loci_ref_sequences/Capsule_locus_12genes_X77921.gb",
    #     genomes_folder="/Users/josediogomoura/Documents/BioFago/BioFago/data/genomes_erwinia_amylovora",
    #     db_folder="/Users/josediogomoura/Documents/BioFago/BioFago/data/flanking_genes/db",
    #     results_folder="/Users/josediogomoura/Documents/BioFago/BioFago/data/flanking_genes/results")
    # # Extract flanking genes and convert them to FASTA
    # locus_identifier.extract_flanking_genes()
    # # Run BLAST for the flanking genes
    # locus_identifier.run_blast_for_flanking_genes()

    processor = BlastResultProcessor('/Users/josediogomoura/Documents/BioFago/BioFago/data/flanking_genes/results/amsG',
                                     '/Users/josediogomoura/Documents/BioFago/BioFago/data/flanking_genes/results/amsL',
                                     '/Users/josediogomoura/Documents/BioFago/BioFago/data/extracted_sequences_erw_amy/flanked/capsule/output.csv_results')
    processor.process_results()


if __name__ == "__main__":
    main()
