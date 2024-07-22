import os
import shutil
from collections import defaultdict
import glob
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from database.BlastRunner import BlastRunner
import logging
from Bio import SeqIO
from pathlib import Path
import csv


# Configure logging at the top of your script
def setup_logging(output_csv):
    log_file = os.path.join(os.path.dirname(output_csv), 'post_blast_output.log')
    logging.basicConfig(level=logging.INFO, filename=log_file, filemode='w',
                        format='%(asctime)s - %(levelname)s - %(message)s')
    return log_file


class FlankGeneExtractor:
    def __init__(self, ref_seq_gb: str):
        self.ref_seq_gb = ref_seq_gb
        self.first_flank_gene = None
        self.last_flank_gene = None
        self.ref_sequence = None  # Attribute to store the reference genome sequence

    def extract_flank_genes(self):
        try:
            with open(self.ref_seq_gb, "r") as input_handle:
                sequences = list(SeqIO.parse(input_handle, "genbank"))
                if not sequences:
                    raise ValueError("No sequences found in the provided GenBank file.")

                # Store the whole sequence of the reference genome
                self.ref_sequence = sequences[0].seq

                # Iterate through features to find the first and last gene
                for feature in sequences[0].features:
                    if feature.type == "gene":
                        if self.first_flank_gene is None:  # The first gene found
                            self.first_flank_gene = feature
                        self.last_flank_gene = feature  # Keep updating last_gene until the end

                # Check if we found any gene at all
                if self.first_flank_gene is None or self.last_flank_gene is None:
                    raise ValueError("The first or last feature is not a gene.")

                logging.info(
                    f"Extracted flank genes: First gene - {self.first_flank_gene.qualifiers['gene'][0]}, Last gene - {self.last_flank_gene.qualifiers['gene'][0]}")
        except Exception as e:
            logging.error(f"Error extracting flank genes: {e}")
            raise

    def write_genes_to_fasta(self, output_file_path):
        try:
            if self.ref_sequence is None:
                raise ValueError("Reference sequence not stored. Cannot extract gene sequences.")

            # Create SeqRecord objects for the first and last genes from the stored attributes
            first_gene_record = SeqRecord(Seq(str(self.first_flank_gene.extract(self.ref_sequence))),
                                          id=self.first_flank_gene.qualifiers['gene'][0],
                                          description="First gene")
            last_gene_record = SeqRecord(Seq(str(self.last_flank_gene.extract(self.ref_sequence))),
                                         id=self.last_flank_gene.qualifiers['gene'][0],
                                         description="Last gene")

            # Write these records to a FASTA file
            with open(output_file_path, "w") as output_handle:
                SeqIO.write([first_gene_record, last_gene_record], output_handle, "fasta")
            logging.info(f"Successfully wrote genes to FASTA file at {output_file_path}")
        except Exception as e:
            logging.error(f"Error writing genes to FASTA: {e}")
            raise


class PostBlastOutput:
    def __init__(self, blast_results_folder, genomes_folder, output_csv, log_file):
        self.blast_results_folder = blast_results_folder
        self.genomes_folder = genomes_folder
        self.output_csv = output_csv
        self.log_file = log_file
        setup_logging(log_file)

    def compile_blast_results_to_csv(self):
        temp_results = defaultdict(list)

        for filename in os.listdir(self.blast_results_folder):
            if filename.endswith("_blast_results.txt"):
                strain_name = filename.split('_blast_results.txt')[0]
                with open(os.path.join(self.blast_results_folder, filename), 'r') as file:
                    for line in file:
                        row = line.strip().split('\t')
                        key = (strain_name, row[0])  # Filename and Query as key
                        temp_results[key].append(row)

        # Filter isolated queries
        filtered_results = {k: v for k, v in temp_results.items() if len(v) > 1}

        with open(self.output_csv, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            header = ['Filename', 'Query', 'Subject', 'Identity', 'Alignment_length', 'Mismatches', 'Gap_opens',
                      'Query_start', 'Query_end', 'Subject_start', 'Subject_end', 'E-value', 'Bit_score']
            csvwriter.writerow(header)

            for key, rows in filtered_results.items():
                strain_name, query = key
                for row in rows:
                    csvwriter.writerow([strain_name] + row)

    def extract_regions_from_genomes(self, extracted_folder):
        strain_query_positions = defaultdict(lambda: {'min_start': float('inf'), 'max_end': -1})

        with open(self.output_csv, 'r') as csvfile:
            csvreader = csv.DictReader(csvfile)
            for row in csvreader:
                filename_query_key = (row['Filename'], row['Query'])
                start = int(row['Query_start'])
                end = int(row['Query_end'])

                strain_query_positions[filename_query_key]['min_start'] = min(
                    strain_query_positions[filename_query_key]['min_start'], start)
                strain_query_positions[filename_query_key]['max_end'] = max(
                    strain_query_positions[filename_query_key]['max_end'], end)

        Path(extracted_folder).mkdir(parents=True, exist_ok=True)

        for (filename, query), positions in strain_query_positions.items():
            min_start = positions['min_start']
            max_end = positions['max_end']

            fasta_filename = f"{filename}.fasta"
            fasta_file_path = Path(self.genomes_folder) / fasta_filename

            if fasta_file_path.is_file():
                for record in SeqIO.parse(fasta_file_path, "fasta"):
                    if record.id == query:
                        extracted_sequence = record.seq[min_start - 1:max_end]
                        output_file = Path(extracted_folder) / f"{filename}_{query}_extracted_sequence.fasta"
                        with open(output_file, "w") as output_handle:
                            SeqIO.write([SeqIO.SeqRecord(extracted_sequence, id=record.id, description="")],
                                        output_handle, "fasta")
                            logging.info(
                                f"Extracted sequence for {query} in {filename}, min start: {min_start}, max end: {max_end}, length: {len(extracted_sequence)}")
            else:
                logging.warning(f"No FASTA file found for {filename}. Expected path: {fasta_file_path}")

    def extract_regions_from_genomes_v2(self, extracted_folder):
        strain_query_positions = defaultdict(lambda: {'min_start': float('inf'), 'max_end': -1})

        with open(self.output_csv, 'r') as csvfile:
            csvreader = csv.DictReader(csvfile)
            for row in csvreader:
                filename_query_key = (row['Filename'], row['Query'])
                start = int(row['Query_start'])
                end = int(row['Query_end'])

                # Update the min_start and max_end to potentially include -20 and +20 bps around the region
                strain_query_positions[filename_query_key]['min_start'] = min(
                    strain_query_positions[filename_query_key]['min_start'], start)
                strain_query_positions[filename_query_key]['max_end'] = max(
                    strain_query_positions[filename_query_key]['max_end'], end)

        Path(extracted_folder).mkdir(parents=True, exist_ok=True)

        for (filename, query), positions in strain_query_positions.items():
            min_start = max(1, positions['min_start'] - 200)  # Ensure we do not go below the start of the sequence
            max_end = positions['max_end'] + 200  #

            fasta_filename = f"{filename}.fasta"
            fasta_file_path = Path(self.genomes_folder) / fasta_filename

            if fasta_file_path.is_file():
                for record in SeqIO.parse(fasta_file_path, "fasta"):
                    if record.id == query:
                        # Adjust max_end to not exceed the sequence length
                        max_end = min(len(record.seq), max_end)

                        extracted_sequence = record.seq[min_start - 1:max_end]  # Adjust for Python's 0-based indexing
                        output_file = Path(extracted_folder) / f"{filename}_{query}_extracted_sequence.fasta"
                        with open(output_file, "w") as output_handle:
                            SeqIO.write([SeqIO.SeqRecord(extracted_sequence, id=record.id, description="")],
                                        output_handle, "fasta")
                            logging.info(
                                f"Extracted sequence for {query} in {filename}, min start: {min_start}, max end: {max_end}, length: {len(extracted_sequence)}")
            else:
                logging.warning(f"No FASTA file found for {filename}. Expected path: {fasta_file_path}")


def move_files(src_folder, dest_folder, file_type=".fna"):
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)
    all_files_paths = glob.glob(src_folder + '/**/*' + file_type, recursive=True)
    for file_path in all_files_paths:
        shutil.move(file_path, dest_folder)


def main_move():
    # to move files!
    move_files("/Users/josediogomoura/Documents/BioFago/BioFago/data/genomesAllErwinia/ncbi_dataset/ncbi_dataset/reference_crispr",
               "/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/genomesAllErwinia/ncbi_dataset/ncbi_dataset/fasta",
               file_type=".fna")


def main():
    # Define paths and parameters
    #genomes_folder_path = "/Volumes/Crucial_X9/BioFago/reference_crispr/ApproachFlankGenes/genomes/ErwiniaAmyl/ncbi_dataset/fasta_files"
    db_folder_path = "/reference_crispr/OtherSpecies/Erwinia_billingiae/srl/db_folder"
    results_folder_path = "/reference_crispr/OtherSpecies/Erwinia_billingiae/srl/blast_results"
    # ref_seq_gb_path = "/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/ref_sequences/Capsule_locus_12genes_X77921.gb"
    extracted_genes_fasta = '/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/OtherSpecies/Erwinia_billingiae/srl/flank_genes/flank_genes.fasta'

    # Ensure output folders exist
    Path(db_folder_path).mkdir(parents=True, exist_ok=True)
    Path(results_folder_path).mkdir(parents=True, exist_ok=True)

    # # Extract flank genes and write to a FASTA file
    # logging.info("Extracting flank genes...")
    # flank_gene_extractor = FlankGeneExtractor(ref_seq_gb_path)
    # flank_gene_extractor.extract_flank_genes()
    # flank_gene_extractor.write_genes_to_fasta(extracted_genes_fasta)

    # Collect all .fna genome files from subdirectories
    # Now you can use it like this
    src_folder = "/reference_crispr/OtherSpecies/Erwinia_billingiae/ncbi_dataset/ncbi_dataset/reference_crispr"
    dest_folder = "/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/OtherSpecies/Erwinia_billingiae/ncbi_dataset/ncbi_dataset/fasta"
    move_files(src_folder, dest_folder)

    # # Assuming BlastRunner can take a list of genome files
    logging.info("Initializing BlastRunner with collected genome files...")
    gene_blast_runner = BlastRunner(extracted_genes_fasta, dest_folder, db_folder_path, results_folder_path)

    # Process genomes (this function must be able to handle a list of genome file paths)
    logging.info("Processing collected genome files...")
    gene_blast_runner.process_genomes()

    # Run BLAST on all genomes (this function must also handle the list)
    logging.info("Running BLAST on all genomes...")
    gene_blast_runner.run_blast_on_all_genomes()

    # Step 3: Compile BLAST results into a CSV
    blast_results_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/OtherSpecies/Erwinia_billingiae/srl/blast_results'
    genomes_folder = "/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/OtherSpecies/Erwinia_billingiae/ncbi_dataset/ncbi_dataset/fasta"
    output_csv = '/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/OtherSpecies/Erwinia_billingiae/srl/blast_csv/compiled_results.csv'
    log_file = os.path.join(os.path.dirname(output_csv), 'post_blast_output.log')
    #
    post_blast = PostBlastOutput(blast_results_folder, genomes_folder, output_csv, log_file)
    post_blast.compile_blast_results_to_csv()  # Compiles BLAST results into a CSV
    #
    extracted_path = "/reference_crispr/OtherSpecies/Erwinia_billingiae/srl/extracted_sequences"
    post_blast.extract_regions_from_genomes_v2(extracted_path)

    logging.info("Workflow completed successfully.")


if __name__ == "__main__":
    main()
