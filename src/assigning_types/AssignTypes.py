import os
import subprocess
import pandas as pd
from Bio import SeqIO
import logging
from src.assigning_types.Database import Database
from src.database.BlastRunner import BlastRunner
from src.database.FlankGenesForBlast import FlankGeneExtractor, PostBlastOutput
from src.utils.extract_with_flank_genes import extract_with_flank_genes

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class AssignTypes:
    def __init__(self, db_file, input_genome):
        if not os.path.exists(db_file):
            raise FileNotFoundError(f"Database file '{db_file}' not found.")
        if not os.path.exists(input_genome):
            raise FileNotFoundError(f"Input genome file '{input_genome}' not found.")

        self.db = Database(db_file)  # Store the Database instance
        self.loci = Database(db_file).get_loci_info()
        self.input_genome = input_genome
        self.blast_db = "input_genome_db"

    def make_blast_db(self):
        cmd = f"makeblastdb -in {self.input_genome} -dbtype nucl -out {self.blast_db}"
        logging.info(f"Running command: {cmd}")
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error creating BLAST database: {e}")
            raise

    def run_blast(self, query_file, output_file):
        cmd = (
            f"blastn -query {query_file} -db {self.blast_db} -out {output_file} "
            f"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' "
            f"-max_hsps 1 -evalue 1e-20 -qcov_hsp_perc 100"
        )
        logging.info(f"Running BLAST with command: {cmd}")
        try:
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error running BLAST: {e}")
            raise

    def parse_blast_results(self, file):
        try:
            columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
                       "send", "evalue", "bitscore", "qlen", "slen"]
            blast_results = pd.read_csv(file, sep="\t", header=None, names=columns)
            return blast_results
        except Exception as e:
            logging.error(f"Error parsing BLAST results from {file}: {e}")
            return pd.DataFrame()

    def generate_fasta_files(self, output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for locus, genes in self.loci.items():
            for gene in genes:
                file_name = f"{locus}_{gene['locus_tag']}.fasta"
                file_path = os.path.join(output_dir, file_name)
                with open(file_path, "w") as f:
                    if gene['nucleotide_sequence'] is not None:
                        f.write(f">{gene['locus_tag']}\n{gene['nucleotide_sequence']}\n")
                logging.info(f"Generated FASTA file: {file_path}")

    def assign_types(self, output_dir):
        self.generate_fasta_files(output_dir)
        results = []
        for locus, genes in self.loci.items():
            for gene in genes:
                file_name = f"{locus}_{gene['locus_tag']}.fasta"
                query_file = os.path.join(output_dir, file_name)
                output_file = os.path.join(output_dir, f"{file_name}_blast.tsv")
                self.run_blast(query_file, output_file)
                blast_results = self.parse_blast_results(output_file)
                if blast_results.empty:
                    logging.info(f"No BLAST results for file: {file_name}")
                    continue

                best_hit = blast_results.iloc[0]
                coverage = (best_hit['length'] / best_hit['qlen']) * 100
                identity = best_hit['pident']

                results.append({
                    'Locus': locus,
                    'Gene': gene['gene'],
                    'Locus Tag': gene['locus_tag'],
                    'Coverage (%)': coverage,
                    'Identity (%)': identity
                })
        return results

    def write_results_to_csv(self, results, output_file):
        df = pd.DataFrame(results)
        df.to_csv(output_file, index=False)
        logging.info(f"Results written to: {output_file}")


def write_genes_to_fasta(output_fasta, first_gene_seq, last_gene_seq, first_gene_tag, last_gene_tag, line_length=60):
    def format_sequence(sequence, line_length):
        return '\n'.join(sequence[i:i + line_length] for i in range(0, len(sequence), line_length))

    try:
        with open(output_fasta, 'w') as f:
            f.write(f">{first_gene_tag}\n")
            f.write(f"{format_sequence(first_gene_seq, line_length)}\n")

            f.write(f">{last_gene_tag}\n")
            f.write(f"{format_sequence(last_gene_seq, line_length)}\n")
    except KeyError as e:
        logging.error(f"Error retrieving sequence for gene: {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")


def main():
    # # Reference database file and genomes genome

    ref_db = "/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/assign_types/reference_types_database/cellulose/types_cellulose.gb"
    input_genome = "/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/results_sequenciados/sequences/genomes.fasta"
    output_dir = '/Users/josediogomoura/Documents/BioFago/BioFago/data/assign_types/cellulose/test_1'

    region = extract_with_flank_genes(source_output_dir=output_dir, ref_db=ref_db, input_genome=input_genome,
                             AssignTypes=AssignTypes, BlastRunner=BlastRunner, PostBlastOutput=PostBlastOutput)

    print(region)



    # ref_db = "/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/assign_types/lps/types_lps.gb"
    # input_genome = "/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/results_sequenciados/sequences/genomes.fasta"
    # output_dir = '/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/test1/assigning_types/first_test'
    # results_file = '/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/test1/assigning_types/first_test/results.csv'
    #
    # # Ensure output directories exist
    # os.makedirs(output_dir, exist_ok=True)
    #
    # # Initialize the AssignTypes class
    # assigner = AssignTypes(ref_db, input_genome)
    # first_gene_seq = assigner.db.first_gene_sequence_nt
    # first_gene_tag = assigner.db.first_gene['locus_tag']
    # last_gene_seq = assigner.db.last_gene_sequence_nt
    # last_gene_tag = assigner.db.last_gene['locus_tag']
    #
    # # Unify those sequences into one fasta file to the path
    # extracted_genes_fasta = '/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/test1/assigning_types/first_test/first_last_genes.fasta'
    # output_first_last = '/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/test1/assigning_types/first_test/first_last_genes/flank_genes.fasta'
    #
    # write_genes_to_fasta(output_fasta=extracted_genes_fasta, first_gene_seq=first_gene_seq, last_gene_seq=last_gene_seq,
    #                      first_gene_tag=first_gene_tag, last_gene_tag=last_gene_tag)
    #
    # db_folder_path = '/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/test1/assigning_types/first_test/blast/blast_results/db'
    # results_folder_path = '/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/test1/assigning_types/first_test/blast'
    # # Ensure directories exist
    # os.makedirs(db_folder_path, exist_ok=True)
    # os.makedirs(results_folder_path, exist_ok=True)
    #
    # # Set genomes_folder to the directory containing input_genome
    # genomes_folder = os.path.dirname(input_genome)
    #
    # logging.info("Initializing BlastRunner with collected genome files...")
    # gene_blast_runner = BlastRunner(extracted_genes_fasta, genomes_folder, db_folder_path, results_folder_path)
    # gene_blast_runner.run_blast_on_all_genomes()
    #
    # blast_results_folder = results_folder_path
    # output_csv = '/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/test1/assigning_types/first_test/blast/csv/results.csv'
    # log_file = os.path.join(os.path.dirname(output_csv), 'post_blast_output.log')
    #
    # post_blast = PostBlastOutput(blast_results_folder, genomes_folder, output_csv, log_file)
    # post_blast.compile_blast_results_to_csv()  # Compiles BLAST results into a CSV
    #
    # extracted_path = "/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/test1/assigning_types/first_test/extracted_region"
    # post_blast.extract_regions_from_genomes_v2(extracted_path)
    #
    # logging.info("Workflow completed successfully.")


if __name__ == "__main__":
    main()
