import os
import logging
from pathlib import Path


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


def extract_with_flank_genes(source_output_dir, ref_db, input_genome, AssignTypes, BlastRunner, PostBlastOutput):
    # Create necessary subdirectories inside the source_output_dir
    output_dir = source_output_dir
    first_last_genes_dir = os.path.join(output_dir, 'first_last_genes')
    extracted_region_dir = os.path.join(output_dir, 'extracted_region')
    blast_results_dir = os.path.join(output_dir, 'blast', 'blast_results')
    blast_csv_dir = os.path.join(output_dir, 'blast', 'csv')

    os.makedirs(first_last_genes_dir, exist_ok=True)
    os.makedirs(extracted_region_dir, exist_ok=True)
    os.makedirs(blast_results_dir, exist_ok=True)
    os.makedirs(blast_csv_dir, exist_ok=True)

    # Initialize the AssignTypes class
    assigner = AssignTypes(ref_db, input_genome)
    first_gene_seq = assigner.db.first_gene_sequence_nt
    first_gene_tag = assigner.db.first_gene['locus_tag']
    last_gene_seq = assigner.db.last_gene_sequence_nt
    last_gene_tag = assigner.db.last_gene['locus_tag']

    # Unify those sequences into one fasta file to the path
    extracted_genes_fasta = os.path.join(first_last_genes_dir, 'first_last_genes.fasta')

    write_genes_to_fasta(output_fasta=extracted_genes_fasta, first_gene_seq=first_gene_seq, last_gene_seq=last_gene_seq,
                         first_gene_tag=first_gene_tag, last_gene_tag=last_gene_tag)

    db_folder_path = os.path.join(blast_results_dir, 'db')
    results_folder_path = blast_results_dir
    # Ensure directories exist
    os.makedirs(db_folder_path, exist_ok=True)
    os.makedirs(results_folder_path, exist_ok=True)

    # Set genomes_folder to the directory containing input_genome
    genomes_folder = os.path.dirname(input_genome)

    logging.info("Initializing BlastRunner with collected genome files...")
    gene_blast_runner = BlastRunner(extracted_genes_fasta, genomes_folder, db_folder_path, results_folder_path)
    gene_blast_runner.run_blast_on_all_genomes()

    output_csv = os.path.join(blast_csv_dir, 'results.csv')
    log_file = os.path.join(blast_csv_dir, 'post_blast_output.log')

    post_blast = PostBlastOutput(results_folder_path, genomes_folder, output_csv, log_file)
    post_blast.compile_blast_results_to_csv()  # Compiles BLAST results into a CSV

    post_blast.extract_regions_from_genomes_v2(extracted_region_dir)

    logging.info("Workflow completed successfully.")

    return extracted_region_dir