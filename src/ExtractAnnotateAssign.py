import os
import shutil
import logging
from pathlib import Path
from typing import List, Optional

import pandas as pd

from src.assigning_types.AssignTypes import AssignTypes
from src.utils.extract_with_flank_genes import extract_with_flank_genes
from src.utils.get_prokka_faa_file import get_prokka_faa_file
from src.utils.prokka_docker import run_prokka_docker
from src.assigning_types.Database import BlastProteinv2, TypeAnalysis
from src.database.BlastRunner import BlastRunner
from src.database.FlankGenesForBlast import PostBlastOutput

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Dynamically determine the path to the reference types directory
SCRIPT_DIR = Path(__file__).resolve().parent
REFERENCE_TYPES = SCRIPT_DIR.parents[0] / 'fully_gb_database'


def move_genomes_to_folders(genomes_folder: str) -> None:
    """Move genome files into their own subdirectories."""
    for genome_file in os.listdir(genomes_folder):
        if genome_file.endswith(".fasta"):
            genome_path = os.path.join(genomes_folder, genome_file)
            genome_dir = os.path.join(genomes_folder, genome_file[:-6])
            os.makedirs(genome_dir, exist_ok=True)
            shutil.move(genome_path, os.path.join(genome_dir, genome_file))
            logging.info(f"Moved {genome_file} to {genome_dir}")


def find_extracted_fasta(directory: str) -> Optional[str]:
    """Find the extracted fasta file in the given directory."""
    for file in os.listdir(directory):
        if file.endswith(".fasta"):
            return os.path.join(directory, file)
    return None


def process_genome_with_reference(
        genome_dir: str, input_genome: str, reference_file: Path, genome_output_folder: str
) -> None:
    """Process a genome with a specific reference file."""
    logging.info(f"Processing reference file: {reference_file}")

    gb_output_folder = os.path.join(genome_output_folder, reference_file.stem)
    os.makedirs(gb_output_folder, exist_ok=True)

    # Extract region
    extracted_region_dir = extract_with_flank_genes(
        source_output_dir=gb_output_folder,
        ref_db=str(reference_file),
        input_genome=input_genome,
        AssignTypes=AssignTypes,
        BlastRunner=BlastRunner,
        PostBlastOutput=PostBlastOutput
    )

    extracted_region_fasta = find_extracted_fasta(extracted_region_dir)
    logging.info(f"Expected extracted fasta file at: {extracted_region_fasta}")

    if not extracted_region_fasta:
        logging.warning(
            f"Extracted region fasta file not found at {extracted_region_fasta}. Skipping this genome/reference combination.")
        return

    # Run Prokka annotation
    prokka_base_output_folder = os.path.join(gb_output_folder, 'prokka')
    locus_tag_prefix = 'PREFIX'
    run_prokka_docker(
        fasta_file=extracted_region_fasta,
        base_output_folder=prokka_base_output_folder,
        custom_db_path=str(reference_file),
        locus_tag_prefix=locus_tag_prefix
    )

    proteins_faa = get_prokka_faa_file(prokka_base_output_folder)
    if not proteins_faa:
        logging.error("Prokka .faa file not found.")
        return

    # Define paths for Blast and TypeAnalysis
    db_folder = os.path.join(gb_output_folder, "db_folder")
    results_folder = os.path.join(gb_output_folder, "results")
    final_csv_file = os.path.join(gb_output_folder, "final.csv")
    output_csv = os.path.join(gb_output_folder, "loci_gene_presence.csv")

    blast = BlastProteinv2(reference_types=str(reference_file), db_folder=db_folder, results_folder=results_folder)
    input_csv = os.path.join(gb_output_folder, "loci.csv")
    blast.create_csv_loci(input_csv)
    blast.process_csv(input_csv, proteins_faa, final_csv_file)

    # Analyze locus information
    gene_presence = TypeAnalysis(input_csv=final_csv_file, output_csv=output_csv, proteins_faa=proteins_faa)
    df_analyzed = gene_presence.analyze_locus()
    type_report, final_type_locus, final_type = gene_presence.assign_type(df_analyzed)

    pd.set_option('display.max_rows', None)
    for idx, row in type_report.iterrows():
        logging.info(row)

    logging.info(
        f"Final Assigned Type for {os.path.basename(genome_dir)} with {reference_file.stem}: {final_type_locus} ({final_type})")


def process_genome(genome_dir: str, base_folder: str) -> None:
    """Process a single genome against all reference files."""
    genome_files = [f for f in os.listdir(genome_dir) if f.endswith(".fasta")]
    if not genome_files:
        logging.warning(f"No fasta files found in {genome_dir}")
        return

    input_genome = os.path.join(genome_dir, genome_files[0])
    genome_name = os.path.basename(genome_dir)
    genome_output_folder = os.path.join(base_folder, genome_name)
    os.makedirs(genome_output_folder, exist_ok=True)

    for gb_folder in REFERENCE_TYPES.iterdir():
        if gb_folder.is_dir():
            for gb_file in gb_folder.glob("*.gb"):
                process_genome_with_reference(genome_dir, input_genome, gb_file, genome_output_folder)


def extract_annotate_assign(
        genomes_folder: Path
) -> None:
    """Main function to extract, annotate, and assign loci types to genomes."""
    if not os.path.exists(genomes_folder):
        raise FileNotFoundError("Input genome folder not found.")

    # Dynamically create base_folder path based on genomes_folder
    genomes_folder_path = Path(genomes_folder).resolve()
    base_folder = genomes_folder_path.parent / 'types_finder'
    os.makedirs(base_folder, exist_ok=True)

    move_genomes_to_folders(genomes_folder)
    genome_dirs = [os.path.join(genomes_folder, d) for d in os.listdir(genomes_folder) if
                   os.path.isdir(os.path.join(genomes_folder, d))]
    logging.info(f"Genome directories: {genome_dirs}")

    for genome_dir in genome_dirs:
        logging.info(f"Processing genome directory: {genome_dir}")
        process_genome(genome_dir, base_folder)


if __name__ == "__main__":
    genomes_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/test_2/genomes'

    extract_annotate_assign(genomes_folder)