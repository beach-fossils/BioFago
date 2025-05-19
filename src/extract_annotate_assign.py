import os
import shutil
import logging
from pathlib import Path
from typing import List, Optional, Dict, Any, Tuple

import pandas as pd

from assigning_types.AssignTypes import AssignTypes
from utils.extract_with_flank_genes import extract_with_flank_genes
from utils.get_prokka_faa_file import get_prokka_faa_file
from utils.prokka_docker import run_prokka_docker
from assigning_types.Database import BlastProteinv2, TypeAnalysis
from database.BlastRunner import BlastRunner
from database.FlankGenesForBlast import PostBlastOutput
from assigning_types.OmpAGenotyper import genotype_ompa_for_genome

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Dynamically determine the path to the reference types directory
SCRIPT_DIR = Path(__file__).resolve().parent
REFERENCE_TYPES = SCRIPT_DIR.parents[0] / 'reference_types_database'
REFERENCE_OMPA = REFERENCE_TYPES / 'ompa' / 'types_ompa.gb'


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
):
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
        return None

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
        return None

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
    type_report, final_type_locus, final_type, flagged_genes, final_type_info = gene_presence.assign_type(df_analyzed)

    pd.set_option('display.max_rows', None)
    for idx, row in type_report.iterrows():
        logging.info(row)

    logging.info(
        f"Final Assigned Type for {os.path.basename(genome_dir)} with {reference_file.stem}: {final_type_info}")

    type_report, final_type_locus, final_type, flagged_genes, final_type_info = gene_presence.assign_type(df_analyzed)

    logging.info(f"Type analysis results for {os.path.basename(genome_dir)}:")
    logging.info(f"  Final type locus: {final_type_locus}")
    logging.info(f"  Final type: {final_type}")
    logging.info(f"  Flagged genes: {flagged_genes}")

    return os.path.basename(
        genome_dir), reference_file.stem, final_type_locus, final_type, flagged_genes, final_type_info


def process_genome(genome_dir: str, base_folder: str) -> List[tuple]:
    """Process a single genome against all reference files."""
    genome_files = [f for f in os.listdir(genome_dir) if f.endswith(".fasta")]
    if not genome_files:
        logging.warning(f"No fasta files found in {genome_dir}")
        return []

    input_genome = os.path.join(genome_dir, genome_files[0])
    # Get full genome name without extension
    genome_dir_path = Path(genome_dir)
    genome_name = genome_dir_path.name
    
    # Preserve the full filename pattern without any tampering
    # Only remove the final extension if present
    if any(genome_name.lower().endswith(ext) for ext in ['.fasta', '.fa', '.fna']):
        genome_name = os.path.splitext(genome_name)[0]
    
    # Special handling for GCF/GCA files to ensure full name preservation
    # For example, GCF_002732285.1_GCF_002732285.1_ASM273228v1_genomic.fna should become
    # GCF_002732285.1_GCF_002732285.1_ASM273228v1_genomic in the output
    
    logging.info(f"Using genome name: {genome_name} for processing")
    
    # Add extra logging for debugging name issues
    logging.debug(f"Original directory name: {genome_dir_path}")
    logging.debug(f"Processed genome name: {genome_name}")
    
    genome_output_folder = os.path.join(base_folder, genome_name)
    os.makedirs(genome_output_folder, exist_ok=True)

    results = []
    
    # Process ompA genotyping as a special case
    if REFERENCE_OMPA.exists():
        try:
            logging.info(f"Processing ompA genotyping for genome {genome_name}")
            
            # Run ompA genotyping
            typing_report, ompa_type = genotype_ompa_for_genome(
                Path(input_genome), 
                REFERENCE_OMPA, 
                Path(genome_output_folder)
            )
            
            # Format result as a tuple similar to other typing results to maintain compatibility
            ompa_result = (
                genome_name,              # genome name
                "types_ompa",             # reference type name
                "",                       # locus information (empty for ompA since it's a single gene)
                ompa_type,                # type designation
                None,                     # flagged genes (not applicable for ompA)
                typing_report             # type info with details
            )
            
            results.append(ompa_result)
            logging.info(f"Added ompA genotyping result for {genome_name}: {typing_report}")
            
        except Exception as e:
            logging.error(f"Error in ompA genotyping for {genome_name}: {e}")
    
    # Process other loci-based references
    for gb_folder in REFERENCE_TYPES.iterdir():
        if gb_folder.is_dir() and gb_folder.name != "ompa":  # Skip ompA folder as we handled it separately
            for gb_file in gb_folder.glob("*.gb"):
                logging.info(f"Processing genome {genome_name} with reference {gb_file}")
                result = process_genome_with_reference(genome_dir, input_genome, gb_file, genome_output_folder)
                if result:
                    logging.info(f"Got result for {genome_name}: {result}")
                    results.append(result)
                else:
                    logging.warning(f"No result for {genome_name} with reference {gb_file}")

    logging.info(f"Total results for genome {genome_name}: {len(results)}")
    return results


def extract_annotate_assign(genomes_folder: Path, batch_size: int = 0, num_workers: int = 4, docker_limit: int = 4) -> List[tuple]:
    """Main function to extract, annotate, and assign loci types to genomes.
    
    Args:
        genomes_folder: Path to folder containing genomes
        batch_size: Number of genomes to process in a batch (0 for all at once)
        num_workers: Number of worker processes for parallelization
        docker_limit: Maximum number of concurrent Docker containers
    """
    if not os.path.exists(genomes_folder):
        raise FileNotFoundError("Input genome folder not found.")

    # Initialize Docker semaphore
    from utils.prokka_docker import init_docker_semaphore
    init_docker_semaphore(docker_limit)
    logging.info(f"Initialized Docker semaphore with limit of {docker_limit} containers")

    # Dynamically create base_folder path based on genomes_folder
    genomes_folder_path = Path(genomes_folder).resolve()
    base_folder = genomes_folder_path / 'types_finder'
    os.makedirs(base_folder, exist_ok=True)

    move_genomes_to_folders(genomes_folder)
    genome_dirs = [os.path.join(genomes_folder, d) for d in os.listdir(genomes_folder) if
                   os.path.isdir(os.path.join(genomes_folder, d))]
    logging.info(f"Found {len(genome_dirs)} genome directories")

    all_results = []
    
    # If batch_size is 0 or greater than the number of genomes, process all at once
    if batch_size <= 0 or batch_size >= len(genome_dirs):
        all_results = _process_genome_batch(genome_dirs, base_folder, num_workers)
    else:
        # Process genomes in batches
        logging.info(f"Processing genomes in batches of {batch_size}")
        for i in range(0, len(genome_dirs), batch_size):
            batch = genome_dirs[i:i+batch_size]
            logging.info(f"Processing batch {i//batch_size + 1} with {len(batch)} genomes")
            batch_results = _process_genome_batch(batch, base_folder, num_workers)
            all_results.extend(batch_results)
            logging.info(f"Completed batch {i//batch_size + 1}, total processed: {len(all_results)}/{len(genome_dirs)}")
    
    logging.info("Final results:")
    for result in all_results:
        logging.info(
            f"Genome: {result[0]}, Reference: {result[1]}, Type: {result[2]} ({result[3]}), Flagged genes: {result[4]}")

    return all_results

def _process_genome_batch(genome_dirs: List[str], base_folder: Path, num_workers: int) -> List[tuple]:
    """Process a batch of genomes in parallel."""
    all_results = []
    
    if num_workers <= 1:
        # Sequential processing
        for genome_dir in genome_dirs:
            logging.info(f"Processing genome directory sequentially: {genome_dir}")
            results = process_genome(genome_dir, base_folder)
            all_results.extend(results)
    else:
        # Parallel processing with multiprocessing
        from concurrent.futures import ProcessPoolExecutor, as_completed
        
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            # Submit all genome directories for processing
            future_to_dir = {
                executor.submit(process_genome, genome_dir, base_folder): genome_dir 
                for genome_dir in genome_dirs
            }
            
            # Process results as they complete
            for future in as_completed(future_to_dir):
                genome_dir = future_to_dir[future]
                try:
                    results = future.result()
                    all_results.extend(results)
                    logging.info(f"Successfully processed genome directory: {genome_dir}")
                except Exception as e:
                    logging.error(f"Error processing genome directory {genome_dir}: {e}")
    
    return all_results


if __name__ == "__main__":
    genomes_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/test-data/rubus2/testing3'

    extract_annotate_assign(genomes_folder)
