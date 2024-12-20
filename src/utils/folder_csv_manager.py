import logging
import os
import shutil
from pathlib import Path
from typing import List, Dict

import pandas as pd
import time

from metrics_species_caller import new_run_species_metrics_finder

def create_species_finder_folder(genomes_folder_path: Path) -> Path:
    species_finder_path = genomes_folder_path / 'species_finder'
    species_finder_path.mkdir(exist_ok=True)
    return species_finder_path


def run_species_metrics_for_all(genomes_folder_path: Path, species_finder_path: Path, threshold_species: float):
    """Run species metrics finder for all genomes in their individual folders."""
    logging.info(f"Starting species metrics analysis in {genomes_folder_path}")

    # Find all FASTA files in both root and subdirectories
    genome_files = []

    # Look for files directly in the genomes folder
    for fasta_file in genomes_folder_path.glob('*.fasta'):
        if fasta_file.is_file():
            genome_files.append(fasta_file)

    # Look in subdirectories
    for subdir in genomes_folder_path.iterdir():
        if subdir.is_dir():
            for fasta_file in subdir.glob('*.fasta'):
                if fasta_file.is_file():
                    genome_files.append(fasta_file)

    if not genome_files:
        logging.error(f"No FASTA files found in {genomes_folder_path} or its subdirectories")
        return

    logging.info(f"Found {len(genome_files)} genome files to process")

    # Process each genome
    for genome_file in genome_files:
        logging.info(f"Processing genome file for species analysis: {genome_file}")
        try:
            # Create output directory if it doesn't exist
            species_finder_path.mkdir(parents=True, exist_ok=True)

            # Run metrics finder for this genome
            new_run_species_metrics_finder(
                single_sequence_path=genome_file,
                species_finder_path=species_finder_path,
                threshold_species=threshold_species
            )

            # Verify results
            output_file = species_finder_path / f"{genome_file.stem}.csv"
            if output_file.exists():
                try:
                    df = pd.read_csv(output_file)
                    species = df['Species'].iloc[0] if 'Species' in df.columns else 'Unknown'
                    ani = df['ANI'].iloc[0] if 'ANI' in df.columns else 0.0
                    logging.info(f"Species assignment for {genome_file.name}: {species} (ANI: {ani})")
                except Exception as e:
                    logging.error(f"Error reading species results for {genome_file.name}: {e}")
            else:
                logging.error(f"No results file generated for {genome_file.name}")

        except Exception as e:
            logging.error(f"Error in species metrics finder for {genome_file}: {e}")
            logging.exception("Exception details:")

    logging.info("Completed species metrics analysis for all genomes")

def create_individual_folders(genomes_folder: Path) -> Path:
    """Move each genome to its individual folder, preserving the original file name."""
    genomes_folder_path = Path(genomes_folder).resolve()
    for genome_file in genomes_folder_path.glob('*.fasta'):
        genome_dir = genomes_folder_path / genome_file.stem
        os.makedirs(genome_dir, exist_ok=True)
        shutil.move(str(genome_file), str(genome_dir / genome_file.name))
        logging.info(f"Moved {genome_file} to {genome_dir}")
    return genomes_folder_path


def update_species_csv(species_finder_path: Path, genome_name: str, locus_name: str, locus_type: str,
                       assign_confidence: str):
    species_csv_path = species_finder_path / f"{genome_name}.csv"
    if species_csv_path.exists():
        df = pd.read_csv(species_csv_path)
        column_name = f"{locus_name}"
        df[column_name] = f"{locus_type} ({assign_confidence})"
        df.to_csv(species_csv_path, index=False)
        logging.info(f"Updated {species_csv_path} with {locus_name}: {locus_type} ({assign_confidence})")
    else:
        logging.warning(f"Species CSV file not found at {species_csv_path}")


def update_species_csv_with_results(species_finder_path: Path, genome_name: str, plasmid_info: str, sorbitol_info: str, crispr_info: str, crr_info: str):
    csv_file = species_finder_path / f"{genome_name}.csv"
    if csv_file.exists():
        df = pd.read_csv(csv_file)
        df['Plasmid'] = plasmid_info
        df['Sorbitol'] = sorbitol_info
        df['CRISPR'] = crispr_info
        df['CRR'] = crr_info
        df.to_csv(csv_file, index=False)
        logging.info(f"Updated {csv_file} with plasmid, sorbitol, CRISPR, and CRR information")
    else:
        logging.warning(f"Species CSV file not found at {csv_file}")


# def create_species_finder_folder(genomes_folder_path: Path) -> Path:
#     """Create the species_finder folder at the correct level."""
#     species_finder_path = genomes_folder_path / 'species_finder'
#     species_finder_path.mkdir(exist_ok=True)
#     return species_finder_path


def cleanup_unwanted_species_folders(base_output_path: Path, main_species_finder_path: Path):
    """
    Remove unwanted 'species_finder' folders from output directories.
    """
    for output_dir in ['CRISPR_finder', 'CRR_finder', 'plasmid_finder', 'resistance_finder']:
        dir_path = base_output_path / output_dir
        if dir_path.exists():
            species_finder = dir_path / 'species_finder'
            if species_finder.exists() and species_finder != main_species_finder_path:
                try:
                    shutil.rmtree(species_finder)
                    logging.info(f"Removed unwanted species_finder folder from {output_dir}")
                except Exception as e:
                    logging.error(f"Error removing species_finder folder from {output_dir}: {e}")


