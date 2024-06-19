import logging
import os
from pathlib import Path
import shutil
from src.MetricsSpeciesCaller import run_species_metrics_finder, setup_logging
from src.ExtractAnnotateAssign import extract_annotate_assign

# Configure logging
LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'
LOG_LEVEL = logging.DEBUG
logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)


def create_individual_folders(genomes_folder: str) -> Path:
    """Move each genome to its individual folder."""
    genomes_folder_path = Path(genomes_folder).resolve()
    for genome_file in genomes_folder_path.glob('*.fasta'):
        genome_dir = genomes_folder_path / genome_file.stem
        os.makedirs(genome_dir, exist_ok=True)
        shutil.move(str(genome_file), str(genome_dir / genome_file.name))
        logging.info(f"Moved {genome_file} to {genome_dir}")
    return genomes_folder_path


def run_species_metrics_for_all(genomes_folder_path: Path, threshold_species: float):
    """Run species metrics finder for all genomes in their individual folders."""
    for genome_dir in genomes_folder_path.iterdir():
        if genome_dir.is_dir():
            for genome_file in genome_dir.glob('*.fasta'):
                logging.info(f"Processing genome file: {genome_file}")
                try:
                    run_species_metrics_finder(str(genome_file), threshold_species)
                except Exception as e:
                    logging.error(f"Error in run_species_metrics_finder for {genome_file}: {e}")


# def move_results_to_genome_folder(genome_dir: Path):
#     """Move the species_finder and types_finder folders to the genome's individual folder."""
#     species_finder_src = genome_dir.parent / 'species_finder'
#     types_finder_src = genome_dir.parent / 'types_finder'
#
#     if species_finder_src.exists():
#         shutil.move(str(species_finder_src), str(genome_dir / 'species_finder'))
#         logging.info(f"Moved species_finder to {genome_dir}")
#     if types_finder_src.exists():
#         shutil.move(str(types_finder_src), str(genome_dir / 'types_finder'))
#         logging.info(f"Moved types_finder to {genome_dir}")


def run_species_and_types_finder(genomes_folder: str, threshold_species: float = 0.95):
    """Run species metrics finder and types finder on all genomes in the folder."""
    # Create individual folders for each genome
    genomes_folder_path = create_individual_folders(genomes_folder)

    # Run species metrics finder for each genome
    run_species_metrics_for_all(genomes_folder_path, threshold_species)
    # after this there will be a csv inside the species_finder folder with the species name

    # Run extract and annotate assign on each genome folder
    for genome_dir in genomes_folder_path.iterdir():
        if genome_dir.is_dir():
            extract_annotate_assign(genome_dir)



if __name__ == '__main__':
    genomes_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/test_2/genomes2'
    run_species_and_types_finder(genomes_folder, threshold_species=0.95)