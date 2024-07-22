import logging
from pathlib import Path
from src.metrics_species_caller import run_species_metrics_finder, setup_logging
from src.extract_annotate_assign import extract_annotate_assign
from src.resistance.str_resistance import StrResistance
from src.utils.folder_csv_manager import create_individual_folders, run_species_metrics_for_all, update_species_csv
from src.plasmids.plasmid_finder import PlasmidFinder
from src.crr_genotypes.crr_genotype_finder import CRRFinder

# Configure logging
LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'
LOG_LEVEL = logging.DEBUG
logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)
logger = logging.getLogger(__name__)


def run_species_and_types_finder(genomes_folder: Path, threshold_species: float = 0.95) -> None:
    """Run species metrics finder and types finder on all genomes in the folder."""
    try:
        # Create individual folders for each genome
        logger.info(f"Creating individual folders for genomes in: {genomes_folder}")
        genomes_folder_path = create_individual_folders(genomes_folder)

        # Run species metrics finder for each genome
        logger.info(f"Running species metrics finder for genomes in: {genomes_folder_path}")
        run_species_metrics_for_all(genomes_folder_path, threshold_species)
        # after this there will be a csv inside the species_finder folder with the species name

        # Plasmid finder runs here
        logger.info(f"Running plasmid finder for genomes in: {genomes_folder_path}")
        plasmid_finder = PlasmidFinder()
        plasmid_finder.process_all_genomes_in_folder(Path(genomes_folder_path))

        # Sorbitol resistance finder runs here
        logger.info(f"Running sorbitol resistance finder for genomes in: {genomes_folder_path}")
        sorbitol_resistance = StrResistance()
        sorbitol_resistance.process_all_genomes_in_folder(Path(genomes_folder_path))

        # CRRFinder runs here
        logger.info(f"Running CRRFinder for genomes in: {genomes_folder_path}")
        for genome_dir in Path(genomes_folder_path).iterdir():
            if genome_dir.is_dir():
                for genome_file in genome_dir.glob('*.fasta'):
                    crr_finder = CRRFinder(genome_file, genome_dir)
                    crr_finder.analyze_genome()
                    crr_finder.process_directory()

        # Run extract and annotate assign on each genome folder
        results = extract_annotate_assign(genomes_folder_path)
        for genome_name, locus_name, locus_type, assign_confidence in results:
            update_species_csv(genomes_folder_path, genome_name, locus_name, locus_type, assign_confidence)
    except Exception as e:
        logger.error(f"Error in run_species_and_types_finder: {e}")


if __name__ == '__main__':
    try:
        genomes_folder = Path('/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/test_2/genomes')
        run_species_and_types_finder(genomes_folder, threshold_species=0.95)
    except Exception as e:
        logger.error(f"Error in __main__: {e}")
