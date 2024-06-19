import logging
from pathlib import Path
import pandas as pd
from typing import List
from src.assigning_types.assembly_statistics import FastaStatistics
from src.assigning_types.SpeciesFinderWPyani import LocalANIExecutor, SpeciesTabModifier


#REFERENCE_GENOMES = '/Users/josediogomoura/Documents/BioFago/BioFago/species_reference_genomes'
# Dynamically determine the path to the reference types directory
SCRIPT_DIR = Path(__file__).resolve().parent
REFERENCE_GENOMES = SCRIPT_DIR.parents[0] / 'species_reference_genomes'

LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'
LOG_LEVEL = logging.DEBUG


def setup_logging():
    logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)


def generate_paths(input_path: Path) -> dict:
    base_dir = input_path.parent.parent
    results_dir = base_dir / 'species_finder' / 'results'
    return {
        "results_dir": results_dir,
        "output_stats_path": results_dir / 'stats' / 'stats_results.csv',
        "ani_tab_file": results_dir / 'ani_tab' / 'ani_results.tab',
        "species_output_file": base_dir / 'species_finder' / (input_path.stem + '.csv')
    }


def ensure_directories(paths: dict):
    for key in ['output_stats_path', 'ani_tab_file', 'species_output_file']:
        paths[key].parent.mkdir(parents=True, exist_ok=True)


def cleanup_files(paths: dict):
    for key in ['output_stats_path', 'ani_tab_file', 'output_final_file']:
        try:
            if key in paths and paths[key].exists() and paths[key].is_file():
                paths[key].unlink()
        except Exception as e:
            logging.warning(f"Error deleting {key}: {e}")


def remove_empty_directories(paths: dict):
    directories_to_check = ['output_stats_path', 'ani_tab_file', 'results_dir']
    for key in directories_to_check:
        try:
            if key in paths and not any(paths[key].parent.iterdir()):
                paths[key].parent.rmdir()
        except Exception as e:
            logging.warning(f"Error deleting directory {paths[key].parent}: {e}")


def run_species_metrics_finder(single_sequence_path: str, threshold_species: float = 0.95):
    setup_logging()

    input_path = Path(single_sequence_path)
    paths = generate_paths(input_path)
    ensure_directories(paths)

    try:
        # Run assembly statistics
        fasta_stats = FastaStatistics(single_sequence_path)
        stats = fasta_stats.generate_assembly_statistics()
        df_stats = pd.DataFrame([stats])
        df_stats.to_csv(paths['output_stats_path'], index=False)

        # Prepare for ANI analysis
        ani_executor = LocalANIExecutor(single_sequence_path, REFERENCE_GENOMES, paths['ani_tab_file'])
        ani_executor.execute()
        logging.info("Local ANI analysis completed.")

        # Assign species to the sequences
        logging.debug("Modifying species tab...")
        mapping = SpeciesTabModifier(paths['ani_tab_file'])
        paths['output_final_file'] = paths['ani_tab_file'].with_name(paths['ani_tab_file'].stem + '_final.tab')
        mapping.modify_tab_file(paths['output_final_file'])
        species_list = mapping.check_species_above_threshold(paths['output_final_file'], threshold=threshold_species)

        # Add species information to the final DataFrame
        final_df = pd.read_csv(paths['output_stats_path'])
        final_df['Species'] = ','.join(species_list)
        final_df.to_csv(paths['species_output_file'], index=False)
        logging.info(f"Species assignment completed. Results saved to {paths['species_output_file']}")
    finally:
        cleanup_files(paths)
        remove_empty_directories(paths)
        # Check if the 'results' directory itself is empty and delete it
        try:
            if not any(paths['results_dir'].iterdir()):
                paths['results_dir'].rmdir()
            logging.info("Empty directories deleted.")
        except Exception as e:
            logging.warning(f"Error deleting results directory: {e}")


if __name__ == '__main__':
    # Example call to the function with parameters
    single_sequence_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/test_2/genomes2/PRR28_CAB/PRR28_CAB.fasta'
    run_species_metrics_finder(single_sequence_path, threshold_species=0.95)
