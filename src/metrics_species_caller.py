import logging
import shutil
from pathlib import Path
import pandas as pd
import time
import subprocess
from typing import List
from assigning_types.assembly_statistics import FastaStatistics
from assigning_types.SpeciesFinderWPyani import NewSpeciesTabModifier, OptimizedLocalANIExecutor


# Dynamically determine the path to the reference types directory
SCRIPT_DIR = Path(__file__).resolve().parent
REFERENCE_GENOMES = SCRIPT_DIR.parents[0] / 'reference_real_species_genomes'

LOG_FORMAT = '%(asctime)s - %(levelname)s - %(message)s'
LOG_LEVEL = logging.DEBUG


def setup_logging():
    logging.basicConfig(level=LOG_LEVEL, format=LOG_FORMAT)


def generate_paths(input_path: Path, species_finder_path: Path) -> dict:
    results_dir = species_finder_path / 'results'
    return {
        "results_dir": results_dir,
        "output_stats_path": results_dir / 'stats' / 'stats_results.csv',
        "ani_tab_file": results_dir / 'ani_tab' / 'ani_results.tab',
        "species_output_file": species_finder_path / f"{input_path.stem}.csv"
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


def new_run_species_metrics_finder(single_sequence_path: Path, species_finder_path: Path,
                                   threshold_species: float = 0.95):
    setup_logging()

    input_path = Path(single_sequence_path)
    species_finder_path = Path(species_finder_path)

    paths = generate_paths(input_path, species_finder_path)
    ensure_directories(paths)

    try:
        # Run assembly statistics
        fasta_stats = FastaStatistics(single_sequence_path)
        stats = fasta_stats.generate_assembly_statistics()
        df_stats = pd.DataFrame([stats])

        # Prepare for ANI analysis
        ani_executor = OptimizedLocalANIExecutor(single_sequence_path, REFERENCE_GENOMES, paths['ani_tab_file'],
                                                 threshold_species)
        ani_success = ani_executor.execute()

        if ani_success:
            logging.info("Local ANI analysis completed.")
            if ani_executor.match_found:
                species = ani_executor.best_match['species']
                ani_value = ani_executor.best_match['ani']
                logging.info(f"Match found: {species} with ANI value {ani_value:.4f}")
            else:
                species = "Unknown"
                ani_value = 0.0
                logging.info(f"No match found above threshold. Best match: {ani_executor.best_match['species']} with ANI {ani_executor.best_match['ani']:.4f}")

            # Add species information to the DataFrame
            df_stats['Species'] = species
            df_stats['ANI'] = ani_value
        else:
            logging.warning("ANI analysis failed. Setting species as 'Unknown'.")
            df_stats['Species'] = 'Unknown'
            df_stats['ANI'] = 0.0

        # Ensure the output file doesn't exist before writing
        if paths['species_output_file'].exists():
            paths['species_output_file'].unlink()

        df_stats.to_csv(paths['species_output_file'], index=False)
        logging.info(f"Results saved to {paths['species_output_file']}")

    except Exception as e:
        logging.error(f"Error in run_species_metrics_finder: {e}")
        df_stats = pd.DataFrame([{
            'Species': 'Unknown',
            'ANI': 0.0,
            'Error': str(e)
        }])
        df_stats.to_csv(paths['species_output_file'], index=False)
        logging.info(f"Results saved to {paths['species_output_file']} with Unknown species due to error")

    finally:
        cleanup_files(paths)
        remove_empty_directories(paths)

# def new_run_species_metrics_finder(single_sequence_path: Path, species_finder_path: Path,
#                                    threshold_species: float = 0.95):
#     setup_logging()
#
#     input_path = Path(single_sequence_path)
#     species_finder_path = Path(species_finder_path)
#
#     paths = generate_paths(input_path, species_finder_path)
#     ensure_directories(paths)
#
#     try:
#         # Run assembly statistics
#         fasta_stats = FastaStatistics(single_sequence_path)
#         stats = fasta_stats.generate_assembly_statistics()
#         df_stats = pd.DataFrame([stats])
#
#         # Add placeholder species information
#         df_stats['Species'] = 'Not Analyzed'
#         df_stats['ANI'] = 0.0
#
#         logging.info("Skipping ANI analysis as requested.")
#
#         # Ensure the output file doesn't exist before writing
#         if paths['species_output_file'].exists():
#             paths['species_output_file'].unlink()
#
#         df_stats.to_csv(paths['species_output_file'], index=False)
#         logging.info(f"Results saved to {paths['species_output_file']}")
#
#     except Exception as e:
#         logging.error(f"Error in run_species_metrics_finder: {e}")
#         df_stats = pd.DataFrame([{
#             'Species': 'Unknown',
#             'ANI': 0.0,
#             'Error': str(e)
#         }])
#         df_stats.to_csv(paths['species_output_file'], index=False)
#         logging.info(f"Results saved to {paths['species_output_file']} with Unknown species due to error")
#
#     finally:
#         cleanup_files(paths)
#         remove_empty_directories(paths)
#
#     return df_stats


if __name__ == '__main__':

    pass


