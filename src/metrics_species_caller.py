# At the top of metrics_species_caller.py
import logging
import shutil
from pathlib import Path
import pandas as pd
import time
import subprocess
from typing import List

from assigning_types.assembly_statistics import FastaStatistics, logger
from assigning_types.SpeciesFinderWPyani import NewSpeciesTabModifier, OptimizedLocalANIExecutor

# Dynamically determine the path to the reference types directory
SCRIPT_DIR = Path(__file__).resolve().parent
REFERENCE_GENOMES = SCRIPT_DIR.parents[0] / 'reference_real_species_genomes'

# Add logging to debug path resolution
logging.info(f"SCRIPT_DIR: {SCRIPT_DIR}")
logging.info(f"REFERENCE_GENOMES: {REFERENCE_GENOMES}")
logging.info(f"REFERENCE_GENOMES exists: {REFERENCE_GENOMES.exists()}")


def setup_logging():
    """Do not initialize logging again if it's already configured"""
    if not logger.handlers:  # Only add handler if none exists
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        if not logger.level:  # Only set level if not already set
            logger.setLevel(logging.INFO)


def generate_paths(input_path: Path, species_finder_path: Path, output_name: str = None) -> dict:
    results_dir = species_finder_path / 'results'
    
    # If output_name is provided, use it instead of input_path filename without extension
    if output_name:
        file_name = output_name
    else:
        # Get the filename without the very last extension only
        filename = input_path.name
        file_name = filename.rsplit('.', 1)[0] if '.' in filename else filename
    
    return {
        "results_dir": results_dir,
        "output_stats_path": results_dir / 'stats' / 'stats_results.csv',
        "ani_tab_file": results_dir / 'ani_tab' / 'ani_results.tab',
        "species_output_file": species_finder_path / f"{file_name}.csv"
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
                                   threshold_species: float = 0.95, skip_species_assignment: bool = False,
                                   output_name: str = None):
    setup_logging()
    logger.info(f"[Species Analysis] Starting analysis for {single_sequence_path}")
    logger.info(f"[Species Analysis] Skip assignment: {skip_species_assignment}")
    if output_name:
        logger.info(f"[Species Analysis] Using custom output name: {output_name}")

    input_path = Path(single_sequence_path).resolve()
    species_finder_path = Path(species_finder_path).resolve()
    logger.info(f"[Species Analysis] Input path resolved to: {input_path}")
    logger.info(f"[Species Analysis] Species finder path resolved to: {species_finder_path}")

    # Use custom output name if provided
    paths = generate_paths(input_path, species_finder_path, output_name)
    logger.info(f"[Species Analysis] Generated paths: {paths}")

    try:
        # Run assembly statistics
        logger.info(f"[Species Analysis] Running assembly statistics")
        fasta_stats = FastaStatistics(single_sequence_path)
        stats = fasta_stats.generate_assembly_statistics()
        df_stats = pd.DataFrame([stats])

        # Create a directory for the genome in the species finder folder
        # Use the output_name if provided, otherwise clean the input path stem
        if output_name:
            dir_name = output_name
        else:
            # Use the full filename without just the final extension
            dir_name = input_path.name
            dir_name = os.path.splitext(dir_name)[0]
        
        genome_dir = species_finder_path / dir_name
        genome_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(input_path, genome_dir / input_path.name)
        logger.info(f"[Species Analysis] Created and copied to {genome_dir} (using name: {dir_name})")

        if skip_species_assignment:
            logger.info(f"[Species Analysis] Species assignment skipped for {input_path.stem}")
            species = "Not Analyzed (Species Assignment Skipped)"
            ani_value = 0.0
        else:
            try:
                logger.info("[Species Analysis] Starting ANI analysis")
                ani_executor = OptimizedLocalANIExecutor(
                    single_sequence_path=genome_dir / input_path.name,
                    genomes_directory=REFERENCE_GENOMES,
                    results_file=paths['ani_tab_file'],
                    threshold=threshold_species
                )
                logger.info(f"[Species Analysis] Reference genomes path: {REFERENCE_GENOMES}")

                ani_success = ani_executor.execute()

                # ----------------------------------------------------
                if ani_success:
                    species = ani_executor.best_match['species']
                    ani_value = ani_executor.best_match['ani']
                else:
                    species = "ANI Analysis Failed"
                    ani_value = 0.0
                # ----------------------------------------------------

                logger.info(f"[Species Analysis] ANI analysis completed for {input_path.stem} with status: {ani_success}")
            except Exception as e:
                logger.error(f"[Species Analysis] Error in ANI: {str(e)}")
                logger.exception("[Species Analysis] Exception details:")
                species = "ANI Analysis Error"
                ani_value = 0.0

        df_stats['Species'] = species
        df_stats['ANI'] = ani_value
        logger.info(f"[Species Analysis] Added species info: {species}, ANI: {ani_value}")

        # Save results
        if paths['species_output_file'].exists():
            paths['species_output_file'].unlink()
            logger.info("[Species Analysis] Removed existing output file")

        df_stats.to_csv(paths['species_output_file'], index=False)
        logger.info(f"[Species Analysis] Results saved to {paths['species_output_file']}")

    except Exception as e:
        logger.error(f"[Species Analysis] Error: {str(e)}")
        logger.exception("[Species Analysis] Full exception details:")
        df_stats = pd.DataFrame([{
            'Species': 'Unknown',
            'ANI': 0.0,
            'Error': str(e)
        }])
        df_stats.to_csv(paths['species_output_file'], index=False)
        logger.info("[Species Analysis] Saved error results")

    finally:
        logger.info("[Species Analysis] Starting cleanup")
        cleanup_files(paths)
        remove_empty_directories(paths)
        logger.info("[Species Analysis] Cleanup completed")

    return df_stats


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



