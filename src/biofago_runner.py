import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from typing import List, Dict

import pandas as pd

from extract_annotate_assign import extract_annotate_assign
from utils.folder_csv_manager import (create_individual_folders, run_species_metrics_for_all,
                                          create_species_finder_folder, cleanup_unwanted_species_folders)
from utils.genome_processing import process_genome, write_results_to_csv, cleanup_analysis_folders, keep_loci_files
from utils.logging_config import setup_logging
from utils.config import Config


def run_species_and_types_finder(genomes_folder: Path, threshold_species: float = 0.95,
                                 keep_sequence_loci: bool = False) -> None:
    log_file = genomes_folder / "process.log"
    setup_logging(log_file)
    logging.info(f"Starting analysis for genomes in {genomes_folder}")

    try:
        genomes_folder_path = create_individual_folders(genomes_folder)
        base_output_path = genomes_folder.parent
        species_finder_path = create_species_finder_folder(base_output_path)

        run_species_metrics_for_all(genomes_folder_path, species_finder_path, threshold_species)

        process_genome_with_path = partial(process_genome)
        with ProcessPoolExecutor() as executor:
            genome_dirs = [d for d in genomes_folder_path.iterdir() if d.is_dir()]
            results = list(executor.map(process_genome_with_path, genome_dirs))

        extract_annotate_results = extract_annotate_assign(genomes_folder_path)
        logging.debug(f"Extract and annotate results: {extract_annotate_results}")

        final_results = process_results(results, species_finder_path, extract_annotate_results)

        if final_results:
            output_csv = species_finder_path / "all_results.csv"
            write_results_to_csv(final_results, output_csv)
            cleanup_individual_csv_files(species_finder_path, output_csv)
            logging.info(f"Final results written to {output_csv}")
        else:
            logging.error("No genomes were processed successfully.")

        cleanup_unwanted_species_folders(base_output_path, species_finder_path)

        if keep_sequence_loci:
            keep_loci_files(genomes_folder_path)
        else:
            cleanup_analysis_folders(genomes_folder_path)

    except Exception as e:
        logging.error(f"Error in run_species_and_types_finder: {str(e)}")
        logging.exception("Exception details:")


def process_results(results: List[Dict], species_finder_path: Path, extract_annotate_results: List) -> List[Dict]:
    final_results = []
    for result in results:
        if result:
            genome_name = result['name']
            species_csv = species_finder_path / f"{genome_name}.csv"
            if species_csv.exists():
                df = pd.read_csv(species_csv)
                result.update({
                    'species': df['Species'].iloc[0] if 'Species' in df.columns else 'Unknown',
                    'ANI_species': round(df['ANI'].iloc[0], 2) if 'ANI' in df.columns else 0.0,
                })

            result.update({
                'capsule_locus': 'Unknown',
                'cellulose_locus': 'Unknown',
                'lps_locus': 'Unknown',
                'sorbitol_locus': 'Unknown'
            })

            update_locus_information(result, extract_annotate_results)
            clean_crispr_info(result)
            final_results.append(result)
    return final_results


def update_locus_information(result: Dict, extract_annotate_results: List) -> None:
    for annotation_result in extract_annotate_results:
        if len(annotation_result) == 4:
            annotation_genome, locus_type, final_type_locus, final_type = annotation_result
            if annotation_genome == result['name']:
                if locus_type == 'types_capsule':
                    result['capsule_locus'] = f"{final_type_locus} ({final_type})"
                elif locus_type == 'types_cellulose':
                    result['cellulose_locus'] = f"{final_type_locus} ({final_type})"
                elif locus_type == 'types_lps':
                    result['lps_locus'] = f"{final_type_locus} ({final_type})"
                elif locus_type == 'types_srl':
                    result['sorbitol_locus'] = f"{final_type_locus} ({final_type})"


def clean_crispr_info(result: Dict) -> None:
    if 'crispr_spacers' in result:
        result['crispr_spacers'] = result['crispr_spacers'].replace('CRISPR: ', '')
    if 'crispr_genotype' in result:
        result['crispr_genotype'] = result['crispr_genotype'].replace('CRR: ', '')


def cleanup_individual_csv_files(species_finder_path: Path, output_csv: Path) -> None:
    for csv_file in species_finder_path.glob('*.csv'):
        if csv_file != output_csv:
            csv_file.unlink()


if __name__ == '__main__':
    try:
        # Dynamically determine the path to the configuration file
        config_path = Path(__file__).resolve().parent.parent / 'config.yaml'
        logging.info(f"Reading configuration from: {config_path}")

        # Create Config instance
        config = Config(str(config_path))

        # Setup logging
        log_file = Path(config.output_folder) / "process.log" if config.output_folder else Path("process.log")
        setup_logging(log_file, config.log_level)

        genomes_folder = Path(config.genomes_folder)
        logging.info(f"Starting analysis with genomes folder: {genomes_folder}")

        run_species_and_types_finder(genomes_folder,
                                     threshold_species=config.threshold_species,
                                     keep_sequence_loci=config.keep_sequence_loci)
    except FileNotFoundError as e:
        print(f"Configuration file error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

