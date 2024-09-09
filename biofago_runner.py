import os
import shutil
import subprocess
import sys
import tempfile

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))

import argparse
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from typing import List, Dict

import pandas as pd

from extract_annotate_assign import extract_annotate_assign
from utils.config import Config
from utils.folder_csv_manager import (create_individual_folders, run_species_metrics_for_all,
                                      create_species_finder_folder, cleanup_unwanted_species_folders)
from utils.genome_processing import process_genome, write_results_to_csv, cleanup_analysis_folders, keep_loci_files
from utils.logging_config import setup_logging


def check_docker():
    try:
        subprocess.run(["docker", "info"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError:
        return False


def copy_types_folders(source_folder: Path, destination_folder: Path):
    types_folders = ['types_capsule', 'types_cellulose', 'types_lps', 'types_srl']
    for genome_folder in source_folder.iterdir():
        if genome_folder.is_dir():
            for types_folder in types_folders:
                source_path = genome_folder / types_folder
                if source_path.exists():
                    dest_path = destination_folder / genome_folder.name / types_folder
                    shutil.copytree(source_path, dest_path, dirs_exist_ok=True)
    logging.info(f"Copied types folders to {destination_folder}")

def run_species_and_types_finder(genomes_folder: Path, threshold_species: float = 0.95,
                                 keep_sequence_loci: bool = False) -> None:
    log_file = genomes_folder / "process.log"
    setup_logging(log_file)
    logging.info(f"Starting analysis for genomes in {genomes_folder}")

    if not check_docker():
        logging.error("Docker is not running or not accessible. Please start Docker and try again.")
        sys.exit(1)

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

            # Copy types folders to the output directory
        copy_types_folders(genomes_folder_path, species_finder_path)

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


def main_yml():
    try:
        # Dynamically determine the path to the configuration file
        config_path = Path(__file__).resolve().parent.parent / 'config.yaml'
        logging.info(f"Reading configuration from: {config_path}")

        # Create Config instance
        config = Config(Path(config_path))

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


def is_fasta(file_path):
    """Check if a file is in FASTA format."""
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            return first_line.startswith('>')
    except:
        return False


def galaxy_runner():
    parser = argparse.ArgumentParser(description='BioFago Erwinia Analysis')
    parser.add_argument('--input', required=True, help='Input genome file or folder containing genome files')
    parser.add_argument('--threshold_species', type=float, default=0.95, help='ANI threshold for species assignment')
    parser.add_argument('--keep_sequence_loci', default=False, action='store_true', help='Keep sequence loci files')
    parser.add_argument('--log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help='Logging level')
    parser.add_argument('--output_dir', required=True, help='Output directory')

    args = parser.parse_args()

    try:

        if not check_docker():
            print("Error: Docker is not running or not accessible. Please start Docker and try again.")
            sys.exit(1)

        input_path = Path(args.input)
        output_dir = Path(args.output_dir)
        output_dir.mkdir(exist_ok=True)
        log_file = output_dir / "process.log"

        setup_logging(log_file, args.log_level)

        # Create a temporary directory to hold input files
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_genomes_folder = Path(temp_dir) / "genomes"
            temp_genomes_folder.mkdir(exist_ok=True)

            if input_path.is_file():
                if is_fasta(input_path):
                    shutil.copy2(input_path, temp_genomes_folder)
                else:
                    logging.error(f"Input file is not a FASTA file: {input_path}")
                    sys.exit(1)
            elif input_path.is_dir():
                fasta_files = [f for f in input_path.glob('*') if f.is_file() and is_fasta(f)]
                if not fasta_files:
                    logging.error(f"No FASTA files found in the input directory: {input_path}")
                    sys.exit(1)
                for fasta_file in fasta_files:
                    shutil.copy2(fasta_file, temp_genomes_folder)
            else:
                logging.error(f"Input is neither a file nor a directory: {input_path}")
                sys.exit(1)

            logging.info(f"Starting analysis with genomes folder: {temp_genomes_folder}")
            logging.info(f"Results will be saved in: {output_dir}")

            # Run the analysis
            run_species_and_types_finder(temp_genomes_folder,
                                         threshold_species=args.threshold_species,
                                         keep_sequence_loci=args.keep_sequence_loci)

            # Copy results from temp directory to output directory
            for item in Path(temp_dir).glob('*'):
                if item.is_dir():
                    shutil.copytree(item, output_dir / item.name, dirs_exist_ok=True)
                else:
                    shutil.copy2(item, output_dir)

    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        logging.exception("Exception details:")
        sys.exit(1)


if __name__ == '__main__':
    galaxy_runner()

