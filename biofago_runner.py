import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from typing import List, Dict, Optional
import multiprocessing
from functools import partial

import argparse
import logging
import pandas as pd

# Add project root and src to sys.path
project_root = Path(__file__).resolve().parent
src_path = project_root / 'src'
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(src_path))

matrix_path = os.path.join(project_root, 'reference_crispr', 'clades_groups_correlation.csv')

from extract_annotate_assign import extract_annotate_assign
from utils.folder_csv_manager import (
    create_individual_folders,
    run_species_metrics_for_all,
    create_species_finder_folder,
    cleanup_unwanted_species_folders
)
from utils.genome_processing import process_single_genome, write_results_to_csv, keep_loci_files
from utils.clade_assigner import CRISPRCladeClassifier, parse_spacer_counts


def setup_logging(log_file: Path, log_level: str) -> None:
    """Set up logging configuration."""
    logging.basicConfig(
        filename=log_file,
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    )


def create_temp_dir():
    temp_dir = tempfile.mkdtemp()
    os.chmod(temp_dir, 0o755)  # Ensures read and execute permissions for all users
    return Path(temp_dir)


def check_docker() -> bool:
    """Check if Docker is running and accessible."""
    try:
        subprocess.run(["docker", "info"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError:
        return False


def check_blast_version() -> bool:
    """Check if BLAST is installed and get its version."""
    try:
        result = subprocess.run(['blastn', '-version'], capture_output=True, text=True, check=True)
        version = result.stdout.strip()
        logging.info(f"BLAST version: {version}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Error checking BLAST version: {e}")
        return False
    except FileNotFoundError:
        logging.error("BLAST is not installed or not in the system PATH.")
        return False


def is_fasta(file_path: Path) -> bool:
    """Check if a file is in FASTA format."""
    try:
        with open(file_path, 'rb') as f:
            first_line = f.readline().decode('utf-8', errors='ignore').strip()
            return first_line.startswith('>')
    except Exception as e:
        logging.error(f"Error checking if file is FASTA: {file_path}: {e}")
        return False


def clean_fasta_name(file_path: Path) -> Path:
    """Clean up FASTA file name, ensuring it has a .fasta extension."""
    path = Path(file_path)
    while path.suffix.lower() in ['.fasta', '.fa', '.fna']:
        path = path.with_suffix('')
    return path.with_suffix('.fasta')


def copy_fasta_files(input_path: Path, temp_genomes_folder: Path) -> List[Path]:
    """Copy FASTA files to a temporary folder."""
    fasta_files = []
    if input_path.is_file():
        if is_fasta(input_path):
            fasta_files = [input_path]
        else:
            raise ValueError(f"Input file is not a FASTA file: {input_path}")
    elif input_path.is_dir():
        try:
            for entry in os.scandir(str(input_path)):
                if entry.is_file():
                    file_path = Path(entry.path)
                    if is_fasta(file_path):
                        fasta_files.append(file_path)
                    else:
                        logging.warning(f"Skipping non-FASTA file: {file_path}")
        except Exception as e:
            logging.error(f"Error reading directory {input_path}: {e}")
            raise
        if not fasta_files:
            raise ValueError(f"No FASTA files found in the input directory: {input_path}")
    else:
        raise ValueError(f"Input is neither a file nor a directory: {input_path}")

    for fasta_file in fasta_files:
        dest_file = clean_fasta_name(temp_genomes_folder / fasta_file.name)
        try:
            shutil.copy2(fasta_file, dest_file)
            logging.info(f"Copied and renamed {fasta_file} to {dest_file}")
        except Exception as e:
            logging.error(f"Error copying {fasta_file} to {dest_file}: {e}")

    return [f for f in temp_genomes_folder.glob('*.fasta')]


def copy_types_folders(source_folder: Path, destination_folder: Path) -> None:
    """Copy type-specific folders to the destination."""
    types_folders = ['types_capsule', 'types_cellulose', 'types_lps', 'types_srl']
    for genome_folder in source_folder.iterdir():
        if genome_folder.is_dir():
            for types_folder in types_folders:
                source_path = genome_folder / types_folder
                if source_path.exists():
                    dest_path = destination_folder / genome_folder.name / types_folder
                    shutil.copytree(source_path, dest_path, dirs_exist_ok=True)
    logging.info(f"Copied types folders to {destination_folder}")


def cleanup_analysis_folders(output_dir: Path, keep_sequence_loci: bool) -> None:
    """Clean up analysis folders and remove empty directories."""
    folders_to_remove = ['CRISPR_finder', 'CRR_finder', 'plasmid_finder', 'resistance_finder']
    if not keep_sequence_loci:
        folders_to_remove.append('types_finder')

    for folder in folders_to_remove:
        folder_path = output_dir / folder
        if folder_path.exists():
            try:
                shutil.rmtree(folder_path)
                logging.info(f"Removed folder: {folder_path}")
            except Exception as e:
                logging.error(f"Error removing folder {folder_path}: {str(e)}")

    # Clean up empty 'results' folder in species_finder
    species_finder_path = output_dir / 'species_finder'
    if species_finder_path.exists():
        results_folder = species_finder_path / 'results'
        if results_folder.exists() and not any(results_folder.iterdir()):
            try:
                results_folder.rmdir()
                logging.info(f"Removed empty folder: {results_folder}")
            except Exception as e:
                logging.error(f"Error removing empty folder {results_folder}: {str(e)}")

    logging.info(f"Cleanup completed. Remaining contents of {output_dir}: {list(output_dir.glob('*'))}")


def process_results(results: List[Dict], species_finder_path: Path, extract_annotate_results: List,
                    clade_classifier: CRISPRCladeClassifier) -> List[Dict]:
    final_results = []
    processed_genomes = set()
    logging.info(f"Starting process_results with {len(results)} results")

    for result in results:
        if result and result['name'] not in processed_genomes:
            genome_name = result['name']
            processed_genomes.add(genome_name)

            # Process species information
            species_csv = species_finder_path / f"{genome_name}.csv"
            if species_csv.exists():
                try:
                    df = pd.read_csv(species_csv)
                    result.update({
                        'species': df['Species'].iloc[0] if 'Species' in df.columns else 'Unknown',
                        'ANI_species': round(df['ANI'].iloc[0], 2) if 'ANI' in df.columns else 0.0,
                    })
                    species_csv.unlink()
                except Exception as e:
                    logging.error(f"Error reading species CSV for {genome_name}: {str(e)}")
                    result.update({'species': 'Unknown', 'ANI_species': 0.0})

            # Initialize all locus types with unknown status
            locus_types = [
                'capsule', 'cellulose', 'lps', 'sorbitol',
                'flag_i', 'flag_ii', 'flag_iii', 'flag_iv',
                't3ss_i', 't3ss_ii', 't6ss_i', 't6ss_ii'
            ]

            for locus_type in locus_types:
                result[f'{locus_type}_locus'] = '(Unknown)'

            # Process locus information from extract_annotate_results
            locus_types_mapping = {
                'types_capsule': 'capsule_locus',
                'types_cellulose': 'cellulose_locus',
                'types_lps': 'lps_locus',
                'types_srl': 'sorbitol_locus',
                'types_flag_I': 'flag_i_locus',
                'types_flag_II': 'flag_ii_locus',
                'types_flag_III': 'flag_iii_locus',
                'types_flag_IV': 'flag_iv_locus',
                'types_T3SS_I': 't3ss_i_locus',
                'types_T3SS_II': 't3ss_ii_locus',
                'types_T6SS_I': 't6ss_i_locus',
                'types_T6SS_II': 't6ss_ii_locus'
            }

            for annotation_result in extract_annotate_results:
                if len(annotation_result) == 6:
                    annotation_genome, reference_type, final_type_locus, final_type, flagged_genes, _ = annotation_result
                    if annotation_genome == result['name']:
                        locus_key = locus_types_mapping.get(reference_type)
                        if locus_key:
                            if final_type_locus:
                                formatted_locus = f"{final_type_locus} ({final_type})"
                            else:
                                formatted_locus = f"({final_type})"
                            if flagged_genes:
                                formatted_locus += f" - Flagged genes: {flagged_genes}"
                            result[locus_key] = formatted_locus.strip()

            # Process other information (plasmids, CRISPR, etc.)
            if 'present_plasmids' not in result:
                result['present_plasmids'] = 'None'

            # Process virulence genes
            virulence_columns = [col for col in result if col.endswith('_genes')]
            for col in virulence_columns:
                if result[col] != 'None':
                    result[col] = ', '.join(sorted(result[col].split(', ')))

            # Process CRISPR clade classification
            genotype = result.get('crispr_genotype', '')
            spacers = result.get('crispr_spacers', '')
            spacer_counts = parse_spacer_counts(spacers)
            clade, score, spacer_score, genotype_score, confidence_level, subgroup = clade_classifier.determine_clade(
                genotype, spacer_counts)

            result.update({
                'clade': f"{clade} {subgroup}".strip() if clade != "Unknown" else clade,
                'clade_confidence_score': round(score, 2),
                'clade_confidence_level': confidence_level
            })

            final_results.append(result)
            logging.info(f"Added result for {genome_name} to final_results")

    logging.info(f"Finished process_results. Final results count: {len(final_results)}")
    return final_results


def update_locus_information(result: Dict, extract_annotate_results: List) -> None:
    """Update locus information in the result dictionary."""
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
    """Clean up CRISPR information in the result dictionary."""
    if 'crispr_spacers' in result:
        result['crispr_spacers'] = result['crispr_spacers'].replace('CRISPR: ', '')
    if 'crispr_genotype' in result:
        result['crispr_genotype'] = result['crispr_genotype'].replace('CRR: ', '')


def run_species_and_types_finder(genomes_folder: Path, output_dir: Path, threshold_species: float,
                                 keep_sequence_loci: bool) -> None:
    """Main function to run species and types finder analysis."""
    logging.info(f"Starting analysis for genomes in {genomes_folder}")
    logging.info(f"Output directory: {output_dir}")

    try:
        species_finder_path = create_species_finder_folder(output_dir)

        genome_files = list(genomes_folder.glob('*.fasta'))
        if not genome_files:
            raise ValueError(f"No FASTA files found in {genomes_folder}")

        logging.info(f"Found genome files: {genome_files}")

        logging.info("Starting species metrics analysis")
        run_species_metrics_for_all(genomes_folder, species_finder_path, threshold_species)
        logging.info("Completed species metrics analysis")

        # Initialize the CRISPRCladeClassifier
        clade_classifier = CRISPRCladeClassifier(matrix_path)

        with ProcessPoolExecutor() as executor:
            # Create a list of tuples (genome_file, clade_classifier)
            genome_classifier_pairs = [(genome_file, clade_classifier) for genome_file in genome_files]
            results = list(executor.map(process_genome_with_classifier, genome_classifier_pairs))

        logging.info(f"Processed {len(results)} genomes")
        logging.info(f"Results preview: {results[:2]}")  # Log first two results for debugging

        extract_annotate_results = extract_annotate_assign(genomes_folder)
        logging.info(f"Extract and annotate results: {extract_annotate_results}")

        final_results = process_results(results, species_finder_path, extract_annotate_results, clade_classifier)

        if final_results:
            output_csv = species_finder_path / "all_results.csv"
            write_results_to_csv(final_results, output_csv)
            logging.info(f"Final results written to {output_csv}")
        else:
            logging.error("No genomes were processed successfully.")

        logging.info("Cleaning up analysis folders...")
        logging.info(f"Contents of output directory before cleanup: {list(output_dir.glob('*'))}")
        if keep_sequence_loci:
            keep_loci_files(output_dir)
            logging.info("Kept loci files as requested.")
        else:
            cleanup_analysis_folders(output_dir, keep_sequence_loci)
        logging.info(f"Contents of output directory after cleanup: {list(output_dir.glob('*'))}")

        logging.info("Analysis completed successfully.")


    except Exception as e:

        logging.error(f"Error in run_species_and_types_finder: {str(e)}")

        logging.exception("Exception details:")

        raise


def copy_final_results(temp_dir: Path, output_dir: Path, keep_sequence_loci: bool) -> None:
    """Copy final results from temporary directory to output directory."""
    species_finder_src = temp_dir / 'species_finder'
    if species_finder_src.exists():
        species_finder_dest = output_dir / 'species_finder'
        species_finder_dest.mkdir(parents=True, exist_ok=True)

        all_results_src = species_finder_src / 'all_results.csv'
        if all_results_src.exists():
            shutil.copy2(all_results_src, species_finder_dest / 'all_results.csv')
            logging.info(f"Copied all_results.csv to {species_finder_dest}")
        else:
            logging.warning("all_results.csv not found in temporary directory")

    if keep_sequence_loci:
        types_finder_src = temp_dir / 'types_finder'
        if types_finder_src.exists():
            shutil.copytree(types_finder_src, output_dir / 'types_finder', dirs_exist_ok=True)
            logging.info(f"Copied types_finder folder to {output_dir / 'types_finder'}")
        else:
            logging.warning(f"types_finder folder not found in {temp_dir}")

    for file in temp_dir.glob('*.log'):
        shutil.copy2(file, output_dir)

    logging.info(f"Copied final results to {output_dir}")
    logging.info(f"Final contents of output directory: {list(output_dir.glob('*'))}")

    types_finder_output = output_dir / 'types_finder'
    if types_finder_output.exists():
        logging.info("Contents of types_finder folder:")
        for root, dirs, files in os.walk(types_finder_output):
            for file in files:
                logging.info(f"  {os.path.join(root, file)}")

def process_genome_with_classifier(args):
    genome_file, clade_classifier = args
    return process_single_genome(genome_file, clade_classifier)


def galaxy_runner():
    parser = argparse.ArgumentParser(description='BioFago Erwinia Analysis')
    parser.add_argument('--input', required=True, help='Input genome file or folder containing genome files')
    parser.add_argument('--threshold_species', type=float, default=0.95, help='ANI threshold for species assignment')
    parser.add_argument('--keep_sequence_loci', action='store_true', help='Keep sequence loci files')
    parser.add_argument('--log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help='Logging level')
    parser.add_argument('--output_dir', required=True, help='Output directory')

    args = parser.parse_args()

    try:
        input_path = Path(args.input).resolve()
        output_dir = Path(args.output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        log_file = output_dir / "process.log"

        setup_logging(log_file, args.log_level)

        logging.info(f"Input path: {input_path}")
        logging.info(f"Input path exists: {input_path.exists()}")
        logging.info(f"Input path is file: {input_path.is_file()}")
        logging.info(f"Input path is dir: {input_path.is_dir()}")

        if not input_path.exists():
            logging.error(f"Input file or directory not found: {input_path}")
            sys.exit(1)

        if not check_docker():
            logging.error("Docker is not running or not accessible. Please start Docker and try again.")
            sys.exit(1)

        if not check_blast_version():
            logging.error("BLAST is not installed or not accessible. Please install BLAST and try again.")
            sys.exit(1)

        with create_temp_dir() as temp_dir:
            temp_dir_path = Path(temp_dir)
            temp_genomes_folder = temp_dir_path / "genomes"
            temp_genomes_folder.mkdir(exist_ok=True)

            try:
                logging.info(f"Attempting to copy FASTA files from {input_path} to {temp_genomes_folder}")
                fasta_files = copy_fasta_files(input_path, temp_genomes_folder)
                logging.info(f"Copied {len(fasta_files)} FASTA files to temporary folder")
            except ValueError as e:
                logging.error(str(e))
                sys.exit(1)
            except Exception as e:
                logging.error(f"Unexpected error in copy_fasta_files: {str(e)}")
                logging.exception("Exception details:")
                sys.exit(1)

            logging.info(f"Starting analysis with genomes folder: {temp_genomes_folder}")
            logging.info(f"Results will be saved in: {output_dir}")

            try:
                run_species_and_types_finder(temp_genomes_folder, temp_dir_path,
                                             threshold_species=args.threshold_species,
                                             keep_sequence_loci=args.keep_sequence_loci)
            except Exception as e:
                logging.error(f"Error in run_species_and_types_finder: {str(e)}")
                logging.exception("Exception details:")
                sys.exit(1)

            logging.info("Copying final results to output directory")
            try:
                copy_final_results(temp_dir_path, output_dir, args.keep_sequence_loci)
            except Exception as e:
                logging.error(f"Error in copy_final_results: {str(e)}")
                logging.exception("Exception details:")
                sys.exit(1)

        # Remove empty 'results' folder in species_finder
        species_finder_results = output_dir / 'species_finder' / 'results'
        if species_finder_results.exists() and not any(species_finder_results.iterdir()):
            try:
                shutil.rmtree(species_finder_results)
                logging.info(f"Removed empty 'results' folder in species_finder")
            except Exception as e:
                logging.error(f"Error removing empty 'results' folder: {str(e)}")

        logging.info(f"Final contents of output directory: {list(output_dir.glob('*'))}")
        logging.info("Analysis completed successfully.")

    except Exception as e:
        logging.error(f"An unexpected error occurred: {str(e)}")
        logging.exception("Exception details:")
        sys.exit(1)


if __name__ == '__main__':
    galaxy_runner()
