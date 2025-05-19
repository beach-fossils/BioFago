import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import List, Dict, Optional
import multiprocessing
from functools import partial
import argparse
# Silence warnings first
import warnings
warnings.filterwarnings("ignore")

# Add project root and src to sys.path before importing our modules
project_root = Path(__file__).resolve().parent
src_path = project_root / 'src'
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(src_path))

# Import quiet_mode early to handle all console output suppression
from src.quiet_mode import QUIET_MODE, ORIGINAL_STDOUT, ORIGINAL_STDERR, NULL_OUTPUT

# Import logging AFTER output redirection
import logging
import pandas as pd
import tqdm
import io
import contextlib
import importlib


from src.assigning_types.assembly_statistics import FastaStatistics
from src.metrics_species_caller import new_run_species_metrics_finder, REFERENCE_GENOMES
from src.extract_annotate_assign import extract_annotate_assign
from src.utils.folder_csv_manager import (
    create_individual_folders,
    run_species_metrics_for_all,
    create_species_finder_folder,
    cleanup_unwanted_species_folders,
    final_enhance_csv
)
from src.utils.genome_processing import process_single_genome, write_results_to_csv, keep_loci_files, process_genomes
from src.utils.clade_assigner import CRISPRCladeClassifier, parse_spacer_counts
from src.resistance.levan_synthesis import (
    run_levan_analysis,
    REFERENCE_LEVAN,
    format_levan_result
)

matrix_path = os.path.join(project_root, 'reference_crispr', 'clades_groups_correlation.csv')


# Define a global guard that will detect quiet mode
QUIET_MODE = False

# Create a global null handler for all loggers when in quiet mode
class QuietHandler(logging.Handler):
    """A handler that does nothing in quiet mode, logs to file in normal mode."""
    def __init__(self, log_file, level=logging.INFO):
        super().__init__(level)
        self.log_file = log_file
        self.file_handler = None
        self.setup_file_handler()
        
    def setup_file_handler(self):
        """Set up a file handler for logging to file."""
        self.file_handler = logging.FileHandler(self.log_file)
        self.file_handler.setLevel(self.level)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        self.file_handler.setFormatter(formatter)
        
    def emit(self, record):
        """Only emit records to file, never to console in quiet mode."""
        # Always emit to file
        if self.file_handler:
            self.file_handler.emit(record)

# Store original stdout and stderr for restoration
ORIGINAL_STDOUT = sys.stdout
ORIGINAL_STDERR = sys.stderr

def setup_logging(log_file: Path, log_level: str, quiet: bool = False) -> None:
    """Set up logging configuration with complete suppression in quiet mode.
    
    Args:
        log_file: Path to the log file
        log_level: Log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        quiet: If True, suppress ALL console output except for progress bars
    """
    global QUIET_MODE
    
    # IMPORTANT: This is already set by quiet_mode.py, but we keep this for clarity
    QUIET_MODE = quiet
    
    # Get numeric log level
    numeric_level = getattr(logging, log_level.upper())
    
    # Reset the root logger (required to avoid duplicate handlers)
    root_logger = logging.getLogger()
    
    # Remove all existing handlers before attaching new ones
    for hdlr in root_logger.handlers[:]:
        hdlr.close()
        root_logger.removeHandler(hdlr)
    
    # Always set up file logging regardless of quiet mode
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(numeric_level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    
    # Set up the root logger with file handler
    root_logger.setLevel(numeric_level)
    root_logger.addHandler(file_handler)
    
    # In quiet mode, we don't add any console handlers
    if not quiet:
        # Add console handler for non-quiet mode
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter('%(levelname)s: %(message)s')
        console_handler.setFormatter(console_formatter)
        root_logger.addHandler(console_handler)
        
    # Monkey-patch subprocesses by redirecting their output in quiet mode
    if quiet:
        # Define a null handler for all loggers
        null_handler = logging.NullHandler()
        
        # Override all module loggers to use NullHandler for console output
        for logger_name in logging.root.manager.loggerDict:
            module_logger = logging.getLogger(logger_name)
            # Remove all console handlers
            for hdlr in module_logger.handlers[:]:
                if isinstance(hdlr, logging.StreamHandler) and not isinstance(hdlr, logging.FileHandler):
                    module_logger.removeHandler(hdlr)
            # Add the null handler
            module_logger.addHandler(null_handler)
            # Ensure propagation to the root logger (with file handler)
            module_logger.propagate = True
        
        # More aggressive monkey-patching of subprocess.run to ensure ALL output is suppressed
        original_run = subprocess.run
        def quiet_run(*args, **kwargs):
            # Handle the case when capture_output is used
            if kwargs.get('capture_output'):
                # Don't modify stdout/stderr if capture_output is True
                pass
            else:
                # Set stdout/stderr to PIPE only if capture_output is not used
                kwargs['stdout'] = kwargs.get('stdout', subprocess.PIPE)
                kwargs['stderr'] = kwargs.get('stderr', subprocess.PIPE)
            return original_run(*args, **kwargs)
        
        # Apply the monkey patch
        subprocess.run = quiet_run
        
        # Log that we're in quiet mode (to file only)
        logging.info("Logging initialized in quiet mode (console output suppressed)")
    else:
        # Log that we're in standard mode (to console and file)
        logging.info("Logging initialized in standard mode")


def create_temp_dir():
    temp_dir = tempfile.mkdtemp()
    os.chmod(temp_dir, 0o755)  # Ensures read and execute permissions for all users
    return Path(temp_dir)


def check_docker() -> bool:
    """Check if Docker is running and accessible."""
    try:
        # Use /dev/null for output redirection in quiet mode
        if QUIET_MODE:
            with open(os.devnull, 'w') as devnull:
                subprocess.run(["docker", "info"], check=True, stdout=devnull, stderr=devnull)
        else:
            subprocess.run(["docker", "info"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError:
        return False


def check_blast_version() -> bool:
    """Check if BLAST is installed and get its version."""
    try:
        # Use /dev/null for output redirection instead of capture_output
        if QUIET_MODE:
            with open(os.devnull, 'w') as devnull:
                result = subprocess.run(['blastn', '-version'], stdout=devnull, stderr=devnull, check=True)
            version = "BLAST is installed (output suppressed in quiet mode)"
        else:
            result = subprocess.run(['blastn', '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
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
    """Clean up FASTA file name, ensuring it has a .fasta extension.
    
    Simply removes the last extension and replaces with .fasta,
    preserving the full filename regardless of its format.
    """
    path = Path(file_path)
    
    # Only remove the last extension (suffix) and replace with .fasta
    # This preserves all parts of complex filenames without any special handling
    # Example: any_genome_name.fna -> any_genome_name.fasta
    base_name = path.stem  # Gets name without last extension
    return path.with_name(f"{base_name}.fasta")


def copy_fasta_files(input_path: Path, temp_genomes_folder: Path) -> List[Path]:
    """Copy FASTA files to a temporary folder with subdirectories."""
    fasta_files = []
    result_files = []

    # Collect FASTA files
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

    # Copy files to subdirectories
    for fasta_file in fasta_files:
        # Get the base name without just the last extension
        # This is important for GCF_ files with complex names
        file_basename = os.path.splitext(fasta_file.name)[0]
        
        # Create a subdirectory with the same name as the file (without extension)
        subdir_name = file_basename
        subdir_path = temp_genomes_folder / subdir_name
        subdir_path.mkdir(exist_ok=True)
        
        logging.info(f"Created subdirectory {subdir_path} for file {fasta_file.name}")

        # Copy the file to its subdirectory
        dest_file = clean_fasta_name(subdir_path / fasta_file.name)
        try:
            shutil.copy2(fasta_file, dest_file)
            logging.info(f"Copied {fasta_file.name} to {dest_file.parent.name}/{dest_file.name}")
            result_files.append(dest_file)
        except Exception as e:
            logging.error(f"Error copying {fasta_file} to {dest_file}: {e}")

    if not result_files:
        logging.error("No files were successfully copied")

    logging.info(f"Successfully copied {len(result_files)} files to subdirectories")

    # Double check the files exist
    for file in result_files:
        if not file.exists():
            logging.error(f"Expected file does not exist after copy: {file}")

    return result_files


def copy_types_folders(source_folder: Path, destination_folder: Path) -> None:
    """Copy type-specific folders to the destination."""
    types_folders = ['types_capsule', 'types_cellulose', 'types_lps', 'types_srl', 'types_ompa']
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
    logging.info(f"Number of extract_annotate_results: {len(extract_annotate_results)}")

    # Enhanced genome name mapping - preserving full names without splitting
    genome_name_mapping = {}
    for result in results:
        if result and 'name' in result:
            # Get the full name without any processing or splitting
            full_name = result['name']
            logging.info(f"Processing genome with name: {full_name}")
            
            # Only remove the final extension if present
            if any(full_name.lower().endswith(ext) for ext in ['.fasta', '.fa', '.fna']):
                full_name = os.path.splitext(full_name)[0]
                
            # Create additional mappings for shortened versions that might be used in other parts of the code
            # These shortened versions should ONLY be used for lookups, never for final output
            short_name_1 = full_name.split('.')[0]  # First part before any dot
            
            # Create a more specific mapping for cases like GCF_002732285.1_GCF_002732285.1_ASM273228v1_genomic
            # where we need to match both the full name and the shortened GCF_* name
            parts = full_name.split('_')
            if len(parts) > 2 and parts[0] in ('GCF', 'GCA'):
                # Create a mapping for GCF_number only (for lookup purposes)
                gcf_prefix = '_'.join(parts[:2])  # e.g., GCF_002732285
                genome_name_mapping[gcf_prefix] = full_name
                logging.debug(f"Added GCF prefix mapping: {gcf_prefix} -> {full_name}")
                
                # Another common pattern is GCF_number.version_GCF_number
                if len(parts) > 3 and '.' in parts[1] and parts[2] == parts[0]:
                    gcf_with_version = f"{parts[0]}_{parts[1]}"  # e.g., GCF_002732285.1
                    gcf_repeated = f"{gcf_with_version}_{parts[0]}_{parts[1]}"  # e.g., GCF_002732285.1_GCF_002732285.1
                    genome_name_mapping[gcf_with_version] = full_name
                    genome_name_mapping[gcf_repeated] = full_name
                    logging.debug(f"Added GCF version mapping: {gcf_with_version} -> {full_name}")
                    logging.debug(f"Added GCF repeated mapping: {gcf_repeated} -> {full_name}")
            
            # Add the main mappings
            genome_name_mapping[short_name_1] = full_name  # For lookups only
            genome_name_mapping[full_name] = full_name     # Keep original full name unchanged
            logging.debug(f"Added genome name mappings: {short_name_1} -> {full_name}, {full_name} -> {full_name}")

    for result in results:
        if not result:
            continue

        # Create a copy of the result to avoid modifying the original
        processed_result = result.copy()
        genome_name = processed_result.get('name')

        if not genome_name or genome_name in processed_genomes:
            continue

        logging.info(f"Processing results for genome: {genome_name}")
        processed_genomes.add(genome_name)

        try:
            # Process species information
            species_csv = species_finder_path / f"{genome_name}.csv"
            if species_csv.exists():
                try:
                    df = pd.read_csv(species_csv)
                    processed_result.update({
                        'species': df['Species'].iloc[0] if 'Species' in df.columns else 'Unknown',
                        'ANI_species': round(float(df['ANI'].iloc[0]), 2) if 'ANI' in df.columns else 0.0,
                    })
                except Exception as e:
                    logging.error(f"Error reading species CSV for {genome_name}: {str(e)}")
                    processed_result.update({'species': 'Unknown', 'ANI_species': 0.0})
            else:
                processed_result.update({'species': 'Unknown', 'ANI_species': 0.0})

            # Initialize locus fields
            locus_types = {
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
                'types_T6SS_II': 't6ss_ii_locus',
                'types_flag3': 'flag3_locus',
                'types_ompa': 'ompa_locus'  # Add ompA locus type
            }

            # Initialize all locus fields with default value
            for locus_key in locus_types.values():
                if locus_key not in processed_result:
                    processed_result[locus_key] = '(Unknown)'

            # Get both full and short names for matching
            full_genome_name = genome_name
            short_genome_name = genome_name.split('.')[0]

            # Process annotation results
            for annotation_result in extract_annotate_results:
                if len(annotation_result) >= 6:
                    ann_genome, ref_type, type_locus, final_type, flagged_genes, final_type_info = annotation_result
                    ann_genome_base = ann_genome.split('.')[0]

                    logging.debug(f"Comparing names - Ann: {ann_genome}, Full: {full_genome_name}, "
                                  f"Short: {short_genome_name}, Ann Base: {ann_genome_base}")

                    if ((ann_genome == full_genome_name or
                         ann_genome == short_genome_name or
                         ann_genome_base == short_genome_name) and
                            ref_type in locus_types):

                        locus_key = locus_types[ref_type]
                        
                        # Special handling for ompA which uses a different format
                        if ref_type == 'types_ompa':
                            # For ompA, final_type_info already contains the formatted string
                            if final_type_info:
                                formatted_output = final_type_info
                            else:
                                formatted_output = f"({final_type})"
                            processed_result[locus_key] = formatted_output
                            logging.info(f"Added ompA information for {genome_name}: {formatted_output}")
                        else:
                            # Regular loci formatting
                            formatted_output = f"({final_type})"
                            if type_locus and type_locus.strip():
                                formatted_output = f"{type_locus} {formatted_output}"

                            if flagged_genes:
                                if isinstance(flagged_genes, list):
                                    flagged_genes = ', '.join(flagged_genes)
                                formatted_output += f" - Flagged genes: {flagged_genes}"

                            processed_result[locus_key] = formatted_output.strip()
                        logging.info(f"Added {ref_type} information for {genome_name}: {formatted_output}")

            final_results.append(processed_result)
            logging.info(f"Successfully processed genome {genome_name} with {len(locus_types)} locus types")
            logging.debug(f"Final result for {genome_name}: {processed_result}")

        except Exception as e:
            logging.error(f"Error processing genome {genome_name}: {str(e)}")
            logging.exception("Exception details:")

    # Write and enhance CSV after all genomes are processed
    if final_results:
        output_csv = species_finder_path / "all_results.csv"
        write_results_to_csv(final_results, output_csv)
        # Move the enhance call to after the write is confirmed
        try:
            final_enhance_csv(str(output_csv))
            logging.info(f"Successfully enhanced CSV at {output_csv}")
        except Exception as e:
            logging.error(f"Error enhancing CSV: {str(e)}")
        logging.info(f"Wrote and enhanced {len(final_results)} results to {output_csv}")


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
                                 keep_sequence_loci: bool, clade_classifier: CRISPRCladeClassifier,
                                 skip_species_assignment: bool = False, batch_size: int = 0,
                                 num_workers: int = 4, docker_limit: int = 4, quiet: bool = False) -> None:
    try:
        logging.info("========= Starting Species and Types Analysis =========")
        logging.info(f"Skip species assignment: {skip_species_assignment}")
        logging.info(f"Genomes folder: {genomes_folder}")
        logging.info(f"Output directory: {output_dir}")
        
        # Create progress bar
        # Simplified context manager for progress bars - no longer needs special handling
        # since we're not redirecting stdout/stderr globally
        @contextlib.contextmanager
        def temporary_stdout(capture_to_log=False):
            """Simple pass-through context manager for progress bars and messages"""
            try:
                yield
            finally:
                # Ensure output is flushed immediately
                sys.stdout.flush()
                sys.stderr.flush()
                
        def create_progress_bar(total, desc="Processing", color="green"):
            """Create a colorful progress bar appropriate for the current mode.
            
            Args:
                total: Total number of items
                desc: Description text
                color: Color of the progress bar ('green', 'blue', 'cyan', etc.)
            """
            # Color formatting using ANSI escape codes
            colors = {
                'green': '\033[92m',
                'blue': '\033[94m',
                'cyan': '\033[96m',
                'magenta': '\033[95m',
                'yellow': '\033[93m',
                'white': '\033[97m',
                'reset': '\033[0m'
            }
            
            selected_color = colors.get(color, colors['green'])
            
            # Use the temporary_stdout context to ensure progress bars show in quiet mode
            with temporary_stdout():
                return tqdm.tqdm(
                    total=total,
                    desc=f"{selected_color}{desc}{colors['reset']}",
                    ncols=80,
                    bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]',
                    # Use stderr to avoid interference with logging 
                    file=sys.stderr
                )

        # Check if reference directory exists
        if not REFERENCE_LEVAN.exists():
            raise FileNotFoundError(f"Levan synthesis reference directory not found: {REFERENCE_LEVAN}")

        # Verify reference files exist
        required_files = ['lsc.fasta', 'rlsA.fasta', 'rlsB.fasta']
        missing_files = [f for f in required_files if not (REFERENCE_LEVAN / f).exists()]
        if missing_files:
            raise FileNotFoundError(f"Missing levan reference files: {', '.join(missing_files)}")

        # Create species finder directory and log its path
        species_finder_path = create_species_finder_folder(output_dir)
        logging.info(f"Created species finder directory at: {species_finder_path}")

        # Get and validate genome files
        with temporary_stdout():
            print("Collecting genome files...")
        
        logging.info("Collecting genome files...")
        genome_files = []
        
        # Just collect all FASTA files without detailed logging of each one
        for fasta_file in genomes_folder.glob('*.fasta'):
            if fasta_file.is_file():
                genome_files.append(fasta_file)
                logging.debug(f"Found root FASTA file: {fasta_file}")

        for subdir in genomes_folder.iterdir():
            if subdir.is_dir():
                logging.debug(f"Checking subdirectory: {subdir}")
                for fasta_file in subdir.glob('*.fasta'):
                    if fasta_file.is_file():
                        genome_files.append(fasta_file)
                        logging.debug(f"Found FASTA file in subdirectory: {fasta_file}")

        if not genome_files:
            logging.error(f"Directory contents of {genomes_folder}:")
            for item in genomes_folder.iterdir():
                logging.error(f"  {item}")
            raise ValueError(f"No FASTA files found in {genomes_folder} or its subdirectories")

        # Display overall progress using progress stages
        total_stages = 4
        stage = 0
        
        with temporary_stdout():
            print(f"\n✺ Found {len(genome_files)} genome files to analyze")
            print(f"✺ Starting analysis with {num_workers} worker processes...")
            overall_progress = create_progress_bar(total_stages, "Overall Progress", color="cyan")

        # Step 1: Process genomes first to get clade classifications
        stage += 1
        with temporary_stdout():
            overall_progress.update(1)
            print(f"Stage {stage}/{total_stages}: Processing genomes for clade classification...")
            
        logging.info("Starting genome processing with classifier...")
        logging.info(f"Using batch_size={batch_size}, num_workers={num_workers}")
        
        # Create a progress bar for genome processing
        genome_progress = create_progress_bar(len(genome_files), "Genome Processing", color="green")
        
        # Custom wrapper to track progress while processing genomes
        def track_genome_progress(genome_files, clade_classifier, max_workers, batch_size):
            results = []
            
            # If batch size is set, process in batches
            if batch_size > 0 and batch_size < len(genome_files):
                for i in range(0, len(genome_files), batch_size):
                    batch = genome_files[i:i+batch_size]
                    batch_results = process_genomes(batch, clade_classifier, max_workers=max_workers, batch_size=0)
                    results.extend(batch_results)
                    genome_progress.update(len(batch))
            else:
                # Otherwise process all at once with custom tracking
                with multiprocessing.Pool(processes=max_workers) as pool:
                    for result in pool.imap_unordered(
                        process_genome_with_classifier, 
                        [(g, clade_classifier) for g in genome_files]
                    ):
                        if result:
                            results.append(result)
                        genome_progress.update(1)
            
            return results
            
        all_results = track_genome_progress(
            genome_files, 
            clade_classifier, 
            max_workers=num_workers, 
            batch_size=batch_size
        )
        
        genome_progress.close()
        logging.info(f"Completed genome processing. Got {len(all_results)} results")

        # Step 2: Run species metrics analysis
        stage += 1
        with temporary_stdout():
            overall_progress.update(1)
            print(f"Stage {stage}/{total_stages}: Running species metrics analysis...")
            
        logging.info("========= Starting Species Metrics Analysis =========")

        # Check reference genome directory
        logging.info(f"Reference genomes directory path: {REFERENCE_GENOMES}")
        if not REFERENCE_GENOMES.exists():
            logging.error(f"Reference genomes directory does not exist: {REFERENCE_GENOMES}")
            raise FileNotFoundError(f"Reference genomes directory not found: {REFERENCE_GENOMES}")
        logging.info("Reference genomes directory exists")

        if skip_species_assignment:
            logging.info("Species assignment is set to be skipped")
        else:
            logging.info("Species assignment will be performed")

        try:
            run_species_metrics_for_all(
                genomes_folder_path=genomes_folder,
                species_finder_path=species_finder_path,
                threshold_species=threshold_species,
                skip_species_assignment=skip_species_assignment
            )
            logging.info("Species metrics analysis completed successfully")
        except Exception as e:
            logging.error(f"Error in species metrics analysis: {str(e)}")
            logging.exception("Species metrics exception details:")
            raise

        # Step 3: Run type analysis with parallelization
        stage += 1
        with temporary_stdout():
            overall_progress.update(1)
            print(f"Stage {stage}/{total_stages}: Running typing analysis...")
            
        logging.info(f"Starting type analysis with batch_size={batch_size}, num_workers={num_workers}, docker_limit={docker_limit}")
        extract_annotate_results = extract_annotate_assign(
            genomes_folder, 
            batch_size=batch_size, 
            num_workers=num_workers, 
            docker_limit=docker_limit
        )
        logging.info(f"Obtained {len(extract_annotate_results)} annotation results")

        # Step 4: Run levan analysis for each genome
        with temporary_stdout():
            print("Running levan synthesis analysis...")
            levan_progress = create_progress_bar(len(genome_files), "Levan Analysis", color="magenta")
            
        logging.info("Starting levan synthesis analysis")
        for genome_file in genome_files:
            logging.info(f"Full path for genome file: {genome_file}")
            try:
                logging.info(f"Running levan analysis on {genome_file}")
                levan_result = run_levan_analysis(
                    genome_path=genome_file,
                    strain_id=genome_file.stem
                )
                if levan_result is not None and not levan_result.empty:
                    for result in all_results:
                        if result['name'] == genome_file.stem:
                            result['levan_synthesis'] = format_levan_result(levan_result.iloc[0].to_dict())
                            logging.info(f"Added levan synthesis results for {genome_file.stem}")
                            break
                else:
                    logging.error(f"Levan analysis returned empty results for {genome_file.stem}")
                    for result in all_results:
                        if result['name'] == genome_file.stem:
                            result['levan_synthesis'] = "No levan synthesis results found"
                            break
            except Exception as e:
                logging.error(f"Error in levan analysis for {genome_file.name}: {e}")
                for result in all_results:
                    if result['name'] == genome_file.stem:
                        result['levan_synthesis'] = f"Analysis failed: {str(e)}"
                        break
                        
            with temporary_stdout():
                levan_progress.update(1)
                
        with temporary_stdout():
            levan_progress.close()

        # Step 5: Process final results
        stage += 1
        with temporary_stdout():
            overall_progress.update(1)
            print(f"Stage {stage}/{total_stages}: Processing final results...")
            
        # Process final results just once, after all analyses are complete
        final_results = process_results(all_results, species_finder_path, extract_annotate_results, clade_classifier)

        if final_results:
            output_csv = species_finder_path / "all_results.csv"
            write_results_to_csv(final_results, output_csv)
            logging.info(f"Wrote {len(final_results)} results to {output_csv}")

        if keep_sequence_loci:
            keep_loci_files(output_dir)
        else:
            cleanup_analysis_folders(output_dir, keep_sequence_loci)
            
        with temporary_stdout():
            overall_progress.close()
            print(f"Analysis completed successfully. Results saved to {output_dir}")

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
        # Copy types_finder folder
        types_finder_src = temp_dir / 'types_finder'
        if types_finder_src.exists():
            shutil.copytree(types_finder_src, output_dir / 'types_finder', dirs_exist_ok=True)
            logging.info(f"Copied types_finder folder to {output_dir / 'types_finder'}")
        else:
            logging.warning(f"types_finder folder not found in {temp_dir}")
            
        # Check for genome-specific types_ompa folders that might not be in types_finder
        for genome_dir in temp_dir.glob('*'):
            if genome_dir.is_dir() and genome_dir.name != 'types_finder' and genome_dir.name != 'species_finder':
                ompa_dir = genome_dir / 'types_ompa'
                if ompa_dir.exists():
                    dest_dir = output_dir / 'types_finder' / genome_dir.name / 'types_ompa'
                    dest_dir.mkdir(parents=True, exist_ok=True)
                    # Only copy .gbk and .fna files as requested
                    files_to_copy = []
                    # Check for GenBank files
                    files_to_copy.extend(ompa_dir.glob('*.gbk'))
                    # Check for FASTA files
                    files_to_copy.extend(ompa_dir.glob('*.fasta'))
                    files_to_copy.extend(ompa_dir.glob('*.fna'))
                    
                    for file in files_to_copy:
                        if file.is_file():
                            shutil.copy2(file, dest_dir)
                            logging.info(f"Copied ompa file from alternate location: {file} to {dest_dir}")

    for file in temp_dir.glob('*.log'):
        shutil.copy2(file, output_dir)

    logging.info(f"Copied final results to {output_dir}")
    logging.info(f"Final contents of output directory: {list(output_dir.glob('*'))}")
    
    # Log contents of types_finder to help with debugging
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
    parser.add_argument('--skip_species_assignment', action='store_true', help='Skip species assignment step')
    parser.add_argument('--log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                      help='Logging level')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--batch_size', type=int, default=0, 
                      help='Number of genomes to process in parallel (0 = process all at once)')
    parser.add_argument('--num_workers', type=int, default=4,
                      help='Number of worker processes for parallelization')
    parser.add_argument('--docker_limit', type=int, default=4,
                      help='Maximum number of concurrent Docker containers')
    parser.add_argument('--quiet', action='store_true', 
                      help='Reduce console output and show only progress bars and essential messages')
    
    # Note: We already handled early patching at the very start of the script
    # so no need to do it again here
    
    args = parser.parse_args()
    
    # Record start time for analysis
    start_time = time.time()
    
    # Create a simplified context manager for printing
    @contextlib.contextmanager
    def safe_print():
        """Simple context manager to ensure output is visible"""
        try:
            yield
        finally:
            # Make sure output is flushed immediately
            sys.stdout.flush()
            sys.stderr.flush()
    
    # Determine terminal capabilities for colored output
    def supports_color():
        """Check if the terminal supports color."""
        plat = sys.platform
        supported_platform = plat != 'win32' or 'ANSICON' in os.environ
        is_a_tty = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()
        return supported_platform and is_a_tty
    
    # Only use colors if the terminal supports them
    use_colors = supports_color()
    
    # ANSI color codes
    if use_colors:
        blue = '\033[94m'
        green = '\033[92m'
        cyan = '\033[96m'
        magenta = '\033[95m'
        yellow = '\033[93m'
        white = '\033[97m'
        reset = '\033[0m'
        bold = '\033[1m'
    else:
        # No colors if not supported
        blue = green = cyan = magenta = yellow = white = reset = bold = ''
    
    # Show a nice banner (always shown even in quiet mode)
    with safe_print():
        banner = f"""
{blue}{bold}╔═══════════════════════════════════════════════════════╗{reset}
{blue}{bold}║                  {white}ErwinATyper v2.0{white}                     ║{reset}
{blue}{bold}║         Comprehensive Erwinia Genome Analysis         ║{reset}
{blue}{bold}╚═══════════════════════════════════════════════════════╝{reset}

{cyan}✪ Processing input:{reset} {args.input}
{cyan}✪ Output directory:{reset} {args.output_dir}
{cyan}✪ Quiet mode:{reset} {"Enabled" if args.quiet else "Disabled"}
"""
        print(banner)

    try:
        input_path = Path(args.input).resolve()
        output_dir = Path(args.output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        log_file = output_dir / "process.log"

        setup_logging(log_file, args.log_level, args.quiet)

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

            # Initialize the clade classifier here
            logging.info(f"Initializing CRISPR clade classifier with matrix: {matrix_path}")
            clade_classifier = CRISPRCladeClassifier(matrix_path)

            try:
                run_species_and_types_finder(
                    temp_genomes_folder,
                    temp_dir_path,
                    threshold_species=args.threshold_species,
                    keep_sequence_loci=args.keep_sequence_loci,
                    clade_classifier=clade_classifier,
                    skip_species_assignment=args.skip_species_assignment,
                    batch_size=args.batch_size,
                    num_workers=args.num_workers,
                    docker_limit=args.docker_limit,
                    quiet=args.quiet
                )
            except Exception as e:
                logging.error(f"Error in run_species_and_types_finder: {str(e)}")
                logging.exception("Exception details:")
                sys.exit(1)

            with safe_print():
                print("Finalizing results...")
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
        
        # Always show the success message regardless of quiet mode
        with safe_print():
            # Print final success message with path to results
            output_csv = output_dir / 'species_finder' / 'all_results.csv'
            if output_csv.exists():
                # Calculate elapsed time (best approximation)
                elapsed_seconds = int(time.time() - start_time)
                minutes, seconds = divmod(elapsed_seconds, 60)
                hours, minutes = divmod(minutes, 60)
                time_str = f"{hours}h {minutes}m {seconds}s" if hours > 0 else f"{minutes}m {seconds}s"
                
                # Count analyzed genomes 
                genome_count = sum(1 for _ in Path(output_dir).glob('*/*/*.fasta'))
                
                # Create a fancy success message box
                success_msg = f"""
{green}{bold}╔═══════════════════════════════════════════════════════╗{reset}
{green}{bold}║ ✬ Analysis completed successfully!                    ║{reset}
{green}{bold}╚═══════════════════════════════════════════════════════╝{reset}

{cyan}▸ Results file:{reset} {output_csv}
{cyan}▸ Output directory:{reset} {output_dir}
{cyan}▸ Time elapsed:{reset} {time_str}
{cyan}▸ Genomes analyzed:{reset} {genome_count if genome_count > 0 else "All input genomes"}

{yellow}To view results:{reset} Open the CSV file in your preferred spreadsheet application
"""
                print(success_msg)


    except Exception as e:

        logging.error(f"An unexpected error occurred: {str(e)}")

        logging.exception("Exception details:")

        sys.exit(1)


if __name__ == '__main__':
    galaxy_runner()
