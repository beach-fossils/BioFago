import fcntl
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import shutil
from dataclasses import dataclass, asdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from Bio import SeqIO

from crr_genotypes.crispr_analyzer import CRISPRAnalyzer
from resistance.str_resistance import StrResistance
from resistance.virulence_genes import VirulenceGenes
from plasmids.plasmid_finder import PlasmidFinder
from crr_genotypes.crr_genotype_finder import CRRFinder
from assigning_types.assembly_statistics import FastaStatistics

from utils.clade_assigner import CRISPRCladeClassifier, parse_spacer_counts
import os
import contextlib

# Assuming the matrix file is in the project root
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
matrix_path = os.path.join(SCRIPT_DIR, '..', 'reference_crispr', 'clades_groups_correlation.csv')

@dataclass
class GenomeAnalysisResult:
    name: str
    contig_count: int
    N50_value: int
    largest_contig: str
    largest_contig_size_bp: int
    total_size_bp: int
    ambiguous_bases: int
    GC_content_percent: float
    present_plasmids: str
    streptomycin: str
    virulence_genes: str
    crispr_spacers: str
    crispr_genotype: str
    clade: str
    clade_confidence_score: float
    clade_confidence_level: str
    #clade_spacer_match_score: float
    #clade_genotype_score: float

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def process_genomes(genome_files: List[Path], clade_classifier: CRISPRCladeClassifier, max_workers: int = 4) -> List[
    Dict]:
    processed = set()
    results = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_genome = {
            executor.submit(process_single_genome, genome, clade_classifier): genome
            for genome in genome_files
            if genome.stem not in processed
        }

        for future in as_completed(future_to_genome):
            genome = future_to_genome[future]
            try:
                result = future.result()
                if result:
                    results.append(result)
                    processed.add(genome.stem)
            except Exception as e:
                logging.error(f"Exception occurred while processing {genome}: {str(e)}")

    return results

@contextlib.contextmanager
def genome_lock(genome_file: Path):
    lock_file = genome_file.with_suffix('.lock')
    try:
        with open(lock_file, 'w') as f:
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            yield
    finally:
        if lock_file.exists():
            lock_file.unlink()


def process_single_genome(genome_file: Path, clade_classifier: CRISPRCladeClassifier) -> Optional[Dict]:
    with genome_lock(genome_file):
        try:
            logging.info(f"Processing genome file: {genome_file}")

            fasta_stats = FastaStatistics(genome_file)
            stats = fasta_stats.generate_assembly_statistics()

            plasmid_finder = PlasmidFinder()
            present_plasmids = plasmid_finder.run_blast_for_genome(genome_file)

            plasmid_info = ', '.join(sorted(present_plasmids)) if present_plasmids else 'None'
            logging.info(f"Plasmids found for {genome_file.name}: {plasmid_info}")

            str_resistance = StrResistance()
            str_resistance.run_blast_for_genome(genome_file)
            str_result = str_resistance.get_results()

            virulence_genes = VirulenceGenes()
            virulence_genes.run_blast_for_genome(genome_file)
            virulence_result = virulence_genes.get_results()

            crispr_output_dir = genome_file.parent / "CRISPR_finder"
            crispr_analyzer = CRISPRAnalyzer(genome_file.parent, crispr_output_dir)
            crispr_analyzer.run_analysis()
            crispr_summaries = crispr_analyzer.get_crispr_summary()
            crispr_info = crispr_summaries.get(genome_file.stem, "CRISPR: No results")

            crr_finder = CRRFinder(genome_file, genome_file.parent)
            crr_finder.analyze_genome()
            crr_finder.process_directory()
            crr_info = crr_finder.get_crr_summary()

            # CRISPR clade classification
            spacer_counts = parse_spacer_counts(crispr_info)
            clade, score, spacer_score, genotype_score, confidence_level, subgroup = clade_classifier.determine_clade(
                crr_info, spacer_counts)

            result = {
                'name': stats['name'],
                'contig_count': stats['contig_count'],
                'N50_value': stats['N50_value'],
                'largest_contig': stats['largest_contig'],
                'largest_contig_size_bp': stats['largest_contig_size_bp'],
                'total_size_bp': stats['total_size_bp'],
                'ambiguous_bases': stats['ambiguous_bases'],
                'GC_content_percent': stats['GC_content_percent'],
                'present_plasmids': plasmid_info,
                'streptomycin': str(str_result),
                'crispr_spacers': str(crispr_info),
                'crispr_genotype': str(crr_info),
                'clade': f"{clade} {subgroup}".strip() if clade != "Unknown" else clade,
                'clade_confidence_score': round(score, 2),
                'clade_confidence_level': confidence_level,
            }

            # Process virulence genes
            for cluster, genes in virulence_result.items():
                result[f'{cluster}_genes'] = ', '.join(sorted(genes)) if genes else 'None'

            return result

        except Exception as e:
            logging.error(f"Error processing genome {genome_file.name}: {str(e)}")
            return None

def get_locus_info(genome_file: Path, locus_type: str) -> Optional[Tuple[str, str, str]]:
    locus_folder = genome_file.parent / 'types_finder' / genome_file.stem / f'types_{locus_type}'
    logging.info(f"Checking locus folder: {locus_folder}")
    if locus_folder.exists():
        final_csv = locus_folder / 'final.csv'
        if final_csv.exists():
            try:
                df = pd.read_csv(final_csv)
                if not df.empty:
                    row = df.iloc[0]
                    locus = row['Locus']
                    type_assigned = row['Type']
                    different_genes = row.get('Different Genes', '')
                    logging.info(f"Found locus info for {locus_type}: {locus}, {type_assigned}, {different_genes}")
                    return locus, type_assigned, different_genes
                else:
                    logging.warning(f"Empty CSV file for {locus_type} locus: {final_csv}")
            except Exception as e:
                logging.error(f"Error reading CSV file for {locus_type} locus: {final_csv}. Error: {str(e)}")
        else:
            logging.warning(f"Final CSV not found for {locus_type} locus: {final_csv}")
    else:
        logging.warning(f"Locus folder not found for {locus_type}: {locus_folder}")
    return None


def write_results_to_csv(results: List[Dict], output_path: Path) -> None:
    try:
        logging.info(f"Starting write_results_to_csv with {len(results)} results")
        df = pd.DataFrame(results)

        base_columns = [
            'name', 'species', 'ANI_species', 'contig_count', 'N50_value',
            'largest_contig', 'largest_contig_size_bp', 'total_size_bp',
            'ambiguous_bases', 'GC_content_percent', 'present_plasmids',
            'streptomycin', 'crispr_spacers', 'crispr_genotype',
            'capsule_locus', 'cellulose_locus', 'lps_locus', 'sorbitol_locus',
            'flag_i_locus', 'flag_ii_locus', 'flag_iii_locus', 'flag_iv_locus',
            't3ss_i_locus', 't3ss_ii_locus', 't6ss_i_locus', 't6ss_ii_locus',
            'clade', 'clade_confidence_score', 'clade_confidence_level'
        ]

        # Define the virulence gene columns - only include 'other_genes'
        virulence_columns = ['other_genes']

        # Combine base columns and remaining virulence columns
        all_columns = base_columns + virulence_columns

        # Reindex the DataFrame with all columns, filling missing values with 'None'
        df = df.reindex(columns=all_columns, fill_value='None')

        logging.info(f"DataFrame shape before writing: {df.shape}")
        logging.info(f"DataFrame columns: {df.columns.tolist()}")
        df.to_csv(output_path, index=False)
        logging.info(f"Final results written to {output_path}")

        verify_csv(output_path)
    except Exception as e:
        logging.error(f"Error writing results to CSV: {str(e)}")
        logging.exception("Exception details:")

def verify_csv(output_path: Path) -> None:
    if output_path.exists():
        written_df = pd.read_csv(output_path)
        logging.info(f"Verified written CSV. Shape: {written_df.shape}")
        logging.info(f"Verified CSV columns: {written_df.columns.tolist()}")
        logging.info(f"Verified CSV preview:\n{written_df.head().to_string()}")
    else:
        logging.error(f"Failed to create CSV file at {output_path}")

def cleanup_analysis_folders(output_dir: Path) -> None:
    folders_to_remove = ['CRISPR_finder', 'CRR_finder', 'plasmid_finder', 'resistance_finder', 'types_finder']
    logging.info(f"Starting cleanup of analysis folders in {output_dir}")

    for folder in folders_to_remove:
        folder_path = output_dir / folder
        remove_folder(folder_path)

    cleanup_species_finder(output_dir)

    logging.info(f"Cleanup completed. Remaining contents of {output_dir}: {list(output_dir.glob('*'))}")

def remove_folder(folder_path: Path) -> None:
    if folder_path.exists():
        try:
            shutil.rmtree(folder_path)
            logging.info(f"Removed folder: {folder_path}")
        except Exception as e:
            logging.error(f"Error removing folder {folder_path}: {str(e)}")
    else:
        logging.info(f"Folder not found (already removed or never created): {folder_path}")

def cleanup_species_finder(output_dir: Path) -> None:
    species_finder_path = output_dir / 'species_finder'
    if species_finder_path.exists():
        results_folder = species_finder_path / 'results'
        if results_folder.exists():
            if not any(results_folder.iterdir()):
                try:
                    results_folder.rmdir()
                    logging.info(f"Removed empty folder: {results_folder}")
                except Exception as e:
                    logging.error(f"Error removing empty folder {results_folder}: {str(e)}")
            else:
                logging.info(f"Folder not empty, skipping removal: {results_folder}")
        else:
            logging.info(f"Results folder not found: {results_folder}")
    else:
        logging.info(f"Species finder folder not found: {species_finder_path}")


def keep_loci_files(output_dir: Path) -> None:
    types_finder_folder = output_dir / 'genomes' / 'types_finder'
    if types_finder_folder.exists():
        new_types_finder_folder = output_dir / 'types_finder'
        new_types_finder_folder.mkdir(parents=True, exist_ok=True)

        # List of all possible type folders
        type_folders = [
            'types_capsule', 'types_cellulose', 'types_lps', 'types_srl',
            'types_flag_I', 'types_flag_II', 'types_flag_III', 'types_flag_IV',
            'types_T3SS_I', 'types_T3SS_II', 'types_T6SS_I', 'types_T6SS_II'
        ]

        def process_types_finder_folders(types_finder_folder: Path, new_types_finder_folder: Path) -> None:
            for genome_folder in types_finder_folder.iterdir():
                if genome_folder.is_dir():
                    for type_folder in type_folders:
                        locus_folder = genome_folder / type_folder
                        if locus_folder.exists():
                            process_locus_folder(locus_folder, new_types_finder_folder, genome_folder.name, type_folder)

        def process_locus_folder(locus_folder: Path, new_types_finder_folder: Path, genome_name: str,
                                 type_folder: str) -> None:
            prokka_folder = locus_folder / 'prokka'
            if prokka_folder.exists():
                for subfolder in prokka_folder.iterdir():
                    if subfolder.is_dir():
                        # Get both .gbk and .fna files
                        keep_files = list(subfolder.glob('*.gbk')) + list(subfolder.glob('*.fna'))
                        for file in keep_files:
                            try:
                                new_location = new_types_finder_folder / genome_name / type_folder
                                new_location.mkdir(parents=True, exist_ok=True)
                                shutil.copy2(str(file), str(new_location / file.name))
                                logging.info(f"Copied {file} to {new_location / file.name}")
                            except Exception as e:
                                logging.error(f"Error copying file {file}: {str(e)}")

        process_types_finder_folders(types_finder_folder, new_types_finder_folder)

    remove_unnecessary_folders(output_dir)

    logging.info(f"Loci files kept. Remaining contents of {output_dir}: {list(output_dir.glob('*'))}")

def process_types_finder_folders(types_finder_folder: Path, new_types_finder_folder: Path) -> None:
    for genome_folder in types_finder_folder.iterdir():
        if genome_folder.is_dir():
            for locus_type in ['types_capsule', 'types_cellulose', 'types_lps', 'types_srl']:
                locus_folder = genome_folder / locus_type
                if locus_folder.exists():
                    copy_loci_files(locus_folder, new_types_finder_folder, genome_folder.name, locus_type)

def copy_loci_files(locus_folder: Path, new_types_finder_folder: Path, genome_name: str, locus_type: str) -> None:
    prokka_folder = locus_folder / 'prokka'
    if prokka_folder.exists():
        for subfolder in prokka_folder.iterdir():
            if subfolder.is_dir():
                keep_files = list(subfolder.glob('*.gbk')) + list(subfolder.glob('*.fna'))
                for file in keep_files:
                    try:
                        new_location = new_types_finder_folder / genome_name / locus_type
                        new_location.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(str(file), str(new_location / file.name))
                        logging.info(f"Copied {file} to {new_location / file.name}")
                    except Exception as e:
                        logging.error(f"Error copying file {file}: {str(e)}")

def remove_unnecessary_folders(output_dir: Path) -> None:
    folders_to_remove = ['CRISPR_finder', 'CRR_finder', 'plasmid_finder', 'resistance_finder']
    for folder in folders_to_remove:
        folder_path = output_dir / folder
        remove_folder(folder_path)

# Setup logging when the module is imported
setup_logging()
