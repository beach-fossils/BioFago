import logging
import os
import re
import shutil
from pathlib import Path
from typing import List, Dict

import pandas as pd
import time

from src.metrics_species_caller import new_run_species_metrics_finder

logger = logging.getLogger(__name__)

def create_species_finder_folder(genomes_folder_path: Path) -> Path:
    species_finder_path = genomes_folder_path / 'species_finder'
    species_finder_path.mkdir(exist_ok=True)
    return species_finder_path

def run_species_metrics_for_all(genomes_folder_path: Path, species_finder_path: Path,
                                threshold_species: float, skip_species_assignment: bool = False):
    logging.info(f"Starting species metrics analysis in {genomes_folder_path}")
    logging.info(f"Skip species assignment: {skip_species_assignment}")
    logging.info(f"Species finder path: {species_finder_path}")
    logging.info(f"Threshold species: {threshold_species}")

    # Find all FASTA files in both root and subdirectories
    genome_files = []

    # Look for files directly in the genomes folder
    logger.info(f"Searching for FASTA files in: {genomes_folder_path}")
    for fasta_file in genomes_folder_path.glob('*.fasta'):
        if fasta_file.is_file():
            genome_files.append(fasta_file)
            logger.info(f"Found FASTA file in root: {fasta_file}")

    # Look in subdirectories
    for subdir in genomes_folder_path.iterdir():
        if subdir.is_dir():
            logger.info(f"Searching subdirectory: {subdir}")
            for fasta_file in subdir.glob('*.fasta'):
                if fasta_file.is_file():
                    genome_files.append(fasta_file)
                    logger.info(f"Found FASTA file in subdirectory: {fasta_file}")

    if not genome_files:
        logger.error(f"No FASTA files found in {genomes_folder_path} or its subdirectories")
        logger.error(f"Directory contents: {list(genomes_folder_path.glob('*'))}")
        return

    logger.info(f"Found {len(genome_files)} genome files to process")
    logger.info(f"Genome files: {[str(f) for f in genome_files]}")

    # Process each genome
    for genome_file in genome_files:
        logger.info(f"Processing genome file for species analysis: {genome_file}")
        try:
            # Create output directory if it doesn't exist
            species_finder_path.mkdir(parents=True, exist_ok=True)
            logger.info(f"Created species finder directory: {species_finder_path}")

            # Get complete filename without extension
            file_stem = genome_file.name
            # Remove only the final extension
            file_stem = os.path.splitext(file_stem)[0]
            logger.info(f"Extracted full filename without extension: {file_stem}")
            
            # For debugging GCF/GCA special cases
            if 'GCF_' in file_stem or 'GCA_' in file_stem:
                parts = file_stem.split('_')
                logger.info(f"GCF/GCA file detected. Parts: {len(parts)}, First parts: {parts[:3]}")
                
            # Run metrics finder for this genome
            logger.info(f"Running metrics finder for {genome_file.name} (output name: {file_stem})")
            new_run_species_metrics_finder(
                single_sequence_path=genome_file,
                species_finder_path=species_finder_path,
                threshold_species=threshold_species,
                skip_species_assignment=skip_species_assignment,
                output_name=file_stem  # Pass clean name without extension
            )

            # Verify results
            output_file = species_finder_path / f"{file_stem}.csv"
            if output_file.exists():
                try:
                    df = pd.read_csv(output_file)
                    species = df['Species'].iloc[0] if 'Species' in df.columns else 'Unknown'
                    ani = df['ANI'].iloc[0] if 'ANI' in df.columns else 0.0
                    if skip_species_assignment:
                        logger.info(f"Species assignment skipped for {genome_file.name} as requested")
                    else:
                        logger.info(f"Species assignment for {genome_file.name}: {species} (ANI: {ani})")
                except Exception as e:
                    logger.error(f"Error reading species results for {genome_file.name}: {e}")
            else:
                logger.error(f"No results file generated for {genome_file.name} at {output_file}")

        except Exception as e:
            logger.error(f"Error in species metrics finder for {genome_file}: {e}")
            logger.exception("Exception details:")

    logger.info("=========== Completed species metrics analysis for all genomes ===========")

def create_individual_folders(genomes_folder: Path) -> Path:
    """Move each genome to its individual folder, preserving the original file name."""
    genomes_folder_path = Path(genomes_folder).resolve()
    for genome_file in genomes_folder_path.glob('*.fasta'):
        genome_dir = genomes_folder_path / genome_file.stem
        os.makedirs(genome_dir, exist_ok=True)
        shutil.move(str(genome_file), str(genome_dir / genome_file.name))
        logging.info(f"Moved {genome_file} to {genome_dir}")
    return genomes_folder_path

def update_species_csv(species_finder_path: Path, genome_name: str, locus_name: str, locus_type: str,
                       assign_confidence: str):
    species_csv_path = species_finder_path / f"{genome_name}.csv"
    if species_csv_path.exists():
        df = pd.read_csv(species_csv_path)
        column_name = f"{locus_name}"
        df[column_name] = f"{locus_type} ({assign_confidence})"
        df.to_csv(species_csv_path, index=False)
        logging.info(f"Updated {species_csv_path} with {locus_name}: {locus_type} ({assign_confidence})")
    else:
        logging.warning(f"Species CSV file not found at {species_csv_path}")

def update_species_csv_with_results(species_finder_path: Path, genome_name: str, plasmid_info: str, sorbitol_info: str, crispr_info: str, crr_info: str):
    csv_file = species_finder_path / f"{genome_name}.csv"
    if csv_file.exists():
        df = pd.read_csv(csv_file)
        df['Plasmid'] = plasmid_info
        df['Sorbitol'] = sorbitol_info
        df['CRISPR'] = crispr_info
        df['CRR'] = crr_info
        df.to_csv(csv_file, index=False)
        logging.info(f"Updated {csv_file} with plasmid, sorbitol, CRISPR, and CRR information")
    else:
        logging.warning(f"Species CSV file not found at {csv_file}")

def cleanup_unwanted_species_folders(base_output_path: Path, main_species_finder_path: Path):
    """
    Remove unwanted 'species_finder' folders from output directories.
    """
    for output_dir in ['CRISPR_finder', 'CRR_finder', 'plasmid_finder', 'resistance_finder']:
        dir_path = base_output_path / output_dir
        if dir_path.exists():
            species_finder = dir_path / 'species_finder'
            if species_finder.exists() and species_finder != main_species_finder_path:
                try:
                    shutil.rmtree(species_finder)
                    logging.info(f"Removed unwanted species_finder folder from {output_dir}")
                except Exception as e:
                    logging.error(f"Error removing species_finder folder from {output_dir}: {e}")

def determine_flag_type(row: pd.Series) -> str:
    """Determine flag type based on combinations of flag loci."""
    def extract_type(value):
        if pd.isna(value) or value == '(Unknown)':
            return 'NA'
        return value.split()[0]

    flag_i = extract_type(row['flag_i_locus'])
    flag_ii = extract_type(row['flag_ii_locus'])
    flag_iii = extract_type(row['flag_iii_locus'])
    flag_iv = extract_type(row['flag_iv_locus'])

    # Combination for FL01
    if all(x.endswith('01') for x in [flag_i, flag_ii, flag_iii, flag_iv]):
        return 'FL01'
    # Combination for FL02
    elif (flag_i == 'FLI01' and flag_ii == 'FLII01' and
          flag_iii == 'FLIII01' and flag_iv == 'FLIV02'):
        return 'FL02'
    # Any other combination
    return 'NA'

def determine_t6ss_type(row: pd.Series) -> str:
    """Determine T6SS type based on combinations of T6SS loci."""
    def extract_type(value):
        if pd.isna(value) or value == '(Unknown)':
            return 'NA'
        return value.split()[0]

    t6ss_i = extract_type(row['t6ss_i_locus'])
    t6ss_ii = extract_type(row['t6ss_ii_locus'])

    combinations = {
        ('TSI01', 'TSII01'): 'TS01',
        ('TSI01', 'TSII05'): 'TS02',
        ('TSI03', 'TSII01'): 'TS03',
        ('TSI01', 'TSII02'): 'TS04',
        ('TSI04', 'TSII01'): 'TS05'
    }

    return combinations.get((t6ss_i, t6ss_ii), 'NA')

def determine_t3ss_type(row: pd.Series) -> str:
    """Determine T3SS type based on combinations of T3SS loci."""
    def extract_type(value):
        if pd.isna(value) or value == '(Unknown)':
            return 'NA'
        return value.split()[0]

    t3ss_i = extract_type(row['t3ss_i_locus'])
    t3ss_ii = extract_type(row['t3ss_ii_locus'])

    combinations = {
        ('TTI01', 'TTII01'): 'TT01',
        ('TTI01', 'TTII02'): 'TT02'
    }

    return combinations.get((t3ss_i, t3ss_ii), 'NA')

def final_enhance_csv(csv_path: str) -> None:
    """
    Enhance the CSV file with column transformations and reordering:
    1. Add new crispr_genotype column
    2. Reformat crispr_genotype_extra_info column
    3. Rename cr4_spacers to cr3_spacers
    4. Reformat streptomycin column
    5. Update crispr_spacers_count format
    6. Add type determination columns (flag1_type, ts_type, tt_type)
    7. Reorder specific columns
    8. Rename certain locus columns to *_type
    9. Drop the 'other_genes' column if it exists
    """
    logging.info(f"Starting CSV enhancement for {csv_path}")

    # Read the CSV file
    df = pd.read_csv(csv_path)

    # Drop the 'other_genes' column if it exists
    if 'other_genes' in df.columns:
        df.drop(columns=['other_genes'], inplace=True)

    def extract_genotype(info: str) -> str:
        """Extract CR1/CR2/CR4 genotype from extra info string."""
        if pd.isna(info) or not isinstance(info, str) or info == 'None':
            return "Unknown/Unknown/Unknown"

        info = re.sub(r'^CRR:\s*', '', info)
        cr1_group = cr2_group = cr4_group = "Unknown"
        pattern = r'CRR(\d+)\s*-\s*Group:\s*([^,>]+)(?=,|\s*>|\s*;|\s*$)'

        for match in re.finditer(pattern, info):
            crr_num = match.group(1)
            group = match.group(2).strip()
            group = re.sub(r'^Group\s+', '', group)
            group = re.sub(r'\s*\([^)]*\)', '', group)

            if crr_num == '1':
                cr1_group = group
            elif crr_num == '2':
                cr2_group = group
            elif crr_num == '4':
                cr4_group = group

        return f"{cr1_group}/{cr2_group}/{cr4_group}"

    def reformat_extra_info(info: str) -> str:
        """Reformat the extra info string to the new structure."""
        if pd.isna(info) or not isinstance(info, str) or info == 'None':
            return "No CRISPR information available"

        info = re.sub(r'^CRR:\s*', '', info)
        cr1_info = cr2_info = cr4_info = "No information"
        crr_pattern = r'(CRR(\d+)[^;]+)(?:\s*;\s*|$)'

        for match in re.finditer(crr_pattern, info):
            full_section = match.group(1)
            crr_num = match.group(2)

            if crr_num == '4':
                section = 'CR3' + full_section[4:]
            else:
                section = 'CR' + full_section[3:]

            if crr_num == '1':
                cr1_info = section
            elif crr_num == '2':
                cr2_info = section
            elif crr_num == '4':
                cr4_info = section

        return f"{cr1_info} | {cr2_info} | {cr4_info}"

    def format_streptomycin(info: str) -> str:
        """Format streptomycin information in a more readable way with sorted genes."""
        if pd.isna(info) or not isinstance(info, str) or info == 'None':
            return "No streptomycin information available"

        try:
            data = eval(info)  # Evaluate the string representation of the dictionary
            present_genes = data.get('present_genes', [])
            missing_genes = sorted(data.get('missing_genes', []))  # Sort missing genes alphabetically
            rpsL_aa = data.get('rpsL_amino_acid_43', 'Unknown')

            # Sort present genes by gene name
            present_genes.sort(key=lambda x: x[0])

            # Format present genes
            present_str = "Present genes: "
            if present_genes:
                present_str += ", ".join([
                    f"{gene} ({round(identity, 2)}, {coverage})"
                    for gene, identity, coverage in present_genes
                ])
            else:
                present_str += "None"

            # Format missing genes (already sorted)
            missing_str = "Missing genes: "
            if missing_genes:
                missing_str += ", ".join(missing_genes)
            else:
                missing_str += "None"

            return f"{present_str} | {missing_str} | Amino acid at 43_rpsL: {rpsL_aa}"
        except:
            return "Error parsing streptomycin information"

    def update_spacers_count(info: str) -> str:
        """Update the format of spacers count information."""
        if pd.isna(info) or not isinstance(info, str) or info == 'None':
            return "No spacers information available"

        # Replace CRR with CR and CRR4 with CR3
        info = info.replace('CRR', 'CR')
        info = info.replace('CR4', 'CR3')

        return info

    # Process columns
    if 'crispr_genotype_extra_info' in df.columns:
        # Add crispr_genotype column
        df.insert(
            df.columns.get_loc('crispr_genotype_extra_info'),
            'crispr_genotype',
            df['crispr_genotype_extra_info'].apply(extract_genotype)
        )
        df['crispr_genotype_extra_info'] = df['crispr_genotype_extra_info'].apply(reformat_extra_info)

    # Rename cr4_spacers to cr3_spacers
    if 'cr4_spacers' in df.columns:
        df = df.rename(columns={'cr4_spacers': 'cr3_spacers'})

    # Format streptomycin column
    if 'streptomycin' in df.columns:
        df['streptomycin'] = df['streptomycin'].apply(format_streptomycin)

    # Update crispr_spacers_count format
    if 'crispr_spacers_count' in df.columns:
        df['crispr_spacers_count'] = df['crispr_spacers_count'].apply(update_spacers_count)

    # Rename locus columns to type
    rename_map = {
        'capsule_locus': 'capsule_type',
        'cellulose_locus': 'cellulose_type',
        'lps_locus': 'lps_type',
        'sorbitol_locus': 'sorbitol_type',
        'flag3_locus': 'flag3_type',
        'ompa_locus': 'ompa_type'  # Add ompA locus to type renaming
    }
    df = df.rename(columns=rename_map)

    # Add new type columns
    df['flag1_type'] = df.apply(determine_flag_type, axis=1)
    df['ts_type'] = df.apply(determine_t6ss_type, axis=1)
    df['tt_type'] = df.apply(determine_t3ss_type, axis=1)

    # Get all column names for reordering
    cols = list(df.columns)

    # Reorder clade columns
    if all(col in df.columns for col in
           ['crispr_genotype_extra_info', 'clade', 'clade_confidence_score', 'clade_confidence_level']):
        for col in ['clade', 'clade_confidence_score', 'clade_confidence_level']:
            cols.remove(col)
        insert_idx = cols.index('crispr_genotype_extra_info') + 1
        cols[insert_idx:insert_idx] = ['clade', 'clade_confidence_score', 'clade_confidence_level']

    # Reorder type columns
    if 'flag_i_locus' in cols:
        insert_idx = cols.index('flag_i_locus')
        cols.remove('flag1_type')
        cols.insert(insert_idx, 'flag1_type')

    if 't6ss_i_locus' in cols:
        insert_idx = cols.index('t6ss_i_locus')
        cols.remove('ts_type')
        cols.insert(insert_idx, 'ts_type')

    if 't3ss_i_locus' in cols:
        insert_idx = cols.index('t3ss_i_locus')
        cols.remove('tt_type')
        cols.insert(insert_idx, 'tt_type')

    # Move flag3_type after flag_iv_locus
    if 'flag3_type' in cols and 'flag_iv_locus' in cols:
        cols.remove('flag3_type')
        insert_idx = cols.index('flag_iv_locus') + 1
        cols.insert(insert_idx, 'flag3_type')
        
    # Move ompa_type after t6ss_ii_locus
    if 'ompa_type' in cols and 't6ss_ii_locus' in cols:
        cols.remove('ompa_type')
        insert_idx = cols.index('t6ss_ii_locus') + 1
        cols.insert(insert_idx, 'ompa_type')

    # Apply final column ordering
    df = df[cols]

    # Save enhanced CSV
    df.to_csv(csv_path, index=False)
    logging.info("Saved enhanced CSV with all modifications including type determinations and dropped 'other_genes' column")




