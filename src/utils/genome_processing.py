import logging
from pathlib import Path
from typing import Dict, List, Optional
import shutil
from dataclasses import dataclass, asdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from Bio import SeqIO

from crr_genotypes.crispr_analyzer import CRISPRAnalyzer
from resistance.str_resistance import StrResistance
from plasmids.plasmid_finder import PlasmidFinder
from crr_genotypes.crr_genotype_finder import CRRFinder
from assigning_types.assembly_statistics import FastaStatistics

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
    crispr_spacers: str
    crispr_genotype: str

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def process_single_genome(genome_file: Path) -> Optional[Dict]:
    try:
        logging.info(f"Processing genome file: {genome_file}")

        fasta_stats = FastaStatistics(genome_file)
        stats = fasta_stats.generate_assembly_statistics()

        plasmid_finder = PlasmidFinder()
        present_plasmids = plasmid_finder.run_blast_for_genome(genome_file)
        plasmid_info = ', '.join(present_plasmids) if present_plasmids else 'None'
        logging.info(f"Plasmids found for {genome_file.name}: {plasmid_info}")

        str_resistance = StrResistance()
        str_resistance.run_blast_for_genome(genome_file)
        str_result = str_resistance.get_results()

        crispr_output_dir = genome_file.parent / "CRISPR_finder"
        crispr_analyzer = CRISPRAnalyzer(genome_file.parent, crispr_output_dir)
        crispr_analyzer.run_analysis()
        crispr_summaries = crispr_analyzer.get_crispr_summary()
        crispr_info = crispr_summaries.get(genome_file.stem, "CRISPR: No results")

        crr_finder = CRRFinder(genome_file, genome_file.parent)
        crr_finder.analyze_genome()
        crr_finder.process_directory()
        crr_info = crr_finder.get_crr_summary()

        result = GenomeAnalysisResult(
            name=stats['name'],
            contig_count=stats['contig_count'],
            N50_value=stats['N50_value'],
            largest_contig=stats['largest_contig'],
            largest_contig_size_bp=stats['largest_contig_size_bp'],
            total_size_bp=stats['total_size_bp'],
            ambiguous_bases=stats['ambiguous_bases'],
            GC_content_percent=stats['GC_content_percent'],
            present_plasmids=plasmid_info,
            streptomycin=str(str_result),
            crispr_spacers=str(crispr_info),
            crispr_genotype=str(crr_info),
        )
        return asdict(result)
    except Exception as e:
        logging.error(f"Error processing genome {genome_file.name}: {str(e)}")
        return None

def process_genomes(genome_files: List[Path], max_workers: int = 4) -> List[Dict]:
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_genome = {executor.submit(process_single_genome, genome): genome for genome in genome_files}
        for future in as_completed(future_to_genome):
            genome = future_to_genome[future]
            try:
                result = future.result()
                if result:
                    results.append(result)
            except Exception as e:
                logging.error(f"Exception occurred while processing {genome}: {str(e)}")
    return results

def write_results_to_csv(results: List[Dict], output_path: Path) -> None:
    try:
        logging.info(f"Starting write_results_to_csv with {len(results)} results")
        df = pd.DataFrame(results)
        columns = [
            'name', 'species', 'ANI_species', 'contig_count', 'N50_value',
            'largest_contig', 'largest_contig_size_bp', 'total_size_bp',
            'ambiguous_bases', 'GC_content_percent', 'present_plasmids',
            'streptomycin', 'crispr_spacers', 'crispr_genotype',
            'capsule_locus', 'cellulose_locus', 'lps_locus', 'sorbitol_locus'
        ]
        df = df.reindex(columns=columns, fill_value='Unknown')
        logging.info(f"DataFrame shape before writing: {df.shape}")
        logging.info(f"DataFrame columns: {df.columns.tolist()}")
        logging.info(f"DataFrame preview:\n{df.head().to_string()}")
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
