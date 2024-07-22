import logging
from pathlib import Path
from typing import Dict, List
import shutil

import pandas as pd
from Bio import SeqIO

from src.crr_genotypes.crispr_analyzer import CRISPRAnalyzer
from src.resistance.str_resistance import StrResistance
from src.plasmids.plasmid_finder import PlasmidFinder
from src.crr_genotypes.crr_genotype_finder import CRRFinder
from src.assigning_types.assembly_statistics import FastaStatistics


def process_genome(genome_dir: Path) -> Dict:
    try:
        genome_file = genome_dir / f"{genome_dir.name}.fasta"
        logging.info(f"Processing genome: {genome_file}")

        fasta_stats = FastaStatistics(genome_file)
        stats = fasta_stats.generate_assembly_statistics()

        plasmid_finder = PlasmidFinder()
        plasmid_finder.run_blast_for_genome(genome_file)
        plasmid_result = plasmid_finder.get_present_plasmids()

        str_resistance = StrResistance()
        str_resistance.run_blast_for_genome(genome_file)
        str_result = str_resistance.get_results()

        crispr_output_dir = Path(genome_dir).parent / "CRISPR_finder"
        crispr_analyzer = CRISPRAnalyzer(Path(genome_dir), crispr_output_dir)
        crispr_analyzer.run_analysis()
        crispr_info = crispr_analyzer.get_crispr_summary()

        crr_finder = CRRFinder(genome_file, genome_dir.parent)
        crr_finder.analyze_genome()
        crr_finder.process_directory()
        crr_info = crr_finder.get_crr_summary()

        return {
            "name": stats['name'],
            "contig_count": stats['contig_count'],
            "N50_value": stats['N50_value'],
            "largest_contig": stats['largest_contig'],
            "largest_contig_size_bp": stats['largest_contig_size_bp'],
            "total_size_bp": stats['total_size_bp'],
            "ambiguous_bases": stats['ambiguous_bases'],
            "GC_content_percent": stats['GC_content_percent'],
            "present_plasmids": ', '.join(plasmid_result) if plasmid_result else 'None',
            "streptomycin": str(str_result),
            "crispr_spacers": str(crispr_info),
            "crispr_genotype": str(crr_info),
            "genome": genome_dir.name,
            "plasmid": plasmid_result,
            "str": str_result,
            "crispr": crispr_info,
            "crr": crr_info
        }
    except Exception as e:
        logging.error(f"Error processing genome {genome_dir.name}: {str(e)}")
        return None


def write_results_to_csv(results: List[Dict], output_path: Path) -> None:
    try:
        df = pd.DataFrame(results)
        columns = [
            'name', 'species', 'ANI_species', 'contig_count', 'N50_value',
            'largest_contig', 'largest_contig_size_bp', 'total_size_bp',
            'ambiguous_bases', 'GC_content_percent', 'present_plasmids',
            'streptomycin', 'crispr_spacers', 'crispr_genotype',
            'capsule_locus', 'cellulose_locus', 'lps_locus', 'sorbitol_locus'
        ]
        df = df.reindex(columns=columns)
        df.to_csv(output_path, index=False)
        logging.info(f"Final results written to {output_path}")
        logging.debug(f"CSV contents:\n{df.to_string()}")
    except Exception as e:
        logging.error(f"Error writing results to CSV: {str(e)}")


def cleanup_analysis_folders(genomes_folder_path: Path) -> None:
    folders_to_remove = ['CRISPR_finder', 'CRR_finder', 'plasmid_finder', 'resistance_finder', 'types_finder']
    for folder in folders_to_remove:
        folder_path = genomes_folder_path / folder
        if folder_path.exists():
            try:
                shutil.rmtree(folder_path)
                logging.info(f"Removed folder: {folder_path}")
            except Exception as e:
                logging.error(f"Error removing folder {folder_path}: {str(e)}")


def keep_loci_files(genomes_folder_path: Path) -> None:
    types_finder_folder = genomes_folder_path / 'types_finder'
    if types_finder_folder.exists():
        for genome_folder in types_finder_folder.iterdir():
            if genome_folder.is_dir():
                for locus_type in ['types_capsule', 'types_cellulose', 'types_lps', 'types_srl']:
                    locus_folder = genome_folder / locus_type
                    if locus_folder.exists():
                        prokka_folder = locus_folder / 'prokka'
                        if prokka_folder.exists():
                            for subfolder in prokka_folder.iterdir():
                                if subfolder.is_dir():
                                    keep_files = list(subfolder.glob('*.gbk')) + list(subfolder.glob('*.fna'))
                                    for file in keep_files:
                                        try:
                                            new_location = genomes_folder_path / genome_folder.name / locus_type
                                            new_location.mkdir(parents=True, exist_ok=True)
                                            shutil.copy(str(file), str(new_location / file.name))
                                            logging.info(f"Copied {file} to {new_location / file.name}")
                                        except Exception as e:
                                            logging.error(f"Error copying file {file}: {str(e)}")

    # Remove analysis folders after copying files
    cleanup_analysis_folders(genomes_folder_path)
