import json
import subprocess
from pathlib import Path
import logging
import pandas as pd
from typing import Dict, List, Set
import tempfile
from Bio import SeqIO
import os
import sys

# Import quiet mode
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from quiet_mode import QUIET_MODE

# Set up logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Dynamically get the reference_resistance_genes folder
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
REFERENCE_GENES = PROJECT_ROOT / "reference_resistance_genes" / "virulence_genes"
GENES_MAP = REFERENCE_GENES / "genes_and_locus_tags_by_accession_and_cluster.json"

class VirulenceGenes:
    def __init__(self, reference_folder: Path = REFERENCE_GENES) -> None:
        self.reference_folder = reference_folder
        self.gene_lengths = {}
        self.results = {}
        self.present_genes = {}
        logger.info(f"Initialized VirulenceGenes with reference folder: {self.reference_folder}")
        self.gene_mapping = self.load_gene_mapping()

    def combine_gene_sequences(self, temp_dir: Path) -> Path:
        combined_gene_fasta = temp_dir / "combined_genes.fasta"
        logger.info(f"Combining gene sequences into {combined_gene_fasta}")
        gene_count = 0
        with open(combined_gene_fasta, 'w') as outfile:
            for strain_folder in self.reference_folder.glob('*'):
                if strain_folder.is_dir():
                    strain_name = strain_folder.name
                    logger.debug(f"Processing strain folder: {strain_name}")
                    for gene_group in strain_folder.glob('*'):
                        if gene_group.is_dir():
                            gene_group_name = gene_group.name
                            logger.debug(f"Processing gene group: {gene_group_name}")
                            for gene_file in gene_group.glob('*.fasta'):
                                logger.debug(f"Processing gene file: {gene_file}")
                                with open(gene_file, 'r') as infile:
                                    seq_records = list(SeqIO.parse(infile, "fasta"))
                                    for record in seq_records:
                                        record.id = f"{strain_name}|{gene_group_name}|{gene_file.stem}"
                                        record.description = ""
                                        self.gene_lengths[record.id] = len(record.seq)
                                        SeqIO.write(record, outfile, "fasta")
                                        gene_count += 1
        logger.info(f"Combined {gene_count} gene sequences into {combined_gene_fasta}")
        return combined_gene_fasta

    def create_blast_db(self, combined_gene_fasta: Path, temp_dir: Path) -> Path:
        blast_db_name = temp_dir / "combined_genes"
        make_db_cmd = [
            'makeblastdb',
            '-in', str(combined_gene_fasta),
            '-dbtype', 'nucl',
            '-out', str(blast_db_name)
        ]
        logger.info(f"Creating BLAST database with command: {' '.join(make_db_cmd)}")
        try:
            if QUIET_MODE:
                with open(os.devnull, 'w') as devnull:
                    subprocess.run(make_db_cmd, check=True, stdout=devnull, stderr=devnull)
            else:
                subprocess.run(make_db_cmd, check=True, capture_output=True, text=True)
            logger.info(f"BLAST database created at {blast_db_name}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error creating BLAST database: {e}")
            logger.error(f"STDOUT: {e.stdout if hasattr(e, 'stdout') else 'N/A'}")
            logger.error(f"STDERR: {e.stderr if hasattr(e, 'stderr') else 'N/A'}")
            raise
        return blast_db_name

    def run_blast_for_genome(self, genome_file: Path) -> None:
        logger.info(f"Running BLAST for genome: {genome_file}")
        self.genome_file = genome_file
        self.output_folder = genome_file.parent.parent / "virulence_finder" / self.genome_file.stem
        self.output_folder.mkdir(parents=True, exist_ok=True)
        self.results_folder = self.output_folder / "results"
        self.results_folder.mkdir(exist_ok=True)

        self.present_genes = {}  # Reset present_genes for each genome

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = Path(temp_dir)
            combined_gene_fasta = self.combine_gene_sequences(temp_dir_path)
            blast_db_name = self.create_blast_db(combined_gene_fasta, temp_dir_path)

            result_file = self.results_folder / (self.genome_file.stem + "_blast_results.csv")
            blast_cmd = [
                'blastn',
                '-query', str(self.genome_file),
                '-db', str(blast_db_name),
                '-outfmt',
                "10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                '-dust', 'no',
                '-perc_identity', '90.0',
                '-out', str(result_file)
            ]

            logger.info(f"Running BLAST command: {' '.join(blast_cmd)}")
            try:
                if QUIET_MODE:
                    with open(os.devnull, 'w') as devnull:
                        subprocess.run(blast_cmd, check=True, stdout=devnull, stderr=devnull)
                else:
                    subprocess.run(blast_cmd, check=True, capture_output=True, text=True)
                logger.info(f"BLAST search completed for {self.genome_file} with results saved to {result_file}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error running BLAST: {e}")
                logger.error(f"STDOUT: {e.stdout if hasattr(e, 'stdout') else 'N/A'}")
                logger.error(f"STDERR: {e.stderr if hasattr(e, 'stderr') else 'N/A'}")
                raise

            self._parse_results(result_file)

        # Store the results for this genome
        self.results[self.genome_file.stem] = self.present_genes

    def _parse_results(self, result_file: Path) -> None:
        logger.info(f"Parsing BLAST results from {result_file}")
        try:
            df = pd.read_csv(result_file, names=[
                'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
            ])
            logger.info(f"Read {len(df)} rows from BLAST results")

            for _, row in df.iterrows():
                _, cluster, gene_name = row['sseqid'].split('|')
                identity = row['pident']
                alignment_length = row['length']
                coverage = alignment_length / self.gene_lengths[row['sseqid']]

                if coverage >= 0.85 and identity >= 90.0:
                    if cluster not in self.results:
                        self.results[cluster] = set()
                    self.results[cluster].add(gene_name)

            self.replace_locus_tags_with_gene_names()
            logger.info(f"Parsed results: {self.results}")
        except Exception as e:
            logger.error(f"Error parsing BLAST results: {e}")
            raise

    def get_results(self) -> Dict[str, Set[str]]:
        return self.results

    def print_results(self, genome_name: str) -> None:
        print(f"\nResults for genome: {genome_name}")
        print("=" * 50)
        for cluster, genes in self.results.items():
            print(f"  {cluster}: {', '.join(sorted(genes))}")
        print("=" * 50)


    def process_all_genomes_in_folder(self, genomes_folder: Path) -> None:
        genomes_folder_path = genomes_folder
        logger.info(f"Processing all genomes in folder: {genomes_folder_path}")
        genome_count = 0

        # Process fastas in subdirectories
        for subfolder in genomes_folder_path.glob('*'):
            if subfolder.is_dir():
                for genome_file in subfolder.glob('*.fasta'):
                    if genome_file.is_file():
                        logger.info(f"Processing genome file from subfolder: {genome_file}")
                        try:
                            self.run_blast_for_genome(genome_file)
                            species_csv = genomes_folder_path / "species_finder" / f"{genome_file.stem}.csv"
                            self.update_species_csv(species_csv, genome_file.stem)
                            self.print_results(genome_file.stem)
                            genome_count += 1
                        except Exception as e:
                            logger.error(f"Error processing genome {genome_file.stem}: {e}")

        # Process fastas directly in the provided folder
        for genome_file in genomes_folder_path.glob('*.fasta'):
            if genome_file.is_file():
                logger.info(f"Processing genome file from main folder: {genome_file}")
                try:
                    self.run_blast_for_genome(genome_file)
                    species_csv = genomes_folder_path.parent / "species_finder" / f"{genome_file.stem}.csv"
                    self.update_species_csv(species_csv, genome_file.stem)
                    self.print_results(genome_file.stem)
                    genome_count += 1
                except Exception as e:
                    logger.error(f"Error processing genome {genome_file.stem}: {e}")

        logger.info(f"Processed {genome_count} genomes in total")

        if genome_count == 0:
            logger.warning("No genome files were found or processed.")


    def update_species_csv(self, csv_file: Path, genome_name: str) -> None:
        logger.info(f"Updating species CSV file: {csv_file}")
        if not csv_file.exists():
            logger.error(f"CSV file not found: {csv_file}")
            return
        try:
            df = pd.read_csv(csv_file)

            for strain in self.present_genes:
                for gene_group, genes in self.present_genes[strain].items():
                    if gene_group in ['other', 'T3SS']:
                        column_name = f"{gene_group}_presence"
                    else:
                        column_name = f"{strain}_{gene_group}_presence"

                    df[column_name] = ', '.join(genes)

            df.to_csv(csv_file, index=False)
            logger.info(f"Updated CSV file {csv_file} for genome {genome_name} with virulence genes")

        except Exception as e:
            logger.error(f"Error updating CSV file {csv_file}: {e}")

    def load_gene_mapping(self) -> Dict:
        if not GENES_MAP.exists():
            logger.error(f"Gene mapping file not found: {GENES_MAP}")
            raise FileNotFoundError(f"Gene mapping file not found: {GENES_MAP}")

        with open(GENES_MAP, 'r') as f:
            return json.load(f)

    def locus_tag_to_gene_name(self, locus_tag: str) -> str:
        for accession, clusters in self.gene_mapping.items():
            for cluster, genes in clusters.items():
                for gene, locus_tags in genes.items():
                    if locus_tag in locus_tags:
                        return gene
        return locus_tag  # Return the original locus_tag if not found

    def replace_locus_tags_with_gene_names(self) -> None:
        for cluster, genes in self.results.items():
            self.results[cluster] = {self.locus_tag_to_gene_name(gene) for gene in genes}


# Example usage
if __name__ == "__main__":
    virulence_genes = VirulenceGenes()
    virulence_genes.process_all_genomes_in_folder(
        Path("/Users/josediogomoura/Documents/BioFago/BioFago/test-data/random_ea_genomes"))

    # Print overall results after processing all genomes
    print("\nOverall Results:")
    print("=" * 50)
    for genome, results in virulence_genes.get_results().items():
        print(f"\nGenome: {genome}")
        for strain, gene_groups in results.items():
            print(f"  Strain: {strain}")
            for gene_group, genes in gene_groups.items():
                print(f"    {gene_group}: {', '.join(genes)}")
    print("=" * 50)