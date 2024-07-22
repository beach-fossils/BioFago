import subprocess
from pathlib import Path
from Bio import SeqIO, SeqRecord
import logging
import pandas as pd
from typing import List, Dict
import tempfile
import os
from src.utils.get_prokka_faa_file import get_prokka_faa_file
from src.utils.prokka_docker import run_prokka_docker

# Set up logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Dynamically get the reference_resistance_genes folder
current_script_path = Path(__file__).resolve()
project_root = current_script_path.parent.parent.parent  # Adjust as necessary based on script location
REFERENCE_GENES = project_root / "reference_resistance_genes"
RPSL_REF = REFERENCE_GENES / "rpsL.fasta"


class StrResistance:
    def __init__(self, reference_folder: Path = REFERENCE_GENES) -> None:
        self.reference_folder = reference_folder
        self.gene_lengths = {}
        self.results = []
        self.present_genes = []
        self.missing_genes = []
        self.rpsL_amino_acid_43 = None


    def combine_gene_sequences(self, temp_dir: Path) -> Path:
        """Combine all gene sequences into a single FASTA file and store lengths."""
        try:
            combined_gene_fasta = temp_dir / "combined_genes.fasta"
            with open(combined_gene_fasta, 'w') as outfile:
                for gene_file in self.reference_folder.glob('*.fasta'):
                    with open(gene_file, 'r') as infile:
                        seq_records = list(SeqIO.parse(infile, "fasta"))
                        for record in seq_records:
                            record.id = gene_file.stem
                            record.description = ""
                            self.gene_lengths[record.id] = len(record.seq)
                            SeqIO.write(record, outfile, "fasta")
            logger.info(f"Combined gene sequences into {combined_gene_fasta}")
            return combined_gene_fasta
        except Exception as e:
            logger.error(f"Error combining gene sequences: {e}")
            raise

    def create_blast_db(self, combined_gene_fasta: Path, temp_dir: Path) -> Path:
        """Create a single BLAST database for all combined gene sequences."""
        blast_db_name = temp_dir / "combined_genes"
        try:
            make_db_cmd = [
                'makeblastdb',
                '-in', str(combined_gene_fasta),
                '-dbtype', 'nucl',
                '-out', str(blast_db_name)
            ]
            subprocess.run(make_db_cmd, check=True)
            logger.info(f"BLAST database created at {blast_db_name}")
            return blast_db_name
        except subprocess.CalledProcessError as e:
            logger.error(f"Error in creating BLAST database: {e}")
            raise

    def run_blast_for_genome(self, genome_file: Path) -> None:
        """Run BLAST for the genome file against the combined gene database."""
        self.genome_file = genome_file
        self.output_folder = genome_file.parent.parent / "resistance_finder" / self.genome_file.stem
        self.output_folder.mkdir(parents=True, exist_ok=True)
        self.results_folder = self.output_folder / "results"
        self.results_folder.mkdir(exist_ok=True)

        self.results = []
        self.present_genes = []
        self.missing_genes = []

        try:
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
                subprocess.run(blast_cmd, check=True)
                logger.info(f"BLAST search completed for {self.genome_file} with results saved to {result_file}")

                # Parse results
                self._parse_results(result_file)
                # Handle rpsL separately
                self.extract_rpsL_from_results(result_file)

        except subprocess.CalledProcessError as e:
            logger.error(f"Error in running BLAST search: {e}")
        except Exception as e:
            logger.error(f"Unexpected error: {e}")


    def _parse_results(self, result_file: Path) -> None:
        try:
            df = pd.read_csv(result_file, names=[
                'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
            ])
            logger.info(f"Number of lines in BLAST results: {len(df)}")

            found_genes = set()
            for _, row in df.iterrows():
                gene_name = row['sseqid']
                identity = row['pident']
                alignment_length = row['length']
                coverage = alignment_length / self.gene_lengths[gene_name]

                if coverage >= 0.8 and identity >= 90.0:
                    self.present_genes.append((gene_name, identity, coverage))
                    found_genes.add(gene_name)

            all_genes = set(self.gene_lengths.keys())
            self.missing_genes = list(all_genes - found_genes)

            self.present_genes = list(set(self.present_genes))  # Remove duplicates
            logger.info(f"Present genes: {self.present_genes}")
            logger.info(f"Missing genes: {self.missing_genes}")

        except FileNotFoundError:
            logger.error(f"Result file not found: {result_file}")
        except Exception as e:
            logger.error(f"Error parsing results from {result_file}: {e}")

    def extract_rpsL_from_results(self, result_file: Path) -> None:
        """Extract the rpsL gene with flanking regions from the BLAST results."""
        try:
            df = pd.read_csv(result_file, names=[
                'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
            ])

            rpsL_row = df[df['sseqid'] == 'rpsL'].iloc[0]
            if not rpsL_row.empty:
                contig = rpsL_row['qseqid']
                start = rpsL_row['qstart']
                end = rpsL_row['qend']
                flanking_start = max(0, start - 50)
                flanking_end = end + 50

                self.extract_gene_sequence_with_flanking(
                    genome_file=self.genome_file,
                    contig=contig,
                    start=flanking_start,
                    end=flanking_end,
                    output_path=self.output_folder / f"{self.genome_file.stem}_rpsL_extracted.fasta"
                )
                # Annotate rpsL after extracting
                self.annotate_rpsL(self.genome_file.stem)
            else:
                logger.error(f"rpsL gene not found in BLAST results for {self.genome_file}")

        except Exception as e:
            logger.error(f"Error extracting rpsL gene from BLAST results: {e}")

    def process_all_genomes_in_folder(self, genomes_folder: Path) -> None:
        """Process all genome files in their respective subfolders."""
        genomes_folder_path = genomes_folder
        for genome_subfolder in genomes_folder_path.glob('*'):
            if genome_subfolder.is_dir():
                genome_file = genome_subfolder / f"{genome_subfolder.stem}.fasta"
                if genome_file.exists():
                    logger.info(f"Processing genome file: {genome_file}")
                    try:
                        self.run_blast_for_genome(genome_file)
                        self.update_species_csv(genomes_folder_path / "species_finder" / f"{genome_subfolder.stem}.csv",
                                                genome_subfolder.stem)
                        present_genes = self.present_genes
                        logger.info(f"Present genes for {genome_subfolder.stem}: {present_genes}")
                        logger.info(f"Missing genes for {genome_subfolder.stem}: {self.missing_genes}")
                    except Exception as e:
                        logger.error(f"Error processing genome {genome_subfolder.stem}: {e}")

    def extract_gene_sequence_with_flanking(self, genome_file: Path, contig: str, start: int, end: int,
                                            output_path: Path) -> None:
        """Extract the gene sequence with flanking regions from the genome."""
        try:
            for record in SeqIO.parse(genome_file, "fasta"):
                if record.id == contig:
                    extracted_seq = record.seq[start:end]
                    extracted_record = SeqRecord.SeqRecord(
                        seq=extracted_seq,
                        id=f"{contig}_{start}_{end}"[:37],  # Ensure ID is <= 37 characters
                        description=f"Extracted sequence from {contig}:{start}-{end}"
                    )

                    with open(output_path, 'w') as out_file:
                        SeqIO.write(extracted_record, out_file, "fasta")
                    logger.info(f"Extracted gene sequence from {genome_file} and saved to {output_path}")
                    return output_path
            logger.error(f"Contig {contig} not found in {genome_file}")
        except Exception as e:
            logger.error(f"Error extracting gene sequence from {genome_file}: {e}")

    def annotate_rpsL(self, genome_name: str) -> None:
        """Annotate the extracted rpsL gene sequence using Prokka."""
        try:
            base_output_folder = self.output_folder.parent / "prokka" / genome_name
            base_output_folder.mkdir(parents=True, exist_ok=True)
            rpsL_fasta = self.output_folder / f"{genome_name}_rpsL_extracted.fasta"

            if rpsL_fasta.exists():
                run_prokka_docker(
                    fasta_file=rpsL_fasta,
                    base_output_folder=base_output_folder,
                    custom_db_path=str(RPSL_REF),
                    locus_tag_prefix='rpsL'
                )

                faa_file = get_prokka_faa_file(base_output_folder)
                if faa_file:
                    logger.info(f"Prokka .faa file for {genome_name}: {faa_file}")
                    self.save_amino_acid_at_position(faa_file, 43, genome_name)
                else:
                    logger.error(f"Prokka .faa file not found for {genome_name}")

        except Exception as e:
            logger.error(f"Error annotating rpsL for {genome_name}: {e}")

    def save_amino_acid_at_position(self, faa_file: Path, position: int, genome_name: str) -> None:
        try:
            for record in SeqIO.parse(faa_file, "fasta"):
                if len(record.seq) >= position:
                    amino_acid = record.seq[position - 1]  # Convert to 0-based index
                    self.rpsL_amino_acid_43 = amino_acid  # Store the amino acid
                    result_file = self.output_folder / f"{genome_name}_rpsL_amino_acid_position_{position}.txt"
                    with open(result_file, 'w') as out_file:
                        out_file.write(f"Amino acid at position {position} in {genome_name}: {amino_acid}\n")
                    logger.info(f"Saved amino acid at position {position} for {genome_name} to {result_file}")
                else:
                    logger.error(f"Sequence in {faa_file} is shorter than {position} amino acids")
        except Exception as e:
            logger.error(f"Error saving amino acid at position {position} for {genome_name}: {e}")

    def update_species_csv(self, csv_file: Path, genome_name: str) -> None:
        """Update the species CSV file with the resistance gene results and amino acid at position 43 of rpsL."""
        if not csv_file.exists():
            logger.error(f"CSV file not found: {csv_file}")
            return
        try:
            df = pd.read_csv(csv_file)
            # Add strA and strB presence/absence
            for gene in ['strA', 'strB']:
                df[f'{gene} presence/absence'] = 'Present' if any(
                    g == gene for g, _, _ in self.present_genes) else 'Absent'

            # Add rpsL amino acid at position 43
            amino_acid_file = self.output_folder / f"{genome_name}_rpsL_amino_acid_position_43.txt"
            if amino_acid_file.exists():
                with open(amino_acid_file, 'r') as f:
                    aa_line = f.readline().strip()
                    aa_at_43 = aa_line.split(': ')[-1]
                df['rpsL_aa_pos43'] = aa_at_43

            df.to_csv(csv_file, index=False)
            logger.info(
                f"Updated CSV file {csv_file} for genome {genome_name} with resistance genes and rpsL amino acid position 43")

        except Exception as e:
            logger.error(f"Error updating CSV file {csv_file}: {e}")

    def get_results(self) -> Dict:
        return {
            'present_genes': self.present_genes,
            'missing_genes': self.missing_genes,
            'rpsL_amino_acid_43': self.rpsL_amino_acid_43
        }


# Example usage
if __name__ == "__main__":
    sorbitol_resistance = StrResistance()
    sorbitol_resistance.process_all_genomes_in_folder(
        Path("/Users/josediogomoura/Documents/BioFago/BioFago/data/test_resistance/genome/"))
