#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Set
import tempfile
import shutil
from Bio import SeqIO
from dataclasses import dataclass
import json

# Dynamic Path Setup
BASE_PATH = Path(__file__).resolve().parent.parent.parent
REFERENCE_MLST_PATH = BASE_PATH / 'reference_mlst'
JSON_PATTERNS = REFERENCE_MLST_PATH / 'mlst_patterns.json'

# Configure logging with more details
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class BlastResult:
    """Data class for storing BLAST results."""
    query_id: str
    subject_id: str
    identity: float
    alignment_length: int
    mismatches: int
    gaps: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float


class MLSTTyper:
    """MLST typing system with dynamic path handling."""

    def __init__(
            self,
            min_identity: float = 99.0
    ):
        """Initialize MLSTTyper."""
        self.reference_dir = REFERENCE_MLST_PATH
        self.min_identity = min_identity

        # Create temporary directory for all files
        self.temp_dir = Path(tempfile.mkdtemp(prefix='mlst_typer_'))
        logger.info(f"Created temporary directory: {self.temp_dir}")

        # Create subdirectory for BLAST databases
        self.blast_db_dir = self.temp_dir / 'blast_dbs'
        self.blast_db_dir.mkdir()
        logger.info(f"Created BLAST databases directory: {self.blast_db_dir}")

        # Gene order for MLST scheme
        self.genes = ['adk', 'arcA', 'mdh', 'recA', 'dnaA',
                      'pgi', 'fusA', 'gyrB', 'infB', 'rpoB']

        logger.info(f"Initializing MLST typer with {len(self.genes)} genes")

        # Validate and setup
        self._validate_reference_directory()
        self._create_blast_databases()

    def _validate_reference_directory(self) -> None:
        """Validate reference directory structure."""
        logger.info("Validating reference directory structure...")

        if not self.reference_dir.is_dir():
            raise ValueError(f"Reference directory not found: {self.reference_dir}")

        for gene in self.genes:
            gene_dir = self.reference_dir / gene
            if not gene_dir.is_dir():
                raise ValueError(f"Missing gene directory: {gene}")

            allele_files = list(gene_dir.glob("allele_*.fasta"))
            if not allele_files:
                raise ValueError(f"No allele files found for gene: {gene}")
            logger.info(f"Found {len(allele_files)} allele files for {gene}")

    def _create_blast_databases(self) -> None:
        """Create temporary BLAST databases for each gene."""
        logger.info("Creating BLAST databases...")

        for gene in self.genes:
            logger.info(f"Processing {gene}...")
            gene_dir = self.reference_dir / gene

            # Create temp directory for this gene
            temp_gene_dir = self.blast_db_dir / gene
            temp_gene_dir.mkdir()

            # Combine alleles in temporary location with proper headers
            combined_fasta = temp_gene_dir / "combined_alleles.fasta"
            with open(combined_fasta, 'w') as outfile:
                for allele_file in sorted(gene_dir.glob("allele_*.fasta")):
                    allele_num = allele_file.stem.split('_')[1]
                    with open(allele_file) as infile:
                        # Write with simplified header containing gene and allele number
                        outfile.write(f">{gene}_{allele_num}\n")
                        # Skip original header and write sequence
                        next(infile)
                        for line in infile:
                            if not line.startswith('>'):
                                outfile.write(line)

            # Create BLAST database in temporary location
            db_path = temp_gene_dir / gene
            cmd = [
                'makeblastdb',
                '-in', str(combined_fasta),
                '-dbtype', 'nucl',
                '-out', str(db_path)
            ]
            try:
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                logger.info(f"Created BLAST database for {gene}")
            except subprocess.CalledProcessError as e:
                logger.error(f"BLAST database creation failed for {gene}")
                logger.error(f"Error output: {e.stderr}")
                raise RuntimeError(f"Failed to create BLAST database for {gene}")

    def _run_blast(self, genome_path: Path, gene: str) -> Optional[BlastResult]:
        """Run BLAST for a specific gene."""
        logger.info(f"Running BLAST for {gene}...")

        # Use temporary directory for database
        db_path = self.blast_db_dir / gene / gene
        output_file = self.temp_dir / f"{gene}_blast.txt"

        cmd = [
            'blastn',
            '-query', str(genome_path),
            '-db', str(db_path),
            '-outfmt', '6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore',
            '-out', str(output_file),
            '-max_target_seqs', '1'
        ]

        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            if not output_file.stat().st_size:
                logger.warning(f"No BLAST hits found for {gene}")
                return None

            with open(output_file) as f:
                line = f.readline().strip()
                if not line:
                    return None

                fields = line.split('\t')
                fields[2] = float(fields[2])  # Convert identity to float

                # Log BLAST result details
                logger.info(f"{gene} BLAST result: Identity={fields[2]}%, Length={fields[3]}, E-value={fields[10]}")

                return BlastResult(*fields)

        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST error for {gene}: {e.stderr}")
            return None
        except Exception as e:
            logger.error(f"Error processing BLAST results for {gene}: {e}")
            return None


    def cleanup(self) -> None:
        """Clean up all temporary files."""
        logger.info("Cleaning up temporary files...")
        try:
            shutil.rmtree(self.temp_dir)
            logger.info("Successfully removed temporary directory")
        except Exception as e:
            logger.warning(f"Error during cleanup: {e}")

    # Add this method to your MLSTTyper class
    def match_mlst_pattern(self, allele_combination: str) -> dict:
        """
        Match allele combination with known patterns from JSON.

        Args:
            allele_combination: String of allele numbers separated by hyphens

        Returns:
            dict: Pattern information including ST and clades, or Novel status
        """
        try:
            with open(JSON_PATTERNS) as f:
                patterns = json.load(f)

            # Search for the combination in patterns
            for st, info in patterns.items():
                if info['allele_combination'] == allele_combination:
                    return {
                        'status': 'Known',
                        'st': st,
                        'clades': info['clades']
                    }

            # If no match found, it's a novel combination
            return {
                'status': 'Novel',
                'st': 'Novel',
                'clades': ['Novel combination']
            }

        except Exception as e:
            logger.error(f"Error matching MLST pattern: {e}")
            return {
                'status': 'Error',
                'st': 'Unknown',
                'clades': ['Error matching pattern']
            }

    # Then modify your type_genome method to include pattern matching:
    def type_genome(self, genome_path: str) -> Optional[dict]:
        """Perform MLST typing on input genome."""
        genome_path = Path(genome_path)
        if not genome_path.exists():
            raise FileNotFoundError(f"Genome file not found: {genome_path}")

        logger.info(f"Starting genome typing for: {genome_path}")
        allele_calls = {}

        for gene in self.genes:
            blast_result = self._run_blast(genome_path, gene)

            if not blast_result:
                logger.error(f"No BLAST result for {gene}")
                return None

            if blast_result.identity < self.min_identity:
                logger.error(
                    f"Low identity match for {gene}: {blast_result.identity}% (threshold: {self.min_identity}%)")
                return None

            try:
                gene_name, allele_num = blast_result.subject_id.split('_')
                allele_num = int(allele_num)
                allele_calls[gene] = allele_num
                logger.info(f"Assigned allele {allele_num} to {gene}")
            except (IndexError, ValueError) as e:
                logger.error(f"Could not determine allele number for {gene}: {e}")
                logger.error(f"BLAST subject_id was: {blast_result.subject_id}")
                return None

        # Create the allele combination string
        allele_combination = "-".join(str(allele_calls[gene]) for gene in self.genes)
        logger.info(f"Final MLST profile: {allele_combination}")

        # Match with known patterns
        pattern_info = self.match_mlst_pattern(allele_combination)

        return {
            'allele_combination': allele_combination,
            'st': pattern_info['st'],
            'status': pattern_info['status'],
            'clades': pattern_info['clades']
        }


def main():
    """Command line interface."""
    try:
        # Initialize the typer
        typer = MLSTTyper()

        # Test genome path
        test_genome = "/Users/josediogomoura/Documents/BioFago/BioFago/tests/genomes/GCA_000027205.1_ASM2720v1_genomic.fasta"

        logger.info(f"Starting MLST typing for genome: {test_genome}")
        result = typer.type_genome(test_genome)

        if result:
            print("\n=== MLST Typing Results ===")
            print(f"Successfully typed genome!")
            print(f"Allele Combination: {result['allele_combination']}")
            print(f"Sequence Type: {result['st']}")
            print(f"Status: {result['status']}")
            print("Clades:")
            for clade in result['clades']:
                print(f"  - {clade}")
            print("=========================\n")
        else:
            print("\n=== MLST Typing Failed ===")
            print("Could not determine complete MLST profile")
            print("Check logs for detailed error information")
            print("=========================\n")

    except Exception as e:
        logger.error(f"Critical error during MLST typing: {e}")
        sys.exit(1)
    finally:
        typer.cleanup()


if __name__ == "__main__":
    main()