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
import argparse

# Import quiet mode module
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from quiet_mode import QUIET_MODE

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

# Allow external configuration of log level
def set_log_level(level):
    """Set the log level for this module."""
    logger.setLevel(level)
    logger.debug(f"MLST typer log level set to {level}")


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
            min_identity: float = 98.5
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
                if QUIET_MODE:
                    with open(os.devnull, 'w') as devnull:
                        subprocess.run(cmd, check=True, stdout=devnull, stderr=devnull)
                else:
                    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                logger.info(f"Created BLAST database for {gene}")
            except subprocess.CalledProcessError as e:
                logger.error(f"BLAST database creation failed for {gene}")
                logger.error(f"Error output: {e.stderr if hasattr(e, 'stderr') else 'N/A'}")
                raise RuntimeError(f"Failed to create BLAST database for {gene}")

    def _run_blast(self, genome_path: Path, gene: str) -> Optional[BlastResult]:
        """
        Run BLAST for a specific gene and select the best hit.
        The best hit is determined first by highest identity, then by coverage relative to the allele.
        """
        logger.info(f"Running BLAST for {gene}...")
        logger.debug(f"Genome path: {genome_path}")

        # Use temporary directory for database
        db_path = self.blast_db_dir / gene / gene
        output_file = self.temp_dir / f"{gene}_blast.txt"
        logger.debug(f"BLAST database path: {db_path}")
        logger.debug(f"BLAST output file: {output_file}")

        cmd = [
            'blastn',
            '-query', str(genome_path),
            '-db', str(db_path),
            '-outfmt', '6 qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore',
            '-out', str(output_file),
            '-max_target_seqs', '10'  # Increased to allow finding multiple hits
        ]
        logger.debug(f"BLAST command: {' '.join(cmd)}")

        try:
            if QUIET_MODE:
                logger.debug("Running in QUIET mode, suppressing subprocess output")
                with open(os.devnull, 'w') as devnull:
                    subprocess.run(cmd, check=True, stdout=devnull, stderr=devnull)
            else:
                logger.debug("Running with output visible")
                result = subprocess.run(cmd, check=True, capture_output=True, text=True)
                logger.debug(f"BLAST stdout: {result.stdout[:200] if result.stdout else 'empty'}")
                logger.debug(f"BLAST stderr: {result.stderr[:200] if result.stderr else 'empty'}")

            if not output_file.stat().st_size:
                logger.warning(f"No BLAST hits found for {gene}")
                return None

            # First group hits by subject_id (allele)
            allele_hits = {}
            with open(output_file) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    fields = line.split('\t')
                    if len(fields) < 12:  # Ensure all fields are present
                        logger.debug(f"Invalid BLAST output line: {line}")
                        continue
                    
                    # Convert numeric fields
                    fields[2] = float(fields[2])  # identity
                    fields[3] = int(fields[3])    # alignment length
                    fields[11] = float(fields[11])  # bitscore
                    
                    # Group by subject_id (allele)
                    subject_id = fields[1]
                    if subject_id not in allele_hits:
                        allele_hits[subject_id] = []
                    allele_hits[subject_id].append(fields)
            
            # For each allele, select the hit with highest bitscore
            best_hits_per_allele = []
            for subject_id, hits in allele_hits.items():
                best_hit = max(hits, key=lambda x: float(x[11]))  # sort by bitscore
                best_hits_per_allele.append(best_hit)
            
            # First prioritize 100% identity matches
            perfect_matches = [hit for hit in best_hits_per_allele if hit[2] == 100.0]
            
            if perfect_matches:
                logger.debug(f"Found {len(perfect_matches)} perfect matches for {gene}")
                # Among perfect matches, select the one with highest coverage
                best_hit = max(perfect_matches, key=lambda x: int(x[3]))
            else:
                # Otherwise, sort by identity first, then by alignment length
                logger.debug(f"No perfect matches, selecting best among {len(best_hits_per_allele)} candidate alleles")
                best_hit = max(best_hits_per_allele, key=lambda x: (float(x[2]), int(x[3])))
            
            # Log the selected best hit
            logger.info(f"{gene} BLAST result: Identity={best_hit[2]}%, Length={best_hit[3]}, E-value={best_hit[10]}")
            logger.debug(f"Selected best hit: Subject={best_hit[1]}, Bitscore={best_hit[11]}")
            
            # Create BlastResult from the best hit
            return BlastResult(*best_hit)

        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST error for {gene}: {e.stderr if hasattr(e, 'stderr') else 'N/A'}")
            return None
        except Exception as e:
            logger.error(f"Error processing BLAST results for {gene}: {e}")
            logger.debug(f"Exception details:", exc_info=True)
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
            dict: Pattern information including ST and status
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
                        'clades': [st]  # Use ST as the clade if no clades key exists
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
        # Allow log level to be set from command line
        parser = argparse.ArgumentParser(description='MLST Typer')
        parser.add_argument('--genome', type=str, help='Path to genome FASTA file',
                           default="/Volumes/Crucial_X9/biofago_2025/c_BioFago_ompA_names_fixed_last_version/tests/GCF_012367585.1_GCF_012367585.1_ASM1236758v1_genomic.fna")
        parser.add_argument('--log_level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], 
                           default='INFO', help='Set logging level')
        args = parser.parse_args()
        
        # Set log level
        log_level = getattr(logging, args.log_level)
        logger.setLevel(log_level)
        logger.debug(f"Log level set to {args.log_level}")
        
        # Initialize the typer
        typer = MLSTTyper()

        # Test genome path (use command line arg if provided)
        test_genome = args.genome

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