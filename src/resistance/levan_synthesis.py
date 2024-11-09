import pandas as pd
from Bio import SeqIO
import subprocess
import os
import tempfile
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union

# Set up logging with more detailed output
logging.basicConfig(
    level=logging.DEBUG,  # Changed to DEBUG for more detail
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('levan_synthesis')

# Dynamically get the reference_plasmids folder
current_script_path = Path(__file__).resolve()
project_root = current_script_path.parent.parent.parent
REFERENCE_LEVAN = project_root / "reference_resistance_genes" / "levan_synthesis"


class LevanGenesClassifier:
    def __init__(self, reference_dir: Path) -> None:
        self.identity_threshold = 75
        self.coverage_threshold = 80
        self.reference_dir = Path(reference_dir)
        self.reference_files = {
            'lsc': self.reference_dir / 'lsc.fasta',
            'rlsA': self.reference_dir / 'rlsA.fasta',
            'rlsB': self.reference_dir / 'rlsB.fasta'
        }

        # Load and verify reference sequences
        self.reference_sequences = self._load_reference_sequences()

    def _load_reference_sequences(self) -> Dict[str, str]:
        """Load and verify reference sequences"""
        sequences = {}
        for gene, file_path in self.reference_files.items():
            if not file_path.exists():
                logger.error(f"Reference file not found for {gene}: {file_path}")
                raise FileNotFoundError(f"Reference file not found: {file_path}")

            try:
                with open(file_path) as f:
                    record = next(SeqIO.parse(f, "fasta"))
                    sequences[gene] = str(record.seq)
                    logger.info(f"Loaded reference sequence for {gene}: {len(sequences[gene])} bp")
            except Exception as e:
                logger.error(f"Error loading reference sequence for {gene}: {e}")
                raise

        return sequences

    def run_blast(self, query_file: Path, reference_file: Path) -> Tuple[float, float]:
        """
        Run BLASTN with proper coverage calculation
        We're looking for the gene in the genome, so we care about coverage of the reference gene
        """
        try:
            if not os.path.getsize(query_file) > 0:
                logger.error(f"Query file is empty: {query_file}")
                return 0.0, 0.0

            if not os.path.getsize(reference_file) > 0:
                logger.error(f"Reference file is empty: {reference_file}")
                return 0.0, 0.0

            cmd = [
                'blastn',
                '-query', str(query_file),
                '-subject', str(reference_file),
                '-outfmt', '6 qseqid sseqid pident length qlen slen',  # Simplified format
                '-max_target_seqs', '1'
            ]

            logger.debug(f"Running BLAST command: {' '.join(cmd)}")

            result = subprocess.run(cmd,
                                 capture_output=True,
                                 text=True,
                                 check=True)

            if result.stdout.strip():
                fields = result.stdout.strip().split('\t')
                logger.debug(f"Raw BLAST output: {fields}")

                if len(fields) >= 6:
                    identity = float(fields[2])
                    alignment_length = int(fields[3])
                    subject_length = int(fields[5])  # Reference gene length

                    # Calculate coverage based on how much of the reference gene is aligned
                    coverage = (alignment_length / subject_length) * 100

                    logger.info(
                        f"BLAST hit found:\n"
                        f"  Identity: {identity:.1f}%\n"
                        f"  Coverage: {coverage:.1f}% ({alignment_length} bp aligned out of {subject_length} bp reference)"
                    )

                    return identity, coverage

            logger.warning("No BLAST hits found")
            return 0.0, 0.0

        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST error: {e}")
            logger.error(f"BLAST stderr: {e.stderr}")
            return 0.0, 0.0
        except Exception as e:
            logger.error(f"Unexpected error running BLAST: {e}")
            return 0.0, 0.0

    def analyze_strain(self, sequences: Dict[str, str], strain_id: str) -> Dict:
        """Analyze strain with correct BLAST search direction"""
        logger.info(f"Analyzing strain: {strain_id}")

        results = {
            'strain_id': strain_id,
            'lsc_identity': 0.0,
            'lsc_coverage': 0.0,
            'rlsA_identity': 0.0,
            'rlsA_coverage': 0.0,
            'rlsB_identity': 0.0,
            'rlsB_coverage': 0.0
        }

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = Path(temp_dir)

            # Write genome sequences to a single file
            genome_file = temp_dir_path / "genome.fasta"
            with open(genome_file, 'w') as f:
                for seq_id, sequence in sequences.items():
                    f.write(f">{seq_id}\n{sequence}\n")

            # Process each reference gene
            for gene in ['lsc', 'rlsA', 'rlsB']:
                try:
                    cmd = [
                        'blastn',
                        '-query', str(self.reference_files[gene]),  # Reference gene as query
                        '-subject', str(genome_file),  # Genome as subject
                        '-outfmt', '6 qseqid sseqid pident length qlen slen sstart send',
                        '-max_target_seqs', '5',
                        '-task', 'blastn',
                        '-dust', 'no',
                        '-evalue', '1e-10',
                        '-perc_identity', '85'
                    ]

                    logger.debug(f"Running BLAST command for {gene}: {' '.join(cmd)}")

                    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

                    if result.stdout.strip():
                        hits = result.stdout.strip().split('\n')
                        best_identity = 0.0
                        best_coverage = 0.0
                        query_length = None

                        for hit in hits:
                            fields = hit.split('\t')
                            if len(fields) >= 8:
                                identity = float(fields[2])
                                alignment_length = int(fields[3])
                                query_length = int(fields[4])  # Reference gene length
                                sstart = min(int(fields[6]), int(fields[7]))
                                send = max(int(fields[6]), int(fields[7]))

                                # Track best identity
                                best_identity = max(best_identity, identity)

                                # Calculate coverage of reference gene
                                coverage = (alignment_length / query_length) * 100
                                best_coverage = max(best_coverage, coverage)

                        if best_coverage > 0:
                            logger.info(
                                f"{gene} hits found:\n"
                                f"  Best Identity: {best_identity:.1f}%\n"
                                f"  Best Coverage: {best_coverage:.1f}%"
                            )
                            results[f'{gene}_identity'] = best_identity
                            results[f'{gene}_coverage'] = best_coverage
                    else:
                        logger.warning(f"No BLAST hits found for {gene}")

                except subprocess.CalledProcessError as e:
                    logger.error(f"BLAST error for {gene}: {e.stderr}")
                except Exception as e:
                    logger.error(f"Error processing {gene}: {e}")

            # Classify strain based on results
            results['strain_type'] = self._classify_strain_type(
                results['lsc_identity'], results['lsc_coverage'],
                results['rlsA_identity'], results['rlsA_coverage'],
                results['rlsB_identity'], results['rlsB_coverage']
            )

        return results

    def _classify_strain_type(self, lsc_id: float, lsc_cov: float,
                            rlsA_id: float, rlsA_cov: float,
                            rlsB_id: float, rlsB_cov: float) -> str:
        """Classify strain type with improved logging"""
        # Check if genes are present (coverage >= threshold)
        lsc_present = lsc_cov >= self.coverage_threshold
        rlsA_present = rlsA_cov >= self.coverage_threshold
        rlsB_present = rlsB_cov >= self.coverage_threshold

        # Check if genes have high identity when present
        rlsA_variant = rlsA_present and rlsA_id < self.identity_threshold

        logger.debug(
            f"\nGene Status:\n"
            f"LSC:  Present={lsc_present}, Identity={lsc_id:.1f}%, Coverage={lsc_cov:.1f}%\n"
            f"rlsA: Present={rlsA_present}, Identity={rlsA_id:.1f}%, Coverage={rlsA_cov:.1f}%\n"
            f"rlsB: Present={rlsB_present}, Identity={rlsB_id:.1f}%, Coverage={rlsB_cov:.1f}%\n"
            f"rlsA variant: {rlsA_variant}"
        )

        # Classification logic
        if all([lsc_present, rlsA_present, rlsB_present]):
            if all([lsc_id >= self.identity_threshold,
                   rlsA_id >= self.identity_threshold,
                   rlsB_id >= self.identity_threshold]):
                return "Similar to Widely-prevalent E. amylovora (Spiraeoideae-infecting)"

        if not lsc_present:
            return "Similar to E. pyrifoliae (limited host range)"

        if (not rlsA_present or rlsA_variant) and lsc_present:
            return "Similar to Rubus-infecting strain"

        if not rlsB_present and lsc_present:
            return "Similar to blossom-limited strain"

        return "Unusual pattern - further investigation needed"

def read_fasta_sequences(fasta_file: Path) -> Dict[str, str]:
    """Read sequences from a FASTA file"""
    logger.info(f"Reading sequences from {fasta_file}")
    try:
        sequences = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences[record.id] = str(record.seq)
            logger.debug(f"Read sequence {record.id}: {len(record.seq)} bp")
        logger.info(f"Read {len(sequences)} sequences from file")
        return sequences
    except Exception as e:
        logger.error(f"Error reading FASTA file: {e}")
        raise

def process_multiple_strains(reference_dir: Path, query_sequences: Dict[str, Dict[str, str]]) -> pd.DataFrame:
    """Process multiple strains and return a formatted pandas DataFrame"""
    logger.info(f"Processing {len(query_sequences)} strains")
    try:
        classifier = LevanGenesClassifier(reference_dir)
        results = []

        for strain_id, sequences in query_sequences.items():
            try:
                result = classifier.analyze_strain(sequences, strain_id)

                # Format numeric values to 2 decimal places
                for key in result.keys():
                    if isinstance(result[key], float):
                        result[key] = round(result[key], 2)

                results.append(result)
            except Exception as e:
                logger.error(f"Error processing strain {strain_id}: {e}")
                raise

        # Create DataFrame and format numeric columns
        df = pd.DataFrame(results)
        numeric_columns = ['lsc_identity', 'lsc_coverage',
                           'rlsA_identity', 'rlsA_coverage',
                           'rlsB_identity', 'rlsB_coverage']

        for col in numeric_columns:
            df[col] = df[col].round(2)

        return df
    except Exception as e:
        logger.error(f"Error in process_multiple_strains: {e}")
        raise


def run_levan_analysis(genome_path: Union[Path, str],
                       reference_dir: Union[Path, str, None] = None,
                       strain_id: Optional[str] = None,
                       output_file: Union[Path, str, None] = None) -> pd.DataFrame:
    try:
        genome_path = Path(genome_path)
        if reference_dir is None:
            reference_dir = REFERENCE_LEVAN
        else:
            reference_dir = Path(reference_dir)

        logging.info(f"Starting levan analysis for genome: {genome_path}")
        logging.info(f"Reference directory: {reference_dir}")

        # Verify input file exists and is readable
        if not genome_path.exists():
            raise FileNotFoundError(f"Genome file not found: {genome_path}")
        if not genome_path.is_file():
            raise ValueError(f"Genome path is not a file: {genome_path}")

        # Check file size
        file_size = genome_path.stat().st_size
        if file_size == 0:
            raise ValueError(f"Genome file is empty: {genome_path}")
        logging.info(f"Genome file size: {file_size} bytes")

        # Verify reference directory and files
        if not reference_dir.exists():
            raise FileNotFoundError(f"Reference directory not found: {reference_dir}")

        # Check each reference file
        required_files = ['lsc.fasta', 'rlsA.fasta', 'rlsB.fasta']
        missing_files = [f for f in required_files if not (reference_dir / f).exists()]
        if missing_files:
            raise FileNotFoundError(f"Missing reference files: {', '.join(missing_files)}")

        # Set strain ID
        if strain_id is None:
            strain_id = genome_path.stem

        # Read and validate sequences
        sequences = read_fasta_sequences(genome_path)
        if not sequences:
            raise ValueError(f"No valid sequences found in genome file: {genome_path}")

        # Process sequences
        results = process_multiple_strains(reference_dir, {strain_id: sequences})

        if output_file is not None:
            results.to_csv(output_file, index=False)
            logging.info(f"Results saved to {output_file}")

        return results

    except Exception as e:
        logging.error(f"Error in levan analysis for {genome_path}: {str(e)}")
        raise


def format_levan_result(levan_result: dict) -> str:
    """Format levan synthesis results into a single string."""
    try:
        strain_type = levan_result.get('strain_type', 'Unknown')
        if strain_type == "Similar to Widely-prevalent E. amylovora (Spiraeoideae-infecting)":
            type_text = "Similar to Widely-prevalent E. amylovora (Spiraeoideae-infecting)"
        else:
            type_text = strain_type

        return (
            f"lsc: {levan_result.get('lsc_identity', 0):.2f}% identity, {levan_result.get('lsc_coverage', 0):.2f}% coverage | "
            f"rlsA: {levan_result.get('rlsA_identity', 0):.2f}% identity, {levan_result.get('rlsA_coverage', 0):.2f}% coverage | "
            f"rlsB: {levan_result.get('rlsB_identity', 0):.2f}% identity, {levan_result.get('rlsB_coverage', 0):.2f}% coverage | "
            f"Note: {type_text}"
        )
    except Exception as e:
        logging.error(f"Error formatting levan result: {e}")
        return "Error formatting levan synthesis results"

if __name__ == "__main__":
    try:
        # Single genome analysis
        genome_path = "/Users/josediogomoura/Documents/BioFago/BioFago/test-data/random_ea_genomes/GCF_000367685.1_ASM36768v2_genomic.fasta"
        results = run_levan_analysis(
            genome_path=genome_path,
            output_file="/Users/josediogomoura/Documents/BioFago/BioFago/test-data/levan/results.csv"
        )
        print(results.to_string())


    except Exception as e:
        logger.error(f"Error in main: {e}")
        raise
