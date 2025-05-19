import os
import logging
import subprocess
import tempfile
import statistics
from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO, AlignIO, pairwise2, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from typing import Dict, Tuple, Optional, List, Any, Set, Union
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from quiet_mode import QUIET_MODE

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class OmpAGenotyper:
    """
    Enhanced OmpA genotyping class that analyzes genome sequences using advanced bioinformatics approaches.
    This approach focuses exclusively on the ompA gene rather than broader loci.
    """
    # Class constants for standardized analysis
    IDENT_PERFECT = 99.95  # Perfect match identity threshold
    IDENT_VERY_HIGH = 99.5  # Very high match identity threshold
    IDENT_HIGH = 99.0  # High match identity threshold
    IDENT_MODERATE = 98.0  # Moderate match identity threshold
    BLAST_EVALUE = 1e-15  # Stricter BLAST e-value threshold
    MIN_COVERAGE = 90.0  # Minimum coverage percentage
    PADDING = 100  # Standard sequence padding
    OMPA_LENGTH = 1068  # Standard ompA gene length
    PARALLEL_PROCESSES = min(4, max(1, multiprocessing.cpu_count() - 1))  # Sensible parallel processing
    
    # Paths to external tools
    # Search for Clustal Omega in common locations
    CLUSTALO_PATHS = [
        "clustalo",                        # System PATH
        "/usr/bin/clustalo",               # Common Linux location
        "/usr/local/bin/clustalo",         # Common macOS Homebrew location
        "/opt/homebrew/bin/clustalo",      # Apple Silicon macOS Homebrew location
        "C:\\Program Files\\Clustal Omega\\clustalo.exe"  # Windows location
    ]
    CLUSTALO_PATH = None  # Will be set during initialization
    
    def __init__(self, genome_path: Path, reference_ompa_db: Path, output_dir: Path):
        """
        Initialize the OmpA genotyper with optimized settings.
        
        Args:
            genome_path: Path to the input genome FASTA file
            reference_ompa_db: Path to the reference ompA database GenBank file
            output_dir: Path to the output directory
        """
        self.genome_path = Path(genome_path)
        self.reference_ompa_db = Path(reference_ompa_db)
        self.output_dir = Path(output_dir)
        self.ompa_output_dir = self.output_dir / "types_ompa"
        self.ompa_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Metadata storage for debugging and result reporting
        self.metadata = {
            "version": "2.0.0",
            "analysis_steps": [],
            "timing": {},
            "errors": [],
            "warnings": []
        }
        
        # Create temporary directory for intermediate files - using a context manager approach
        self.temp_dir = Path(tempfile.mkdtemp(prefix="ompa_genotyper_"))
        self.reference_cache = {}  # Cache for reference data to minimize file operations
        
        # Extract the genome name for reporting
        self.genome_name = self.genome_path.stem
        
        # Find Clustal Omega binary
        self._find_clustalo()
        
        # Parse reference database once during initialization
        logger.info(f"Initializing OmpA genotyper for {self.genome_name}")
        self.reference_types = self._parse_reference_db()
        
        # Pre-compute reference alignments for faster comparisons
        self._precompute_reference_alignments()
        
    def _find_clustalo(self) -> None:
        """Find Clustal Omega binary in common locations."""
        # Start by assuming clustalo is not available
        self.CLUSTALO_PATH = None
        
        # Try each potential path
        for path in self.CLUSTALO_PATHS:
            try:
                # Run a simple version check
                result = subprocess.run(
                    [path, "--version"], 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE, 
                    timeout=5
                )
                if result.returncode == 0:
                    self.CLUSTALO_PATH = path
                    logger.info(f"Found Clustal Omega at: {path}")
                    break
            except (subprocess.SubprocessError, FileNotFoundError):
                continue
        
        if self.CLUSTALO_PATH:
            logger.info("Clustal Omega is available for multiple sequence alignment.")
        else:
            logger.warning("Clustal Omega not found. Will use Biopython for alignments instead.")
            self.metadata["warnings"].append("Clustal Omega not available - using Biopython for alignments")
        
    def _parse_reference_db(self) -> Dict[str, Dict[str, Any]]:
        """
        Parse the reference GenBank file containing ompA types with enhanced error handling.
        
        Returns:
            Dictionary of ompA types with their properties
        """
        reference_types = {}
        
        try:
            # Load the reference database
            records = list(SeqIO.parse(self.reference_ompa_db, "genbank"))
            if not records:
                raise ValueError(f"No valid records found in reference database: {self.reference_ompa_db}")
                
            logger.info(f"Found {len(records)} reference records to parse")
            
            for record in records:
                locus_tag = None
                amino_acid_seq = None
                gene_coords = None
                
                # Check for the target gene in features
                ompa_features = []
                for feature in record.features:
                    if feature.type == "CDS" and "gene" in feature.qualifiers:
                        if feature.qualifiers["gene"][0].lower() == "ompa":
                            ompa_features.append(feature)
                
                if not ompa_features:
                    logger.warning(f"No ompA gene found in record {record.id}")
                    continue
                
                # Use the first ompA feature (should only be one, but just in case)
                feature = ompa_features[0]
                
                # Extract information from feature
                if "locus_tag" in feature.qualifiers:
                    locus_tag = feature.qualifiers["locus_tag"][0]
                else:
                    # Generate a locus tag if not present
                    locus_tag = f"OM_{record.id}"
                    
                if "translation" in feature.qualifiers:
                    amino_acid_seq = feature.qualifiers["translation"][0]
                else:
                    # Translate the feature sequence if translation not provided
                    try:
                        seq = feature.extract(record.seq)
                        amino_acid_seq = str(seq.translate(table=11))
                        logger.info(f"Generated translation for {locus_tag}")
                    except Exception as e:
                        logger.error(f"Failed to translate feature: {e}")
                
                # Extract coordinates for the gene
                if hasattr(feature, 'location'):
                    gene_coords = (int(feature.location.start), int(feature.location.end))
                
                # Extract additional information
                notes = []
                if "note" in feature.qualifiers:
                    notes = feature.qualifiers["note"]
                    
                cluster_info = next((note for note in notes if "cluster" in note.lower()), None)
                
                # Store all useful information about this reference
                reference_types[locus_tag] = {
                    "nucleotide_seq": str(record.seq),
                    "amino_acid_seq": amino_acid_seq,
                    "gene_seq": str(feature.extract(record.seq)) if hasattr(feature, 'location') else None,
                    "locus_tag": locus_tag,
                    "cluster_info": cluster_info,
                    "record_id": record.id,
                    "description": record.description,
                    "gene_coords": gene_coords,
                    "strand": feature.location.strand if hasattr(feature, 'location') else None
                }
                
                logger.info(f"Parsed reference type: {locus_tag} - {cluster_info}")
                
        except Exception as e:
            logger.error(f"Error parsing reference database: {e}")
            self.metadata["errors"].append(f"Reference DB parsing error: {str(e)}")
            
        # Validate the parsed data
        if not reference_types:
            logger.error("No valid reference types found in the database")
            self.metadata["errors"].append("No valid reference types found")
            
        self.metadata["reference_count"] = len(reference_types)
        return reference_types
    
    def _precompute_reference_alignments(self) -> None:
        """
        Precompute sequence alignments between reference sequences for faster comparisons.
        This creates a multiple sequence alignment of all reference protein sequences.
        """
        try:
            # Extract protein sequences from references
            aa_sequences = []
            for locus_tag, data in self.reference_types.items():
                if data.get("amino_acid_seq"):
                    aa_sequences.append(SeqRecord(
                        Seq(data["amino_acid_seq"]),
                        id=locus_tag,
                        description=""
                    ))
            
            if len(aa_sequences) <= 1:
                logger.warning("Not enough reference sequences for multiple alignment")
                return
                
            # Write sequences to a temporary file
            aa_fasta = self.temp_dir / "reference_aa.fasta"
            SeqIO.write(aa_sequences, aa_fasta, "fasta")
            
            # Run multiple sequence alignment using Clustal Omega if available
            try:
                aligned_fasta = self.temp_dir / "reference_aligned.fasta"
                
                # Skip Clustal Omega if not found during initialization
                if not self.CLUSTALO_PATH:
                    logger.info("Skipping Clustal Omega (not available). Using Biopython alignment instead.")
                    self._align_with_biopython(aa_sequences)
                    return
                
                cmd = [
                    self.CLUSTALO_PATH,
                    "-i", str(aa_fasta),
                    "-o", str(aligned_fasta),
                    "--force",
                    "--threads", str(self.PARALLEL_PROCESSES)
                ]
                
                result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120)
                
                if aligned_fasta.exists():
                    # Store the alignment in the reference cache
                    self.reference_cache["msa"] = AlignIO.read(aligned_fasta, "fasta")
                    logger.info(f"Precomputed multiple sequence alignment of {len(aa_sequences)} references")
                    
                    # Calculate a distance matrix for the alignment
                    self._calculate_distance_matrix()
                    
            except (subprocess.SubprocessError, FileNotFoundError) as e:
                logger.info(f"Using Biopython alignment as fallback (Clustal Omega error: {e}).")
                self._align_with_biopython(aa_sequences)
                
        except Exception as e:
            logger.error(f"Error in reference alignment preprocessing: {e}")
            
    def _align_with_biopython(self, sequences: List[SeqRecord]) -> None:
        """
        Fallback method to create a pairwise alignment matrix using Biopython.
        
        Args:
            sequences: List of protein sequence records
        """
        try:
            # Only process if we have sequences
            if len(sequences) <= 1:
                logger.info("Not enough sequences for alignment (need at least 2)")
                return
                
            # Try to use the newer Align module first
            try:
                from Bio import Align
                
                logger.info("Using Biopython's Align module for pairwise alignment")
                # Use a simple progressive alignment strategy
                aligner = Align.PairwiseAligner()
                aligner.mode = 'global'
                aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
                aligner.open_gap_score = -10
                aligner.extend_gap_score = -0.5
                
                # Calculate all pairwise alignments
                alignment_scores = {}
                    
                for i, seq1 in enumerate(sequences[:-1]):
                    for seq2 in sequences[i+1:]:
                        try:
                            alignments = aligner.align(str(seq1.seq), str(seq2.seq))
                            if alignments:
                                score = alignments[0].score
                                alignment_scores[(seq1.id, seq2.id)] = score
                                alignment_scores[(seq2.id, seq1.id)] = score
                        except Exception as align_err:
                            logger.warning(f"Could not align sequences {seq1.id} and {seq2.id}: {align_err}")
                
                self.reference_cache["pairwise_scores"] = alignment_scores
                logger.info(f"Computed {len(alignment_scores)//2} pairwise alignments for {len(sequences)} sequences")
                
            except (ImportError, AttributeError):
                # Fallback to older pairwise2 module
                logger.info("Using Biopython's legacy pairwise2 module for alignment")
                from Bio import pairwise2
                
                alignment_scores = {}
                for i, seq1 in enumerate(sequences[:-1]):
                    for seq2 in sequences[i+1:]:
                        try:
                            # Use a simpler and more reliable alignment algorithm
                            alignments = pairwise2.align.globalxx(str(seq1.seq), str(seq2.seq))
                            if alignments:
                                score = alignments[0].score
                                alignment_scores[(seq1.id, seq2.id)] = score
                                alignment_scores[(seq2.id, seq1.id)] = score
                        except Exception as align_err:
                            logger.warning(f"Could not align sequences {seq1.id} and {seq2.id}: {align_err}")
                
                self.reference_cache["pairwise_scores"] = alignment_scores
                logger.info(f"Computed {len(alignment_scores)//2} pairwise alignments using legacy method")
            
        except Exception as e:
            logger.warning(f"Biopython alignment failed: {e}. Will proceed without reference alignment preprocessing.")
            # Create an empty cache to prevent further errors
            self.reference_cache["pairwise_scores"] = {}
    
    def _calculate_distance_matrix(self) -> None:
        """
        Calculate a distance matrix from the multiple sequence alignment.
        """
        if "msa" not in self.reference_cache:
            return
            
        try:
            msa = self.reference_cache["msa"]
            n_seqs = len(msa)
            
            # Initialize distance matrix
            dist_matrix = np.zeros((n_seqs, n_seqs))
            
            # Calculate pairwise distances
            for i in range(n_seqs):
                for j in range(i+1, n_seqs):
                    seq1 = str(msa[i].seq)
                    seq2 = str(msa[j].seq)
                    
                    # Calculate p-distance (proportion of differing sites)
                    matches = sum(a == b for a, b in zip(seq1, seq2) if a != '-' and b != '-')
                    aligned_positions = sum(1 for a, b in zip(seq1, seq2) if a != '-' or b != '-')
                    
                    if aligned_positions > 0:
                        distance = 1.0 - (matches / aligned_positions)
                    else:
                        distance = 1.0
                        
                    dist_matrix[i, j] = distance
                    dist_matrix[j, i] = distance  # Matrix is symmetric
            
            # Store the distance matrix with sequence IDs
            self.reference_cache["distance_matrix"] = {
                "matrix": dist_matrix,
                "ids": [seq.id for seq in msa]
            }
            
            logger.info(f"Calculated distance matrix for {n_seqs} sequences")
            
        except Exception as e:
            logger.error(f"Error calculating distance matrix: {e}")
    
    def _find_closest_references(self, query_seq: str, top_n: int = 3) -> List[str]:
        """
        Find the closest reference sequences to the query sequence.
        Uses precomputed distance matrix or pairwise alignments if available.
        
        Args:
            query_seq: Query protein sequence
            top_n: Number of closest references to return
            
        Returns:
            List of locus tags of the closest references
        """
        # If we don't have precomputed data, compare directly
        if not self.reference_cache.get("msa") and not self.reference_cache.get("pairwise_scores"):
            # Compare directly using sequence identity
            identities = []
            for locus_tag, data in self.reference_types.items():
                if data.get("amino_acid_seq"):
                    identity = self._calculate_identity(query_seq, data["amino_acid_seq"])
                    identities.append((locus_tag, identity))
                    
            # Sort by identity (highest first) and return top N
            identities.sort(key=lambda x: x[1], reverse=True)
            return [locus_tag for locus_tag, _ in identities[:top_n]]
        
        # If we have a multiple sequence alignment, add the query to it
        if "msa" in self.reference_cache:
            try:
                from Bio.Align import MultipleSeqAlignment
                from Bio.SeqRecord import SeqRecord
                
                # Create a temporary copy of the MSA
                msa = self.reference_cache["msa"]
                temp_msa = MultipleSeqAlignment([])
                for record in msa:
                    temp_msa.append(record)
                
                # Add query to MSA using a profile alignment method
                # This is a simplified approach - would be better with a proper profile aligner
                #  will just use Biopython's pairwise alignment
                aligner = pairwise2.PairwiseAligner()
                aligner.mode = 'global'
                aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
                aligner.open_gap_score = -10
                aligner.extend_gap_score = -0.5
                
                # Find best match to align against
                distances = []
                for ref_locus, ref_data in self.reference_types.items():
                    if ref_data.get("amino_acid_seq"):
                        identity = self._calculate_identity(query_seq, ref_data["amino_acid_seq"])
                        distances.append((ref_locus, identity))
                
                if not distances:
                    return []
                    
                # Sort by identity (highest first)
                distances.sort(key=lambda x: x[1], reverse=True)
                best_match_locus = distances[0][0]
                
                # Find best match in MSA
                best_match_idx = None
                for i, record in enumerate(msa):
                    if record.id == best_match_locus:
                        best_match_idx = i
                        break
                        
                if best_match_idx is None:
                    return [locus for locus, _ in distances[:top_n]]
                
                # Calculate distances to all sequences in MSA
                distances = []
                for i, record in enumerate(msa):
                    # Extract sequence without gaps
                    ref_seq = str(record.seq).replace('-', '')
                    identity = self._calculate_identity(query_seq, ref_seq)
                    distances.append((record.id, identity))
                
                # Sort by identity and return top N
                distances.sort(key=lambda x: x[1], reverse=True)
                return [locus for locus, _ in distances[:top_n]]
                
            except Exception as e:
                logger.error(f"Error finding closest references using MSA: {e}")
        
        # Fallback to direct comparison
        identities = []
        for locus_tag, data in self.reference_types.items():
            if data.get("amino_acid_seq"):
                identity = self._calculate_identity(query_seq, data["amino_acid_seq"])
                identities.append((locus_tag, identity))
                
        # Sort by identity (highest first) and return top N
        identities.sort(key=lambda x: x[1], reverse=True)
        return [locus_tag for locus_tag, _ in identities[:top_n]]
    
    def run_blast_search(self) -> Optional[Path]:
        """
        Run optimized BLAST search to find ompA gene in the input genome.
        Uses parallel processing for efficiency with large genomes.
        
        Returns:
            Path to the BLAST results file or None if failed
        """
        self.metadata["analysis_steps"].append("run_blast_search")
        
        try:
            # Extract reference sequences to a FASTA file
            reference_fasta = self.temp_dir / "ompa_references.fasta"
            with open(reference_fasta, "w") as f:
                for locus_tag, data in self.reference_types.items():
                    f.write(f">{locus_tag}\n{data['gene_seq'] if data.get('gene_seq') else data['nucleotide_seq']}\n")
            
            # Create a BLAST database from the genome
            genome_db = self.temp_dir / "genome_db"
            cmd = [
                "makeblastdb",
                "-in", str(self.genome_path),
                "-dbtype", "nucl",
                "-out", str(genome_db)
            ]
            
            # Redirect output to /dev/null in quiet mode
            if QUIET_MODE:
                with open(os.devnull, 'w') as devnull:
                    subprocess.run(cmd, check=True, stdout=devnull, stderr=devnull)
            else:
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Run BLAST search with optimized parameters
            blast_output = self.temp_dir / "ompa_blast_results.tsv"
            cmd = [
                "blastn",
                "-query", str(reference_fasta),
                "-db", str(genome_db),
                "-out", str(blast_output),
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
                "-max_target_seqs", "5",
                "-evalue", str(self.BLAST_EVALUE),
                "-num_threads", str(self.PARALLEL_PROCESSES),
                "-task", "blastn",  # Use standard BLASTN for better sensitivity
                "-word_size", "11"  # Smaller word size for better sensitivity
            ]
            
            # Redirect output to /dev/null in quiet mode
            if QUIET_MODE:
                with open(os.devnull, 'w') as devnull:
                    subprocess.run(cmd, check=True, stdout=devnull, stderr=devnull)
            else:
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Check if results exist and have content
            if not blast_output.exists() or blast_output.stat().st_size == 0:
                logger.warning(f"BLAST search produced no results for {self.genome_name}")
                self.metadata["warnings"].append("BLAST search produced no results")
                return None
                
            logger.info(f"BLAST search completed. Results saved to {blast_output}")
            return blast_output
        
        except Exception as e:
            logger.error(f"Error running BLAST search: {e}")
            self.metadata["errors"].append(f"BLAST search error: {str(e)}")
            return None
    
    def _validate_blast_results(self, blast_results_file: Path) -> Tuple[Optional[pd.DataFrame], List[str]]:
        """
        Validate and filter BLAST results to ensure high-quality hits.
        
        Args:
            blast_results_file: Path to the BLAST results file
            
        Returns:
            Tuple of (filtered DataFrame or None, list of warnings)
        """
        warnings = []
        
        try:
            # Define columns for the BLAST output
            columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                      "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]
            
            # Read the blast results
            blast_df = pd.read_csv(blast_results_file, sep="\t", names=columns)
            
            if blast_df.empty:
                warnings.append("BLAST results file is empty")
                return None, warnings
            
            # Calculate coverage percentages
            blast_df['qcoverage'] = (blast_df['length'] / blast_df['qlen']) * 100
            
            # Filter for good hits based on coverage and identity
            good_hits = blast_df[(blast_df['qcoverage'] >= self.MIN_COVERAGE) & 
                                (blast_df['pident'] >= self.IDENT_MODERATE)]
            
            if good_hits.empty:
                filtered_hits = blast_df.sort_values(by=['bitscore', 'length'], ascending=[False, False])
                warnings.append(f"No high-quality matches found. Best hit: {filtered_hits.iloc[0]['qcoverage']:.1f}% coverage, {filtered_hits.iloc[0]['pident']:.1f}% identity")
                return filtered_hits, warnings
            
            return good_hits, warnings
            
        except Exception as e:
            warnings.append(f"Error validating BLAST results: {e}")
            return None, warnings
    
    def extract_ompa_sequence(self, blast_results_file: Path) -> Optional[Path]:
        """
        Extract the ompA gene sequence based on BLAST results with improved accuracy and robustness.
        Uses advanced bioinformatics approaches to ensure complete gene extraction.
        
        Args:
            blast_results_file: Path to the BLAST results file
            
        Returns:
            Path to the extracted ompA sequence file or None if failed
        """
        self.metadata["analysis_steps"].append("extract_ompa_sequence")
        
        try:
            # Validate BLAST results
            filtered_hits, warnings = self._validate_blast_results(blast_results_file)
            for warning in warnings:
                self.metadata["warnings"].append(warning)
                logger.warning(warning)
                
            if filtered_hits is None:
                logger.warning("No valid BLAST hits found for ompA sequence extraction")
                return None
            
            # Get the best hit
            best_hit = filtered_hits.sort_values(by=['bitscore', 'length', 'pident'], 
                                               ascending=[False, False, False]).iloc[0]
            
            # Load genome sequence
            genome_records = list(SeqIO.parse(self.genome_path, "fasta"))
            
            if not genome_records:
                logger.error(f"Could not parse genome file: {self.genome_path}")
                return None
                
            # Find the genome record with the matching ID
            genome_record = None
            for record in genome_records:
                if record.id == best_hit["sseqid"]:
                    genome_record = record
                    break
                    
            if genome_record is None:
                # If no exact match, use the first record
                genome_record = genome_records[0]
                logger.warning(f"Could not find contig {best_hit['sseqid']} in genome, using first record")
            
            # Determine strand and coordinates
            is_reverse = best_hit["sstart"] > best_hit["send"]
            start = min(best_hit["sstart"], best_hit["send"])
            end = max(best_hit["sstart"], best_hit["send"])
            
            # Look up reference sequence info for the best hit
            ref_id = best_hit["qseqid"]
            ref_data = self.reference_types.get(ref_id, {})
            
            # Get reference sequence length
            ref_seq_length = len(ref_data.get("gene_seq", "")) if ref_data.get("gene_seq") else best_hit["qlen"]
            
            # Calculate coverage percentage
            match_coverage = best_hit["length"] / ref_seq_length * 100
            
            # Determine if the match is fragmentary
            is_fragmentary = match_coverage < 95.0
            
            # Calculate padding based on match characteristics
            if is_fragmentary:
                logger.warning(f"Fragmentary match detected: {match_coverage:.1f}% coverage")
                
                # Calculate query relative start/end positions
                qstart_rel = best_hit["qstart"] / ref_seq_length
                qend_rel = best_hit["qend"] / ref_seq_length
                
                # Add more padding to the side(s) where we're missing sequence
                missing_start = qstart_rel > 0.05  # Missing >5% from start
                missing_end = qend_rel < 0.95  # Missing >5% from end
                
                padding_start = max(300, int(ref_seq_length * qstart_rel * 1.2)) if missing_start else self.PADDING
                padding_end = max(300, int(ref_seq_length * (1 - qend_rel) * 1.2)) if missing_end else self.PADDING
                
                # Apply strand-specific padding
                if is_reverse:
                    start = max(1, start - padding_end)
                    end = min(len(genome_record.seq), end + padding_start)
                else:
                    start = max(1, start - padding_start)
                    end = min(len(genome_record.seq), end + padding_end)
                    
                logger.info(f"Applied adaptive padding: start={padding_start}bp, end={padding_end}bp")
            else:
                # Standard padding for complete matches
                start = max(1, start - self.PADDING)
                end = min(len(genome_record.seq), end + self.PADDING)
            
            # Extract sequence with the specified coordinates
            extracted_seq = genome_record.seq[start-1:end]
            
            # If the BLAST hit is on the reverse strand, reverse complement
            if is_reverse:
                extracted_seq = extracted_seq.reverse_complement()
            
            # Create a detailed sequence record
            extracted_record = SeqRecord(
                seq=extracted_seq,
                id=f"{genome_record.id}_ompA",
                description=(
                    f"ompA gene extracted from {genome_record.id} positions {start}-{end} "
                    f"(coverage: {match_coverage:.1f}%, identity: {best_hit['pident']:.1f}%, "
                    f"strand: {'-' if is_reverse else '+'}, "
                    f"reference: {ref_id})"
                )
            )
            
            # Save both FASTA and GenBank format
            extracted_fasta = self.ompa_output_dir / "extracted_ompa.fasta"
            SeqIO.write(extracted_record, extracted_fasta, "fasta")
            
            # Attempt to save in GenBank format with annotation
            try:
                from Bio.SeqFeature import SeqFeature, FeatureLocation
                
                # Create a feature for the ompA gene
                feature = SeqFeature(
                    FeatureLocation(self.PADDING, len(extracted_seq) - self.PADDING),
                    type="gene",
                    strand=1,  # Always set to forward strand since we've already reverse-complemented if needed
                    qualifiers={
                        "gene": ["ompA"],
                        "locus_tag": [f"OMPA_{self.genome_name}"],
                        "note": [
                            f"Extracted from {genome_record.id} positions {start}-{end}",
                            f"Coverage: {match_coverage:.1f}%, Identity: {best_hit['pident']:.1f}%",
                            f"Reference: {ref_id}"
                        ]
                    }
                )
                
                # Add the feature to the record
                extracted_record.features = [feature]
                
                # Save GenBank file
                extracted_gbk = self.ompa_output_dir / "extracted_ompa.gbk"
                SeqIO.write(extracted_record, extracted_gbk, "genbank")
                
            except Exception as e:
                logger.warning(f"Could not save GenBank format: {e}")
            
            # Store metadata about the extraction
            self.metadata["extraction"] = {
                "start": int(start),
                "end": int(end),
                "strand": "-" if is_reverse else "+",
                "length": len(extracted_seq),
                "match_coverage": float(match_coverage),
                "identity": float(best_hit["pident"]),
                "reference": ref_id,
                "is_fragmentary": is_fragmentary
            }
            
            logger.info(
                f"Extracted ompA sequence from {genome_record.id} positions {start}-{end}, "
                f"length: {len(extracted_seq)}bp, coverage: {match_coverage:.1f}%, "
                f"identity: {best_hit['pident']:.1f}%"
            )
            
            return extracted_fasta
        
        except Exception as e:
            logger.error(f"Error extracting ompA sequence: {e}")
            self.metadata["errors"].append(f"Sequence extraction error: {str(e)}")
            return None

    def run_prokka_annotation(self, extracted_fasta: Path) -> Tuple[Optional[Path], Optional[Path]]:
        """
        Run Prokka annotation on the extracted ompA sequence with improved settings.
        Uses more precise gene calling parameters for better annotation of short sequences.
        
        Args:
            extracted_fasta: Path to the extracted ompA sequence
            
        Returns:
            Tuple of (GenBank file path, FAA file path) or (None, None) if failed
        """
        self.metadata["analysis_steps"].append("run_prokka_annotation")
        
        try:
            prokka_outdir = self.ompa_output_dir / "prokka_output"
            prokka_outdir.mkdir(exist_ok=True)
            
            # Prepare Prokka command with optimized parameters for gene finding
            cmd = [
                "docker", "run", "--rm",
                "-v", f"{self.ompa_output_dir.absolute()}:/data",
                "staphb/prokka:latest", "prokka",
                "--outdir", "/data/prokka_output",
                "--prefix", "OMPA",
                "--locustag", "OMPA",
                "--kingdom", "Bacteria",
                "--genus", "Erwinia",
                "--species", "amylovora",
                "--proteins", "/data/prokka_proteins.faa",  # Use reference proteins as a guide
                "--mincontiglen", "500",  # Allow smaller contigs
                "--compliant",  # NCBI compliant output
                "--force",
                "--cpus", str(self.PARALLEL_PROCESSES)
            ]
            
            # In quiet mode, add flags to silence Prokka's output
            if QUIET_MODE:
                # Add parameters to silence Prokka's output
                cmd.extend(["--quiet", "--norrna", "--notrna"])
                
            # Add the input file to the command
            cmd.append("/data/extracted_ompa.fasta")
            
            # Create a protein FASTA file with reference ompA sequences to guide annotation
            ref_proteins_file = self.ompa_output_dir / "prokka_proteins.faa"
            with open(ref_proteins_file, "w") as f:
                for locus_tag, data in self.reference_types.items():
                    if data.get("amino_acid_seq"):
                        f.write(f">{locus_tag}|ompA\n{data['amino_acid_seq']}\n")
            
            # Check if Docker is available
            try:
                subprocess_kwargs = {'check': True, 'timeout': 10}
                if QUIET_MODE:
                    subprocess_kwargs.update({'stdout': subprocess.PIPE, 'stderr': subprocess.PIPE})
                docker_check = subprocess.run(["docker", "version"], **subprocess_kwargs)
            except (subprocess.SubprocessError, FileNotFoundError):
                logger.warning("Docker not available. Using fallback annotation method.")
                return self._try_manual_annotation(extracted_fasta)
                
            # Run Prokka
            try:
                subprocess_kwargs = {'check': True, 'timeout': 300}
                if QUIET_MODE:
                    subprocess_kwargs.update({'stdout': subprocess.PIPE, 'stderr': subprocess.PIPE})
                result = subprocess.run(cmd, **subprocess_kwargs)
                
                # Check if output files exist
                gbk_file = prokka_outdir / "OMPA.gbk"
                faa_file = prokka_outdir / "OMPA.faa"
                
                if gbk_file.exists() and faa_file.exists():
                    logger.info(f"Prokka annotation completed. Output saved to {prokka_outdir}")
                    
                    # Copy files to the main output directory for easier access
                    final_gbk = self.ompa_output_dir / "OMPA.gbk"
                    final_faa = self.ompa_output_dir / "OMPA.faa"
                    
                    with open(gbk_file, "rb") as src, open(final_gbk, "wb") as dst:
                        dst.write(src.read())
                        
                    with open(faa_file, "rb") as src, open(final_faa, "wb") as dst:
                        dst.write(src.read())
                    
                    return final_gbk, final_faa
                else:
                    logger.warning(f"Prokka annotation failed to produce expected output files")
                    
                    # Fallback: Try a manual annotation approach if Prokka fails
                    return self._try_manual_annotation(extracted_fasta)
            
            except subprocess.TimeoutExpired:
                logger.warning("Prokka annotation timed out, using fallback method")
                return self._try_manual_annotation(extracted_fasta)
                
            except subprocess.CalledProcessError as e:
                logger.warning(f"Prokka annotation failed: {e}, using fallback method")
                return self._try_manual_annotation(extracted_fasta)
            
        except Exception as e:
            logger.error(f"Error running Prokka annotation: {e}")
            self.metadata["errors"].append(f"Annotation error: {str(e)}")
            return None, None
    
    def _try_manual_annotation(self, extracted_fasta: Path) -> Tuple[Optional[Path], Optional[Path]]:
        """
        Fallback method when Prokka fails - attempt to manually identify and extract the ompA gene.
        
        Args:
            extracted_fasta: Path to the extracted sequence region
            
        Returns:
            Tuple of (GenBank file path, FAA file path) or (None, None) if failed
        """
        try:
            logger.info("Attempting manual annotation as fallback")
            
            # Load the extracted sequence
            extracted_record = next(SeqIO.parse(extracted_fasta, "fasta"))
            
            # Find ORFs that might be ompA
            orfs = self._find_orfs(str(extracted_record.seq))
            
            if not orfs:
                logger.warning("No suitable ORFs found in the extracted sequence")
                return None, None
                
            # Score ORFs based on similarity to reference sequences
            scored_orfs = []
            for start, end, frame in orfs:
                # Extract and translate the ORF
                if frame > 0:
                    orf_seq = extracted_record.seq[start:end]
                else:
                    orf_seq = extracted_record.seq[start:end].reverse_complement()
                    
                aa_seq = str(orf_seq.translate(table=11))
                
                # Compare with reference sequences
                best_score = 0
                best_ref = None
                
                for locus_tag, data in self.reference_types.items():
                    if data.get("amino_acid_seq"):
                        identity = self._calculate_identity(aa_seq, data["amino_acid_seq"])
                        if identity > best_score:
                            best_score = identity
                            best_ref = locus_tag
                
                scored_orfs.append((start, end, frame, aa_seq, best_score, best_ref))
            
            # Sort ORFs by score (highest first)
            scored_orfs.sort(key=lambda x: x[4], reverse=True)
            
            if not scored_orfs or scored_orfs[0][4] < 75.0:
                logger.warning(f"No ORFs with good match to ompA references (best: {scored_orfs[0][4] if scored_orfs else 0}%)")
                return None, None
                
            # Use the best scoring ORF
            start, end, frame, aa_seq, score, ref_id = scored_orfs[0]
            
            # Create GenBank and FAA files
            manual_gbk = self.ompa_output_dir / "OMPA.gbk"
            manual_faa = self.ompa_output_dir / "OMPA.faa"
            
            # Create a new record for the GenBank file
            from Bio.SeqFeature import SeqFeature, FeatureLocation
            
            gbk_record = SeqRecord(
                seq=extracted_record.seq,
                id=f"{self.genome_name}_ompA",
                name=f"{self.genome_name}_ompA",
                description=f"Manually annotated ompA gene from {self.genome_name}"
            )
            
            # Add the ompA feature
            feature = SeqFeature(
                FeatureLocation(start, end),
                type="CDS",
                strand=(1 if frame > 0 else -1),
                qualifiers={
                    "translation": [aa_seq],
                    "product": ["Outer membrane protein A"],
                    "gene": ["ompA"],
                    "locus_tag": [f"OMPA_{self.genome_name}"],
                    "note": [
                        f"Manually annotated ompA gene",
                        f"Best match: {ref_id} ({score:.1f}% identity)"
                    ]
                }
            )
            
            gbk_record.features = [feature]
            
            # Write the GenBank file
            SeqIO.write(gbk_record, manual_gbk, "genbank")
            
            # Create the FAA file
            faa_record = SeqRecord(
                seq=Seq(aa_seq),
                id=f"OMPA_{self.genome_name}",
                name=f"OMPA_{self.genome_name}",
                description=f"Outer membrane protein A [Erwinia amylovora {self.genome_name}]"
            )
            
            SeqIO.write(faa_record, manual_faa, "fasta")
            
            logger.info(f"Manual annotation completed. Found ompA gene with {score:.1f}% identity to {ref_id}")
            return manual_gbk, manual_faa
            
        except Exception as e:
            logger.error(f"Error in manual annotation: {e}")
            return None, None
    
    def _find_orfs(self, sequence: str, min_length: int = 300) -> List[Tuple[int, int, int]]:
        """
        Find open reading frames in a DNA sequence.
        
        Args:
            sequence: DNA sequence to search
            min_length: Minimum ORF length in nucleotides
            
        Returns:
            List of tuples (start, end, frame) for each ORF
        """
        # Standard genetic code
        start_codon = "ATG"
        stop_codons = {"TAA", "TAG", "TGA"}
        
        orfs = []
        seq_length = len(sequence)
        
        # Check all six reading frames
        for frame in range(3):
            # Forward strand
            i = frame
            while i < seq_length - 2:
                # Look for start codon
                if sequence[i:i+3] == start_codon:
                    j = i + 3
                    # Look for stop codon
                    while j < seq_length - 2:
                        codon = sequence[j:j+3]
                        if codon in stop_codons:
                            # Found a complete ORF
                            if j - i >= min_length:
                                orfs.append((i, j+3, 1))  # +1 for forward frame
                            break
                        j += 3
                    i = j + 3 if j < seq_length - 2 else seq_length
                else:
                    i += 3
                    
            # Reverse strand
            rev_comp = str(Seq(sequence).reverse_complement())
            i = frame
            while i < seq_length - 2:
                # Look for start codon
                if rev_comp[i:i+3] == start_codon:
                    j = i + 3
                    # Look for stop codon
                    while j < seq_length - 2:
                        codon = rev_comp[j:j+3]
                        if codon in stop_codons:
                            # Found a complete ORF
                            if j - i >= min_length:
                                # Convert to original strand coordinates
                                rev_start = seq_length - (j + 3)
                                rev_end = seq_length - i
                                orfs.append((rev_start, rev_end, -1))  # -1 for reverse frame
                            break
                        j += 3
                    i = j + 3 if j < seq_length - 2 else seq_length
                else:
                    i += 3
        
        return orfs
    
    def compare_amino_acid_sequences(self, faa_file: Path) -> Dict[str, Any]:
        """
        Advanced method to compare the amino acid sequence against reference types.
        Uses multiple alignment strategies and sensitive comparison methods.
        
        Args:
            faa_file: Path to the FAA file with proteins
            
        Returns:
            Dictionary with detailed typing results
        """
        self.metadata["analysis_steps"].append("compare_amino_acid_sequences")
        
        try:
            # Extract amino acid sequence from the FAA file
            prokka_proteins = list(SeqIO.parse(faa_file, "fasta"))
            
            if not prokka_proteins:
                logger.warning(f"No proteins found in {faa_file}")
                return {
                    "ompa_type": "Unknown",
                    "match_quality": "No proteins found",
                    "identity": 0,
                    "best_match": None,
                    "all_matches": [],
                    "confidence": 0.0
                }
            
            # Find the ompA protein using multiple identification strategies
            ompa_protein = None
            
            # Strategy 1: Look for ompA in description
            for protein in prokka_proteins:
                desc_lower = protein.description.lower()
                if ("ompa" in desc_lower or 
                    "outer membrane protein a" in desc_lower or 
                    "outer membrane protein porin" in desc_lower):
                    ompa_protein = protein
                    logger.info(f"Found ompA protein by description: {protein.description}")
                    break
            
            # Strategy 2: Look for the protein with the right length
            if ompa_protein is None:
                expected_length = 355  # Typical length of ompA
                length_tolerance = 20  # Allow some variation
                
                for protein in prokka_proteins:
                    protein_length = len(protein.seq)
                    if abs(protein_length - expected_length) <= length_tolerance:
                        ompa_protein = protein
                        logger.info(f"Found potential ompA protein by length: {protein_length}aa")
                        break
            
            # Strategy 3: Compare all proteins to reference sequences
            if ompa_protein is None and len(prokka_proteins) > 0:
                best_match_score = 0
                for protein in prokka_proteins:
                    # Skip very short proteins
                    if len(protein.seq) < 200:
                        continue
                        
                    # Compare with reference sequences
                    for locus_tag, data in self.reference_types.items():
                        if data.get("amino_acid_seq"):
                            identity = self._calculate_identity(str(protein.seq), data["amino_acid_seq"])
                            if identity > best_match_score:
                                best_match_score = identity
                                ompa_protein = protein
                
                if ompa_protein:
                    logger.info(f"Found potential ompA protein by sequence comparison: {best_match_score:.1f}% match")
            
            if ompa_protein is None:
                logger.warning("No ompA protein found in the annotation")
                return {
                    "ompa_type": "Unknown",
                    "match_quality": "ompA protein not found in annotation",
                    "identity": 0,
                    "best_match": None,
                    "all_matches": [],
                    "confidence": 0.0
                }
            
            # Compare with all reference types using multiple methods
            matches = []
            query_seq = str(ompa_protein.seq)
            
            # Method 1: Direct sequence identity
            for locus_tag, data in self.reference_types.items():
                if not data.get("amino_acid_seq"):
                    continue
                
                ref_seq = data["amino_acid_seq"]
                identity = self._calculate_identity(query_seq, ref_seq)
                
                # Calculate additional metrics for very similar sequences
                similarity = self._calculate_similarity(query_seq, ref_seq)
                
                matches.append({
                    "type_id": locus_tag,
                    "identity": identity,
                    "similarity": similarity,
                    "cluster_info": data.get("cluster_info"),
                    "length_diff": abs(len(query_seq) - len(ref_seq))
                })
            
            # Sort matches by identity
            matches.sort(key=lambda x: (x["identity"], -x["length_diff"]), reverse=True)
            
            # If no matches, return unknown
            if not matches:
                return {
                    "ompa_type": "Novel",
                    "match_quality": "No reference matches found",
                    "identity": 0,
                    "best_match": None,
                    "all_matches": [],
                    "confidence": 0.0
                }
            
            # Get the best match
            best_match = matches[0]
            
            # Calculate confidence score
            # Higher when the best match is significantly better than the second best
            confidence = 1.0
            if len(matches) > 1:
                # Calculate drop in identity to second-best match
                identity_drop = best_match["identity"] - matches[1]["identity"]
                # Normalize to a 0-1 scale (larger drops give higher confidence)
                confidence = min(1.0, identity_drop / 1.0)  # 1% difference gives full confidence
            
            # Determine match quality with adjusted thresholds for very similar sequences
            match_quality = "Perfect" if best_match["identity"] >= self.IDENT_PERFECT else \
                           "Very High" if best_match["identity"] >= self.IDENT_VERY_HIGH else \
                           "High" if best_match["identity"] >= self.IDENT_HIGH else \
                           "Moderate" if best_match["identity"] >= self.IDENT_MODERATE else \
                           "Low"
            
            # Find distinguishing mutations if sequences are very similar
            mutations = []
            if best_match["identity"] > 95.0:
                ref_id = best_match["type_id"]
                ref_seq = self.reference_types[ref_id]["amino_acid_seq"]
                mutations = self._find_mutations(query_seq, ref_seq)
            
            # Return comprehensive results
            result = {
                "ompa_type": best_match["type_id"],
                "match_quality": match_quality,
                "identity": best_match["identity"],
                "similarity": best_match["similarity"],
                "confidence": confidence,
                "best_match": best_match,
                "all_matches": matches[:5],  # Return top 5 matches
                "mutations": mutations,
                "query_length": len(query_seq)
            }
            
            logger.info(
                f"OmpA typing result: {result['ompa_type']} "
                f"({result['match_quality']}, {result['identity']:.2f}% identity, "
                f"{confidence:.2f} confidence)"
            )
            
            if mutations:
                logger.info(f"Distinguishing mutations: {', '.join(mutations[:5])}")
                if len(mutations) > 5:
                    logger.info(f"...and {len(mutations) - 5} more mutations")
            
            return result
            
        except Exception as e:
            logger.error(f"Error comparing amino acid sequences: {e}")
            self.metadata["errors"].append(f"Sequence comparison error: {str(e)}")
            return {
                "ompa_type": "Error",
                "match_quality": f"Error: {str(e)}",
                "identity": 0,
                "best_match": None,
                "all_matches": [],
                "confidence": 0.0
            }
    
    def _calculate_identity(self, seq1: str, seq2: str) -> float:
        """
        Calculate the percentage identity between two sequences using optimized methods.
        
        Args:
            seq1: First amino acid sequence
            seq2: Second amino acid sequence
            
        Returns:
            Percentage identity (0-100)
        """
        if not seq1 or not seq2:
            return 0
        
        try:
            # For very similar sequences, use Smith-Waterman alignment for accurate comparison
            # This is the most accurate method for calculating identity
            from Bio import Align
            
            aligner = Align.PairwiseAligner()
            aligner.mode = 'global'
            aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
            aligner.open_gap_score = -10
            aligner.extend_gap_score = -0.5
            
            # Get the best alignment
            alignment = aligner.align(seq1, seq2)[0]
            
            # Calculate identity over aligned positions
            matches = 0
            aligned_positions = 0
            
            for a, b in zip(alignment.aligned[0], alignment.aligned[1]):
                for i in range(a[0], a[1]):
                    if seq1[i] == seq2[b[0] + i - a[0]]:
                        matches += 1
                    aligned_positions += 1
            
            if aligned_positions > 0:
                return (matches / aligned_positions) * 100
                
        except (ImportError, IndexError) as e:
            logger.warning(f"Advanced alignment failed: {e}. Using fallback method.")
        
        # Fallback method: Simple sequence comparison
        min_len = min(len(seq1), len(seq2))
        matches = sum(a == b for a, b in zip(seq1[:min_len], seq2[:min_len]))
        
        # Account for length differences
        length_diff = abs(len(seq1) - len(seq2))
        effective_len = min_len + length_diff  # Count length differences as mismatches
        
        return (matches / effective_len) * 100
    
    def _calculate_similarity(self, seq1: str, seq2: str) -> float:
        """
        Calculate amino acid similarity (accounting for conservative substitutions).
        
        Args:
            seq1: First amino acid sequence
            seq2: Second amino acid sequence
            
        Returns:
            Percentage similarity (0-100)
        """
        if not seq1 or not seq2:
            return 0
            
        try:
            from Bio.SubsMat import MatrixInfo
            
            # Use BLOSUM62 matrix to score similarity
            matrix = MatrixInfo.blosum62
            
            # Align sequences first
            min_len = min(len(seq1), len(seq2))
            similar_positions = 0
            
            for i in range(min_len):
                a, b = seq1[i], seq2[i]
                if a == b:
                    similar_positions += 1
                else:
                    # Check for positive score in BLOSUM matrix (similar amino acids)
                    pair = (a, b) if (a, b) in matrix else (b, a)
                    if pair in matrix and matrix[pair] > 0:
                        similar_positions += 1
            
            # Account for length differences
            length_diff = abs(len(seq1) - len(seq2))
            effective_len = min_len + length_diff
            
            return (similar_positions / effective_len) * 100
            
        except ImportError:
            # If Bio.SubsMat not available, use identity as fallback
            return self._calculate_identity(seq1, seq2)
    
    def _find_mutations(self, query_seq: str, ref_seq: str) -> List[str]:
        """
        Find amino acid mutations between query and reference sequence.
        
        Args:
            query_seq: Query protein sequence
            ref_seq: Reference protein sequence
            
        Returns:
            List of mutation strings in the format "A123B" (ref_aa, position, query_aa)
        """
        mutations = []
        
        try:
            # First align the sequences if necessary
            from Bio import pairwise2
            
            alignments = pairwise2.align.globalxx(ref_seq, query_seq)
            if not alignments:
                return []
                
            alignment = alignments[0]
            ref_aligned, query_aligned = alignment[0], alignment[1]
            
            # Find differences, accounting for gaps
            ref_pos = 0
            for i in range(len(ref_aligned)):
                if ref_aligned[i] != '-':
                    ref_pos += 1
                    
                if ref_aligned[i] != query_aligned[i] and ref_aligned[i] != '-' and query_aligned[i] != '-':
                    mutations.append(f"{ref_aligned[i]}{ref_pos}{query_aligned[i]}")
            
            return mutations
            
        except ImportError:
            # Fallback to simple comparison for equal length sequences
            if len(query_seq) != len(ref_seq):
                return []
                
            for i in range(len(ref_seq)):
                if ref_seq[i] != query_seq[i]:
                    mutations.append(f"{ref_seq[i]}{i+1}{query_seq[i]}")
                    
            return mutations
    
    def generate_typing_report(self, typing_result: Dict[str, Any]) -> str:
        """
        Generate a simplified report string for the ompA typing result.
        Only shows type and match quality for CSV display purposes.
        
        Args:
            typing_result: Dictionary with typing results
            
        Returns:
            Formatted typing report string
        """
        if typing_result["ompa_type"] == "Unknown" or typing_result["ompa_type"] == "Error":
            return f"(Unknown)"
        
        if typing_result["ompa_type"] == "Novel":
            return f"(Novel)"
        
        # Build a simplified report - just type and quality
        result_str = f"{typing_result['ompa_type']} ({typing_result['match_quality']})"
        
        # Save detailed information to metadata for debugging
        detailed_info = {
            "type": typing_result['ompa_type'],
            "match_quality": typing_result['match_quality'],
            "identity": f"{typing_result['identity']:.2f}%",
            "confidence": typing_result.get('confidence', 0)
        }
        
        if typing_result.get('mutations'):
            detailed_info["mutations"] = typing_result['mutations']
            
        if typing_result.get("all_matches") and len(typing_result["all_matches"]) > 1:
            detailed_info["other_matches"] = [
                {"type": m['type_id'], "identity": m['identity']} 
                for m in typing_result["all_matches"][1:3]
            ]
            
        self.metadata["detailed_result"] = detailed_info
        logger.debug(f"Detailed OmpA typing result for {self.genome_name}: {detailed_info}")
        
        return result_str
    
    def cleanup(self):
        """Clean up temporary files and resources."""
        try:
            if hasattr(self, 'temp_dir') and self.temp_dir.exists():
                import shutil
                shutil.rmtree(self.temp_dir)
                logger.info(f"Removed temporary directory: {self.temp_dir}")
                
            # Clear caches to free memory
            if hasattr(self, 'reference_cache'):
                self.reference_cache.clear()
                
        except Exception as e:
            logger.warning(f"Error cleaning up temporary files: {e}")
    
    def genotype_ompa(self) -> Tuple[str, str]:
        """
        Run the full ompA genotyping pipeline with robust error handling and fallbacks.
        
        Returns:
            Tuple of (typing report string, ompA type)
        """
        import time
        start_time = time.time()
        
        try:
            logger.info(f"Starting ompA genotyping pipeline for {self.genome_name}")
            
            # Step 1: Run BLAST search
            self.metadata["timing"]["blast_start"] = time.time()
            blast_results = self.run_blast_search()
            self.metadata["timing"]["blast_end"] = time.time()
            
            if not blast_results:
                return "(Unknown) - BLAST search failed", "Unknown"
            
            # Step 2: Extract ompA sequence
            self.metadata["timing"]["extract_start"] = time.time()
            extracted_fasta = self.extract_ompa_sequence(blast_results)
            self.metadata["timing"]["extract_end"] = time.time()
            
            if not extracted_fasta:
                return "(Unknown) - Failed to extract ompA sequence", "Unknown"
            
            # Step 3: Run Prokka annotation
            self.metadata["timing"]["prokka_start"] = time.time()
            gbk_file, faa_file = self.run_prokka_annotation(extracted_fasta)
            self.metadata["timing"]["prokka_end"] = time.time()
            
            if not gbk_file or not faa_file:
                return "(Unknown) - Gene annotation failed", "Unknown"
            
            # Step 4: Compare amino acid sequences
            self.metadata["timing"]["compare_start"] = time.time()
            typing_result = self.compare_amino_acid_sequences(faa_file)
            self.metadata["timing"]["compare_end"] = time.time()
            
            # Step 5: Generate typing report
            typing_report = self.generate_typing_report(typing_result)
            
            # Add timing information
            self.metadata["timing"]["total"] = time.time() - start_time
            
            # Write metadata to a JSON file for debugging and reproducibility
            self._write_metadata()
            
            return typing_report, typing_result["ompa_type"]
            
        except Exception as e:
            logger.error(f"Error in ompA genotyping pipeline: {e}")
            self.metadata["errors"].append(f"Pipeline error: {str(e)}")
            self._write_metadata()  # Write metadata even on error
            return f"(Error) - {str(e)}", "Error"
        finally:
            self.metadata["timing"]["cleanup_start"] = time.time()
            self.cleanup()
            self.metadata["timing"]["cleanup_end"] = time.time()
    
    def _write_metadata(self) -> None:
        """Write metadata to a JSON file for debugging and reproducibility."""
        try:
            import json
            metadata_file = self.ompa_output_dir / "metadata.json"
            with open(metadata_file, "w") as f:
                json.dump(self.metadata, f, indent=2, default=str)
        except Exception as e:
            logger.warning(f"Error writing metadata: {e}")

# Function to run ompA genotyping on a single genome
def genotype_ompa_for_genome(genome_path: Path, reference_db: Path, output_dir: Path) -> Tuple[str, str]:
    """
    Run ompA genotyping for a single genome with robust error handling.
    
    Args:
        genome_path: Path to the genome FASTA file
        reference_db: Path to the reference ompA database GenBank file
        output_dir: Path to the output directory
        
    Returns:
        Tuple of (typing report string, ompA type)
    """
    try:
        genotyper = OmpAGenotyper(genome_path, reference_db, output_dir)
        return genotyper.genotype_ompa()
    except Exception as e:
        logger.error(f"Error in genotype_ompa_for_genome: {e}")
        return f"(Error) - {str(e)}", "Error"

# Function to process multiple genomes in parallel
def process_genomes_in_folder(genomes_folder: Path, reference_db: Path, 
                             output_dir: Optional[Path] = None,
                             num_processes: int = None) -> List[Dict[str, Any]]:
    """
    Process all genomes in a folder for ompA genotyping with parallel execution.
    
    Args:
        genomes_folder: Path to the folder containing genome FASTA files
        reference_db: Path to the reference ompA database GenBank file
        output_dir: Path to the output directory (if None, uses genomes_folder/types_finder)
        num_processes: Number of processes to use (defaults to CPU count - 1)
        
    Returns:
        List of dictionaries with typing results for each genome
    """
    results = []
    output_dir = output_dir or (genomes_folder / 'types_finder')
    num_processes = num_processes or max(1, multiprocessing.cpu_count() - 1)
    
    try:
        logger.info(f"Processing genomes in folder: {genomes_folder} using {num_processes} processes")
        genomes_folder = Path(genomes_folder)
        
        # Gather all genome files first
        genome_files = []
        genome_subdirs = []
        
        for item in genomes_folder.iterdir():
            if item.is_dir():
                for fasta_file in item.glob("*.fasta"):
                    genome_files.append(fasta_file)
                    genome_subdirs.append(item)
                    break  # Just take the first FASTA file in each directory
        
        logger.info(f"Found {len(genome_files)} genome files to process")
        
        if not genome_files:
            logger.warning(f"No genome files found in {genomes_folder}")
            return []
        
        # Function to process a single genome
        def process_genome(args):
            genome_path, genome_subdir = args
            genome_name = genome_subdir.name
            
            try:
                # Create output directory
                genome_output_dir = output_dir / genome_name
                genome_output_dir.mkdir(parents=True, exist_ok=True)
                
                logger.info(f"Processing genome: {genome_name}")
                typing_report, ompa_type = genotype_ompa_for_genome(genome_path, reference_db, genome_output_dir)
                
                return {
                    "genome_name": genome_name,
                    "typing_report": typing_report,
                    "ompa_type": ompa_type,
                    "output_dir": str(genome_output_dir),
                    "success": True
                }
            except Exception as e:
                logger.error(f"Error processing genome {genome_name}: {e}")
                return {
                    "genome_name": genome_name,
                    "typing_report": f"(Error) - {str(e)}",
                    "ompa_type": "Error",
                    "output_dir": str(output_dir / genome_name),
                    "success": False,
                    "error": str(e)
                }
        
        # Process genomes in parallel
        with ThreadPoolExecutor(max_workers=num_processes) as executor:
            future_to_genome = {
                executor.submit(process_genome, (genome_file, genome_subdir)): (genome_file, genome_subdir)
                for genome_file, genome_subdir in zip(genome_files, genome_subdirs)
            }
            
            for future in future_to_genome:
                try:
                    result = future.result()
                    results.append(result)
                    genome_file, _ = future_to_genome[future]
                    logger.info(f"Completed processing {genome_file.name}: {result['typing_report']}")
                except Exception as e:
                    genome_file, genome_subdir = future_to_genome[future]
                    logger.error(f"Failed to process {genome_file.name}: {e}")
                    results.append({
                        "genome_name": genome_subdir.name,
                        "typing_report": f"(Error) - {str(e)}",
                        "ompa_type": "Error",
                        "output_dir": str(output_dir / genome_subdir.name),
                        "success": False,
                        "error": str(e)
                    })
        
        # Summarize results
        success_count = sum(1 for result in results if result.get("success", False))
        logger.info(f"Successfully processed {success_count} out of {len(results)} genomes")
        
        # Write overall results to CSV
        summary_csv = output_dir / "ompa_typing_results.csv"
        try:
            import csv
            with open(summary_csv, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Genome', 'OmpA Type', 'Typing Report', 'Success'])
                for result in results:
                    writer.writerow([
                        result['genome_name'],
                        result['ompa_type'],
                        result['typing_report'],
                        result.get('success', False)
                    ])
            logger.info(f"Wrote summary results to {summary_csv}")
        except Exception as e:
            logger.error(f"Error writing summary CSV: {e}")
        
        return results
        
    except Exception as e:
        logger.error(f"Error in process_genomes_in_folder: {e}")
        return results

if __name__ == "__main__":
    # Example usage
    reference_db = Path(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "reference_types_database/ompa/types_ompa.gb"))
    genome_path = Path(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "test-data/GCF_012367655.1_ASM1236765v1_genomic_Ea01-03.fasta"))
    output_dir = Path(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "test-data/output"))
    
    # Configure logging
    logger.info("Starting OmpA genotyping example")
    
    # Run genotyping
    genotyper = OmpAGenotyper(genome_path, reference_db, output_dir)
    typing_report, ompa_type = genotyper.genotype_ompa()
    
    print(f"\nOmpA Typing Result: {typing_report}")
    print(f"OmpA Type: {ompa_type}")