import pandas as pd
from pathlib import Path
import logging
from typing import List, Tuple, Dict, Any
from collections import Counter

# Set up logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class FastaStatistics:
    def __init__(self, file_path: Path) -> None:
        self.file_path = Path(file_path)
        self.contigs: List[Tuple[str, str]] = []
        self._read_fasta()

    def _read_fasta(self) -> None:
        """Reads a FASTA file and stores contigs in a list along with their names."""
        try:
            with self.file_path.open('r') as file:
                contig_name = ''
                contig_sequence = []
                for line in file:
                    line = line.strip()
                    if line.startswith('>'):
                        if contig_sequence:
                            self.contigs.append((contig_name, ''.join(contig_sequence)))
                            contig_sequence = []
                        contig_name = line[1:]  # Remove '>' from the name
                    else:
                        contig_sequence.append(line)
                if contig_sequence:  # Add the last contig if exists
                    self.contigs.append((contig_name, ''.join(contig_sequence)))
        except IOError:
            logger.error(f"Error: The file {self.file_path} could not be opened.")

    @property
    def number_of_contigs(self) -> int:
        """Returns the number of contigs in the FASTA file."""
        return len(self.contigs)

    @property
    def largest_contig(self) -> Tuple[int, str]:
        """Returns the length of the largest contig and its name."""
        if not self.contigs:
            return (0, None)
        return max((len(seq), name) for name, seq in self.contigs)

    @property
    def total_length_of_contigs(self) -> int:
        """Returns the total length of all contigs."""
        return sum(len(seq) for _, seq in self.contigs)

    @property
    def average_contig_length(self) -> float:
        """Returns the average length of contigs."""
        return round(self.total_length_of_contigs / self.number_of_contigs, 2) if self.number_of_contigs else 0

    @property
    def gc_content(self) -> float:
        """Calculates the overall GC content of the contigs."""
        gc_count = sum(seq.count('G') + seq.count('C') for _, seq in self.contigs)
        return round(gc_count / self.total_length_of_contigs * 100, 2) if self.total_length_of_contigs > 0 else 0

    @property
    def n50(self) -> int:
        """Calculates the N50 value of the contigs."""
        lengths = sorted((len(seq) for _, seq in self.contigs), reverse=True)
        total_length = sum(lengths)
        cumulative_sum = 0
        for length in lengths:
            cumulative_sum += length
            if cumulative_sum >= total_length / 2:
                return length
        return 0

    @property
    def non_standard_bases_count(self) -> int:
        """Counts non-standard nucleotide bases (anything other than A, C, G, T) in all contigs."""
        standard_bases = set('ACGT')
        return sum(1 for _, seq in self.contigs for base in seq.upper() if base not in standard_bases)

    def generate_assembly_statistics(self) -> dict:
        """Generates a dictionary with assembly statistics."""
        largest_contig_length, largest_contig_name = self.largest_contig
        
        # We need to preserve the full filename except for the very last extension
        # This is critical for files like GCF_002732285.1_GCF_002732285.1_ASM273228v1_genomic.fna
        # where all components need to be preserved
        filename = self.file_path.name
        name_without_ext = filename.rsplit('.', 1)[0] if '.' in filename else filename
        
        # Add extra logging to debug filename issues
        logging.debug(f"Original filename: {self.file_path}")
        logging.debug(f"Name without extension: {name_without_ext}")
        
        return {
            'name': name_without_ext,
            'contig_count': self.number_of_contigs,
            'N50_value': self.n50,
            'largest_contig': largest_contig_name,
            'largest_contig_size_bp': largest_contig_length,
            'total_size_bp': self.total_length_of_contigs,
            'ambiguous_bases': self.non_standard_bases_count,
            'GC_content_percent': self.gc_content
        }

    def output_to_csv(self, output_path: Path) -> None:
        """Outputs the assembly statistics to a CSV file."""
        stats = self.generate_assembly_statistics()
        df = pd.DataFrame([stats])
        try:
            df.to_csv(output_path, index=False)
            logger.info(f"Statistics saved to {output_path}")
        except IOError:
            logger.error(f"Error: Could not write to file {output_path}.")



if __name__ == '__main__':
    file_path = Path('/Users/josediogomoura/Documents/BioFago/BioFago/data/recent_outputs/stats/genomes/PRR_75_modified.fasta')
    output = '/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/recent_outputs/stats/output/stats.csv'
    main(file_path, output)

