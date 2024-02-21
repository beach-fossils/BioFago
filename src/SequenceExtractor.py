from Bio import SeqIO
import pandas as pd
import os
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class SequenceExtractor:
    def __init__(self, fasta_folder, csv_file, output_folder):
        self.fasta_folder = fasta_folder
        self.csv_file = csv_file
        self.output_folder = output_folder
        self.fasta_files = self._load_fasta_files()

    def _load_fasta_files(self):
        fasta_files = {}
        for filename in os.listdir(self.fasta_folder):
            if filename.endswith('.fasta'):
                path = os.path.join(self.fasta_folder, filename)
                file_id = filename.replace('_genomic.fasta', '')  # Remove '_genomic.fasta' to match the CSV identifier
                fasta_files[file_id] = SeqIO.to_dict(SeqIO.parse(path, "fasta"))
                logging.info(f"Loaded {len(fasta_files[file_id])} records from {filename}")
        return fasta_files

    def extract_sequences(self):
        df = pd.read_csv(self.csv_file)
        for file_name, group in df.groupby('file_name'):
            output_path = os.path.join(self.output_folder, f"{file_name}.fasta")
            with open(output_path, 'w') as output_fasta:
                for index, row in group.iterrows():
                    seq_id = row['qseqid']
                    start = int(row['qstart']) - 1  # Convert to 0-based index
                    end = int(row['qend'])
                    # The seq_id must match an ID within the FASTA file
                    for fasta_id, fasta_seqs in self.fasta_files.items():
                        if seq_id in fasta_seqs:
                            sequence = fasta_seqs[seq_id].seq[start:end]
                            output_fasta.write(f">{seq_id}\n{sequence}\n")
                            logging.info(f"Wrote sequence {seq_id} from {file_name} to {output_path}")
                            break
                    else:
                        logging.warning(f"Sequence {seq_id} not found in any FASTA files for {file_name}.")


def main():
    fasta_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/genomes_erwinia_amylovora'
    csv_file = '/Users/josediogomoura/Documents/BioFago/BioFago/data/BLAST/LPS_locus_13genes_FN434113/results/csv_folder/compiled_results.csv'
    output_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/extracted_sequences_erw_amy/LPS_locus_13genes_FN434113'

    # Create an instance of the class
    extractor = SequenceExtractor(fasta_folder, csv_file, output_folder)

    # Extract the sequences
    extractor.extract_sequences()


if __name__ == "__main__":
    main()


