from Bio import SeqIO
import os
import logging

# set logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


def gb_to_fasta(input_path, output_path):
    for subdir, dirs, files in os.walk(input_path):
        for file in files:
            if file.endswith((".gbff", ".gbk", ".gb")):
                gb_file_path = os.path.join(subdir, file)

                # Extract the strain name from the path (name of the folder containing the .gbff file)
                strain_name = os.path.basename(subdir)

                # Adjust the output path to use the strain name for the fasta file, preserving a common structure
                fasta_file_name = strain_name + ".fasta"
                fasta_file_path = os.path.join(output_path, fasta_file_name)

                # Extract all SeqRecord objects from genbank file
                records = list(SeqIO.parse(gb_file_path, "genbank"))

                # Write all records to fasta file once
                with open(fasta_file_path, "w") as output_handle:
                    SeqIO.write(records, output_handle, "fasta")
                logging.info(f"Converted {gb_file_path} to {fasta_file_path}")


def main():
    input_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/genomes_erwinia_amylovora/ncbi_dataset/data'
    output_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/genomes_erwinia_amylovora/ncbi_dataset/fasta'
    gb_to_fasta(input_path, output_path)


if __name__ == "__main__":
    main()
