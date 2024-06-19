import csv
import subprocess
import pandas as pd
from pathlib import Path


def parse_protein_faa(file_path):
    data = {'Protein ID': [], 'EC Number': [], 'Gene Name': [], 'Description': [], 'Sequence': []}
    with open(file_path, 'r') as file:
        content = file.read().strip().split('>')
        for entry in content:
            if entry:
                header, sequence = entry.split('\n', 1)
                sequence = sequence.replace('\n', '')
                parts = header.split('~~~')
                protein_id = parts[0].split()[0]
                ec_number = parts[0].split()[1] if len(parts[0].split()) > 1 else 'N/A'
                gene_name = parts[1]
                description = parts[2] if len(parts) > 2 else 'N/A'
                data['Protein ID'].append(protein_id)
                data['EC Number'].append(ec_number)
                data['Gene Name'].append(gene_name)
                data['Description'].append(description)
                data['Sequence'].append(sequence)

    return pd.DataFrame(data)

def prepare_blast_db(loci_info, output_fasta):
    with open(output_fasta, "w") as fasta_file:
        for locus_id, genes in loci_info.items():
            for gene in genes:
                if gene["translation"]:
                    fasta_file.write(f">{locus_id}~~~{gene['gene']}~~~{gene['product']}\n{gene['translation']}\n")

class BlastProtein:
    def __init__(self, reference_types, db_folder, results_folder, Database):
        self.db_folder = Path(db_folder)
        self.results_folder = Path(results_folder)
        self.db_name = "protein_db"
        self.blast_results = {}
        self.db = Database(reference_types)  # Store the Database instance
        self.loci = self.db.loci  # Store the loci information

        # Ensure results folder exists
        self.results_folder.mkdir(parents=True, exist_ok=True)

    def create_blast_db(self, protein_fasta):
        try:
            make_db_cmd = [
                'makeblastdb',
                '-in', protein_fasta,
                '-dbtype', 'prot',
                '-out', str(self.db_folder / self.db_name)
            ]
            subprocess.run(make_db_cmd, check=True)
        except subprocess.CalledProcessError as e:
            raise Exception(f"Error in creating BLAST database: {e}")

    def blast_protein(self, protein_sequence, query_id):
        try:
            query_fasta = self.db_folder / f"{query_id}.fasta"
            with open(query_fasta, "w") as query_file:
                query_file.write(f">{query_id}\n{protein_sequence}\n")

            blast_output_file = self.results_folder / f"{query_id}_blast_output.txt"
            blast_cmd = [
                'blastp',
                '-query', str(query_fasta),
                '-db', str(self.db_folder / self.db_name),
                '-out', str(blast_output_file),
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
            ]
            subprocess.run(blast_cmd, check=True)
            return blast_output_file
        except subprocess.CalledProcessError as e:
            raise Exception(f"Error in BLAST search: {e}")

    def blast_all_proteins(self, protein_df, output_csv):
        with open(output_csv, "w", newline="") as csvfile:
            fieldnames = ["Protein ID", "EC Number", "Gene Name", "Description", "subject_id", "percent_identity",
                          "alignment_length", "e_value", "bit_score"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for _, row in protein_df.iterrows():
                protein_id = row["Protein ID"]
                protein_sequence = row["Sequence"]
                blast_output_file = self.blast_protein(protein_sequence, protein_id)

                with open(blast_output_file, "r") as output:
                    for line in output:
                        columns = line.strip().split("\t")
                        writer.writerow({
                            "Protein ID": protein_id,
                            "EC Number": row["EC Number"],
                            "Gene Name": row["Gene Name"],
                            "Description": row["Description"],
                            "subject_id": columns[1],
                            "percent_identity": columns[2],
                            "alignment_length": columns[3],
                            "e_value": columns[10],
                            "bit_score": columns[11]
                        })