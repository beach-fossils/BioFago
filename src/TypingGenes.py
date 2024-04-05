import csv
import glob
import subprocess
from collections import defaultdict
import os
import re
import pandas as pd


class BlastTyping:
    def __init__(self, csv_path, output_dir):
        self.csv_path = csv_path
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def read_csv(self):
        # Read sequences from CSV
        return pd.read_csv(self.csv_path)

    def prepare_query_files(self, df):
        # Prepare individual query files for each sequence
        for index, row in df.iterrows():
            gene_type = row['Type']
            gene_name = row['Gene/Locus Tag'].split("_")[0]  # Only take the gene name before "_"
            sequence = row['Sequence']
            gene_folder = os.path.join(self.output_dir, gene_type)
            os.makedirs(gene_folder, exist_ok=True)
            file_path = os.path.join(gene_folder, f"{gene_name}.fasta")
            with open(file_path, 'a') as file:  # Append to the file if it already exists
                file.write(f">{gene_type}_{row['Gene/Locus Tag']}\n{sequence}\n")

    def create_blast_db(self, df):
        # Create a separate BLAST database for each gene
        for gene in df['Gene/Locus Tag'].str.split("_").str[0].unique():  # Only take the gene name before "_"
            db_path = os.path.join(self.output_dir, f"{gene}_db.fasta")
            with open(db_path, 'w') as db_file:
                for index, row in df[df['Gene/Locus Tag'].str.contains(
                        gene)].iterrows():  # Check if the gene name is contained in 'Gene/Locus Tag'
                    db_file.write(f">{row['Type']}_{row['Gene/Locus Tag']}\n{row['Sequence']}\n")
            subprocess.run(['makeblastdb', '-in', db_path, '-dbtype', 'prot', '-out', db_path.replace('.fasta', '')])

    def run_blast(self, df):
        # Run BLAST for each sequence against its corresponding gene database

        # for each of the sequence (gene) in the folder of each type we run the blast against the database of the gene (retrieve the name of it)
        for type_dir in os.listdir(self.output_dir):
            type_path = os.path.join(self.output_dir, type_dir)
            if os.path.isdir(type_path):
                for file in os.listdir(type_path):
                    if file.endswith(".fasta"):
                        gene_name = file.split(".")[0]
                        query_file_path = os.path.join(type_path, file)
                        db_path = os.path.join(self.output_dir, f"{gene_name}_db")  # Remove '.fasta'
                        output_file = os.path.join(type_path, f"blast_{gene_name}.txt")
                        subprocess.run([
                            'blastp', '-query', query_file_path, '-db', db_path,
                            '-outfmt', "6 qseqid sseqid pident length mismatches gapopen evalue bitscore",
                            '-qcov_hsp_perc', '70',  # Minimum query coverage per hsp
                            '-out', output_file
                        ], check=True)

    def process(self):
        df = self.read_csv()
        self.prepare_query_files(df)
        self.create_blast_db(df)
        self.run_blast(df)


class BlastNomenclature:
    def __init__(self, base_folder, reference_folder, output_path, identity_threshold):
        self.base_folder = base_folder
        self.reference_folder = reference_folder
        self.output_path = output_path
        self.identity_threshold = identity_threshold
        self.typing_dict = {}

    def parse_reference_folder(self):
        for blast_file in glob.glob(os.path.join(self.reference_folder, 'blast_*.txt')):
            gene = os.path.basename(blast_file).split('_')[1].split('.')[0]
            print(f"Processing {blast_file} for gene {gene}")
            with open(blast_file, 'r') as file:
                for line in file:
                    parts = line.strip().split('\t')
                    query, subject, identity = parts[0], parts[1], float(parts[2])
                    self.add_subject(gene, subject, identity)
            print(f"Completed processing {blast_file}")

    def add_subject(self, gene, subject, identity):
        if gene not in self.typing_dict:
            self.typing_dict[gene] = {'main': [], 'variants': []}

        if identity >= self.identity_threshold:
            self.typing_dict[gene]['main'].append(subject)
        else:
            # Ensure uniqueness in variants
            added = False
            for variant in self.typing_dict[gene]['variants']:
                if identity == variant['identity']:
                    variant['subjects'].append(subject)
                    added = True
                    break
            if not added:
                self.typing_dict[gene]['variants'].append({'identity': identity, 'subjects': [subject]})

    def verify_variants(self):
        # Iterate over each gene and its variants
        for gene, data in self.typing_dict.items():
            # Temporary dictionary to hold updated groups after verification
            updated_variants = {'main': data['main'], 'variants': []}
            variant_identities = {}

            # Verify each variant group
            for variant in data['variants']:
                variant_identity = variant['identity']
                variant_key = f"{gene}{variant_identity}"

                # Check if BLAST results for this variant are in line with others
                for subject in variant['subjects']:
                    subject_prefix = subject.split('_')[0]  # e.g., KL02 from KL02_amsE
                    blast_results_file = os.path.join(self.base_folder, subject_prefix, f'blast_{gene}.txt')

                    # Read the BLAST results for this subject
                    if os.path.exists(blast_results_file):
                        with open(blast_results_file, 'r') as file:
                            for line in file:
                                parts = line.strip().split('\t')
                                query, target, identity = parts[0], parts[1], float(parts[2])

                                # Check identity against other subjects in the same gene group
                                if target in variant_identities:
                                    if identity >= self.identity_threshold:
                                        # If similar, merge this subject with the previously identified group
                                        updated_variants['variants'][variant_identities[target]]['subjects'].append(
                                            subject)
                                        break
                                else:
                                    # If distinct, keep this subject in its current variant group
                                    if identity < self.identity_threshold:
                                        variant_identities[subject] = len(updated_variants['variants'])
                                        updated_variants['variants'].append(
                                            {'identity': variant_identity, 'subjects': [subject]})
                                        break

            # Update the typing dictionary with verified variant groups
            self.typing_dict[gene] = updated_variants

    def print_dict(self):
        for gene, data in self.typing_dict.items():
            print(f"{gene}: {', '.join(data['main'])}")
            variant_num = 2
            for variant in sorted(data['variants'], key=lambda x: x['identity'], reverse=True):
                print(f"{gene}{variant_num}: {', '.join(variant['subjects'])}")
                variant_num += 1

    def export_to_csv(self):
        rows = []  # Collect all rows here before sorting

        for gene, data in self.typing_dict.items():
            # Main group
            for subject in data['main']:
                rows.append([subject, f"{gene}"])  # Collect row

            # Variants
            variant_num = 2  # Start naming variants from 2
            for variant in sorted(data['variants'], key=lambda x: x['identity'], reverse=True):
                for subject in variant['subjects']:
                    rows.append([subject, f"{gene}{variant_num}"])  # Collect row
                variant_num += 1

        # Now sort the rows based on the "Gene_Name" which is now the first element of each row
        sorted_rows = sorted(rows, key=lambda x: x[0])

        # Write sorted rows to CSV
        with open(self.output_path, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(["Gene_Name", "Type_Tag"])  # Header row
            csvwriter.writerows(sorted_rows)  # Write all sorted rows at once


class CSVManipulator:
    def __init__(self, division_genes_path, nomenclature_path):
        self.division_genes_df = pd.read_csv(division_genes_path)
        self.nomenclature_df = pd.read_csv(nomenclature_path)

    def process_csvs(self):
        # Add a unified column based on Type and Gene/Locus Tag
        self.division_genes_df['Unified'] = self.division_genes_df['Type'] + '_' + self.division_genes_df[
            'Gene/Locus Tag']

        # Merge the two dataframes on the unified column
        merged_df = pd.merge(self.division_genes_df, self.nomenclature_df, left_on='Unified', right_on='Gene_Name',
                             how='left')

        # Fill missing Type_Tag values with Gene/Locus Tag from the first dataframe
        merged_df['Type_Tag'] = merged_df['Type_Tag'].fillna(merged_df['Gene/Locus Tag'])

        # Reorder columns and drop unnecessary ones
        final_df = merged_df[['Type', 'Unified', 'Type_Tag', 'Sequence']]

        return final_df


def main():
    # part one for typing w BLAST:
    csv_path = "/Users/josediogomoura/Documents/BioFago/BioFago/data/AApurify/LPS/90_blastp/division_genes.csv"
    output_dir = "/Users/josediogomoura/Documents/BioFago/BioFago/data/AApurify/LPS/90_blastp/Typing"
    blast_typing = BlastTyping(csv_path, output_dir)
    blast_typing.process()

    # part two for nomenclature:
    output_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/AApurify/LPS/90_blastp/Nomenclature/nomenclature.csv'
    base_dir = "/Users/josediogomoura/Documents/BioFago/BioFago/data/AApurify/LPS/90_blastp/Typing"
    identity_threshold = 90
    reference_dir = "/Users/josediogomoura/Documents/BioFago/BioFago/data/AApurify/LPS/90_blastp/Typing/OL01"

    nomenclature = BlastNomenclature(base_dir, reference_dir, output_path, identity_threshold)
    nomenclature.parse_reference_folder()
    nomenclature.verify_variants()
    nomenclature.print_dict()
    nomenclature.export_to_csv()

    # part three for linking the csvs:
    # Example usage
    gene_info_csv = '/Users/josediogomoura/Documents/BioFago/BioFago/data/AApurify/LPS/90_blastp/division_genes.csv'
    gene_name_type_tag_csv = '/Users/josediogomoura/Documents/BioFago/BioFago/data/AApurify/LPS/90_blastp/Nomenclature/nomenclature.csv'
    output_csv = '/Users/josediogomoura/Documents/BioFago/BioFago/data/AApurify/LPS/90_blastp/Nomenclature/final_curation.csv'

    csv_linker = CSVManipulator(gene_info_csv, gene_name_type_tag_csv)
    final_df = csv_linker.process_csvs()
    final_df.to_csv(output_csv, index=False)


if __name__ == "__main__":
    main()
