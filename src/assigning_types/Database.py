import subprocess
import logging
from pathlib import Path
import pandas as pd
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio import SeqIO, SeqRecord, Seq
from src.utils.get_prokka_faa_file import get_prokka_faa_file
import csv
import tempfile
import os

_PROTEIN_ALIGNER = PairwiseAligner(scoring='blastp', mode='local')

# Setup logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


class Database:
    def __init__(self, gb_file):
        """
        Initialize the Database with a GenBank file.

        :param gb_file: Path to the GenBank file.
        """
        self.gb_file = gb_file
        self.loci = self._parse_gb_file()

    def _parse_gb_file(self):
        """
        Parse the GenBank file to extract loci information.

        :return: Dictionary of loci, each containing a list of genes.
        """
        loci_dict = {}

        with open(self.gb_file, "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                locus_id = record.name.split('_')[0]
                if locus_id not in loci_dict:
                    loci_dict[locus_id] = []

                for feature in record.features:
                    if feature.type == "CDS":
                        locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                        if locus_tag:
                            gene = feature.qualifiers.get("gene", ["unknown_gene"])[0]
                            product = feature.qualifiers.get("product", ["unknown_product"])[0]
                            protein_id = feature.qualifiers.get("protein_id", [""])[0]
                            translation = feature.qualifiers.get("translation", [""])[0]
                            nucleotide_sequence = self._find_nt_sequence(locus_tag)

                            loci_dict[locus_id].append({
                                "locus_tag": locus_tag,
                                "gene": gene,
                                "product": product,
                                "protein_id": protein_id,
                                "translation": translation,
                                "nucleotide_sequence": nucleotide_sequence  # Fetch during parsing
                            })
        logging.info(f"Parsed {len(loci_dict)} loci from GenBank file.")
        for locus_id, genes in loci_dict.items():
            logging.debug(f"Locus {locus_id} has {len(genes)} genes.")
        return loci_dict

    def _find_nt_sequence(self, locus_tag):
        """
        Find the nucleotide sequence for a given locus tag.

        :param locus_tag: The locus tag to find the sequence for.
        :return: Nucleotide sequence as a string, or None if not found.
        """
        with open(self.gb_file, "r") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if (feature.type == "CDS" and "locus_tag" in feature.qualifiers
                            and feature.qualifiers["locus_tag"][0] == locus_tag):
                        return str(feature.location.extract(record).seq)
        return None

    @property
    def first_gene(self):
        """
        Get the first gene in the first locus.

        :return: Dictionary containing first gene information.
        """
        first_locus_key = next(iter(self.loci))
        return self.loci[first_locus_key][0]

    @property
    def last_gene(self):
        """
        Get the last gene in the first locus.

        :return: Dictionary containing last gene information.
        """
        first_locus_key = next(iter(self.loci))
        return self.loci[first_locus_key][-1]

    @property
    def first_gene_sequence_nt(self):
        """
        Get the nucleotide sequence of the first gene.

        :return: Nucleotide sequence as a string.
        """
        return self.first_gene['nucleotide_sequence']

    @property
    def last_gene_sequence_nt(self):
        """
        Get the nucleotide sequence of the last gene.

        :return: Nucleotide sequence as a string.
        """
        return self.last_gene['nucleotide_sequence']

    def create_faa_file(self, output_faa):
        """
        Create a FASTA file with protein sequences from the loci reference_crispr.

        :param output_faa: Path to the output FASTA file.
        """
        records = []
        for locus, genes in self.loci.items():
            for gene in genes:
                if gene["translation"]:
                    record = SeqRecord.SeqRecord(
                        Seq.Seq(gene["translation"]),
                        id=f"{gene['locus_tag']} {gene['product']}~~~{gene['gene']}",
                        description=""
                    )
                    records.append(record)

        with open(output_faa, "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")
        logging.info(f"Created FASTA file with {len(records)} protein sequences.")

    def get_loci_info(self):
        """
        Get the loci information.

        :return: Dictionary of loci with gene information.
        """
        return self.loci

    def write_loci_to_csv(self, output_csv):
        """
        Write the loci information to a CSV file.

        :param output_csv: Path to the output CSV file.
        """
        with open(output_csv, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(
                ["Locus", "Gene", "Locus Tag", "Product", "Protein ID", "Translation", "Nucleotide Sequence"])
            for locus, genes in self.loci.items():
                for gene in genes:
                    writer.writerow([locus, gene["gene"], gene["locus_tag"], gene["product"], gene["protein_id"],
                                     gene["translation"], gene["nucleotide_sequence"]])
        logging.info(f"Created CSV file for loci information: {output_csv}")

    def _get_first_gene(self):
        """
        Get the first gene in the loci.

        :return: Dictionary of the first gene's information.
        """
        for locus, genes in self.loci.items():
            if genes:
                return genes[0]
        return None

    def _get_last_gene(self):
        """
        Get the last gene in the loci.

        :return: Dictionary of the last gene's information.
        """
        for locus, genes in self.loci.items():
            if genes:
                return genes[-1]
        return None

    def get_flank_genes(self):
        """
        Get the first and last genes of the first locus.

        :return: Tuple containing dictionaries of the first and last genes.
        """
        first_locus_key = next(iter(self.loci))
        first_gene = self.loci[first_locus_key][0]
        last_gene = self.loci[first_locus_key][-1]
        return first_gene, last_gene


class BlastProteinv2:
    def __init__(self, reference_types, db_folder, results_folder):
        self.db_folder = Path(db_folder) if db_folder else None
        self.results_folder = Path(results_folder) if results_folder else None
        self.db_name = "protein_db"
        self.blast_results = {}
        self.db = Database(reference_types)
        self.loci = self.db.loci

        # Ensure results folder exists if provided
        if self.results_folder:
            self.results_folder.mkdir(parents=True, exist_ok=True)
        logging.info("BlastProteinv2 initialized.")

    def create_csv_loci(self, output_file):
        with open(output_file, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(
                ["Locus", "Gene", "Locus Tag", "Product", "Protein ID", "Translation", "Nucleotide Sequence"])
            for locus, genes in self.loci.items():
                for gene in genes:
                    writer.writerow([locus, gene["gene"], gene["locus_tag"], gene["product"], gene["protein_id"],
                                     gene["translation"], gene["nucleotide_sequence"]])
        logging.info(f"Created CSV file for loci information: {output_file}")

    def make_blast_db(self, input_fasta):
        try:
            makeblastdb_cmd = [
                'makeblastdb',
                '-in', input_fasta,
                '-dbtype', 'prot',
                '-out', str(self.db_folder / self.db_name)
            ]
            logging.debug(f"Running makeblastdb command: {' '.join(makeblastdb_cmd)}")
            subprocess.run(makeblastdb_cmd, check=True)
            logging.info("BLAST database created successfully.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error in making BLAST database: {e}")
            raise Exception(f"Error in making BLAST database: {e}")

    def run_blast(self, query_fasta, output_file):
        try:
            blast_cmd = [
                'blastp',
                '-query', str(query_fasta),
                '-db', str(self.db_folder / self.db_name),
                '-outfmt', '6 qseqid sseqid pident length qlen slen',
                '-out', str(output_file)
            ]
            logging.debug(f"Running blastp command: {' '.join(blast_cmd)}")
            subprocess.run(blast_cmd, check=True)
            logging.info(f"BLAST results written to {output_file}.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error in running BLAST: {e}")
            raise Exception(f"Error in running BLAST: {e}")

    def analyze_blast_results(self, blast_output_file):
        best_hit = None

        try:
            with open(blast_output_file, 'r') as file:
                for line in file:
                    parts = line.strip().split('\t')
                    qseqid, sseqid, pident, length, qlen, slen = parts
                    pident = float(pident)
                    length = int(length)
                    qlen = int(qlen)
                    slen = int(slen)

                    if not best_hit or best_hit['pident'] < pident:
                        best_hit = {'sseqid': sseqid, 'pident': pident, 'length': length, 'qlen': qlen, 'slen': slen}
                        logging.debug(f"Best hit updated: {best_hit}")

        except Exception as e:
            logging.error(f"Error in analyzing BLAST results from {blast_output_file}: {e}")

        return best_hit

    def process_csv(self, input_csv, proteins_faa, final_csv_file):
        # Make sure the BLAST database is created
        self.make_blast_db(proteins_faa)

        results = []
        with open(input_csv, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if row['Translation']:
                    # Create a temporary FASTA file for the query sequence
                    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".fasta") as temp_fasta:
                        temp_fasta.write(
                            f">{row['Locus Tag']} {row['Product']}~~~{row['Gene']}\n{row['Translation']}\n")
                        temp_fasta_name = temp_fasta.name

                    output_filename = f"{row['Locus']}_{row['Gene']}_{row['Locus Tag']}_blast_output.txt"
                    temp_output_file = self.results_folder / output_filename
                    self.run_blast(temp_fasta_name, temp_output_file)
                    best_hit = self.analyze_blast_results(temp_output_file)

                    # Add locus and gene information to the BLAST results file
                    with open(temp_output_file, 'a') as blast_file:
                        blast_file.write(f"\nLocus: {row['Locus']}, Gene: {row['Gene']}\n")

                    row['Best BLAST Identity'] = best_hit['pident'] if best_hit else ''
                    row['Alignment Length'] = best_hit['length'] if best_hit else ''
                    row['Translation Length'] = len(row['Translation']) if row['Translation'] else ''
                    row['Reference Length'] = best_hit['slen'] if best_hit else ''

                    # Clean up the temporary file
                    Path(temp_fasta_name).unlink()

                results.append(row)

        with open(final_csv_file, 'w', newline='') as csvfile:
            fieldnames = ["Locus", "Gene", "Locus Tag", "Product", "Protein ID", "Translation", "Nucleotide Sequence",
                          "Best BLAST Identity", "Alignment Length", "Translation Length", "Reference Length"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for result in results:
                writer.writerow(result)
        logging.info(f"Created final CSV file: {final_csv_file}")


class TypeAnalysis:
    def __init__(self, input_csv, output_csv, proteins_faa):
        self.input_csv = input_csv
        self.output_csv = output_csv
        self.proteins_faa = proteins_faa
        self.number_of_genes = self.calculate_number_of_genes()
        self.total_protein_length = self.calculate_total_protein_length()
        self.type_hierarchy = ['None', 'Low', 'Good', 'High', 'Very high', 'Perfect']

    def calculate_number_of_genes(self):
        """
        Calculate the total number of genes from the protein FASTA file.
        """
        with open(self.proteins_faa, 'r') as file:
            return sum(1 for line in file if line.startswith('>'))

    def calculate_total_protein_length(self):
        """
        Calculate the total length of all protein sequences from the protein FASTA file.
        """
        total_length = 0
        with open(self.proteins_faa, 'r') as file:
            for line in file:
                if not line.startswith('>'):
                    total_length += len(line.strip())
        return total_length

    def analyze_locus(self):
        """
        Analyze the loci by calculating the translation and reference coverage and determining presence/absence.
        """
        df = pd.read_csv(self.input_csv)
        df['Translation Coverage'] = df['Alignment Length'] / df['Translation Length']
        df['Reference Coverage'] = df['Alignment Length'] / df['Reference Length']
        df['Presence/Absence'] = ((df['Best BLAST Identity'] >= 95) &
                                  (df['Translation Coverage'] >= 0.9) &
                                  (df['Reference Coverage'] >= 0.9)).astype(int)
        df['Truncated'] = ((df['Best BLAST Identity'] >= 95) &
                           (df['Translation Coverage'] < 0.9) &
                           ((df['Translation Coverage'] + df['Reference Coverage']) > 0.05)).astype(int)
        df.to_csv(self.output_csv, index=False)
        return df

    def count_truncated_genes(self, group):
        """
        Count the number of truncated genes in the group.
        """
        return group['Truncated'].sum()

    def count_extra_genes(self, group):
        """
        Count the number of extra genes in the group based on proteins.faa file.
        """
        present_genes = group['Presence/Absence'].sum()
        truncated_genes = self.count_truncated_genes(group)
        return self.number_of_genes - (present_genes + truncated_genes)

    def calculate_locus_coverage(self, group):
        """
        Calculate the coverage of the locus.
        """
        total_alignment_length = group['Alignment Length'].sum()
        return total_alignment_length / group['Translation Length'].sum()

    def assign_type(self, df):
        """
        Assign a type to each locus based on the analysis.
        """
        report = []
        grouped = df.groupby('Locus')

        for locus, group in grouped:
            total_genes = group.shape[0]
            present_genes = group['Presence/Absence'].sum()
            truncated_genes = self.count_truncated_genes(group)
            extra_genes = self.count_extra_genes(group)
            locus_coverage = self.calculate_locus_coverage(group)

            type_assigned = self.determine_type(present_genes, truncated_genes, extra_genes, locus_coverage)

            report.append({
                'Locus': locus,
                'Total Genes': total_genes,
                'Present Genes': present_genes,
                'Truncated Genes': truncated_genes,
                'Extra Genes': extra_genes,
                'Locus Coverage': round(locus_coverage, 2),
                'Type': type_assigned
            })

        type_report = pd.DataFrame(report)
        final_type_locus, final_type = self.select_best_type(type_report)
        return type_report, final_type_locus, final_type

    def determine_type(self, present_genes, truncated_genes, extra_genes, locus_coverage):
        """
        Determine the type of a locus based on various criteria.
        """
        if present_genes + truncated_genes == self.number_of_genes:
            if present_genes == self.number_of_genes:
                if truncated_genes == 0 and extra_genes == 0:
                    if locus_coverage >= 0.99:
                        return 'Perfect'
            if truncated_genes == 0 and extra_genes == 0:
                if locus_coverage >= 0.99:
                    return 'Very high'
            if truncated_genes <= 3 and extra_genes == 0:
                if locus_coverage >= 0.99:
                    return 'High'
            if truncated_genes <= 3 and extra_genes <= 1:
                if locus_coverage >= 0.95:
                    return 'Good'
            if truncated_genes <= 3 and extra_genes <= 2:
                if locus_coverage >= 0.90:
                    return 'Low'
        else:
            if locus_coverage >= 0.95 and truncated_genes <= 3 and extra_genes <= 1:
                return 'Good'
            elif locus_coverage >= 0.90 and truncated_genes <= 3 and extra_genes <= 2:
                return 'Low'
        return 'None'

    def select_best_type(self, type_report):
        """
        Select the best type from the type report based on the hierarchy.
        """
        highest_priority_type = 'None'
        highest_priority_locus = ''
        for idx, row in type_report.iterrows():
            if self.type_hierarchy.index(row['Type']) > self.type_hierarchy.index(highest_priority_type):
                highest_priority_type = row['Type']
                highest_priority_locus = row['Locus']
        return highest_priority_locus, highest_priority_type


if __name__ == "__main__":
    base_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/assign_types/cellulose/test_1'
    reference_types = "/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/assign_types/reference_types_database/cellulose/types_cellulose.gb"
    prokka_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/reference_crispr/assign_types/cellulose/test_1/prokka'
    proteins_faa = get_prokka_faa_file(prokka_folder)

    db_folder = os.path.join(base_folder, "db_folder")
    results_folder = os.path.join(base_folder, "results")
    final_csv_file = os.path.join(base_folder, "final.csv")
    output_csv = os.path.join(base_folder, "loci_gene_presence.csv")

    # Initialize BlastProteinv2
    blast = BlastProteinv2(reference_types=reference_types, db_folder=db_folder, results_folder=results_folder)

    # Create CSV of loci information
    # Assuming input_csv is created as part of the process, you might need to define input_csv path
    input_csv = os.path.join(base_folder, "loci.csv")
    blast.create_csv_loci(input_csv)

    # Process CSV and perform BLAST searches
    blast.process_csv(input_csv, proteins_faa, final_csv_file)

    # Analyze locus information
    gene_presence = TypeAnalysis(input_csv=final_csv_file, output_csv=output_csv, proteins_faa=proteins_faa)
    df_analyzed = gene_presence.analyze_locus()
    type_report, final_type_locus, final_type = gene_presence.assign_type(df_analyzed)

    # Print all values in the report
    pd.set_option('display.max_rows', None)  # Ensure all rows are printed
    for idx, row in type_report.iterrows():
        print(row)

    print(f"Final Assigned Type: {final_type_locus} ({final_type})")