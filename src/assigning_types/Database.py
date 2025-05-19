import multiprocessing
import subprocess
import logging
from functools import partial
from pathlib import Path
import pandas as pd
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio import SeqIO, SeqRecord, Seq
from utils.get_prokka_faa_file import get_prokka_faa_file
import csv
import tempfile
import os
import sys
from typing import Dict, List, Tuple

# Import quiet mode module
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from quiet_mode import QUIET_MODE

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
            
            # Respect quiet mode when running subprocess
            if QUIET_MODE:
                with open(os.devnull, 'w') as devnull:
                    subprocess.run(makeblastdb_cmd, check=True, stdout=devnull, stderr=devnull)
            else:
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
                '-outfmt', '6 qseqid sseqid pident length qlen slen qstart qend sstart send',
                '-out', str(output_file)
            ]
            logging.debug(f"Running blastp command: {' '.join(blast_cmd)}")
            
            # Respect quiet mode when running subprocess
            if QUIET_MODE:
                with open(os.devnull, 'w') as devnull:
                    subprocess.run(blast_cmd, check=True, stdout=devnull, stderr=devnull)
            else:
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
                    qseqid, sseqid, pident, length, qlen, slen, qstart, qend, sstart, send = parts
                    pident = float(pident)
                    length = int(length)
                    qlen = int(qlen)
                    slen = int(slen)
                    qstart = int(qstart)
                    qend = int(qend)
                    sstart = int(sstart)
                    send = int(send)

                    if not best_hit or best_hit['pident'] < pident:
                        best_hit = {
                            'sseqid': sseqid,
                            'pident': pident,
                            'length': length,
                            'qlen': qlen,
                            'slen': slen,
                            'qstart': qstart,
                            'qend': qend,
                            'sstart': sstart,
                            'send': send
                        }
                        logging.debug(f"Best hit updated: {best_hit}")

        except Exception as e:
            logging.error(f"Error in analyzing BLAST results from {blast_output_file}: {e}")

        return best_hit

    def process_csv(self, input_csv, proteins_faa, final_csv_file):
        # Make sure the BLAST database is created
        self.make_blast_db(proteins_faa)

        # Read CSV using pandas
        df = pd.read_csv(input_csv)

        # Use multiprocessing to run BLAST searches in parallel
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            results = pool.map(partial(self.process_row, proteins_faa=proteins_faa), df.to_dict('records'))

        # Convert results back to DataFrame and save
        result_df = pd.DataFrame(results)
        result_df.to_csv(final_csv_file, index=False)
        logging.info(f"Created final CSV file: {final_csv_file}")

    def process_row(self, row, proteins_faa):
        if row['Translation']:
            with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix=".fasta") as temp_fasta:
                temp_fasta.write(f">{row['Locus Tag']} {row['Product']}~~~{row['Gene']}\n{row['Translation']}\n")
                temp_fasta_name = temp_fasta.name

            output_filename = f"{row['Locus']}_{row['Gene']}_{row['Locus Tag']}_blast_output.txt"
            temp_output_file = self.results_folder / output_filename
            self.run_blast(temp_fasta_name, temp_output_file)
            best_hit = self.analyze_blast_results(temp_output_file)

            row['Best BLAST Identity'] = best_hit['pident'] if best_hit else ''
            row['Alignment Length'] = best_hit['length'] if best_hit else ''
            row['Query Length'] = best_hit['qlen'] if best_hit else ''
            row['Subject Length'] = best_hit['slen'] if best_hit else ''
            row['Query Start'] = best_hit['qstart'] if best_hit else ''
            row['Query End'] = best_hit['qend'] if best_hit else ''
            row['Subject Start'] = best_hit['sstart'] if best_hit else ''
            row['Subject End'] = best_hit['send'] if best_hit else ''

            # Clean up the temporary file
            Path(temp_fasta_name).unlink()

        return row


class TypeAnalysis:
    def __init__(self, input_csv, output_csv, proteins_faa):
        self.input_csv = input_csv
        self.output_csv = output_csv
        self.proteins_faa = proteins_faa
        self.number_of_genes = self.calculate_number_of_genes()
        self.total_protein_length = self.calculate_total_protein_length()
        self.type_hierarchy = ['Very low', 'Low', 'Good', 'High', 'Very High', 'Perfect']
        # Updated locus types to include all systems
        self.locus_types = [
            'capsule', 'cellulose', 'lps', 'sorbitol',
            'flag_i', 'flag_ii', 'flag_iii', 'flag_iv',
            't3ss_i', 't3ss_ii', 't6ss_i', 't6ss_ii'
        ]

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
        df = pd.read_csv(self.input_csv)
        df['Query Coverage'] = df['Alignment Length'] / df['Query Length']
        df['Subject Coverage'] = df['Alignment Length'] / df['Subject Length']
        df['Length Ratio'] = df['Query Length'] / df['Subject Length']

        df['Presence/Absence'] = ((df['Best BLAST Identity'] >= 95) &
                                  (df['Query Coverage'] >= 0.9) &
                                  (df['Subject Coverage'] >= 0.9) &
                                  (df['Length Ratio'] >= 0.9) &
                                  (df['Length Ratio'] <= 1.1)).astype(int)

        # Flag ANY difference from perfect identity
        df['Significantly Different'] = ((df['Best BLAST Identity'] < 100) |  # Any difference in identity
                                         (df['Query Coverage'] < 1.0) |  # Any difference in coverage
                                         (df['Subject Coverage'] < 1.0) |  # Any difference in coverage
                                         (df['Length Ratio'] != 1.0)).astype(int)  # Any difference in length

        df['Truncated'] = ((df['Best BLAST Identity'] >= 95) &
                           (df['Query Coverage'] < 0.9) &
                           ((df['Query Coverage'] + df['Subject Coverage']) > 0.05) &
                           (df['Length Ratio'] < 0.9)).astype(int)

        df['Flagged Genes'] = df.apply(
            lambda row: row['Gene'] if row['Significantly Different'] else '',
            axis=1
        )

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
        total_query_coverage = group['Query Coverage'].sum()
        total_subject_coverage = group['Subject Coverage'].sum()
        return (total_query_coverage + total_subject_coverage) / 2

    def assign_type(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, str, str, str, Dict[str, str]]:
        try:
            if df.empty:
                logging.warning("Input DataFrame is empty")
                return pd.DataFrame(), '', '', '', {}

            report = []
            # Initialize results dictionary with all locus types
            locus_results = {locus_type: {'type': 'Very low', 'locus': '', 'flagged_genes': []}
                           for locus_type in self.locus_types}

            for locus, group in df.groupby('Locus'):
                total_genes = group.shape[0]
                present_genes = group['Presence/Absence'].sum()
                truncated_genes = self.count_truncated_genes(group)
                extra_genes = self.count_extra_genes(group)
                locus_coverage = self.calculate_locus_coverage(group)

                flagged_genes = group[group['Flagged Genes'] != '']['Flagged Genes'].tolist()
                flagged_genes = list(
                    set([gene.strip() for genes in flagged_genes for gene in genes.split(',') if gene.strip()]))

                type_assigned = self.determine_type(present_genes, truncated_genes, extra_genes, locus_coverage,
                                                    flagged_genes)

                report.append({
                    'Locus': locus,
                    'Total Genes': total_genes,
                    'Present Genes': present_genes,
                    'Truncated Genes': truncated_genes,
                    'Extra Genes': extra_genes,
                    'Locus Coverage': round(locus_coverage, 2),
                    'Type': type_assigned,
                    'Flagged Genes': ', '.join(flagged_genes) if flagged_genes else ''
                })

                locus_type = self.determine_locus_type(str(locus))
                if locus_type:
                    if self.type_hierarchy.index(type_assigned) > self.type_hierarchy.index(
                            locus_results[locus_type]['type']):
                        locus_results[locus_type] = {
                            'type': type_assigned,
                            'locus': locus,
                            'flagged_genes': flagged_genes
                        }

            type_report = pd.DataFrame(report)

            # Select the best locus based on type
            best_locus = max(report, key=lambda x: self.type_hierarchy.index(x['Type']))
            final_type_locus = best_locus['Locus']
            final_type = best_locus['Type']
            flagged_genes_str = best_locus['Flagged Genes']

            # Format the locus information for each locus type
            formatted_locus_info = {}
            for locus_type in self.locus_types:
                result = locus_results[locus_type]
                locus = result['locus']
                type_assigned = result['type']
                flagged_genes = result['flagged_genes']

                if locus:
                    formatted_info = f"{locus} ({type_assigned})"
                    if flagged_genes:
                        formatted_info += f" - Flagged genes: {', '.join(sorted(flagged_genes))}"
                else:
                    formatted_info = f"({type_assigned})"

                formatted_locus_info[f"{locus_type}_locus"] = formatted_info

            logging.info(f"Locus results: {locus_results}")
            logging.info(
                f"Returning from assign_type: {type_report}, {final_type_locus}, {final_type}, {flagged_genes_str}, {formatted_locus_info}")

            return type_report, final_type_locus, final_type, flagged_genes_str, formatted_locus_info

        except Exception as e:
            logging.error(f"Error in assign_type: {str(e)}")
            logging.error(f"DataFrame info:\n{df.info()}")
            logging.error(f"DataFrame head:\n{df.head()}")
            logging.error(f"Type hierarchy: {self.type_hierarchy}")

            # Create a complete default result dictionary with all locus types
            default_results = {f"{locus_type}_locus": "None" for locus_type in self.locus_types}
            return pd.DataFrame(), '', '', '', default_results

    def determine_locus_type(self, locus: str) -> str:
        """
        Determine the type of locus based on folder name or locus prefix.
        """
        # First check for direct folder name matches
        folder_mappings = {
            'types_flag_I': 'flag_i',
            'types_flag_II': 'flag_ii',
            'types_flag_III': 'flag_iii',
            'types_flag_IV': 'flag_iv',
            'types_T3SS_I': 't3ss_i',
            'types_T3SS_II': 't3ss_ii',
            'types_T6SS_I': 't6ss_i',
            'types_T6SS_II': 't6ss_ii',
            'types_capsule': 'capsule',
            'types_cellulose': 'cellulose',
            'types_lps': 'lps',
            'types_srl': 'sorbitol'
        }

        # Check for folder name match first
        if str(locus) in folder_mappings:
            return folder_mappings[str(locus)]

        # Then check for locus ID prefixes
        locus_prefixes = {
            'FLI': 'flag_i',
            'FLII': 'flag_ii',
            'FLIII': 'flag_iii',
            'FLIV': 'flag_iv',
            'TTI': 't3ss_i',
            'TTII': 't3ss_ii',
            'TFI': 't6ss_i',
            'TFII': 't6ss_ii',
            'KL': 'capsule',
            'CL': 'cellulose',
            'OL': 'lps',
            'SR': 'sorbitol'
        }

        locus_str = str(locus).upper()
        for prefix, locus_type in locus_prefixes.items():
            if locus_str.startswith(prefix):
                return locus_type
        return ''

    def determine_type(self, present_genes, truncated_genes, extra_genes, locus_coverage, different_genes):
        # Perfect - Everything matches exactly
        if not different_genes and present_genes == self.number_of_genes and truncated_genes == 0 and extra_genes == 0 and locus_coverage >= 0.99:
            return 'Perfect'

        num_different_genes = len(different_genes)

        # Very High - Almost perfect with minimal differences
        if num_different_genes <= 2 and locus_coverage >= 0.97 and extra_genes <= 1:
            return 'Very High'

        # High - Good similarity with some variations
        if num_different_genes <= 4 and locus_coverage >= 0.95 and extra_genes <= 2:
            return 'High'

        # Good - Recognizable system with acceptable variations
        if num_different_genes <= 6 and locus_coverage >= 0.90 and extra_genes <= 3:
            return 'Good'

        # Low - System present but with significant variations
        if num_different_genes <= 8 and locus_coverage >= 0.80 and extra_genes <= 4:
            return 'Low'

        # Very low - System barely recognizable
        if locus_coverage >= 0.70 or num_different_genes <= 10:  # This condition is fine as is
            return 'Very low'

        # Default case
        return 'Very low'

    def select_best_type(self, type_report):
        type_hierarchy = ['None', 'Low', 'Good', 'High', 'Very High', 'Perfect']
        best_locus = None
        best_metrics = {'coverage': 0, 'present_genes': 0}

        for idx, row in type_report.iterrows():
            if best_locus is None or (
                    type_hierarchy.index(row['Type']) == type_hierarchy.index(best_locus['Type']) and
                    (row['Locus Coverage'] > best_metrics['coverage'] or
                     row['Present Genes'] > best_metrics['present_genes'])
            ):
                best_locus = row
                best_metrics = {
                    'coverage': row['Locus Coverage'],
                    'present_genes': row['Present Genes']
                }

        return str(best_locus['Locus']), str(best_locus['Type'])


if __name__ == "__main__":
    import glob  # Add this import at the top of the file

    # Define paths
    base_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/test-data/rubus2/more_testing'
    reference_types = "/Users/josediogomoura/Documents/BioFago/BioFago/reference_types_database/T3SS_II/types_T3SS_II.gb"
    prokka_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/test-data/rubus2/more_testing/prokka'

    # Debug: Check if prokka folder exists and list its contents
    print(f"Checking prokka folder: {prokka_folder}")
    if os.path.exists(prokka_folder):
        print("Prokka folder exists")
        print("Contents of prokka folder:")
        print(os.listdir(prokka_folder))
    else:
        print("Prokka folder does not exist!")
        exit(1)

    # Create required folders
    db_folder = os.path.join(base_folder, "db_folder")
    results_folder = os.path.join(base_folder, "results")
    os.makedirs(db_folder, exist_ok=True)
    os.makedirs(results_folder, exist_ok=True)

    # Try to find proteins.faa file
    print("\nLooking for proteins.faa file...")
    faa_pattern = os.path.join(prokka_folder, '*.faa')
    matching_files = glob.glob(faa_pattern)
    print(f"Found .faa files: {matching_files}")

    proteins_faa = get_prokka_faa_file(prokka_folder)
    if not proteins_faa or not os.path.exists(proteins_faa):
        raise FileNotFoundError(f"Could not find proteins.faa file in {prokka_folder}")
    print(f"Using proteins.faa from: {proteins_faa}")

    # Define output files
    final_csv_file = os.path.join(base_folder, "final.csv")
    output_csv = os.path.join(base_folder, "loci_gene_presence.csv")
    input_csv = os.path.join(base_folder, "loci.csv")

    # Initialize BlastProteinv2
    blast = BlastProteinv2(reference_types=reference_types, db_folder=db_folder, results_folder=results_folder)

    # Create CSV of loci information
    blast.create_csv_loci(input_csv)

    # Process CSV and perform BLAST searches
    blast.process_csv(input_csv, proteins_faa, final_csv_file)

    # Analyze locus information
    gene_presence = TypeAnalysis(input_csv=final_csv_file, output_csv=output_csv, proteins_faa=proteins_faa)
    df_analyzed = gene_presence.analyze_locus()
    type_report, final_type_locus, final_type, flagged_genes, locus_info = gene_presence.assign_type(df_analyzed)

    # Print results
    print("\nType Report:")
    print(type_report)
    print(f"\nFinal Type Locus: {final_type_locus}")
    print(f"Final Type: {final_type}")
    print(f"Flagged Genes: {flagged_genes}")
    print("\nLocus Information:")
    for locus_type, info in locus_info.items():
        print(f"{locus_type}: {info}")