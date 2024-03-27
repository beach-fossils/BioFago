import csv
import subprocess
from pathlib import Path
from Bio import SeqIO
import tempfile
from src.PostRoaryPlotter import GenePresenceAbsencePlotter
import pandas as pd
import logging


# setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SpeciesTypePlotter:
    def __init__(self, data_path, json_path, output_csv_path):
        self.data_path = data_path
        self.json_path = json_path
        self.output_csv_path = output_csv_path

    def load_and_process_data(self):
        plotter = GenePresenceAbsencePlotter(self.data_path)
        plotter.load_and_process_data()
        plotter.process_strain_names()  # Simplify the strain names for JSON search
        plotter.identify_unique_patterns()

        # Generate the heatmap
        plotter.plot_heatmap()

        plotter.make_final_csv(self.output_csv_path, self.json_path)
        logger.info(f"Final CSV created at {self.output_csv_path}")

    @staticmethod
    def append_pattern_column_to_csv_sorted_by_frequency(input_csv_path, output_csv_path, prefix):
        df = pd.read_csv(input_csv_path)
        patterns = df.iloc[:, 2:].apply(lambda row: ''.join(row.astype(str)), axis=1)
        pattern_counts = patterns.value_counts()
        sorted_patterns = pattern_counts.index.tolist()
        pattern_to_id = {pattern: f"{prefix}{str(index + 1).zfill(2)}" for index, pattern in enumerate(sorted_patterns)}
        df['type'] = patterns.map(pattern_to_id)
        df.to_csv(output_csv_path, index=False)
        logger.info(f"Modified CSV saved to {output_csv_path}")


class ProkkaCDSExtractor:
    def __init__(self, input_csv, prokka_base_path, output_csv):
        self.input_csv = input_csv
        self.prokka_base_path = Path(prokka_base_path)
        self.output_csv = output_csv

    def extract_strain_types(self):
        """Extract unique strain types and their corresponding accession number."""
        strain_types = {}
        with open(self.input_csv, mode='r') as infile:
            reader = csv.DictReader(infile)
            for row in reader:
                strain_type = row['type']
                if strain_type not in strain_types:
                    strain_types[strain_type] = row['accession']
        return strain_types

    def find_genbank_path(self, accession):
        """Find the GenBank file path for the given accession."""
        for path in self.prokka_base_path.glob(f"**/*{accession}*/*.gbk"):
            return path
        return None

    def extract_cds_sequences(self, filepath):
        """Extract all CDS sequences from a GenBank file."""
        cds_info = []
        with open(filepath, 'r') as genbank_file:
            for record in SeqIO.parse(genbank_file, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        gene_or_locus = \
                            feature.qualifiers.get('gene', [feature.qualifiers.get('locus_tag', ['NA'])[0]])[0]
                        sequence = feature.qualifiers.get('translation', ['NA'])[0]
                        cds_info.append((gene_or_locus, sequence))
        return cds_info

    def write_to_csv(self, data):
        """Write extracted data to CSV, sorted by type."""
        sorted_data = sorted(data, key=lambda x: x[0])  # Sort by strain type
        with open(self.output_csv, 'w', newline='') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['Type', 'Gene/Locus Tag', 'Sequence'])
            writer.writerows(sorted_data)

    def run(self):
        all_cds = []
        strain_types = self.extract_strain_types()
        for strain_type, accession in strain_types.items():
            genbank_path = self.find_genbank_path(accession)
            if genbank_path:
                cds_sequences = self.extract_cds_sequences(genbank_path)
                for gene_or_locus, sequence in cds_sequences:
                    all_cds.append((strain_type, gene_or_locus, sequence))
        self.write_to_csv(all_cds)


class DiamondRunner:
    def __init__(self, csv_path, diamond_db, output_path, diamond_executable='diamond', fasta_db_path=None,
                 base_output_dir="/desired/path/for/results"):
        self.csv_path = csv_path
        self.diamond_db = Path(base_output_dir) / diamond_db  # Path object to handle file paths
        self.output_path = Path(base_output_dir) / output_path
        self.diamond_executable = diamond_executable
        self.fasta_db_path = fasta_db_path
        self.diamond_output_dir = Path(base_output_dir) / "diamond_outputs"
        self.diamond_output_dir.mkdir(parents=True, exist_ok=True)  # Ensure the directory exists, create any necessary parents

    def create_diamond_db(self):
        if self.fasta_db_path is None:
            raise ValueError("Path to fasta database file must be provided for database creation.")
        subprocess.run([self.diamond_executable, 'makedb', '--in', self.fasta_db_path, '-d', str(self.diamond_db)],
                       check=True)
        logger.info(f"DIAMOND database created at {self.diamond_db}.dmnd")

    def read_sequences_from_csv(self, prefix):
        with open(self.csv_path, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if prefix in row['Gene/Locus Tag']:  # Use the user-defined prefix
                    yield row['Type'], row['Gene/Locus Tag'], row['Sequence']

    def run_diamond(self, sequence, type_, gene):
        diamond_output_file = self.diamond_output_dir / f"{type_}_{gene}_diamond_output.txt"
        with tempfile.NamedTemporaryFile(mode='w+', delete=True) as temp_fasta:
            temp_fasta.write(f">temp_sequence\n{sequence}")
            temp_fasta.flush()
            subprocess.run([
                'diamond', 'blastp', '-d', self.diamond_db, '-q', temp_fasta.name, '-o', diamond_output_file,
                '-f', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart',
                'send', 'evalue', 'bitscore',
                '-p', '8', '--fast'
            ], check=True)

        return str(diamond_output_file)

    def process_sequences(self, prefix='KL_'):  # Default prefix set to 'KL_', can be overridden by user
        results = []
        for type_, gene, sequence in self.read_sequences_from_csv(prefix):  # Use the prefix in the call
            diamond_output_file = self.run_diamond(sequence, type_, gene)
            with open(diamond_output_file, 'r') as file:
                diamond_output = file.read()
            best_match = next(csv.reader(diamond_output.splitlines(), delimiter='\t'), None)
            if best_match:
                results.append([
                    type_, gene, best_match[1], best_match[2], best_match[3], best_match[4], best_match[5],
                    best_match[6], best_match[7], best_match[8], best_match[9], best_match[10], best_match[11],
                    str(diamond_output_file)  # Including DIAMOND output file path if needed
                ])
            else:
                results.append([type_, gene, 'No match', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', ''])

        # Write the results to the output CSV, including headers for the new metrics
        with open(self.output_path, 'w', newline='') as outcsv:
            writer = csv.writer(outcsv)
            writer.writerow([
                'Type', 'Gene/Locus Tag', 'Best Match', 'Percentage Identity', 'Alignment Length', 'Mismatches',
                'Gap Openings', 'Query Start', 'Query End', 'Subject Start', 'Subject End', 'E-value', 'Bit Score',
                'DIAMOND Output File'
            ])
            writer.writerows(results)



def DIAMONDrunner():
    protein_csv_path = "/Users/josediogomoura/Documents/BioFago/BioFago/data/tests/protein_sequences.csv"
    diamond_db = "erwinia_loci_db"
    output_path = "/Users/josediogomoura/Documents/BioFago/BioFago/data/tests/output.csv"  # Adjusted to be relative to the base_output_dir
    fasta_db_path = "/Users/josediogomoura/Documents/BioFago/BioFago/data/input/loci_ref_sequences/compiled_fasta/compiled_loci.fasta"
    base_output_dir = "/Users/josediogomoura/Documents/BioFago/BioFago/data/tests"  # Base directory for all outputs

    # Initialize DiamondRunner with the base_output_dir parameter
    diamond_runner = DiamondRunner(protein_csv_path, diamond_db, output_path, fasta_db_path=fasta_db_path, base_output_dir=base_output_dir)
    diamond_runner.create_diamond_db()  # Create the database before processing sequences
    diamond_runner.process_sequences(prefix='KL_')

def main_workflow():
    # Paths configuration
    gene_presence_absence_csv = '/Volumes/Crucial_X9/BioFago/data/ApproachFlankGenes/capsule/roary_90_1710520453/gene_presence_absence.csv'
    json_path = '/Volumes/Crucial_X9/BioFago/data/genomesAllErwinia/ncbi_dataset/ncbi_dataset/data/assembly_data_report.jsonl'
    species_to_defined_type_csv_initial = '/Users/josediogomoura/Documents/BioFago/BioFago/data/tests/accession_to_species_to_pattern_90.csv'
    species_to_defined_type_csv_final = '/Users/josediogomoura/Documents/BioFago/BioFago/data/tests/species_to_defined_type_90.csv'
    prokka_base_path = '/Volumes/Crucial_X9/BioFago/data/ALLApproachFlankGenes/capsule/prokka'
    output_csv_for_prokka = '/Users/josediogomoura/Documents/BioFago/BioFago/data/tests/protein_sequences.csv'


    # Step 1: Load, process data, and create initial CSV
    species_plotter = SpeciesTypePlotter(gene_presence_absence_csv, json_path, species_to_defined_type_csv_initial)
    species_plotter.load_and_process_data()

    # Step 1.5: Append pattern column to CSV, sorted by frequency (Mandatory step)
    SpeciesTypePlotter.append_pattern_column_to_csv_sorted_by_frequency(species_to_defined_type_csv_initial, species_to_defined_type_csv_final, prefix='CL')

    # Step 2: Extract sequences with ProkkaCDSExtractor
    extractor = ProkkaCDSExtractor(species_to_defined_type_csv_final, prokka_base_path, output_csv_for_prokka)
    extractor.run()


if __name__ == "__main__":
    main_workflow()

    #DIAMONDrunner()
