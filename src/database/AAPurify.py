import csv
import json
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class LocusTagUpdater:
    def __init__(self, diamond_csv_path, target_csv_path, prefix="KL_", output_dir="output"):
        self.diamond_csv_path = diamond_csv_path
        self.target_csv_path = target_csv_path
        self.prefix = prefix
        self.output_dir = Path(output_dir)
        self.updated_csv_path = self.output_dir / "updated_target.csv_output"
        self.json_path = self.output_dir / "curated_genes.json"
        self.gene_updates = {
            "Curated genes": {},
            "Not possible to fix": []  # Initialized as a list
        }

        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def read_diamond_csv(self):
        with open(self.diamond_csv_path, mode='r', encoding='utf-8') as file:
            reader = csv.DictReader(file)
            for row in reader:
                if self.prefix in row['Gene/Locus Tag']:
                    best_match_details = row['Best Match'].split('|') if row['Best Match'] != "No match" else [
                        "No match"]
                    gene_name = best_match_details[1] if len(best_match_details) > 1 and best_match_details[
                        1] != "N/A" else None
                    ncbi_ref = best_match_details[2] if len(best_match_details) > 2 else None

                    if gene_name:
                        if row['Type'] not in self.gene_updates["Curated genes"]:
                            self.gene_updates["Curated genes"][row['Type']] = []
                        self.gene_updates["Curated genes"][row['Type']].append({
                            "Locus Tag": row['Gene/Locus Tag'],
                            "Gene name": gene_name
                        })
                    else:
                        # Adjust here to include entries with "No match" explicitly
                        no_match_info = {"Type": row['Type'], "Locus Tag": row['Gene/Locus Tag']}
                        if ncbi_ref:  # If there's an NCBI reference, include it
                            no_match_info["NCBI reference"] = ncbi_ref
                        else:  # Otherwise, explicitly note the absence of a match
                            no_match_info["Note"] = "No match found"

                        self.gene_updates["Not possible to fix"].append(no_match_info)

                    logger.info(f"Processed {row['Gene/Locus Tag']} for update.")
        self.export_for_curation()

    def export_for_curation(self):
        """Export gene updates for curation to a JSON file."""
        # Dynamically set JSON path here for export
        json_path = self.output_dir / "curated_genes.json"
        with open(json_path, 'w', encoding='utf-8') as json_file:
            json.dump(self.gene_updates, json_file, ensure_ascii=False, indent=4)
        logger.info(f"Exported gene updates for curation to {json_path}.")

    def apply_updates_to_csv(self, curated_json_path):
        # Load the curated updates from the JSON file
        with open(curated_json_path, 'r', encoding='utf-8') as json_file:
            curated_updates = json.load(json_file)

        # Open the genomes CSV for reading and the updated CSV for writing
        with open(self.target_csv_path, mode='r', encoding='utf-8') as infile, \
                open(self.updated_csv_path, mode='w', encoding='utf-8', newline='') as outfile:
            reader = csv.DictReader(infile)
            writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
            writer.writeheader()

            for row in reader:
                # Flag to check if an update has been applied
                update_applied = False

                # Iterate over each type and its updates in the curated JSON
                for type_, updates in curated_updates["Curated genes"].items():
                    # Check if the row's type matches and if there is an update for its locus tag
                    if row['Type'] == type_:
                        for update in updates:
                            if row['Gene/Locus Tag'] == update["Locus Tag"]:
                                # Update the Gene/Locus Tag column with the curated gene name
                                row['Gene/Locus Tag'] = update["Gene name"]
                                update_applied = True
                                break

                    # Break the outer loop if an update has been applied to avoid unnecessary iterations
                    if update_applied:
                        break

                # Write the potentially updated row to the new CSV
                writer.writerow(row)

        logger.info(f"Updated target CSV based on curated JSON. Output file: {self.updated_csv_path}")


import csv


class FindGeneDivision:
    def __init__(self, csv_path, output_path=None):
        self.csv_path = csv_path
        self.output_path = output_path if output_path is not None else 'output.csv_output'
        # This structure now will be a dictionary of dictionaries of lists to handle types, gene names, and sequences
        self.data = {}

    def process_csv(self):
        with open(self.csv_path, mode='r', encoding='utf-8') as file:
            reader = csv.DictReader(file)
            for row in reader:
                type = row['Type']
                gene_name = row['Gene/Locus Tag']
                sequence = row['Sequence']

                if type not in self.data:
                    self.data[type] = {}
                if gene_name not in self.data[type]:
                    self.data[type][gene_name] = []

                # Add sequence only if it's not already listed for this gene to avoid duplicates
                if sequence not in self.data[type][gene_name]:
                    self.data[type][gene_name].append(sequence)

    def rename_divisions(self):
        for type in self.data:
            for gene_name in list(self.data[type].keys()):  # list() to duplicate keys for safe iteration
                sequences = self.data[type][gene_name]
                if len(sequences) > 1:
                    # More than one sequence for this gene name, so we rename
                    for i, sequence in enumerate(sequences):
                        new_gene_name = f"{gene_name}_{i + 1}"
                        if new_gene_name not in self.data[type]:
                            self.data[type][new_gene_name] = [sequence]
                        else:
                            self.data[type][new_gene_name].append(sequence)
                    del self.data[type][gene_name]  # Remove the original entry after renaming

    def export_divisions(self):
        with open(self.output_path, mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerow(["Type", "Gene/Locus Tag", "Sequence"])
            for type, genes in self.data.items():
                for gene_name, sequences in genes.items():
                    for sequence in sequences:
                        writer.writerow([type, gene_name, sequence])

    def run(self):
        self.process_csv()
        self.rename_divisions()
        self.export_divisions()



def main():
    # part 1
    diamond_csv_path = "/data/PostRoaryForALLErwinia/cellulose/99_blastp/3_output_diamond/protein_diamond_output.csv_output"
    target_csv_path = "/data/PostRoaryForALLErwinia/cellulose/99_blastp/2_csv_with_aasequences/protein_sequences.csv_output"
    output_dir = "/data/AApurify/cellulose/99_blastp"
    prefix = "CL_"

    updater = LocusTagUpdater(diamond_csv_path, target_csv_path, prefix=prefix, output_dir=output_dir)
    # updater.read_diamond_csv()  # Process and optionally export for curation
    # Once curation is done (manually updating the JSON), call to finalize updates: // maybe use genomes for this
    updater.apply_updates_to_csv(
         curated_json_path='/data/AApurify/cellulose/99_blastp/curated_genes.json')

    # part 2
    # now we check for gene division with the updated csv_output
    csv_path = "/data/AApurify/cellulose/99_blastp/updated_target.csv_output"
    output_path = "/data/AApurify/cellulose/99_blastp/division_genes.csv_output"
    gene_division = FindGeneDivision(csv_path, output_path)
    gene_division.run()


if __name__ == '__main__':
    main()
