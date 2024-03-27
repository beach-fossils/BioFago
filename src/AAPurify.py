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
        self.updated_csv_path = self.output_dir / "updated_target.csv"
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

        # Open the target CSV for reading and the updated CSV for writing
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


def main():
    diamond_csv_path = "/Users/josediogomoura/Documents/BioFago/BioFago/data/PostRoaryForALLErwinia/capsule/90_blastp/3_output_diamond/protein_diamond_output.csv"
    target_csv_path = "/Users/josediogomoura/Documents/BioFago/BioFago/data/PostRoaryForALLErwinia/capsule/90_blastp/2_csv_with_aasequences/protein_sequences.csv"
    output_dir = "/Users/josediogomoura/Documents/BioFago/BioFago/data/tests/Purify"

    updater = LocusTagUpdater(diamond_csv_path, target_csv_path, output_dir=output_dir)
    updater.read_diamond_csv()  # Process and optionally export for curation
    # Once curation is done (manually updating the JSON), call to finalize updates:
    updater.apply_updates_to_csv(
        curated_json_path='/Users/josediogomoura/Documents/BioFago/BioFago/data/tests/Purify/curated_genes.json')


if __name__ == '__main__':
    main()
