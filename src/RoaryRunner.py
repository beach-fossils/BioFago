import subprocess
import os
import logging
from pathlib import Path
import shutil


# Set up basic logging
# Configure logging to output to a file
def configure_logging(log_file):
    logging.basicConfig(
        filename=log_file,
        filemode='a',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO
    )


class RoaryRunner:
    def __init__(self, prokka_folder, roary_output_folder, roary_path='roary', threads=8):
        self.prokka_folder = Path(prokka_folder)
        self.roary_output_folder = Path(roary_output_folder)
        self.roary_path = roary_path
        self.threads = threads
        self.gff_files = []

        # Set up the log file
        log_file = self.roary_output_folder / 'roary_runner.log'
        configure_logging(log_file)
        logging.info(f"RoaryRunner initialized. Log file created at {log_file}")

    def find_gff_files(self):

        # Walk through the Prokka output folders and collect all GFF files
        for folder in self.prokka_folder.iterdir():
            if folder.is_dir():
                gff_files = list(folder.glob("*.gff"))
                for gff_file in gff_files:
                    # Create a new filename by combining the folder name and file stem, with version
                    folder_name_with_version = folder.name.replace('_', '.')
                    new_filename = f"{folder_name_with_version}_{gff_file.stem}.gff"
                    new_filepath = self.roary_output_folder / new_filename

                    # Copy the GFF file to the new folder
                    shutil.copy(gff_file, new_filepath)
                    self.gff_files.append(new_filepath)
                    logging.info(f"Copied {gff_file.name} to {new_filename}")
        logging.info(f"Found and copied a total of {len(self.gff_files)} GFF files.")

    def run_roary(self):
        # Ensure the Roary output folder exists
        self.roary_output_folder.mkdir(parents=True, exist_ok=True)

        # Run Roary if more than one GFF file is found
        if len(self.gff_files) > 1:
            gff_paths = ' '.join(str(gff) for gff in self.gff_files)
            roary_cmd = f"{self.roary_path} -f {self.roary_output_folder} -e -n -p {self.threads} -v {gff_paths}"
            logging.info(f"Running Roary: {roary_cmd}")
            try:
                subprocess.run(roary_cmd, shell=True, check=True, text=True)
                logging.info("Roary analysis completed successfully.")
            except subprocess.CalledProcessError as e:
                logging.error(f"Roary failed with error: {e.stderr}")
        else:
            logging.error("Not enough GFF files to run Roary. At least 2 are required.")


def main():
    prokka_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/prokka/LPS_locus_13genes_FN434113'
    roary_output_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/roary/LPS_locus_13genes_FN434113'

    roary_runner = RoaryRunner(prokka_folder, roary_output_folder)
    roary_runner.find_gff_files()
    roary_runner.run_roary()


if __name__ == "__main__":
    main()