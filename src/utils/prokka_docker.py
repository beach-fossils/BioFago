import os
import subprocess
from pathlib import Path
import shutil

def run_prokka_docker(fasta_file, base_output_folder, custom_db_path, locus_tag_prefix):
    # Additional Prokka options
    prokka_options = {
        '--kingdom': 'Bacteria',
        '--genus': 'Erwinia',
        '--species': 'amylovora',
        '--force': 'true'
    }

    # Check if Docker is running
    try:
        # Use 'sudo docker ps' if you truly need sudo for Docker listing.
        subprocess.run(['docker', 'ps'], shell=False, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        print("Error: Docker is not running, or 'sudo' permission is required. Please start Docker and try again.")
        return
    except FileNotFoundError:
        print("Error: Docker is not installed or not found in the system PATH.")
        return

    # Get the filename without extension for output folder name
    fasta_file_name = Path(fasta_file).stem
    output_path = os.path.join(base_output_folder, fasta_file_name)
    Path(output_path).mkdir(parents=True, exist_ok=True)  # Create the output folder

    # Configure the Prokka command for the specific fasta file
    # IMPORTANT: Prepend 'sudo' to the Docker command if your environment requires it.
    prokka_command = [
        'sudo', 'docker', 'run', '--rm',
        '--platform', 'linux/amd64',
        '-v', f'{Path(fasta_file).parent}:/reference_crispr',
        '-v', f'{output_path}:/output',
        '-v', f'{custom_db_path}:/custom_db/sequence.gb',
        'staphb/prokka:latest',
        'prokka', f'/reference_crispr/{Path(fasta_file).name}', '-o', '/output', '--locustag', locus_tag_prefix
    ]

    # Add custom database option if provided
    if custom_db_path:
        prokka_command.extend(['--proteins', '/custom_db/sequence.gb'])

    # Add additional options
    for option, value in prokka_options.items():
        prokka_command.extend([option, value])

    # Execute Prokka command
    try:
        subprocess.run(prokka_command, check=True)
        print(f"Prokka annotation for {fasta_file} completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: Prokka annotation for {fasta_file} failed with exit code {e.returncode}.")

def rename_gbk_files(base_folder):
    renamed_count = 0
    errors = []

    for root, dirs, files in os.walk(base_folder):
        for file in files:
            if file.endswith('.gbk'):
                old_path = os.path.join(root, file)
                parent_folder_name = os.path.basename(root)
                new_filename = f"{parent_folder_name}.gbk"
                new_path = os.path.join(root, new_filename)

                try:
                    shutil.move(old_path, new_path)
                    print(f"Renamed: {old_path} -> {new_path}")
                    renamed_count += 1
                except Exception as e:
                    errors.append(f"Error renaming {old_path}: {str(e)}")

    print(f"\nTotal files renamed: {renamed_count}")
    if errors:
        print("\nErrors encountered:")
        for error in errors:
            print(error)

def copy_gbk_files(base_folder, destination_folder):
    copied_count = 0
    errors = []

    # Ensure the destination folder exists
    os.makedirs(destination_folder, exist_ok=True)

    for root, dirs, files in os.walk(base_folder):
        for file in files:
            if file.endswith('.gbk'):
                source_path = os.path.join(root, file)
                destination_path = os.path.join(destination_folder, file)

                try:
                    shutil.copy2(source_path, destination_path)
                    print(f"Copied: {source_path} -> {destination_path}")
                    copied_count += 1
                except Exception as e:
                    errors.append(f"Error copying {source_path}: {str(e)}")

    print(f"\nTotal files copied: {copied_count}")
    if errors:
        print("\nErrors encountered during copying:")
        for error in errors:
            print(error)

def main():
    # Define the genomes and output directories
    folder_with_fastas = '/Users/josediogomoura/Documents/BioFago/BioFago/test-data/flag3/extracted'
    base_output_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/test-data/flag3/annotated'
    locus_tag_prefix = 'PREFIX'
    custom_db_path = '/Users/josediogomoura/Documents/BioFago/BioFago/test-data/flag3/ref_region/PROKKA_11052024.gbk'

    # Get all .fasta files in the folder
    fasta_files = [f for f in os.listdir(folder_with_fastas) if f.endswith('.fasta')]

    # Iterate over each .fasta file and run Prokka
    for fasta_file in fasta_files:
        full_fasta_path = os.path.join(folder_with_fastas, fasta_file)
        print(f"Processing file: {fasta_file}")
        run_prokka_docker(
            fasta_file=full_fasta_path,
            base_output_folder=base_output_folder,
            locus_tag_prefix=locus_tag_prefix,
            custom_db_path=custom_db_path
        )

if __name__ == "__main__":
    # main()
    base_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/test-data/flag3/annotated'
    destination_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/test-data/flag3/annotated/all_gbk'

    rename_gbk_files(base_folder)
    copy_gbk_files(base_folder, destination_folder)
