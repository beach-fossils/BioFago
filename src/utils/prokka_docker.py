import os
import subprocess
from pathlib import Path


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
        subprocess.run(['docker', 'ps'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        print("Error: Docker is not running. Please start Docker and try again.")
        return
    except FileNotFoundError:
        print("Error: Docker is not installed or not found in the system PATH.")
        return

    # Get the filename without extension for output folder name
    fasta_file_name = Path(fasta_file).stem
    output_path = os.path.join(base_output_folder, fasta_file_name)
    Path(output_path).mkdir(parents=True, exist_ok=True)  # Create the output folder

    # Configure the Prokka command for the specific fasta file
    prokka_command = [
        'docker', 'run', '--rm',
        '--platform', 'linux/amd64',
        '-v', f'{Path(fasta_file).parent}:/data',
        '-v', f'{output_path}:/output',
        '-v', f'{custom_db_path}:/custom_db/sequence.gb',
        'staphb/prokka:latest',
        'prokka', f'/data/{Path(fasta_file).name}', '-o', '/output', '--locustag', locus_tag_prefix
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
        print(e.stderr)


def main():
    # Define the genomes and output directories
    fasta_file = '/Users/josediogomoura/Documents/BioFago/BioFago/data/assign_types/cellulose/test_1/extracted_region/PRR1_INIAV.fasta'

    base_output_folder = '/Users/josediogomoura/Documents/BioFago/BioFago/data/assign_types/cellulose/test_1/prokka'
    custom_db_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/assign_types/fully_gb_database/cellulose/curated_cellulose.gb'
    locus_tag_prefix = 'PREFIX'

    # Run Prokka annotation using Docker
    run_prokka_docker(fasta_file, base_output_folder, custom_db_path, locus_tag_prefix)


if __name__ == "__main__":
    main()
