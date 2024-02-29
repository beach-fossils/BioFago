import subprocess
import logging


def setup_logging(file_path):
    """
    Sets up logging to file, ensuring it saves in the results folder.
    """
    log_format = '%(asctime)s %(message)s'
    logging.basicConfig(filename=file_path, level=logging.DEBUG, format=log_format, filemode='w')

    # Check if there are handlers already to avoid adding a new handler multiple times
    if not logging.getLogger().hasHandlers():
        logging.getLogger().addHandler(logging.FileHandler(file_path))


def download_genome_dataset(keyword, output_directory):
    command = [
        "datasets",
        "download",
        "genome",
        "taxon",
        keyword,
        "--filename",
        f"{output_directory}/ncbi_dataset.zip"
    ]
    try:
        subprocess.run(command, check=True)
        print("Genome dataset downloaded successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading dataset: {e}")


def unzip_genome_dataset(output_directory):
    # Unzip the downloaded file
    command = [
        "unzip",
        f"{output_directory}/ncbi_dataset.zip",
        "-d",
        output_directory
    ]

    try:
        subprocess.run(command, check=True)
        print("Genome dataset unzipped successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Error unzipping dataset: {e}")


def main():
    # Example usage:
    keyword = "Erwinia"
    output_directory = "/Users/josediogomoura/Documents/BioFago/BioFago/data/all_genomes_erwinia/fasta"
    log_file = "/Users/josediogomoura/Documents/BioFago/BioFago/data/all_genomes_erwinia/fasta/download.log"

    setup_logging(log_file)
    download_genome_dataset(keyword, output_directory)
    unzip_genome_dataset(output_directory)


if __name__ == "__main__":
    main()
