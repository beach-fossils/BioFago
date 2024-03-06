import subprocess
import logging
import time


def setup_logging(file_path):
    """
    Sets up logging to file, ensuring it saves in the results folder.
    """
    log_format = '%(asctime)s %(message)s'
    logging.basicConfig(filename=file_path, level=logging.DEBUG, format=log_format, filemode='w')

    # Check if there are handlers already to avoid adding a new handler multiple times
    if not logging.getLogger().hasHandlers():
        logging.getLogger().addHandler(logging.FileHandler(file_path))


def download_genome_dataset_with_retry(keyword, output_directory, attempts=3, delay=5, annotated=True, assembly_source='refseq'):
    """
    Attempts to download the dataset with retries.
    :param keyword: Search keyword for the dataset.
    :param output_directory: Directory to save the downloaded file.
    :param attempts: Number of retry attempts.
    :param delay: Delay between attempts in seconds.
    :param annotated: Boolean indicating whether to limit to annotated genomes.
    :param assembly_source: String indicating whether to limit to RefSeq or GenBank genomes ('refseq', 'genbank', or 'all').
    """
    command = [
        "datasets",
        "download",
        "genome",
        "taxon",
        keyword,
        # "--assembly-level",
        # "complete",  # Filter for complete genome assemblies
        "--assembly-source",
        assembly_source,  # Filter for either 'RefSeq' (GCF_) or 'GenBank' (GCA_) genomes
    ]

    if annotated:
        command.append("--annotated")

    command.extend([
        "--filename",
        f"{output_directory}/ncbi_dataset.zip"
    ])

    for attempt in range(attempts):
        try:
            print(f"Attempt {attempt + 1} of {attempts}")
            subprocess.run(command, check=True)
            print("Genome dataset downloaded successfully!")
            break  # Break out of the loop if successful
        except subprocess.CalledProcessError as e:
            print(f"Error downloading dataset: {e}")
            if attempt < attempts - 1:
                print(f"Retrying in {delay} seconds...")
                time.sleep(delay)
            else:
                print("Maximum attempts reached. Download failed.")


# def unzip_genome_dataset(output_directory):
#     # Unzip the downloaded file
#     command = [
#         "unzip",
#         f"{output_directory}/ncbi_dataset.zip",
#         "-d",
#         output_directory
#     ]
#
#     try:
#         subprocess.run(command, check=True)
#         print("Genome dataset unzipped successfully!")
#     except subprocess.CalledProcessError as e:
#         print(f"Error unzipping dataset: {e}")


def main():
    # Example usage:
    keyword = "552"
    output_directory = "/Users/josediogomoura/Documents/BioFago/BioFago/data/ApproachFlankGenes/genomes/ErwiniaAmyl"
    log_file = "/Users/josediogomoura/Documents/BioFago/BioFago/data/ApproachFlankGenes/genomes/ErwiniaAmyl/download.log"
    annotated = True
    assembly_source = 'all'  # Use 'genbank' for GenBank or 'all' for both, 'refseq'

    setup_logging(log_file)
    download_genome_dataset_with_retry(keyword, output_directory, annotated=annotated, assembly_source=assembly_source)


if __name__ == "__main__":
    main()
