import logging
from Bio import Entrez
import os

# Configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler = logging.StreamHandler()
handler.setFormatter(formatter)
logger.addHandler(handler)


def setup_entrez(email: str, api_key: str):
    """
    Configure the Entrez API with the provided email and API key.
    """
    Entrez.email = email
    Entrez.api_key = api_key


def fetch_nucleotide_sequences(sequence_ids: list, output_dir: str, format: str = 'gb'):
    """
    Fetch nucleotide sequences from NCBI and write them to files.
    """
    for sequence_id in sequence_ids:
        try:
            handle = Entrez.efetch(db="nucleotide", id=sequence_id, rettype=format, retmode="text")
            sequence_data = handle.read()
            handle.close()

            # Define the output file path
            output_path = os.path.join(output_dir, f"{sequence_id}.{format}")

            # Write the sequence data to a file
            with open(output_path, 'w') as file:
                file.write(sequence_data)
            logger.info(f"Downloaded: {output_path}")
        except Exception as e:
            logger.error(f"Failed to download sequence {sequence_id}: {e}")


def get_nucleotide_ids(search_term: str) -> list:
    """
    Retrieve a list of nucleotide sequence IDs from NCBI based on a search term.
    """
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax='500')
    search_results = Entrez.read(handle)
    handle.close()
    return search_results['IdList']


def main():
    """
    Main function to execute the script.
    """
    # Constants
    EMAIL = 'example@gmail.com'
    API_KEY = 'None'
    OUTPUT_DIR = '/Users/josediogomoura/Desktop/BioFago/github/data/output/phages'
    FORMAT = 'gb'
    SEARCH_TERM = '(("Erwinia amylovora"[Organism] OR Erwinia Amylovora[All Fields]) AND phage[All Fields]) AND viruses[filter]'

    # Set up Entrez
    setup_entrez(EMAIL, API_KEY)

    # Make sure output directory exists
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # Get nucleotide sequence IDs
    sequence_ids = get_nucleotide_ids(SEARCH_TERM)
    logger.info(f"Found {len(sequence_ids)} sequences")

    # Fetch and save the nucleotide sequences
    fetch_nucleotide_sequences(sequence_ids, OUTPUT_DIR, FORMAT)


if __name__ == "__main__":
    main()
