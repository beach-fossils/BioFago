import gzip
import logging
import shutil
import ssl
from ftplib import FTP, error_perm
import certifi
from Bio import Entrez
import os

# Setup SSL context
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())

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

    :param email: str - The email address to be used for Entrez API requests.
    :param api_key: str - The API key for accessing the Entrez API.
    """

    if email is not None:
        Entrez.email = email
    if api_key is not None:
        Entrez.api_key = api_key


def download_ftp_file(ftp_server: str, ftp_dir: str, file_name: str, output_path: str):
    """
    Download a file from an FTP server.

    :param ftp_server: str - Address of the FTP server.
    :param ftp_dir: str - Directory on the FTP server where the file is located.
    :param file_name: str - Name of the file to download.
    :param output_path: str - Local path to save the downloaded file.
    """

    try:
        with FTP(ftp_server) as ftp:
            ftp.login()
            ftp.cwd(ftp_dir)
            with open(output_path, 'wb') as f:
                ftp.retrbinary(f'RETR {file_name}', f.write)
        logger.info(f"Downloaded {output_path}")
    except error_perm as e:
        logger.error(f"Error: {e}")
        if os.path.exists(output_path):
            os.remove(output_path)


def list_ftp_directory(ftp_server: str, ftp_dir: str) -> list:
    """
    List files in a specific directory on an FTP server.

    :param ftp_server: str - Address of the FTP server.
    :param ftp_dir: str - Directory on the FTP server to list.
    :return: list - List of file names in the given directory.
    """

    with FTP(ftp_server) as ftp:
        ftp.login()
        ftp.cwd(ftp_dir)
        files = ftp.nlst()
    return files


def unzip_file(file_path: str) -> str:
    """
    Unzip a gzip file.

    :param file_path: str - Path to the gzip file.
    :return: str - Path to the unzipped file.
    """

    output_path = file_path.replace('.gz', '')
    with gzip.open(file_path, 'rb') as f_in:
        with open(output_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(file_path)
    return output_path


def download_and_unzip_file(ftp_server: str, ftp_dir: str, file_name: str, output_path: str, max_retries: int = 3):
    """
    Download and unzip a file from an FTP server with retry logic.

    :param ftp_server: str - Address of the FTP server.
    :param ftp_dir: str - Directory on the FTP server where the file is located.
    :param file_name: str - Name of the file to download and unzip.
    :param output_path: str - Local path to save the downloaded and unzipped file.
    :param max_retries: int - Maximum number of retries for downloading the file.
    """

    for retry in range(max_retries):
        try:
            download_ftp_file(ftp_server, ftp_dir, file_name, output_path)
            unzip_file(output_path)
            return
        except Exception as e:
            logger.error(f"Error downloading or unzipping {file_name}. Retry {retry + 1}/{max_retries}: {e}")

    logger.error(f"Failed to download and unzip {file_name} after {max_retries} retries.")


def get_assembly_ids(search_term: str, email: str = None, api_key: str = None) -> list:
    """
    Retrieve a list of assembly IDs from the NCBI database based on a search term.

    :param search_term: str - Search term for querying the NCBI database.
    :param email: str - Email for NCBI Entrez.
    :param api_key: str - API key for NCBI Entrez.
    :return: list - List of assembly IDs.
    """

    setup_entrez(email, api_key)
    search_handle = Entrez.esearch(db="assembly", term=search_term, retmax='500')
    search_results = Entrez.read(search_handle)
    search_handle.close()
    assembly_ids = search_results['IdList']
    logger.info(f"Found {len(assembly_ids)} assemblies for '{search_term}'.")
    return assembly_ids


def download_genomes(assembly_ids: list, output_dir: str, formats: list) -> int:
    """
    Download genomes based on a list of assembly IDs from the NCBI Assembly database.

    :param assembly_ids: List[str] - A list of assembly IDs for which the genome data will be downloaded.
    :param output_dir: str - The directory path where the downloaded genome files will be stored.
    :param formats: List[str] - A list of file formats to be downloaded for each genome assembly.
                               Common formats include 'gbff' (GenBank flat file format), 'fasta', etc.
    :return: int - The number of successfully downloaded files. This count helps in tracking the number
                   of genomes successfully downloaded and identifying any issues in the download process.
    """

    download_count = 0
    for assembly_id in assembly_ids:
        summary_handle = Entrez.esummary(db="assembly", id=assembly_id)
        summary = Entrez.read(summary_handle)
        summary_handle.close()

        docsum = summary['DocumentSummarySet']['DocumentSummary'][0]
        ftp_path = docsum['FtpPath_GenBank']
        if ftp_path == "":
            logger.warning(f"No download link for assembly ID {assembly_id}.")
            continue

        ftp_server = ftp_path.split('/')[2]
        ftp_dir = '/'.join(ftp_path.split('/')[3:])
        files = list_ftp_directory(ftp_server, ftp_dir)

        for file in files:
            if "wgsmaster" not in file and any(file.endswith(f".{format}.gz") for format in formats):
                file_name = file
                output_path = os.path.join(output_dir, file_name)
                logger.info(f"Downloading {file_name}...")
                try:
                    download_and_unzip_file(ftp_server, ftp_dir, file_name, output_path)
                    download_count += 1
                except Exception as e:
                    logger.error(f"Failed to download {file_name}: {e}")

    return download_count


def search_and_download_genomes(search_term: str, output_dir: str, formats: list, email: str = None,
                                api_key: str = None):
    """
    Search for and download genome data from the NCBI Assembly database based on a given search term.

    :param search_term: str - A query term for searching the NCBI Assembly database. This could include specific
                              organism names, genome types, etc., formatted as a query string.
    :param output_dir: str - Directory path where the downloaded genome files will be stored.
    :param formats: list - A list of file formats to be downloaded for each genome assembly. Common formats include
                           'gbff' (GenBank flat file format), 'fna', etc. It is required to specify the formats.
    :param email: str, optional - Email address associated with the NCBI account. This is required for utilizing the
                                  Entrez programming utilities.
    :param api_key: str, optional - API key for the NCBI Entrez system. Using an API key increases the permitted rate of
                                    requests to the NCBI servers.
    :return: int - The number of successfully downloaded files. This can be used to check how many genome files were
                   downloaded and potentially identify if any assemblies were missed.
    """

    assembly_ids = get_assembly_ids(search_term, email, api_key)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if formats is None:
        formats = ['gbff']
    download_count = download_genomes(assembly_ids, output_dir, formats)
    logger.info(f"Download and unzip process completed. Total successful downloads: {download_count}")


def main():
    """
    Main function to execute the script.
    It parses command-line arguments and triggers the genome downloading process.
    """

    # Constants
    DEFAULT_EMAIL = 'email@example.com'
    #DEFAULT_OUTPUT_DIR = 'your_output_path'
    #DEFAULT_API_KEY = 'your_api_key'
    DEFAULT_TERM = '"Erwinia amylovora"[Organism] AND (latest[filter] AND all[filter] NOT anomalous[filter])'

    import argparse
    parser = argparse.ArgumentParser(description='Download genomes from NCBI Assembly.')
    parser.add_argument('--search_term', type=str, help='Search term for NCBI Assembly', default=DEFAULT_TERM)
    parser.add_argument('--email', type=str, help='Email for NCBI Entrez', default=DEFAULT_EMAIL)
    parser.add_argument('--output_dir', type=str, help='Output directory', default=DEFAULT_OUTPUT_DIR)
    parser.add_argument('--formats', nargs='*', help='Formats to download', default=['gbff'], type=list)
    parser.add_argument('--api_key', type=str, help='Entrez API key', default=DEFAULT_API_KEY)

    args = parser.parse_args()

    search_and_download_genomes(args.search_term, args.output_dir, args.formats, args.email, args.api_key)


if __name__ == "__main__":
    main()

# Run the script with the following command:
# python get_genomes.py --search_term '"Erwinia amylovora"[Organism] AND (latest[filter] AND all[filter] NOT anomalous[filter])'
# --output_dir '/Users/josediogomoura/Desktop/BioFago/github/data/output/erwinia_gbff_genomes'
