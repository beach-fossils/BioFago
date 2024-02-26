import os
import subprocess

def download_genomes(search_term: str, output_path: str, format: str):
    search_term_clean = search_term.replace(' ', '_')
    output_dir = os.path.join(output_path, 'genomes_' + search_term_clean.lower())
    os.makedirs(output_dir, exist_ok=True)

    datasets_command = ["/Users/josediogomoura/bin/datasets",
                        "download",
                        "genome",
                        "taxon",
                        search_term,
                        "--include",
                        format,
                        "--dehydrated",
                        "--filename",
                        os.path.join(output_dir, search_term_clean + "_dehydrated.zip")]

    subprocess.run(datasets_command, check=True)


def unzip_genomes(output_path: str, search_term: str):
    search_term_clean = search_term.replace(' ', '_')
    output_dir = os.path.join(output_path, 'genomes_' + search_term_clean.lower())

    unzip_command = ["unzip",
                     os.path.join(output_dir, search_term_clean + "_dehydrated.zip"),
                     "-d",
                     output_dir]

    subprocess.run(unzip_command, check=True)


def rehydrate_data(output_path: str, search_term: str):
    search_term_clean = search_term.replace(' ', '_')
    output_dir = os.path.join(output_path, 'genomes_' + search_term_clean.lower())

    rehydrate_command = ["/Users/josediogomoura/bin/datasets",
                         "rehydrate",
                         "--directory",
                         output_dir]

    subprocess.run(rehydrate_command, check=True)


def main():
    search_term = "Erwinia amylovora"
    output_path = '/Users/josediogomoura/Documents/BioFago/BioFago/data/all_genomes_erwinia'
    format = "gbff"
    download_genomes(search_term, output_path, format)
    unzip_genomes(output_path, search_term)
    rehydrate_data(output_path, search_term)

if __name__ == "__main__":
    main()