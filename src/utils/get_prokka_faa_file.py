import os


def get_prokka_faa_file(prokka_output_folder):
    """
    Function to get the .faa file from the Prokka output folder.
    :param prokka_output_folder: The folder where Prokka output files are stored.
    :return: Path to the .faa file.
    """
    for root, dirs, files in os.walk(prokka_output_folder):
        for file in files:
            if file.endswith(".faa") and file.startswith("PROKKA_") and "tmp" not in file:
                return os.path.join(root, file)
    return None

if __name__ == "__main__":
    path = '/Users/josediogomoura/Documents/BioFago/BioFago/test-data/T6SS/cluster_1/extracted_seq'

    print(get_prokka_faa_file(path)) # Expected: /Users/josediogomoura/Documents/BioFago/BioFago/data/assign_types/cellulose/test_1/prokka/PRR1_INIAV_Contig_2_consensus_sequence_extracted_sequence/PROKKA_09232021.faa