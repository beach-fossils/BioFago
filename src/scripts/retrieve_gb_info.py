import csv
import os
from Bio import SeqIO


def safe_gc_content(file_path):
    """Calculate GC content by reading the file, looking for the ORIGIN line, and counting from there."""
    with open(file_path, 'r') as file:
        origin_found = False
        gc_count = 0
        atgc_count = 0
        for line in file:
            if origin_found:
                # Remove whitespace and digits
                line = ''.join(filter(str.isalpha, line))
                gc_count += line.upper().count('G') + line.upper().count('C')
                atgc_count += len(line)
            elif line.startswith('ORIGIN'):
                origin_found = True
        if atgc_count > 0:
            return (gc_count / atgc_count) * 100
        else:
            return 'NA'


def extract_info_from_gb(file_path):
    """Extract required information from a GenBank file."""
    with open(file_path, 'r') as gb_file:
        for record in SeqIO.parse(gb_file, 'genbank'):
            gc_content = safe_gc_content(file_path)  # Calculate GC content from file
            size_bp = 'NA'  # Default to 'NA' if sequence is undefined
            if record.seq:  # If sequence is defined, update size_bp
                size_bp = len(record.seq)
            phage_name = record.description
            accession_number = record.id
            num_tRNA = sum(1 for feature in record.features if feature.type == 'tRNA')
            num_cds = sum(1 for feature in record.features if feature.type == 'CDS')
            family = 'NA'
            subfamily = 'NA'
            doi = 'NA'

            # Extract family and subfamily from taxonomy, if available
            for feature in record.features:
                if feature.type == 'source':
                    organism = feature.qualifiers.get('organism', [''])[0]
                    taxonomy = organism.split('; ')
                    if len(taxonomy) > 2:
                        family = taxonomy[-2]
                    if len(taxonomy) > 3:
                        subfamily = taxonomy[-1]

            # Extract DOI from references
            for reference in record.annotations.get('references', []):
                if 'DOI' in reference.journal:
                    doi = reference.journal.split('DOI:')[-1].strip()
                    break  # Only take the first DOI found

            return {
                'Accession Number': accession_number,
                'Phage Name': phage_name,
                'Size (bp)': size_bp,
                '% GC Content': gc_content,
                'Number of tRNA': num_tRNA,
                'Number of CDS': num_cds,
                'Family': family,
                'Subfamily': subfamily,
                'DOI': doi
            }


def write_to_csv(data, csv_file_path):
    """Write the extracted information to a CSV file."""
    with open(csv_file_path, 'w', newline='') as csvfile:
        fieldnames = ['Accession Number', 'Phage Name', 'Size (bp)', '% GC Content',
                      'Number of tRNA', 'Number of CDS', 'Family', 'Subfamily', 'DOI']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for entry in data:
            writer.writerow(entry)


def main(genbank_dir, output_csv):
    """Main function to orchestrate the GenBank file parsing and CSV writing."""
    genbank_files = [os.path.join(genbank_dir, file) for file in os.listdir(genbank_dir) if file.endswith('.gb')]
    extracted_data = []

    for gb_file in genbank_files:
        info = extract_info_from_gb(gb_file)
        if info:
            extracted_data.append(info)

    write_to_csv(extracted_data, output_csv)

# Usage example
genbank_directory = '/Users/josediogomoura/Desktop/BioFago/github/data/output/phages'
output_csv_file = '/Users/josediogomoura/Desktop/BioFago/github/data/output/phages/phages_info.csv'

if __name__ == '__main__':
    main(genbank_directory, output_csv_file)
