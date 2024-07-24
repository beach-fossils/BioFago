# BioFago_Erwinia Tool

BioFago_Erwinia is a comprehensive genomic analysis tool designed for processing and analyzing bacterial genomes, specifically for *Erwinia amylovora*, the causative agent of fire blight in rosaceous plants. This tool focuses on:

1. Relevant locus typing (capsule, cellulose, LPS, and sorbitol loci)
2. Plasmid detection
3. Streptomycin resistance gene identification
4. CRISPR genotype analysis

BioFago provides researchers and plant pathologists with a powerful platform for in-depth genomic characterisation of *E. amylovora* strains, facilitating studies on population dynamics, virulence factors, and antibiotic resistance in the context of fire blight management.

## Prerequisites

- Python 3.8+
- Docker
- BLAST 2.15.0+

## Installation

1. Clone the repository (dev branch):
   ```bash
   git clone -b dev https://github.com/beach-fossils/BioFago.git
   cd BioFago
   ```

2. Set up a Python virtual environment (optional but recommended):

   #### Create a virtual environment named 'biofago_env'
   ```bash
   python -m venv biofago_env
   ```

   #### Activate the virtual environment
   
   ##### On Windows:
   ```bash
   biofago_env\Scripts\activate
   ```

   ##### On Unix or MacOS:
   ```
   source biofago_env/bin/activate
   ```

   ##### Your command prompt should now show (biofago_env), indicating it's active
   
4. Install the required Python packages:
    
    ```bash
   pip install -r requirements.txt
    ```
   
5. Install BLAST:

   
   This script will check if BLAST is already installed and at the correct version. If not, it will attempt to install or update BLAST.

   Note: The script requires sudo privileges on Linux systems and Homebrew on macOS. If you encounter any issues, please refer to the [BLAST manual installation instructions](https://www.ncbi.nlm.nih.gov/books/NBK279671/).
   
   ```bash
     chmod +x external/blast/install_blast.sh
     ./external/blast/install_blast.sh
   ```        


6. Install Docker:

Follow the [official Docker installation guide](https://docs.docker.com/get-docker/) for your operating system.


## Configuration

1. Edit `config.yaml` to match your environment and requirements:

   ```yaml
   genomes_folder: "/path/to/your/genomes/folder"  # Mandatory: specify the folder containing your genome files
   keep_sequence_loci: false  # Optional: set to true to retain sequences for each analyzed locus
   threshold_species: 0.95  # Optional: ANI threshold for species assignment using Pyani
   log_level: "INFO"  # Optional: set logging verbosity (DEBUG, INFO, WARNING, ERROR, CRITICAL)
   

*Note: The output folder is automatically created at the same level as the `genomes_folder`. It is named `species_finder`.*


## Usage

   ```bash
  python src/biofago_runner.py
   ```

This command will process all genomes in the specified genomes_folder and output results to the automatically created `species_finder` folder.

### Output Structure
After running the tool, you can expect the following output structure:
```
parent_folder/
│
├── genomes_folder/
│   ├── genome1.fasta
│   ├── genome2.fasta
    ...
│
└── species_finder/
    ├── all_results.csv
```

### Example Output

For a comprehensive example of the analysis results, you can view a sample `all_results.csv` file in our GitHub repository:

[Example all_results.csv](https://github.com/beach-fossils/BioFago/blob/dev/examples/all_results.csv)

This example file demonstrates the structure and content of the BioFago analysis output, including:

- Genome statistics (e.g., contig count, N50, GC content)
- Species identification and ANI scores
- Plasmid detection results
- Streptomycin resistance gene information
- CRISPR spacer counts and genotypes
- Locus typing results for capsule, cellulose, LPS, and sorbitol loci

Reviewing this example can help you understand the type and format of data produced by the BioFago tool.
