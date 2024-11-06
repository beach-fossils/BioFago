# BioFago_Erwinia Tool

BioFago_Erwinia is a comprehensive genomic analysis tool designed for processing and analyzing bacterial genomes, specifically for *Erwinia amylovora*, the causative agent of fire blight in rosaceous plants. This tool focuses on:

1. Relevant locus typing (capsule, cellulose, LPS, and sorbitol loci)
2. Plasmid detection
3. Streptomycin resistance gene identification
4. CRISPR genotype analysis

BioFago provides researchers and plant pathologists with a powerful platform for in-depth genomic characterisation of *E. amylovora* strains, facilitating studies on population dynamics, virulence factors, and antibiotic resistance in the context of fire blight management.

## Prerequisites

- Python 3.9+
- Docker 26.0.0 (with the required images pulled)
- BLAST 2.15.0+

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/beach-fossils/BioFago.git
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
   

3. Using a Conda environment (alternative method)


   #### Create a Conda environment named 'biofago_env'
   ```bash
   conda create -n biofago_env python=3.9
   ```

   #### Activate the Conda environment
   ```bash
   conda activate biofago_env
   ```
   
   ##### Your command prompt should now show (biofago_env), indicating it's active


   
4. Install the required Python packages:
    
    ```bash
   pip install -r requirements.txt
    ```
   
5. Install BLAST:

   
   This script will check if BLAST is already installed and at the correct version. If not, it will attempt to install or update BLAST.

   Note: The script requires sudo privileges on Linux systems and Homebrew on macOS. Also, if you are running this on Windows you should use Git Bash or WSL to execute the following commands since Windows CMD or PowerShell does not support `chmod`. If you encounter any issues, please refer to the [BLAST manual installation instructions](https://www.ncbi.nlm.nih.gov/books/NBK279671/).
   
   ```bash
     chmod +x external/blast/install_blast.sh
     ./external/blast/install_blast.sh
   ```        


6. Install Docker:

Follow the [official Docker installation guide](https://docs.docker.com/get-docker/) for your operating system.



## Docker Images

The development that has been made until now relies on two Docker images for some of its functionality. Before running the tool, make sure to pull these images:

1. Prokka (for genome annotation):
   ```bash
   docker pull staphb/prokka:latest
   ```

2. Average Nucleotide Identity (ANI) calculator:
   ```bash
   docker pull leightonpritchard/average_nucleotide_identity:v0.2.9
   ```



## Usage

You can run BioFago_Erwinia using command-line arguments.

### Using Command-Line Arguments

Run the tool with the following command-line arguments:


```bash
python biofago_runner.py --input <input_path> --output_dir <output_directory> [options]
```


Available options:

```
--input: Specify the genome file (.fasta) to be processed as input (mandatory)
--output_dir: Specify the output directory for results (mandatory)
--keep_sequence_loci: Flag to retain sequences for each analyzed locus (optional)
--threshold_species: Set the ANI threshold for species assignment (optional, default: 0.95)
--log_level: Set logging verbosity (optional, choices: DEBUG, INFO, WARNING, ERROR, CRITICAL, default: INFO)
```


*Example:*

```bash
python biofago_runner.py --input /path/to/genome1.fasta --output_dir /path/to/output --keep_sequence_loci
```


*Note: The output folders named `species_finder` and `types_finder` (if --keep_sequence_loci is used) are automatically created in the specified output directory.*



### Output Structure
After running the tool, you can expect the following output structure:
```
output_dir/
│
├── species_finder/
│   └── all_results.csv
│
├── types_finder/ (if --keep_sequence_loci is used)
│   └── genome1/
│       ├── types_capsule/
│       ├── types_cellulose/
│       ├── types_lps/
│       └── types_srl/
│           ├── PROKKA_[DATE].fna
│           └── PROKKA_[DATE].gbk


```

### Example `.csv` output

For a comprehensive example of the analysis results, you can view a sample `all_results.csv` file in our GitHub repository:


[Example all_results.csv](https://github.com/beach-fossils/BioFago/blob/main/test-data/csv/PRR1_INIAV.csv)


This example file demonstrates the structure and content of the BioFago analysis output, including:

- Genome statistics (e.g., contig count, N50, GC content)
- Species identification and ANI scores
- Plasmid detection results
- Streptomycin resistance gene information
- CRISPR spacer counts and genotypes
- Locus typing results for capsule, cellulose, LPS, and sorbitol loci

Reviewing this example can help you understand the type and format of data produced by the BioFago tool.
