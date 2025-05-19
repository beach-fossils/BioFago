# ErwinATyper Tool

A comprehensive genomic analysis tool for *Erwinia amylovora* providing:

- Multi-locus sequence typing (MLST) analysis
- Locus typing (capsule, cellulose, LPS, sorbitol)
- Plasmid detection
- Streptomycin resistance gene identification
- CRISPR genotype analysis
- Type III/VI secretion system variant identification
- Identification of flagellar systems


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

You can run ErwinATyper using command-line arguments.

### Using Command-Line Arguments

Run the tool with the following command-line arguments:


```bash
python biofago_runner.py --input <input_path> --output_dir <output_directory> [options]
```


Available options:

```
--input: Specify the genome file (.fasta) or directory containing multiple genome files to be processed as input (mandatory)
--output_dir: Specify the output directory for results (mandatory)
--keep_sequence_loci: Flag to retain sequences for each analyzed locus (optional)
--threshold_species: Set the ANI threshold for species assignment (optional, default: 0.95)
--skip_species_assignment: Flag to skip the module to identify the species (optional)
--log_level: Set logging verbosity (optional, choices: DEBUG, INFO, WARNING, ERROR, CRITICAL, default: INFO)
--batch_size: Number of genomes to process in parallel (optional, default: 0 = process all at once)
--num_workers: Number of worker processes for parallelization (optional, default: 4)
--docker_limit: Maximum number of concurrent Docker containers (optional, default: 4)
--quiet: Reduce console output and show only progress bars and essential messages (optional)
```


*Examples:*

```bash
# Process a single genome
python biofago_runner.py --input /path/to/genome1.fasta --output_dir /path/to/output --keep_sequence_loci

# Process a directory with multiple genomes
python biofago_runner.py --input /path/to/genomes_folder --output_dir /path/to/output --keep_sequence_loci

# Process a large batch of genomes with controlled parallelism
python biofago_runner.py --input /path/to/many_genomes --output_dir /path/to/output --batch_size 10 --num_workers 8 --docker_limit 4

# Process genomes with minimal console output (just progress bars and final result)
python biofago_runner.py --input /path/to/genomes --output_dir /path/to/output --quiet
```

### Batch Processing and Parallelization

ErwinATyper supports processing large batches of genomes efficiently:

- **Batch Size (`--batch_size`)**: Process genomes in smaller batches to control memory usage. A value of 0 means process all genomes at once.
- **Worker Processes (`--num_workers`)**: Control how many genomes are processed in parallel.
- **Docker Container Limit (`--docker_limit`)**: Limit the number of concurrent Docker containers to prevent system overload.

For best performance on large datasets (200+ genomes):

1. Use a batch size of 10-20 genomes
2. Set worker processes based on your CPU cores (typically 4-8)
3. Limit Docker containers to 4-6 to avoid excessive resource usage

### File Naming

ErwinATyper preserves full genome filenames (without extension) in all results. For example:
- Input file: `GCA_023183245.1_GCA_023183245.1_ASM2318324v1_genomic.fna`
- Name in results: `GCA_023183245.1_GCA_023183245.1_ASM2318324v1_genomic`

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

For a comprehensive example of the analysis results, you can view a sample `example.csv` output file in this GitHub repository:


[example.csv](https://github.com/beach-fossils/BioFago/blob/main/examples/example.csv)



## About the Project

**Funding & Timeline**
- Project Reference: PRR-C05-i03-I-000179
- Funded by: Plano de Recuperação e Resiliência (PRR)
 - República Portuguesa
 - União Europeia - NextGenerationEU
- Duration: Jan 2023 - Sep 2025

**Partners**
- Universities: UM (Coordinator), FCUP, IPVC
- Research Centers: INIAV
- Industry: ANP, COTHN, Asfertglobal, Frutus, Granfer, CAB, Coopval, Fruoeste, Cooperfrutas
  

## Scientific Documentation

For a detailed scientific description of the project and ErwinATyper's implementation and methodology, please refer to the [scientific document](https://github.com/beach-fossils/BioFago/blob/main/ErwinATyper_sci_doc.pdf).


## Contact

Bioinformatician/ Software dev: José Diogo Moura (UM)
email: jddmoura@gmail.com
