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

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/biofago_tool.git
   cd biofago_tool
   ```

2. Set up a Python virtual environment (optional but recommended):

   # Create a virtual environment named 'biofago_env'
   python -m venv biofago_env

   # Activate the virtual environment
   
   ## On Windows:
   biofago_env\Scripts\activate

   ## On Unix or MacOS:
   source biofago_env/bin/activate

   # Your command prompt should now show (biofago_env), indicating it's active
   
3. Install the required Python packages:
    
    ```bash
   pip install -r requirements.txt
    ```
   
4. Install BLAST:

   
   This script will check if BLAST is already installed and at the correct version. If not, it will attempt to install or update BLAST.

   Note: The script requires sudo privileges on Linux systems and Homebrew on macOS. If you encounter any issues, please refer to the [BLAST manual installation instructions](https://www.ncbi.nlm.nih.gov/books/NBK279671/).
   
   ```bash
     chmod +x external/blast/install_blast.sh
     ./external/blast/install_blast.sh
   ```        


5. Install Docker:

Follow the [official Docker installation guide](https://docs.docker.com/get-docker/) for your operating system.

## Configuration

1. Copy the example configuration file:

   ```bash
   cp config.yaml.example config.yaml
   ```

2. Edit config.yaml to match your environment and requirements:
   ```yaml
   genomes_folder: "/path/to/your/genomes/folder"
   keep_sequence_loci: false
   threshold_species: 0.95
   log_level: "INFO"
   ```

## Usage

   ```bash
  python src/biofago_runner.py
   ```
