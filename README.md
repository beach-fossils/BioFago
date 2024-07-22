# BioFago Tool

BioFago is a comprehensive genomic analysis tool designed for processing and analyzing bacterial genomes, with a focus on plasmid detection, resistance gene identification, and CRISPR array analysis.

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
    
    ```bash
   python -m venv venv
   
    # On Windows
    source venv/bin/activate  # On Windows, use venv\Scripts\activate
   
    # On Unix or MacOS
    source venv/bin/activate
    ```
   
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
