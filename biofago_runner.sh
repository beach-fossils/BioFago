#!/bin/bash

# Activate conda environment
conda activate biofago_env

# Create logs directory if it doesn't exist
mkdir -p /Users/josediogomoura/Documents/BioFago/BioFago_to_Server/BioFagov4/logs

# Create log filename with timestamp
LOG_FILE="/Users/josediogomoura/Documents/BioFago/BioFago_to_Server/BioFagov4/logs$(date +%Y%m%d_%H%M%S).log"

# Run BioFago in detached mode with logging
nohup python biofago_runner.py \
    --input /Users/josediogomoura/Documents/BioFago/BioFago_to_Server/BioFagov4/test-data/genomes \
    --output_dir /Users/josediogomoura/Documents/BioFago/BioFago_to_Server/BioFagov4/test-data/genomes/results \
    --keep_sequence_loci \
    > "$LOG_FILE" 2>&1 &

# Print the process ID and log file location
echo "BioFago process started with PID: $!"
echo "Log file location: $LOG_FILE"