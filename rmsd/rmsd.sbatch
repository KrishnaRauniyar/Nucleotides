#!/bin/bash
#SBATCH -p workq
#SBATCH -n 64
#SBATCH -t 01:00:00
#SBATCH -A loni_tsr_4
#SBATCH -N 1
#SBATCH -J rmsd
#SBATCH -o output_rmsd.out
#SBATCH -e error_rmsd.err

# parameters \
#   (--csv for input csv file)
#   (--residue for each selection of residues)
#   (--subgroup for all cross residue rmsd calculation with certain subgroup selection like backbone, sugar, phosphate and base atoms)
#   (--processes for CPU utilizaton)

# Printing environment info for debugging
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "Running on host: $(hostname)"
echo "Using Python: $(which python3)"
echo "Python version: $(python3 --version)"

# Define paths
WORK_DIR="/ddnB/work/rauniyar/Spring_2024/rmsd"
VENV_DIR="$WORK_DIR/myenv"
REQUIREMENTS_FILE="$WORK_DIR/requirements.txt"
SCRIPT="$WORK_DIR/rmsd.py"
CSV_FILE="/ddnB/work/rauniyar/Spring_2024/drugs/key_generator/drug_csv/proteinNumKeysDist.csv"

# Set up virtual environment
echo "Creating virtual environment at $VENV_DIR..."
python3 -m venv "$VENV_DIR"

# Activate the virtual environment
source "$VENV_DIR/bin/activate"

# Upgrade pip and install dependencies
echo "Upgrading pip..."
pip install --upgrade pip
echo "Installing dependencies from $REQUIREMENTS_FILE..."
pip install -r "$REQUIREMENTS_FILE"

# Run the script (uncomment the desired command)
echo "Running the RMSD script..."

# For specific residue (e.g., DA, A, ...)
python3 "$SCRIPT" --csv "$CSV_FILE" --residue DA

# For subgroup across all residues (e.g., backbone, sugar, phosphate and base atoms with 4 processes)
# python3 "$SCRIPT" --csv "$CSV_FILE" --subgroup backbone --processes 4

echo "Job completed"
