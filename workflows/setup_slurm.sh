#!/bin/bash

### This file sets up the slurm with snakemake with some predefined parameters

cd ~/cs4641-team24/workflows/ # Force workflow directory

# Constants
LOG_DIR="results/logs"
JOBS=4

# Default parameters
TIME="1:00:00"
CPUS=4
MEM="8G"
JOBS=4
NAME="cs4641-team24"

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --time) TIME="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem) MEM="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    --help)
      echo "Usage: $0 [--time hh:mm:ss] [--cpus N] [--mem SIZE] [--jobs N] [--name JOBNAME] [snakemake args]"
      exit 0
      ;;
    *) break ;;
  esac
done

# Logging setup
LOG_DIR="results/logs"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOG_FILE="$LOG_DIR/${NAME}_${TIMESTAMP}.log"

# Run Snakemake
echo "[$(date)] Starting Snakemake with:"
echo "  TIME=$TIME, CPUS=$CPUS, MEM=$MEM, JOBS=$JOBS, NAME=$NAME"
echo "  Logs -> $LOG_FILE"

snakemake \
  --jobs "$JOBS" \
  --use-conda \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds \
  --latency-wait 60 \
  --cluster "sbatch --job-name=$NAME --time=$TIME --cpus-per-task=$CPUS --mem=$MEM --output=$LOG_DIR/%x_%j.out --error=$LOG_DIR/%x_%j.err" \
  "$@" > "$LOG_FILE" 2>&1
