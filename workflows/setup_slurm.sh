#!/bin/bash
### SLURM + Snakemake submission helper
cd ~/cs4641-team24/workflows/ # Force workflow directory

# Constants
LOG_DIR="results/logs"

# Default parameters
TIME="1:00:00"
CPUS=4
MEM="8G"
JOBS=4
NAME="cs4641-team24"
GPUS=0 # default no GPU
GPU_TYPE="" # optional GPU type specification

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --time) TIME="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem) MEM="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    --part) PARTITION="$2"; shift 2 ;;
    --gpu-type) GPU_TYPE="$2"; shift 2 ;;
    -g|--gpu)
      # Check if next arg is a number (e.g., --gpu 2)
      if [[ "$2" =~ ^[0-9]+$ ]]; then
        GPUS="$2"
        shift 2
      else
        GPUS=1
        shift 1
      fi
      ;;
    --help)
      echo "Usage: $0 [--time hh:mm:ss] [--cpus N] [--mem SIZE] [--jobs N] [--name JOBNAME] [--part PARTITION] [--gpu [N]] [--gpu-type TYPE] [snakemake args]"
      echo "Example: $0 --time 2:00:00 --mem 16G --gpu 2 --gpu-type A100 --jobs 8 target_file"
      echo "GPU types: V100, RTX_6000, A40, A100, H100, H200, L40S, MI210"
      exit 0
      ;;
    *) break ;;
  esac
done

# Logging setup
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOG_DIR="$LOG_DIR/${TIMESTAMP}"

# Build resource flags
RESOURCE_FLAGS=(
  walltime=$TIME
  cpus_per_task=$CPUS
  mem_mb=${MEM%G}000
)

# Add GPU resources if requested
if [[ "$GPUS" -gt 0 ]]; then
  RESOURCE_FLAGS+=(gpus_per_task=$GPUS)
fi

# Build SLURM extra arguments for GPU
SLURM_EXTRA=""
if [[ "$GPUS" -gt 0 ]]; then
  if [[ -n "$GPU_TYPE" ]]; then
    SLURM_EXTRA="--gres=gpu:${GPU_TYPE}:$GPUS"
  else
    SLURM_EXTRA="--gres=gpu:$GPUS"
  fi
fi

# Print config
echo "[$(date)] Starting Snakemake with:"
echo "  TIME=$TIME, CPUS=$CPUS, MEM=$MEM, JOBS=$JOBS, NAME=$NAME"
if [[ "$GPUS" -gt 0 ]]; then
  echo "  GPUS=$GPUS${GPU_TYPE:+, TYPE=$GPU_TYPE}"
  echo "  SLURM_EXTRA: $SLURM_EXTRA"
fi
echo "  Logs -> $LOG_DIR"
echo "  Targets: $@"

# Run Snakemake
if [[ "$GPUS" -gt 0 ]]; then
  snakemake \
    --executor slurm \
    -j "$JOBS" \
    --default-resources "${RESOURCE_FLAGS[@]}" \
    --default-resources "slurm_extra='$SLURM_EXTRA'" \
    --retries 3 \
    --latency-wait 30 \
    --rerun-incomplete \
    "$@"
else
  snakemake \
    --executor slurm \
    -j "$JOBS" \
    --default-resources "${RESOURCE_FLAGS[@]}" \
    --retries 3 \
    --latency-wait 30 \
    --rerun-incomplete \
    "$@"
fi
