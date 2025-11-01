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
QOS="coc-ice"
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
    --qos) QOS="$2"; shift 2 ;;
    --gpu-type) GPU_TYPE="$2"; shift 2 ;;
    -g|--gpu)
      if [[ "$2" =~ ^[0-9]+$ ]]; then
        GPUS="$2"
        shift 2
      else
        GPUS=1
        shift 1
      fi
      ;;
    --help)
      echo "Usage: $0 [--time hh:mm:ss] [--cpus N] [--mem SIZE] [--jobs N] [--name JOBNAME] [--qos QOS] [--gpu [N]] [--gpu-type TYPE] [snakemake args]"
      echo "Example: $0 --time 2:00:00 --mem 16G --gpu 2 --gpu-type A100 --qos inferno --jobs 8 target_file"
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

# Adjust defaults for GPU jobs
if [[ "$GPUS" -gt 0 ]]; then
  if [[ "$MEM" == "8G" ]]; then
    MEM="16G"
    echo "Note: Increased memory to 16G for GPU job"
  fi
  
  if [[ "$CPUS" -lt 4 ]]; then
    CPUS=4
  fi
fi

# Build resource flags - ALL in one array
RESOURCE_FLAGS=(
  walltime=$TIME
  cpus_per_task=$CPUS
  mem_mb=${MEM%G}000
  slurm_qos=$QOS
)

# Add GPU resources if requested
if [[ "$GPUS" -gt 0 ]]; then
  RESOURCE_FLAGS+=(gpus_per_task=$GPUS)
  
  # Build and add slurm_extra for GPU
  if [[ -n "$GPU_TYPE" ]]; then
    RESOURCE_FLAGS+=(slurm_extra="--gres=gpu:${GPU_TYPE}:$GPUS")
  else
    RESOURCE_FLAGS+=(slurm_extra="--gres=gpu:$GPUS")
  fi
fi

# Print config
echo "[$(date)] Starting Snakemake with:"
echo "  TIME=$TIME, CPUS=$CPUS, MEM=$MEM, JOBS=$JOBS, NAME=$NAME"
if [[ -n "$QOS" ]]; then
  echo "  QOS=$QOS"
fi
if [[ "$GPUS" -gt 0 ]]; then
  echo "  GPUS=$GPUS${GPU_TYPE:+, TYPE=$GPU_TYPE}"
fi
echo "  Logs -> $LOG_DIR"
echo "  Targets: $@"
echo "  Resources: ${RESOURCE_FLAGS[@]}"

# Run Snakemake - single --default-resources call
snakemake \
  --executor slurm \
  -j "$JOBS" \
  --default-resources "${RESOURCE_FLAGS[@]}" \
  --retries 3 \
  --latency-wait 30 \
  --rerun-incomplete \
  "$@"
