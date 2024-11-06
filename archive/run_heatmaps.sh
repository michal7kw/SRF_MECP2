#!/bin/bash
#SBATCH --job-name=heatmaps
#SBATCH --account=kubacki.michal
#SBATCH --mem=86GB
#SBATCH --time=48:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/logs/heatmaps.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/logs/heatmaps.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb 

python heatmaps.py \
    --matrix-dir metaprofile_results_opt \
    --output-dir heatmaps \
    --tissues Neuron NSC \
    --categories non up down \
    --max-cores 32

# python heatmaps.py \
#     --matrix-dir metaprofile_results_opt \
#     --output-dir heatmaps \
#     --tissues Neuron \
#     --categories up \
#     --debug