#!/bin/bash
#SBATCH --job-name=heatmaps_conditions_combined_replicates
#SBATCH --account=kubacki.michal
#SBATCH --mem=86GB
#SBATCH --time=48:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/logs/heatmaps_conditions_combined_replicates.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/logs/heatmaps_conditions_combined_replicates.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb 

# Print working directory and list contents of matrix directory
echo "Current working directory: $(pwd)"
echo "Contents of matrix directory:"
ls -R metaprofile_results_conditions

python heatmaps_conditions_combined_replicates.py \
    --matrix-dir metaprofile_results_conditions \
    --bigwig-dir metaprofile_results_conditions/bigwig \
    --output-dir heatmaps_combined \
    --tissues Neuron NSC \
    --categories non up down \
    --conditions Endogenous Exogenous \
    --max-cores 32