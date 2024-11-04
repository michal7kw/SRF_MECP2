#!/bin/bash
#SBATCH --job-name=heatmaps_conditions_combined_replicates_reuse
#SBATCH --account=kubacki.michal
#SBATCH --mem=86GB
#SBATCH --time=48:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/logs/heatmaps_conditions_combined_replicates_reuse.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/logs/heatmaps_conditions_combined_replicates_reuse.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb 

# Print working directory and list contents of matrix directory
echo "Current working directory: $(pwd)"
echo "Contents of matrix directory:"
ls -R /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/metaprofile_results_conditions_combined_replicates/matrix_dir

python heatmaps_conditions_combined_replicates_replicates_reuse.py \
    --matrix-dir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/metaprofile_results_conditions_combined_replicates/matrix_dir \
    --output-dir heatmaps_combined_reuse \
    --tissues Neuron NSC \
    --categories up down non \
    --conditions Endogenous Exogenous \
    --max-cores 32