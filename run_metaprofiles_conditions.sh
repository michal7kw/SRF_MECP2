#!/bin/bash
#SBATCH --job-name=metaprofiles_conditions
#SBATCH --account=kubacki.michal
#SBATCH --mem=86GB
#SBATCH --time=48:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/logs/metaprofiles_conditions.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/logs/metaprofiles_conditions.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb 

# Create a flag file to indicate completion
DONE_FILE="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/metaprofiles_conditions.done"

# Run the first script
python metaprofiles_opt_conditions.py \
    --gtf-file ./DATA/gencode.vM10.annotation.gtf.gz \
    --genome-size-file ./DATA/mm10.chrom.sizes \
    --output-dir metaprofile_results_conditions \
    --max-cores 32

# If the script completed successfully, create the done file
if [ $? -eq 0 ]; then
    touch "${DONE_FILE}"
    echo "Individual replicates processing completed successfully"
else
    echo "Individual replicates processing failed"
    exit 1
fi

# --sample-sheet EXOGENOUS_sample_sheet.csv \