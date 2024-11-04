#!/bin/bash
#SBATCH --job-name=metaprofiles_conditions_combined_replicates
#SBATCH --account=kubacki.michal
#SBATCH --mem=86GB
#SBATCH --time=48:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/logs/metaprofiles_conditions_combined_replicates.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/logs/metaprofiles_conditions_combined_replicates.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb 

# Check if the first script completed
DONE_FILE="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/metaprofiles_conditions.done"
if [ ! -f "${DONE_FILE}" ]; then
    echo "Error: Individual replicates processing must complete first"
    exit 1
fi

# Run the combined replicates script
python metaprofiles_opt_conditions_combined_replicates2.py \
    --gtf-file gencode.vM10.annotation.gtf.gz \
    --genome-size-file mm10.chrom.sizes \
    --output-dir /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/metaprofile_results_conditions_combined_replicates2 \
    --max-cores 32 \
    --reuse-bigwigs /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/metaprofile_results_conditions/bigwig \
    --endo-sample-sheet ENDOGENOUS_sample_sheet.csv \
    --exo-sample-sheet EXOGENOUS_sample_sheet.csv

# --sample-sheet EXOGENOUS_sample_sheet.csv \