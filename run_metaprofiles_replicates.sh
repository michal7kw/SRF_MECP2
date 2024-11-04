#!/bin/bash
#SBATCH --job-name=metaprofiles_replicates
#SBATCH --account=kubacki.michal
#SBATCH --mem=86GB
#SBATCH --time=48:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/logs/metaprofiles_replicates.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MECP2/logs/metaprofiles_replicates.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb 

Full run with 32 cores
python metaprofiles_opt_replicates.py \
    --gtf-file gencode.vM10.annotation.gtf.gz \
    --sample-sheet EXOGENOUS_sample_sheet.csv \
    --genome-size-file mm10.chrom.sizes \
    --output-dir metaprofile_results_opt_replicates \
    --max-cores 32

# Test run
# python metaprofiles_opt.py \
#     --gtf-file gencode.vM10.annotation.gtf.gz \
#     --sample-sheet EXOGENOUS_sample_sheet.csv \
#     --genome-size-file mm10.chrom.sizes \
#     --output-dir metaprofile_results_opt \
#     --max-cores 32 \
#     --test-run