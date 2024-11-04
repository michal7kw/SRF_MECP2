"""
This is a bioinformatics pipeline script designed to generate and visualize gene expression metaprofiles by comparing endogenous and exogenous conditions across different tissues (Neurons and NSCs - Neural Stem Cells). Here's a breakdown of its main functions:

Data Processing:

Takes input from BAM files (genomic alignment data)
Uses sample sheets for both endogenous and exogenous conditions
Processes DESeq2 results to identify up-regulated, down-regulated, and non-regulated genes
Converts BAM files to BigWig format with RPKM normalization


File Handling:

Works with compressed GTF (Gene Transfer Format) files
Creates BED files for different gene categories (up/down/non-regulated)
Manages file conversions and temporary files
Reuses existing BigWig files when available


Metaprofile Generation:

Creates separate profiles for each tissue type (Neuron and NSC)
Analyzes three categories of genes: up-regulated, down-regulated, and non-regulated
Generates matrices using computeMatrix from deepTools
Creates visualization plots showing gene expression profiles
Includes regions 5000 base pairs before and after genes


Visualization Features:

Generates plots with different colors for endogenous (red) and exogenous (blue) conditions
Labels samples appropriately based on their replicate numbers
Shows gene regions from TSS (Transcription Start Site) to TES (Transcription End Site)
Includes proper legends and axis labels


Performance Optimization:

Uses parallel processing with both ProcessPoolExecutor and ThreadPoolExecutor
Automatically manages computational resources
Includes checkpoints to avoid regenerating existing files
Handles multiple samples concurrently with configurable core usage
"""

import argparse
import pandas as pd
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import tempfile
import logging
from functools import partial
import gzip
import shutil
import re

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def decompress_gtf(gtf_gz_path):
    """Decompress GTF file if not already decompressed"""
    gtf_path = str(gtf_gz_path).replace('.gz', '')
    if not os.path.exists(gtf_path):
        logging.info(f"Decompressing {gtf_gz_path}")
        with gzip.open(gtf_gz_path, 'rb') as f_in:
            with open(gtf_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    return gtf_path

def process_deseq_results(deseq_file, padj_threshold=0.05, lfc_threshold=1):
    """Process DESeq2 results using vectorized operations"""
    df = pd.read_csv(deseq_file)
    up_mask = (df['padj'] < padj_threshold) & (df['log2FoldChange'] > lfc_threshold)
    down_mask = (df['padj'] < padj_threshold) & (df['log2FoldChange'] < -lfc_threshold)
    non_mask = ~(up_mask | down_mask)
    
    return {
        'up': df.loc[up_mask, 'gene'].tolist(),
        'down': df.loc[down_mask, 'gene'].tolist(),
        'non': df.loc[non_mask, 'gene'].tolist()
    }

def create_bed_file(genes, gtf_file, output_file):
    """Create BED file for given genes"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        temp_file.write('\n'.join(genes))
        temp_file_path = temp_file.name
    
    try:
        cmd = f"""grep -F -f {temp_file_path} {gtf_file} | 
                 awk -F'\t' '$3=="gene" {{
                     if (match($9, /gene_name "([^"]+)"/, n))
                         print $1"\t"$4-1"\t"$5"\t"n[1]"\t.\t"$7
                 }}' > {output_file}"""
        subprocess.run(cmd, shell=True, check=True)
    finally:
        os.unlink(temp_file_path)
    return output_file

def bam_to_bigwig(args):
    """Convert BAM to BigWig with RPKM normalization"""
    bam_file, output_dir, genome_size_file, threads_per_job = args
    basename = Path(bam_file).stem
    output_bw = f"{output_dir}/{basename}.bw"
    
    if os.path.exists(output_bw):
        logging.info(f"BigWig file exists for {basename}, skipping")
        return output_bw
        
    logging.info(f"Converting {basename} BAM to BigWig")
    with tempfile.TemporaryDirectory(prefix=f"tmp_{basename}_") as temp_dir:
        sorted_bam = f"{temp_dir}/{basename}.sorted.bam"
        subprocess.run(f"samtools sort -@ {threads_per_job} -T {temp_dir}/sort_tmp {bam_file} -o {sorted_bam}",
                      shell=True, check=True)
        subprocess.run(f"samtools index -@ {threads_per_job} {sorted_bam}", 
                      shell=True, check=True)
        
        subprocess.run(f"""
            bamCoverage -b {sorted_bam} \
                -o {output_bw} \
                --binSize 10 \
                --normalizeUsing RPKM \
                --numberOfProcessors {threads_per_job}
        """, shell=True, check=True)
    
    return output_bw

def generate_metaprofile(args):
    """Generate separate metaprofiles for endogenous and exogenous conditions"""
    tissue, category, bigwigs_endo, bigwigs_exo, bed_file, output_dir, threads_per_job = args
    
    # Use absolute paths
    output_dir = os.path.abspath(output_dir)
    matrix_dir = os.path.join(output_dir, "matrix_dir")
    os.makedirs(matrix_dir, exist_ok=True)
    
    # Define final matrix and profile files
    final_matrix = os.path.join(matrix_dir, f"{tissue}_{category}_combined_matrix.gz")
    profile_file = os.path.join(output_dir, f"{tissue}_{category}_profile.png")
    
    # Generate matrix if it doesn't exist
    if not os.path.exists(final_matrix):
        logging.info(f"Generating matrix for {tissue}_{category}")
        
        # Generate matrix for tissue-specific samples only
        all_bigwigs = bigwigs_endo + bigwigs_exo
        
        # Create temporary file with list of bigwig files
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as tmp:
            tmp.write('\n'.join(all_bigwigs))
            bigwig_list = tmp.name
        
        try:
            cmd = f"""
            computeMatrix scale-regions \
                --scoreFileName {bigwig_list} \
                -R {bed_file} \
                -b 5000 -a 5000 \
                --regionBodyLength 5000 \
                --skipZeros \
                --numberOfProcessors {threads_per_job} \
                --missingDataAsZero \
                -o {final_matrix}
            """
            subprocess.run(cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"Error generating matrix for {tissue}_{category}: {str(e)}")
            logging.error(f"Command was: {cmd}")
            return None
        finally:
            os.unlink(bigwig_list)
    
    # Generate plot if needed
    if not os.path.exists(profile_file):
        logging.info(f"Generating profile plot for {tissue}_{category}")
        
        try:
            # Create temporary files for labels and colors
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_labels.txt') as labels_file:
                # Generate labels
                endo_labels = [f"Endogenous_{Path(bw).stem.split('M')[1].split('_')[0]}" 
                             for bw in bigwigs_endo]
                exo_labels = [f"Exogenous_{Path(bw).stem.split('V')[1].split('_')[0]}" 
                             for bw in bigwigs_exo]
                labels_file.write('\n'.join(endo_labels + exo_labels))
                labels_path = labels_file.name
            
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='_colors.txt') as colors_file:
                colors = ['red'] * len(bigwigs_endo) + ['blue'] * len(bigwigs_exo)
                colors_file.write('\n'.join(colors))
                colors_path = colors_file.name
            # Create plotProfile command
            cmd = f"""
            plotProfile \
                -m {final_matrix} \
                --plotTitle "{tissue} {category}-regulated genes" \
                --legendLocation "upper-right" \
                --samplesLabel $(cat {labels_path}) \
                --colors $(cat {colors_path}) \
                --regionsLabel "Genes" \
                --startLabel "TSS" \
                --endLabel "TES" \
                --yAxisLabel "Average signal" \
                --averageType mean \
                -o {profile_file}
            """
            
            logging.info(f"Using labels: {endo_labels + exo_labels}")
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')
            
        except subprocess.CalledProcessError as e:
            logging.error(f"Error generating plot for {tissue}_{category}: {str(e)}")
            logging.error(f"Command was: {cmd}")
            return None
        finally:
            # Clean up temporary files
            for f in [labels_path, colors_path]:
                if os.path.exists(f):
                    os.unlink(f)

    return final_matrix, profile_file

def main(args):
    setup_logging()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define source and destination paths for BED files
    source_dir = "./metaprofile_results_conditions"
    
    # Copy BED files from source to destination
    for tissue in ['Neuron', 'NSC']:
        for category in ['non', 'up', 'down']:
            src_bed = os.path.join(source_dir, f"{tissue}_{category}_regulated_genes.bed")
            dst_bed = os.path.join(args.output_dir, f"{tissue}_{category}_regulated_genes.bed")
            
            if not os.path.exists(src_bed):
                logging.error(f"Source BED file not found: {src_bed}")
                return
                
            if not os.path.exists(dst_bed):
                logging.info(f"Copying BED file from {src_bed} to {dst_bed}")
                try:
                    shutil.copy2(src_bed, dst_bed)
                except Exception as e:
                    logging.error(f"Error copying BED file: {e}")
                    return
    
    # Rest of your existing code for loading sample sheets and processing
    try:
        # Load and filter sample sheets
        samples_endo = pd.read_csv(os.path.join("./DATA", args.endo_sample_sheet))
        samples_exo = pd.read_csv(os.path.join("./DATA", args.exo_sample_sheet))
        
        # Filter out IgM control and keep only relevant samples
        samples_endo = samples_endo[samples_endo['SampleID'] != 'IgM']
        
        # Debug information
        logging.info(f"Loaded {len(samples_endo)} endogenous samples and {len(samples_exo)} exogenous samples")
        logging.info(f"Endogenous samples: {samples_endo['bamReads'].tolist()}")
        logging.info(f"Exogenous samples: {samples_exo['bamReads'].tolist()}")
    except Exception as e:
        logging.error(f"Error reading sample sheets: {e}")
        return

    # Map bigwig files
    bigwig_files = {}
    bigwig_dir = os.path.abspath(os.path.join("./DATA", args.reuse_bigwigs))
    logging.info(f"Looking for bigwig files in: {bigwig_dir}")
    
    # Get all available bigwig files
    available_bigwigs = {os.path.splitext(f)[0]: os.path.join(bigwig_dir, f)
                        for f in os.listdir(bigwig_dir)
                        if f.endswith('.bw')}
    logging.info(f"Available bigwig files: {list(available_bigwigs.keys())}")

    # Map sample names to bigwig files
    for sample in pd.concat([samples_endo, samples_exo])['bamReads']:
        sample_name = Path(sample).stem
        if sample_name in available_bigwigs:
            bigwig_files[sample_name] = available_bigwigs[sample_name]
            logging.info(f"Mapped sample {sample_name} to {bigwig_files[sample_name]}")
        else:
            logging.error(f"Missing bigwig file for sample: {sample_name}")
            logging.error(f"Available files: {list(available_bigwigs.keys())}")
            return

    # Process each tissue and category
    max_concurrent_jobs = min(6, args.max_cores // 5)
    threads_per_job = args.max_cores // max_concurrent_jobs
    logging.info(f"Using {max_concurrent_jobs} concurrent jobs with {threads_per_job} threads each")

    metaprofile_args = []
    for tissue in ['Neuron', 'NSC']:
        # Filter samples by tissue
        tissue_samples_endo = samples_endo[samples_endo['Tissue'] == tissue]
        tissue_samples_exo = samples_exo[samples_exo['Tissue'] == tissue]
        
        # Get bigwig paths with verification
        tissue_bigwigs_endo = []
        tissue_bigwigs_exo = []
        
        for bam in tissue_samples_endo['bamReads']:
            sample_name = Path(bam).stem
            if sample_name in bigwig_files:
                tissue_bigwigs_endo.append(bigwig_files[sample_name])
                logging.info(f"Added endogenous bigwig for {tissue}: {bigwig_files[sample_name]}")
        
        for bam in tissue_samples_exo['bamReads']:
            sample_name = Path(bam).stem
            if sample_name in bigwig_files:
                tissue_bigwigs_exo.append(bigwig_files[sample_name])
                logging.info(f"Added exogenous bigwig for {tissue}: {bigwig_files[sample_name]}")

        if not tissue_bigwigs_endo or not tissue_bigwigs_exo:
            logging.error(f"Missing bigwig files for tissue {tissue}")
            logging.error(f"Endo bigwigs: {tissue_bigwigs_endo}")
            logging.error(f"Exo bigwigs: {tissue_bigwigs_exo}")
            continue

        # Sort bigwig files to ensure consistent ordering
        tissue_bigwigs_endo.sort()
        tissue_bigwigs_exo.sort()

        for category in ['non', 'up', 'down']:
            bed_file = os.path.join(args.output_dir, f"{tissue}_{category}_regulated_genes.bed")
            if not os.path.exists(bed_file):
                logging.error(f"BED file not found: {bed_file}")
                continue
            
            metaprofile_args.append((
                tissue, category, tissue_bigwigs_endo, tissue_bigwigs_exo,
                bed_file, args.output_dir, threads_per_job
            ))

    if not metaprofile_args:
        logging.error("No valid combinations of files found to generate metaprofiles")
        return

    with ProcessPoolExecutor(max_workers=max_concurrent_jobs) as executor:
        results = list(executor.map(generate_metaprofile, metaprofile_args))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate metaprofiles using multiple cores")
    parser.add_argument("--gtf-file", required=True, help="Path to GTF file (gzipped)")
    parser.add_argument("--genome-size-file", required=True, help="Path to genome size file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--max-cores", type=int, default=32, help="Maximum number of cores to use")
    parser.add_argument("--reuse-bigwigs", required=True, help="Path to existing bigwig directory to reuse")
    parser.add_argument("--endo-sample-sheet", required=True, help="Path to endogenous sample sheet")
    parser.add_argument("--exo-sample-sheet", required=True, help="Path to exogenous sample sheet")
    parser.add_argument("--test-run", action="store_true", help="Run in test mode (single tissue, up-regulated only)")
    
    args = parser.parse_args()
    main(args)