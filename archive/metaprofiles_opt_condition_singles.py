import argparse
import pandas as pd
import numpy as np
import subprocess
import matplotlib.pyplot as plt
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import tempfile
import logging
from functools import partial
import gzip
import shutil
import time

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def run_command(cmd, desc=None):
    """Run command with proper logging and error handling"""
    if desc:
        logging.info(desc)
    try:
        result = subprocess.run(cmd, shell=True, check=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {e.stderr}")
        return False

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
        run_command(cmd, f"Creating BED file: {output_file}")
    finally:
        os.unlink(temp_file_path)
    return output_file

def bam_to_bigwig(bam_file, output_dir, genome_size_file, threads):
    """Convert BAM to BigWig with RPKM normalization"""
    basename = Path(bam_file).stem
    output_bw = f"{output_dir}/{basename}.bw"
    
    if os.path.exists(output_bw):
        logging.info(f"BigWig file exists for {basename}, skipping")
        return output_bw
        
    logging.info(f"Starting conversion of {basename} BAM to BigWig")
    
    with tempfile.TemporaryDirectory(prefix=f"tmp_{basename}_") as temp_dir:
        sorted_bam = f"{temp_dir}/{basename}.sorted.bam"
        
        # Sort BAM
        sort_cmd = f"samtools sort -@ {threads} -m 2G -T {temp_dir}/sort_tmp {bam_file} -o {sorted_bam}"
        if not run_command(sort_cmd, f"Sorting {basename} BAM"):
            return None
            
        # Index BAM
        index_cmd = f"samtools index -@ {threads} {sorted_bam}"
        if not run_command(index_cmd, f"Indexing {basename} BAM"):
            return None
            
        # Convert to BigWig
        coverage_cmd = f"""
        bamCoverage -b {sorted_bam} \
            -o {output_bw} \
            --binSize 10 \
            --normalizeUsing RPKM \
            --numberOfProcessors {threads} \
            --verbose
        """
        if not run_command(coverage_cmd, f"Converting {basename} to BigWig"):
            return None
    
    logging.info(f"Successfully converted {basename} to BigWig")
    return output_bw

def generate_metaprofile(args):
    """Generate single metaprofile"""
    tissue, category, bigwigs, bed_file, output_dir, threads = args
    
    matrix_dir = f"{output_dir}/matrix_{tissue}_{category}"
    os.makedirs(matrix_dir, exist_ok=True)
    
    matrix_file = f"{matrix_dir}/{tissue}_{category}_matrix.gz"
    profile_file = f"{output_dir}/{tissue}_{category}_profile.png"
    
    cmd_matrix = f"""
    computeMatrix scale-regions -S {' '.join(bigwigs)} \
        -R {bed_file} \
        -b 5000 -a 5000 \
        --regionBodyLength 5000 \
        --skipZeros \
        --numberOfProcessors {threads} \
        -o {matrix_file}
    """
    
    cmd_plot = f"""
    plotProfile -m {matrix_file} \
        --perGroup \
        --colors red blue \
        --plotTitle "{tissue} {category}-regulated genes" \
        -out {profile_file}
    """
    
    if not run_command(cmd_matrix, f"Generating matrix for {tissue} {category}"):
        return None
    if not run_command(cmd_plot, f"Creating profile plot for {tissue} {category}"):
        return None
    
    return matrix_file, profile_file

def main(args):
    start_time = time.time()
    setup_logging()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Decompress GTF if needed
    gtf_file = decompress_gtf(args.gtf_file)
    
    # Process tissues
    tissues = ['Neuron', 'NSC'] if not args.test_run else ['Neuron']
    deseq_files = {
        'Neuron': 'DEA_NEU.csv',
        'NSC': 'DEA_NSC.csv'
    }
    
    # Calculate resource allocation
    threads_for_bam = min(args.max_cores, 16)  # Cap threads for BAM processing
    threads_for_meta = min(args.max_cores // 2, 8)  # Threads for metaprofile generation
    
    logging.info(f"Using {threads_for_bam} threads for BAM processing")
    logging.info(f"Using {threads_for_meta} threads per metaprofile generation")
    
    # Process DESeq results and create BED files
    tissue_data = {}
    for tissue in tissues:
        logging.info(f"Processing DESeq results for {tissue}")
        gene_sets = process_deseq_results(deseq_files[tissue])
        categories = ['non', 'up', 'down'] if not args.test_run else ['up']
        
        bed_files = {}
        for category in categories:
            output_file = f"{args.output_dir}/{tissue}_{category}_regulated_genes.bed"
            bed_files[category] = create_bed_file(gene_sets[category], gtf_file, output_file)
        
        tissue_data[tissue] = {'gene_sets': gene_sets, 'bed_files': bed_files}
    
    # Convert BAM to BigWig sequentially
    samples = pd.read_csv(args.sample_sheet)
    bigwig_dir = f"{args.output_dir}/bigwig"
    os.makedirs(bigwig_dir, exist_ok=True)
    
    bam_files = samples['bamReads'].tolist()
    if args.test_run:
        bam_files = bam_files[:2]
    
    bigwig_files = []
    for bam_file in bam_files:
        bigwig = bam_to_bigwig(bam_file, bigwig_dir, args.genome_size_file, threads_for_bam)
        if bigwig:
            bigwig_files.append(bigwig)
        else:
            logging.error(f"Failed to process {bam_file}")
            return
    
    # Generate metaprofiles with parallel processing
    metaprofile_args = []
    for tissue in tissues:
        tissue_samples = samples[samples['Tissue'] == tissue]
        tissue_bigwigs = [bw for bw, sample in zip(bigwig_files, bam_files)
                         if sample in tissue_samples['bamReads'].tolist()]
        
        for category, bed_file in tissue_data[tissue]['bed_files'].items():
            metaprofile_args.append((
                tissue, category, tissue_bigwigs, bed_file,
                args.output_dir, threads_for_meta
            ))
    
    with ProcessPoolExecutor(max_workers=min(len(metaprofile_args), args.max_cores // threads_for_meta)) as executor:
        results = list(executor.map(generate_metaprofile, metaprofile_args))
    
    elapsed_time = time.time() - start_time
    logging.info(f"Processing completed in {elapsed_time:.2f} seconds!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate metaprofiles using multiple cores")
    parser.add_argument("--gtf-file", required=True, help="Path to GTF file (gzipped)")
    parser.add_argument("--sample-sheet", required=True, help="Path to sample sheet CSV")
    parser.add_argument("--genome-size-file", required=True, help="Path to genome size file")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--max-cores", type=int, default=32, help="Maximum number of cores to use")
    parser.add_argument("--test-run", action="store_true", help="Run in test mode (single tissue, up-regulated only)")
    
    args = parser.parse_args()
    main(args)