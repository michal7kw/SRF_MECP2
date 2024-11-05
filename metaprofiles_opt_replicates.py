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
                --binSize 50 \
                --normalizeUsing CPM \
                --centerReads \
                --extendReads \
                --numberOfProcessors {threads_per_job}
        """, shell=True, check=True)
    
    return output_bw

def generate_metaprofile(args):
    """Generate single metaprofile"""
    tissue, category, bigwigs, bed_file, output_dir, threads_per_job = args
    
    matrix_dir = f"{output_dir}/matrix_{tissue}_{category}"
    os.makedirs(matrix_dir, exist_ok=True)
    
    matrix_file = f"{matrix_dir}/{tissue}_{category}_matrix.gz"
    profile_file = f"{output_dir}/{tissue}_{category}_profile.png"
    
    # Skip if final profile already exists
    if os.path.exists(profile_file):
        logging.info(f"Profile file exists for {tissue}_{category}, skipping")
        return matrix_file, profile_file
    
    # Generate matrix only if it doesn't exist
    if not os.path.exists(matrix_file):
        cmd_matrix = f"""
        computeMatrix reference-point \
            -S {' '.join(bigwigs)} \
            -R {bed_file} \
            -b 5000 -a 5000 \
            --referencePoint TSS \
            --skipZeros \
            --numberOfProcessors {threads_per_job} \
            -o {matrix_file}
        """
        subprocess.run(cmd_matrix, shell=True, check=True)
    else:
        logging.info(f"Matrix file exists for {tissue}_{category}, using existing")
    
    # Generate plot
    cmd_plot = f"""
    plotProfile -m {matrix_file} \
        --perGroup \
        --colors red blue green \
        --plotTitle "{tissue} {category}-regulated genes" \
        -out {profile_file}
    """
    
    subprocess.run(cmd_plot, shell=True, check=True)
    
    return matrix_file, profile_file

def main(args):
    setup_logging()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Decompress GTF if needed
    gtf_file = decompress_gtf(args.gtf_file)
    
    # Process tissues
    tissues = ['Neuron', 'NSC'] if not args.test_run else ['Neuron']
    deseq_files = {
        'Neuron': './DATA/DEA_NEU.csv',
        'NSC': './DATA/DEA_NSC.csv'
    }
    
    # Calculate threads per job based on available cores and concurrent jobs
    max_concurrent_jobs = min(len(tissues) * 3, args.max_cores // 2)  # 3 categories per tissue
    threads_per_job = max(1, args.max_cores // max_concurrent_jobs)
    logging.info(f"Using {max_concurrent_jobs} concurrent jobs with {threads_per_job} threads each")
    
    # Process DESeq results and create BED files
    tissue_data = {}
    for tissue in tissues:
        gene_sets = process_deseq_results(deseq_files[tissue])
        categories = ['non', 'up', 'down'] if not args.test_run else ['up']
        
        bed_files = {}
        for category in categories:
            output_file = f"{args.output_dir}/{tissue}_{category}_regulated_genes.bed"
            # Skip if BED file already exists
            if os.path.exists(output_file):
                logging.info(f"BED file exists for {tissue}_{category}, skipping")
                bed_files[category] = output_file
            else:
                bed_files[category] = create_bed_file(gene_sets[category], gtf_file, output_file)
        
        tissue_data[tissue] = {'gene_sets': gene_sets, 'bed_files': bed_files}
    
    # Convert BAM to BigWig in parallel (only for missing files)
    samples = pd.read_csv(args.sample_sheet)
    bigwig_dir = f"{args.output_dir}/bigwig"
    os.makedirs(bigwig_dir, exist_ok=True)
    
    bam_files = samples['bamReads'].tolist()
    if args.test_run:
        bam_files = bam_files[:2]
    
    # Only process BAM files that don't have corresponding BigWig files
    bam_conversion_args = [
        (bam, bigwig_dir, args.genome_size_file, threads_per_job)
        for bam in bam_files
        if not os.path.exists(f"{bigwig_dir}/{Path(bam).stem}.bw")
    ]
    
    if bam_conversion_args:
        with ProcessPoolExecutor(max_workers=max_concurrent_jobs) as executor:
            new_bigwig_files = list(executor.map(bam_to_bigwig, bam_conversion_args))
    
    # Get all bigwig files (both existing and newly created)
    bigwig_files = [f"{bigwig_dir}/{Path(bam).stem}.bw" for bam in bam_files]
    
    # Generate metaprofiles in parallel (only for missing files)
    metaprofile_args = []
    for tissue in tissues:
        tissue_samples = samples[samples['Tissue'] == tissue]
        tissue_bigwigs = [bw for bw, sample in zip(bigwig_files, bam_files)
                         if sample in tissue_samples['bamReads'].tolist()]
        
        for category, bed_file in tissue_data[tissue]['bed_files'].items():
            profile_file = f"{args.output_dir}/{tissue}_{category}_profile.png"
            if not os.path.exists(profile_file):
                metaprofile_args.append((
                    tissue, category, tissue_bigwigs, bed_file,
                    args.output_dir, threads_per_job
                ))
            else:
                logging.info(f"Profile already exists for {tissue}_{category}, skipping")
    
    if metaprofile_args:
        with ProcessPoolExecutor(max_workers=max_concurrent_jobs) as executor:
            results = list(executor.map(generate_metaprofile, metaprofile_args))
    
    logging.info("Processing completed successfully!")

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