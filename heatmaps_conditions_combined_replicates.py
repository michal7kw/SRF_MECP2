import argparse
import subprocess
import os
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def generate_single_heatmap(args):
    """Generate a heatmap for a single matrix file"""
    matrix_file, output_dir, sample_name = args
    
    heatmap_file = os.path.join(output_dir, f"{sample_name}_heatmap.png")
    
    # Create a title from sample name by replacing underscores with spaces and capitalizing
    title = sample_name.replace('_', ' ').title()
    
    cmd = [
        "plotHeatmap",
        "-m", matrix_file,
        "--colorMap", "RdBu_r",
        "--whatToShow", "heatmap and colorbar",
        "--zMin", "0",
        "--zMax", "5",
        "--sortRegions", "descend",
        "--heatmapHeight", "15",
        "--heatmapWidth", "8",
        "--xAxisLabel", "",
        "--startLabel", "5kb upstream",
        "--endLabel", "5kb downstream",
        "--legendLocation", "center-right",
        "--dpi", "300",
        "--plotTitle", title,
        "--outFileName", heatmap_file
    ]
    
    logging.info(f"Generating heatmap for {sample_name}")
    subprocess.run(cmd, check=True)
    return heatmap_file

def combine_bigwig_signals(bigwig_files, bed_file, output_dir, sample_prefix):
    """Combine multiple bigwig files into a single matrix"""
    # First, ensure all bigwig files are converted to strings
    bigwig_files = [str(bw) for bw in bigwig_files]
    
    # Extract the bed file from the matrix directory
    bed_path = os.path.join(os.path.dirname(bed_file), "regions.bed")
    
    cmd = [
        "computeMatrix", "scale-regions",
        "-S"] + bigwig_files + [
        "-R", bed_path,
        "--beforeRegionStartLength", "5000",
        "--regionBodyLength", "5000",
        "--afterRegionStartLength", "5000",
        "--skipZeros",
        "-o", os.path.join(output_dir, f"{sample_prefix}_combined_matrix.gz")
    ]
    
    logging.info(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return os.path.join(output_dir, f"{sample_prefix}_combined_matrix.gz")

def main():
    parser = argparse.ArgumentParser(description="Generate heatmaps from Cut&Tag matrix files")
    parser.add_argument("--matrix-dir", required=True, help="Directory containing matrix files")
    parser.add_argument("--output-dir", required=True, help="Output directory for heatmaps")
    parser.add_argument("--tissues", nargs="+", default=["Neuron", "NSC"], 
                        help="List of tissues to process")
    parser.add_argument("--categories", nargs="+", default=["non", "up", "down"],
                        help="List of regulation categories")
    parser.add_argument("--conditions", nargs="+", default=["Endogenous", "Exogenous"],
                        help="List of conditions to process")
    parser.add_argument("--max-cores", type=int, default=4,
                        help="Maximum number of cores to use")
    parser.add_argument("--bigwig-dir", required=True, help="Directory containing bigwig files")
    
    args = parser.parse_args()
    setup_logging()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define replicate patterns
    replicate_patterns = {
        'Neuron': {
            'Endogenous': 'NeuV[123].bw',
            'Exogenous': 'NeuM[23].bw'
        },
        'NSC': {
            'Endogenous': 'NSCv[123].bw',
            'Exogenous': 'NSCM[123].bw'
        }
    }
    
    # Collect and process files
    matrix_files = []
    for tissue in args.tissues:
        for category in args.categories:
            matrix_dir = os.path.join(args.matrix_dir, f"matrix_{tissue}_{category}")
            bed_path = os.path.join(matrix_dir, "regions.bed")
            
            if not os.path.exists(bed_path):
                logging.warning(f"BED file not found: {bed_path}")
                continue
                
            for condition in args.conditions:
                pattern = replicate_patterns[tissue][condition]
                bigwig_files = sorted(Path(args.bigwig_dir).glob(pattern))
                
                if not bigwig_files:
                    logging.warning(f"No bigwig files found for pattern: {pattern}")
                    continue
                
                sample_name = f"{tissue}_{condition}_{category}"
                logging.info(f"Processing {sample_name} with {len(bigwig_files)} replicates")
                
                # Combine replicates
                combined_matrix = combine_bigwig_signals(
                    bigwig_files, bed_path, args.output_dir, sample_name
                )
                matrix_files.append((combined_matrix, args.output_dir, sample_name))
    
    # Generate heatmaps for combined data
    with ProcessPoolExecutor(max_workers=args.max_cores) as executor:
        heatmaps = list(executor.map(generate_single_heatmap, matrix_files))
    
    logging.info("Combined heatmap generation completed successfully!")
    logging.info(f"Generated heatmaps: {len(heatmaps)}")

if __name__ == "__main__":
    main()