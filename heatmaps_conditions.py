import argparse
import subprocess
import os
import logging
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

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
    
    args = parser.parse_args()
    setup_logging()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Collect all matrix files with conditions
    matrix_files = []
    for tissue in args.tissues:
        for category in args.categories:
            matrix_path = os.path.join(args.matrix_dir, 
                                     f"matrix_{tissue}_{category}",
                                     f"{tissue}_{category}_matrix.gz")
            if os.path.exists(matrix_path):
                matrix_files.append((matrix_path, tissue, category))
            else:
                logging.warning(f"Matrix file not found: {matrix_path}")
    
    # Generate individual heatmaps
    heatmap_args = [
        (matrix_path, args.output_dir, f"{tissue}_{category}")
        for matrix_path, tissue, category in matrix_files
    ]
    
    with ProcessPoolExecutor(max_workers=args.max_cores) as executor:
        individual_heatmaps = list(executor.map(generate_single_heatmap, heatmap_args))
    
    logging.info("Heatmap generation completed successfully!")
    logging.info(f"Individual heatmaps: {len(individual_heatmaps)}")

if __name__ == "__main__":
    main()