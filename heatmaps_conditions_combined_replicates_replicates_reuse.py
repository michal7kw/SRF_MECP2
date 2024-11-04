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
    
    try:
        heatmap_file = os.path.join(output_dir, f"{sample_name}_heatmap.png")
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
            "--outFileName", heatmap_file,
            "--averageType", "mean",      # Average across samples
            "--perGroup"                  # Show one heatmap per group instead of per sample
        ]
        
        logging.info(f"Generating heatmap for {sample_name}")
        result = subprocess.run(cmd, check=False, capture_output=True, text=True)
        
        if result.returncode != 0:
            logging.error(f"Failed to generate heatmap for {sample_name}")
            logging.error(f"Error output: {result.stderr}")
            return None
        
        return heatmap_file
        
    except Exception as e:
        logging.error(f"Error processing {sample_name}: {str(e)}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Generate heatmaps from combined replicate matrices")
    parser.add_argument("--matrix-dir", required=True, help="Directory containing combined matrix files")
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
    
    # Collect combined matrix files
    matrix_files = []
    for tissue in args.tissues:
        for category in args.categories:
            for condition in args.conditions:
                # Use the naming convention from the combined matrices
                matrix_name = f"{tissue}_{condition}_{category}_combined_matrix.gz"
                matrix_path = os.path.join(args.matrix_dir, matrix_name)
                
                if os.path.exists(matrix_path):
                    logging.info(f"Found matrix file: {matrix_path}")
                    matrix_files.append((
                        matrix_path,
                        args.output_dir,
                        f"{tissue}_{condition}_{category}"
                    ))
                else:
                    logging.warning(f"Combined matrix file not found: {matrix_path}")
    
    if not matrix_files:
        logging.error("No combined matrix files found!")
        return
    
    # Generate heatmaps for combined data
    with ProcessPoolExecutor(max_workers=args.max_cores) as executor:
        heatmaps = list(executor.map(generate_single_heatmap, matrix_files))
    
    # Filter out None values from failed heatmap generations
    successful_heatmaps = [h for h in heatmaps if h is not None]
    
    logging.info("Heatmap generation completed!")
    logging.info(f"Successfully generated heatmaps: {len(successful_heatmaps)} out of {len(matrix_files)}")
    
    # List all generated heatmaps
    for heatmap in successful_heatmaps:
        logging.info(f"Generated: {os.path.basename(heatmap)}")

if __name__ == "__main__":
    main()