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

def generate_condition_heatmap(args):
    """Generate a heatmap for a condition (combining replicates)"""
    matrix_file, output_dir, sample_name = args
    
    heatmap_file = os.path.join(output_dir, f"{sample_name}_heatmap.png")
    
    cmd = f"""
    plotHeatmap -m {matrix_file} \
        --colorMap RdBu_r \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 \
        --zMax 0.3 \
        --kmeans 1 \
        --sortRegions descend \
        --heatmapHeight 15 \
        --heatmapWidth 8 \
        --xAxisLabel "" \
        --refPointLabel "TSS" \
        --startLabel "" \
        --endLabel "TES" \
        --legendLocation none \
        --dpi 300 \
        --outFileName {heatmap_file}
    """
    
    logging.info(f"Generating heatmap for {sample_name}")
    subprocess.run(cmd, shell=True, check=True)
    return heatmap_file

def generate_combined_heatmap(args):
    """Generate a combined heatmap for multiple samples of the same category"""
    matrix_files, output_dir, category_name = args
    
    combined_heatmap = os.path.join(output_dir, f"{category_name}_combined_heatmap.png")
    
    matrix_files_str = " ".join(matrix_files)
    cmd = f"""
    plotHeatmap \
        -m {matrix_files_str} \
        --colorMap RdBu_r \
        --whatToShow 'heatmap and colorbar' \
        --zMin 0 \
        --zMax 0.3 \
        --kmeans 1 \
        --sortRegions descend \
        --heatmapHeight 15 \
        --heatmapWidth 8 \
        --xAxisLabel "" \
        --refPointLabel "TSS" \
        --startLabel "" \
        --endLabel "TES" \
        --legendLocation none \
        --dpi 300 \
        --outFileName {combined_heatmap}
    """
    
    logging.info(f"Generating combined heatmap for {category_name}")
    subprocess.run(cmd, shell=True, check=True)
    return combined_heatmap

def main():
    parser = argparse.ArgumentParser(description="Generate heatmaps from precomputed matrix files")
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
    
    # Process matrix files by condition
    matrix_files = []
    for tissue in args.tissues:
        for category in args.categories:
            for condition in args.conditions:
                # Use the condition-specific matrix files
                matrix_path = os.path.join(args.matrix_dir, 
                                         f"matrix_{tissue}_{category}",
                                         f"{tissue}_{category}_{condition}_matrix.gz")
                
                if os.path.exists(matrix_path):
                    matrix_files.append((matrix_path, tissue, category, condition))
                else:
                    logging.warning(f"Matrix file not found: {matrix_path}")
    
    # Generate condition-specific heatmaps
    heatmap_args = [
        (matrix_path, args.output_dir, f"{tissue}_{category}_{condition}")
        for matrix_path, tissue, category, condition in matrix_files
    ]
    
    with ProcessPoolExecutor(max_workers=args.max_cores) as executor:
        individual_heatmaps = list(executor.map(generate_condition_heatmap, heatmap_args))
    
    # Generate combined heatmaps per tissue and category
    combined_args = []
    
    # Combine by tissue (all categories for each tissue)
    for tissue in args.tissues:
        tissue_matrices = [m[0] for m in matrix_files if m[1] == tissue]
        if tissue_matrices:
            combined_args.append((
                tissue_matrices,
                args.output_dir,
                f"{tissue}_all_categories"
            ))
    
    # Combine by category (all tissues for each category)
    for category in args.categories:
        category_matrices = [m[0] for m in matrix_files if m[2] == category]
        if category_matrices:
            combined_args.append((
                category_matrices,
                args.output_dir,
                f"all_tissues_{category}"
            ))
    
    # Generate all combined heatmaps
    with ProcessPoolExecutor(max_workers=args.max_cores) as executor:
        combined_heatmaps = list(executor.map(generate_combined_heatmap, combined_args))
    
    logging.info("Heatmap generation completed successfully!")
    logging.info(f"Individual heatmaps: {len(individual_heatmaps)}")
    logging.info(f"Combined heatmaps: {len(combined_heatmaps)}")

if __name__ == "__main__":
    main()