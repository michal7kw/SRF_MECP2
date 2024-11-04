#!/bin/bash

# Function to rename files
rename_files() {
    # Endogenous files
    for file in *Endogenous*matrix.gz; do
        if [[ -f "$file" ]]; then
            # Extract components
            if [[ $file =~ (Neuron|NSC)_(up|down)_Endogenous_matrix\.gz ]]; then
                cell_type=${BASH_REMATCH[1]}
                direction=${BASH_REMATCH[2]}
                new_name="${cell_type}_Endogenous_${direction}_combined_matrix.gz"
                mv -v "$file" "$new_name"
            fi
        fi
    done

    # Exogenous files
    for file in *Exogenous*matrix.gz; do
        if [[ -f "$file" ]]; then
            # Extract components
            if [[ $file =~ (Neuron|NSC)_(up|down)_Exogenous_matrix\.gz ]]; then
                cell_type=${BASH_REMATCH[1]}
                direction=${BASH_REMATCH[2]}
                new_name="${cell_type}_Exogenous_${direction}_combined_matrix.gz"
                mv -v "$file" "$new_name"
            fi
        fi
    done

    # Handle shortened versions (endo/exo)
    for file in *endo*matrix.gz; do
        if [[ -f "$file" ]]; then
            if [[ $file =~ (Neuron|NSC)_(up|down)_endo_matrix\.gz ]]; then
                cell_type=${BASH_REMATCH[1]}
                direction=${BASH_REMATCH[2]}
                new_name="${cell_type}_Endogenous_${direction}_combined_matrix.gz"
                mv -v "$file" "$new_name"
            fi
        fi
    done

    for file in *exo*matrix.gz; do
        if [[ -f "$file" ]]; then
            if [[ $file =~ (Neuron|NSC)_(up|down)_exo_matrix\.gz ]]; then
                cell_type=${BASH_REMATCH[1]}
                direction=${BASH_REMATCH[2]}
                new_name="${cell_type}_Exogenous_${direction}_combined_matrix.gz"
                mv -v "$file" "$new_name"
            fi
        fi
    done
}

# Create a backup of files before renaming
mkdir -p backup
cp *matrix.gz backup/

# Execute the renaming function
rename_files