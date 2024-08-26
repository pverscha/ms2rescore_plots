import os
import pandas as pd
import re

def extract_peptide_sequence(peptide):
    """
    Extracts the sequence of letters between two dots, ignoring text between square brackets.
    """
    # Use a regular expression to find the sequence between dots, excluding content within brackets
    match = re.search(r'\.([^.[]*?(?:\[[^]]*\][^.[]*?)*)\.', peptide)
    if match:
        # Remove content within square brackets
        cleaned_sequence = re.sub(r'\[.*?\]', '', match.group(1))
        return cleaned_sequence
    return ''

def filter_and_save_peptides(input_directory, output_directory):
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)
    
    # List all txt files in the input directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".txt"):
            filepath = os.path.join(input_directory, filename)
            
            # Read the tab-separated file into a DataFrame
            df = pd.read_csv(filepath, sep='\t')
            
            # Ensure the necessary columns are present
            if 'peptide' in df.columns and 'q-value' in df.columns:
                # Define q-value thresholds and corresponding suffixes
                thresholds = [0.001, 0.01, 0.05]
                suffixes = ['0_001', '0_01', '0_05']
                
                for threshold, suffix in zip(thresholds, suffixes):
                    # Filter rows where q-value is <= threshold
                    filtered_df = df[df['q-value'] <= threshold].copy()
                    
                    # Extract and clean the peptide sequences
                    filtered_df['peptide'] = filtered_df['peptide'].apply(extract_peptide_sequence)
                    
                    # Create output filename
                    output_filename = f"{os.path.splitext(filename)[0]}.{suffix}.txt"
                    output_filepath = os.path.join(output_directory, output_filename)
                    
                    # Save the filtered peptides to a new file, ignoring empty sequences
                    filtered_df[filtered_df['peptide'] != '']['peptide'].to_csv(output_filepath, index=False, header=False)
                    
                print(f"Processed {filename} and created filtered files in {output_directory}.")
            else:
                print(f"Skipped {filename}: Required columns not found.")

# Specify the directory containing your input txt files
input_directory = './../ms2rescore_files'

# Specify the directory where the output files should be saved
output_directory = './../ms2rescore_peptides'

# Call the function to process the files
filter_and_save_peptides(input_directory, output_directory)
