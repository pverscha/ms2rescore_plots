import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

# Define function to compute relative frequencies
def compute_relative_frequencies(data):
    genera_of_interest = ["Salmonella", "Tequatrovirus", "Bacillus", "Escherichia"]
    genus_counts = data['genus'].value_counts()
    
    # Calculate frequencies for specified genera and others
    frequencies = {genus: genus_counts.get(genus, 0) for genus in genera_of_interest}
    frequencies['other'] = genus_counts.drop(genera_of_interest).sum()

    # Normalize frequencies to relative values
    total = sum(frequencies.values())
    relative_frequencies = {genus: count / total for genus, count in frequencies.items()}
    return relative_frequencies

# Define function to plot the data
def plot_relative_frequencies(data_dict, output_dir):
    for mix, software_data in data_dict.items():
        fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
        fig.suptitle(f'Relative Frequencies for {mix}')

        for ax, (software, fdr_data) in zip(axes, software_data.items()):
            df = pd.DataFrame(fdr_data).T
            df.plot(kind='barh', stacked=True, ax=ax)  # Create horizontal bar plot
            ax.set_title(software)
            ax.set_ylabel('FDR Level')
            ax.set_xlabel('Relative Frequency')
            ax.legend(title='Genus')

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(os.path.join(output_dir, f'{mix}_relative_frequencies.png'))
        plt.close()

# Main script
def main(input_directory, output_directory):
    data_dict = defaultdict(lambda: defaultdict(dict))
    
    for filename in os.listdir(input_directory):
        if filename.endswith('.tsv'):
            # Extract mix, software, and FDR level from filename
            parts = filename.split('-')
            mix = parts[0]
            software_fdr = parts[1].split('.')[0]
            fdr_level = parts[1].split('_')[1]  # Extract FDR part after '_'

            # Read file and compute relative frequencies
            filepath = os.path.join(input_directory, filename)
            data = pd.read_csv(filepath, sep='\t')
            relative_frequencies = compute_relative_frequencies(data)

            # Store the results
            data_dict[mix][software_fdr][fdr_level] = relative_frequencies

    # Plot the results
    plot_relative_frequencies(data_dict, output_directory)

# Example usage
input_directory = './example_files'
output_directory = './plots'
main(input_directory, output_directory)
