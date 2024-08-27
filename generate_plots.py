import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

from collections import defaultdict


GENERA_OF_INTEREST = ["Salmonella", "Tequatrovirus", "Bacillus", "Escherichia"]
AXES_LABEL_SIZE = 14
FONT_PATH = "./fonts/Roboto-Regular.ttf"

# Define function to compute relative frequencies
def compute_relative_frequencies(data):
    genus_counts = data['genus'].value_counts()
    
    # Calculate frequencies for specified genera and others
    frequencies = {genus: genus_counts.get(genus, 0) for genus in GENERA_OF_INTEREST}
    # We've added the "errors -> ignore" parameter here to make sure that Pandas continues when a genus is not present
    # in the frequencies. Otherwise it does raise an error.
    frequencies['other'] = genus_counts.drop(GENERA_OF_INTEREST, errors="ignore").sum()

    # Normalize frequencies to relative values
    total = sum(frequencies.values())
    relative_frequencies = {genus: count / total if total > 0 else 0 for genus, count in frequencies.items()}
    print(relative_frequencies)
    return relative_frequencies

bar_positions = [0, 1, 2, 4, 5, 6, 8, 9, 10]  # Adding extra spacing after bars 3 and 6

# Define function to plot the data
def plot_relative_frequencies(data_dict, output_dir):
    # Create a font object
    font_prop = fm.FontProperties(fname=FONT_PATH)
    fm.fontManager.addfont(FONT_PATH)

    # Use the font in Matplotlib
    plt.rc('font', family=font_prop.get_name())

    # We generate a new plot file for every mix of files
    for mix, software_data in data_dict.items():

        fig, axes = plt.subplots(nrows=3, ncols=1)
        # fig.suptitle(f'Relative Frequencies for {mix}')

        # We generate every bar that we want to be displayed in the final plot separately
        # for (software, fdr_data) in software_data.items():
        #     df = pd.DataFrame(fdr_data).T
        #     print(f"{df}")
        #     df.plot(kind='barh', stacked=True, ax=ax)
            # print(f"Added barchart for {software}, {fdr_data}")

        for i, (ax, (software, fdr_data)) in enumerate(zip(axes, software_data.items())):
            df = pd.DataFrame(fdr_data).T
            df = df.rename(index={"05": "5%", "01": "1%", "001": "0.1%"})
            df = df.sort_index(ascending=False)
            df.plot(kind='barh', stacked=True, ax=ax, cmap='Accent')  # Create horizontal bar plot
            ax.set_title(software)

            if (i + 1) == len(software_data):
                ax.set_xlabel('Relative Frequency', fontsize=AXES_LABEL_SIZE)

            ax.get_legend().remove()

        fig.supylabel('FDR Level', fontsize=AXES_LABEL_SIZE)

        plt.legend(
            GENERA_OF_INTEREST + ["Other"],
            bbox_to_anchor=(0.85, -0.6),
            ncol=3
        )

        plt.subplots_adjust(hspace=0.75)

        # plt.tight_layout(rect=(0, 0, 1, 1))
        plt.savefig(os.path.join(output_dir, f'{mix}_relative_frequencies.png'), bbox_inches='tight', dpi=300)
        plt.savefig(os.path.join(output_dir, f'{mix}_relative_frequencies.eps'), bbox_inches='tight', dpi=300)
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
input_directory = './lcas'
output_directory = './plots'
main(input_directory, output_directory)
