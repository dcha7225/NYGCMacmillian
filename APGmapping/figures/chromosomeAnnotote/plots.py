import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap

chromosome_df = pd.read_csv('chromosome_file.txt', sep='\t', header=None, names=['Chromosome', 'Start', 'End'])
alignment_df = pd.read_excel('/grid/chambwe/macmillan/DanielCha/APGmapping/alignmentTables/alignments_with_names.xlsx')

# Calculate overlap density in bins
def calculate_overlap_density_in_bins(df, chrom_length, num_bins):
    bin_size = chrom_length / num_bins
    density = np.zeros(num_bins)
    unique_density = np.zeros(num_bins)
    for _, row in df.iterrows():
        start = int(row['start_pos'] // bin_size)
        end = int(row['end_pos'] // bin_size)
        density[start:end + 1] += 1
        if row['unique_region']:
            unique_density[start:end + 1] += 1
    return density, unique_density

num_bins = 50

# Prepare density data
density_data = []
unique_density_data = []
for _, row in chromosome_df.iterrows():
    chrom_length = row['End']
    chrom_name = row['Chromosome']
    chrom_density, unique_density = calculate_overlap_density_in_bins(alignment_df[alignment_df['mapped_chromosome'] == chrom_name], chrom_length, num_bins)
    density_data.append((chrom_name, chrom_density))
    unique_density_data.append((chrom_name, unique_density))

y_positions = {chromosome: i for i, chromosome in enumerate(chromosome_df['Chromosome'].unique())}

sns.set(style="darkgrid")

# Create the figure and axes
fig, ax = plt.subplots(figsize=(15, 10))

# Choose the 'hot' colormap and reverse it
cmap = plt.get_cmap('gist_heat_r')  # '_r' reverses the colormap

# Plotting each chromosome
for _, row in chromosome_df.iterrows():
    chrom_length = row['End']
    chrom_name = row['Chromosome']
    y_pos = y_positions[chrom_name]

    chrom_density = [d for c, d in density_data if c == chrom_name][0]
    unique_density = [d for c, d in unique_density_data if c == chrom_name][0]

    bin_size = chrom_length / num_bins

    norm = plt.Normalize(0, np.max(chrom_density))

    # Draw the chromosome as a white rectangle with a black outline
    ax.add_patch(patches.Rectangle((0, y_pos - 0.3), chrom_length, 0.6, edgecolor='black', facecolor='white'))

    # Plot the density as colored bars
    for i in range(len(chrom_density)):
        color = cmap(norm(chrom_density[i]))
        start_pos = i * bin_size
        end_pos = (i + 1) * bin_size
        ax.add_patch(patches.Rectangle((start_pos, y_pos - 0.3), end_pos - start_pos, 0.6, edgecolor=(0, 0, 0, 0), facecolor=color))

    # Highlight unique regions with yellow lines
    for i in range(len(unique_density)):
        if unique_density[i] > 0:
            start_pos = i * bin_size
            end_pos = (i + 1) * bin_size
            ax.plot([start_pos, end_pos], [y_pos + 0.3, y_pos + 0.3], color='#39FF14', linewidth=2)

# Setting x-axis limits
max_chrom_length = chromosome_df['End'].max()
ax.set_xlim(0, max_chrom_length*1.1)

# Color mapping
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax)
cbar.set_label('Number of Alignments Mapping to Bin')

# Legend for unique regions
unique_patch = patches.Patch(facecolor='#39FF14', edgecolor='#39FF14', label='Unique Region')
ax.legend(handles=[unique_patch], loc='upper right')

ax.set_yticks(list(y_positions.values()))
ax.set_yticklabels(list(y_positions.keys()))
ax.set_xlabel('Position')
ax.set_ylabel('Chromosome')
ax.set_title('T2T-CHM13 Annotations with APG Contig Alignments and Unique Regions')

plt.tight_layout()

plt.savefig('chromosome_plot_with_density_heatmap_unique_regions.png', dpi=300)

plt.close()
