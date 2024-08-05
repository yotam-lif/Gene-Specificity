import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mpl_scatter_density  # Import the scatter density projection
from matplotlib.colors import LinearSegmentedColormap

data_tissue_name = "combinedMDLfeatures_v3.xlsx"

# Define the paths
base_dir = os.path.dirname(os.path.abspath(__file__))
data_tissue_path = os.path.join(base_dir, 'Data_Files', data_tissue_name)
output_dir = os.path.join(base_dir, 'Outputs')

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Load the Excel file
df = pd.read_excel(data_tissue_path)

# Drop the 'GeneStableID' and 'GeneName' columns as they contain string values which cannot be converted to float
df_numeric = df.drop(columns=['GeneStableID', 'GeneName'])

# Extract columns excluding TAU score
columns = df_numeric.columns.difference(['TAUscore_SingleCell'])

horizontal_data = df_numeric['TAUscore_SingleCell'].values

# "Viridis-like" colormap with white background
white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
    (0, '#ffffff'),
    (1e-20, '#440053'),
    (0.2, '#404388'),
    (0.4, '#2a788e'),
    (0.6, '#21a784'),
    (0.8, '#78d151'),
    (1, '#fde624'),
], N=256)


def using_mpl_scatter_density(fig, x, y):
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    density = ax.scatter_density(x, y, cmap=white_viridis)
    fig.colorbar(density, label='Number of points per pixel')


# Generate and save a 2D scatter density plot for each feature
for column in columns:
    fig = plt.figure(figsize=(14, 8))
    vert_data = df_numeric[column].values

    # Remove the maximum value from vert_data
    max_value_index = np.argmax(vert_data)
    vert_data = np.delete(vert_data, max_value_index)
    horizontal_data_adj = np.delete(horizontal_data, max_value_index)

    using_mpl_scatter_density(fig, horizontal_data_adj, vert_data)
    plt.title(f'2D Scatter Density Plot of {column} vs TAU Score')
    plt.xlabel('TAU Score')
    plt.ylabel(f'{column}')
    plt.tight_layout()

    # Save the figure
    output_path = os.path.join(output_dir, f'scatter_density_{column}_vs_tau.png')
    plt.savefig(output_path)
    plt.close()