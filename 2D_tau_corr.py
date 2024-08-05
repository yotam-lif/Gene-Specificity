import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Assuming the existence of your functions in Funcs as F
import Funcs as F

# Define the parameters
data_tissue_name = "rna_tissue_hpa.tsv"
yotau = True
row_key = 'Gene name'
column_key = 'Tissue'
nbins = 20

# Define the paths
base_dir = os.path.dirname(os.path.abspath(__file__))
data_tissue_path = os.path.join(base_dir, 'Data_Files', data_tissue_name)

# Create expression matrix and calculate Pearson matrix
Xi_matrix, gene_names = F.create_expression_matrix(data_tissue_path, row_key=row_key, column_key=column_key)
C = F.Pearson_Matrix(Xi_matrix)

# Calculate tau values
tau_values = {}
for i, gene in enumerate(gene_names):
    Xi = Xi_matrix[i]
    tau_values[gene] = F.first_order_tau(Xi, yotau=yotau)

# Calculate the second order tau values
tau2_values = {}
for i, gene in enumerate(gene_names):
    Xi = Xi_matrix[i]
    tau2_values[gene] = F.second_order_tau(Xi, C)

# Extract tau and tau2 values
tau_list = np.array(list(tau_values.values()))
tau2_list = np.array(list(tau2_values.values()))

# Normalize tau and tau2 values
tau_list /= np.max(tau_list)
tau2_list /= np.max(tau2_list)

# Create a 2D histogram for tau vs tau2
heatmap_data, xedges, yedges = np.histogram2d(tau_list, tau2_list, bins=nbins, range=[[0, 1], [0, 1]], density=True)

# Plot the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(heatmap_data.T, cmap="viridis", cbar=True)

# Adjusting the ticks to ensure (0,0) is at the bottom left corner
ticks = np.linspace(0, 1, num=nbins + 1)
plt.xticks(np.arange(len(ticks)), np.around(ticks, 2))
plt.yticks(np.arange(len(ticks)), np.flip(np.around(ticks, 2)))

plt.xlabel(r'$\tau$')
plt.ylabel(r'$\tau_{2}$')
plt.title(r'Heatmap of $\tau$ vs $\tau_{2}$')
plt.tight_layout()

# Save the figure
figure_path = os.path.join(base_dir, 'tau_vs_tau2_heatmap_refined.png')
plt.savefig(figure_path, dpi=300)

# Show the plot
plt.show()