import pandas as pd
import os
import json
import matplotlib.pyplot as plt
import Funcs as F


data_tissue_name = "rna_tissue_hpa.tsv"
yotau = False
row_key = 'Gene name'
column_key = 'Tissue'

# Define the paths
base_dir = os.path.dirname(os.path.abspath(__file__))
data_tissue_path = os.path.join(base_dir, 'Data_Files', data_tissue_name)

# Create the Xi matrix
if yotau:
    tau_name = "Yotau"
else:
    tau_name = "Tau"

Xi_matrix, gene_names = F.create_expression_matrix(data_tissue_path, row_key=row_key, column_key=column_key)
Xi_nonlog_mat, _ = F.create_expression_matrix(data_tissue_path, row_key=row_key, column_key=column_key, log_transform=False)
C = F.Pearson_Matrix(Xi_nonlog_mat)

# Calculate tau values
tau_values = {}

for i, gene in enumerate(gene_names):
    Xi = Xi_matrix[i]
    tau_values[gene] = F.second_order_tau(Xi, C)

max_val = max(tau_values.values())
# Normalize the tau values
tau_values = {gene: tau / max_val for gene, tau in tau_values.items()}
# Save the tau values to a JSON file
output_file_path = os.path.join(base_dir, 'Outputs', f'tau2_values.json')

with open(output_file_path, 'w') as f:
    json.dump(tau_values, f)

print(f"Values saved to {output_file_path}")

# Plot the histogram of tau values
plt.figure(figsize=(10, 6))
plt.hist(tau_values.values(), bins=100, edgecolor='black')
plt.title('$\\tau_2$')
plt.xlabel('Score')
plt.ylabel('Frequency')

# Save the plot
histogram_path = os.path.join(base_dir, 'Outputs', f'tau2_values_histogram.png')
plt.savefig(histogram_path)
