import pandas as pd
import os
import json
import matplotlib.pyplot as plt
import Funcs as F


data_tissue_name = "rna_tissue_hpa.tsv"
yotau = True
row_key = 'Gene name'
column_key = 'Tissue'

# Define the paths
base_dir = os.path.dirname(os.path.abspath(__file__))
data_tissue_path = os.path.join(base_dir, 'Data_Files', data_tissue_name)

# Create the Xi matrix
# Xi_matrix, gene_names = create_Xi(data_tissue_path)
if yotau:
    tau_name = "Yotau"
else:
    tau_name = "Tau"

Xi_matrix, gene_names = F.create_expression_matrix(data_tissue_path, row_key=row_key, column_key=column_key)

# Calculate tau values
tau_values = {}

for i, gene in enumerate(gene_names):
    Xi = Xi_matrix[i]
    tau_values[gene] = F.first_order_tau(Xi, yotau=yotau)

# Save the tau values to a JSON file
output_file_path = os.path.join(base_dir, 'Outputs', f'{tau_name}_values.json')

with open(output_file_path, 'w') as f:
    json.dump(tau_values, f)

print(f"Values saved to {output_file_path}")

# Plot the histogram of tau values
plt.figure(figsize=(10, 6))
plt.hist(tau_values.values(), bins=100, edgecolor='black')
plt.title(f'{tau_name} Values')
plt.xlabel('Score')
plt.ylabel('Frequency')

# Save the plot
histogram_path = os.path.join(base_dir, 'Outputs', f'{tau_name}_values_histogram.png')
plt.savefig(histogram_path)
