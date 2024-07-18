# Gene-Specificity

This project performs analysis on gene expression data, specifically calculating and visualizing various specificity scores from an RNA tissue expression dataset.

## Project Structure

- `create_Xi(file_path)`: Reads the data, performs log transformation, normalizes each row, and returns a matrix where each row corresponds to a gene and columns correspond to tissues.
- `calculate_pearson_correlation_matrix(X_matrix)`: Calculates the Pearson correlation matrix.
- Calculation of different Tau values:
  - `first_order_tau`
  - `first_order_yotau`
  - `second_order_tau`
- Saving Tau values to JSON files.
- Plotting histograms and scatter plots for the Tau values.

## Requirements

- Python 3
- numpy
- pandas
- matplotlib

## Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/yotam-lif/Gene-Specificity.git
    cd gene-expression-analysis
    ```

2. Install the required Python packages:
    ```sh
    pip install numpy pandas matplotlib
    ```

## Usage

1. Place your RNA data file in the 'Data_Files' directory.

2. Run an analysis script:
    ```sh
    python calc_reg_tau
    ```

3. The script will generate and save Json files with gene name keys and your choice of tau score as values. It will also generate histograms of the different tau scores for all genes in the dataset. These files will be save in 'Outputs' directory.
   
