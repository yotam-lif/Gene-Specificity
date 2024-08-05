import numpy as np
import pandas as pd


def first_order_tau(X: np.array, yotau=False):
    """
    Calculates the first order tau value for a given gene expression vector X
    :param X: Expression vector
    :param yotau: Boolean; Can calculate the modified yotau score
    :return: Tau score
    """
    N = len(X)
    X = 1 - X
    if yotau:
        return np.sqrt(np.sum(X ** 2) / (N - 1))
    else:
        return sum(X) / (N - 1)


def second_order_tau(X: np.array, C: np.array):
    """
    Calculates the second order tau value for a given gene expression vector X
    :param X: Expression vector
    :param C: Matrix of distances. Can be correlation matrix or other.
    :return: Second order tau score
    """
    # return X @ C @ X
    norm = np.dot(X, X)
    if norm == 0:
        return 0
    else:
        return X @ C @ X / norm


def create_expression_matrix(file_path, row_key: str, column_key: str, log_transform=True):
    """
    Create an expression matrix from a file
    :param file_path: Path to the data file
    :param row_key: Key to use for rows
    :param column_key: Key to use for columns
    :return: Expression matrix and list of row keys
    """

    data = pd.read_csv(file_path, sep='\t')

    # Pivot the table to create a matrix
    # each row corresponds to 'row_key' and each column corresponds to 'column_key
    # For instance, row_key can be 'Gene Name' and 'column_key' can be 'Tissue'
    pivot_data = data.pivot_table(index=row_key, columns=column_key, values='nTPM', aggfunc='mean').fillna(0)

    # Convert nTPM values to log10(x+1)
    Xi = pivot_data.values
    if log_transform:
        Xi = np.log10(Xi + 1)

    # Normalize each row by the maximum value in that row, avoiding division by zero
    max_values = Xi.max(axis=1, keepdims=True)
    max_values[max_values == 0] = 1  # Replace zeros with ones to avoid NaN
    Xi = Xi / max_values

    return Xi, pivot_data.index.tolist()


def Pearson_Matrix(X_matrix):
    # Calculate the Pearson correlation matrix between the columns of X_matrix
    C = np.corrcoef(X_matrix, rowvar=False)

    # Make it zero diagonal
    np.fill_diagonal(C, 0)

    return C
