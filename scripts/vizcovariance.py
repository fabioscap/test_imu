import numpy as np
import matplotlib.pyplot as plt

import sys

# Function to read the matrix from a text file
def read_matrix_from_file(filename):
    matrix = []
    with open(filename, 'r') as file:
        for line in file:
            row = line.strip().split()
            matrix.append([float(value) for value in row])
    return np.array(matrix)

# Main function
if __name__ == '__main__':
    input_file = '../examples/covariance_matrix.txt'  # Replace with your input text file name
    matrix = read_matrix_from_file(input_file)
    if len(sys.argv) == 1:

        plt.matshow(matrix)

    elif sys.argv[1] == "normalize":
        max_value = np.max(np.abs(matrix))
        plt.matshow((matrix)/(max_value))
    elif sys.argv[1] == "binary" :
        plt.matshow(matrix!= 0)

    plt.show()
