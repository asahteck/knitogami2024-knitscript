import numpy as np


def parse_knit_file(filename):
    '''Take in .txt file from pattern output, parses K/P as 0/255,
    copies for dimension size
    and outputs as a matrix'''
    with open(filename, 'r') as file:
        dim_info = file.readlines()[0]
        x_dim = dim_info[2]
        y_dim = dim_info[4]
    with open(filename, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line


    matrix = []
    for j in range(0, int(y_dim)):
        for line in lines:
            row = [255 if char == 'P' else 0 for char in line.split()]
            row = row*int(x_dim)
            matrix.append(row)

    return np.array(matrix)

# Example usage
#filename = 'patterns/knit_3a_miura.txt'
#array = parse_knit_file(filename)
#print(array)
#print(len(array))
#print(len(array[0]))