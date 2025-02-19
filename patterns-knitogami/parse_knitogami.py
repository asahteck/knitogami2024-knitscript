import numpy as np


def parse_knit_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line

    matrix = []
    for line in lines:
        row = [255 if char == 'P' else 0 for char in line.split()]
        matrix.append(row)

    return np.array(matrix)


# Example usage
#filename = 'patterns-knitogami/knit_3a_miura.txt'
#array = parse_knit_file(filename)
#print(array)