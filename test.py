import numpy as np 
import pandas as pd



    
#    print(f"printing matrix, {rows}, {cols}")
def print_matrix(mat, rows, cols):

    for i in range(rows): 
        for j in range(cols):
            x = "%.4f" % round(mat[i*cols + j], 4)
            if j==cols-1:
                print(x)
            else: 
                print(f"{x}," , end="")
                

# using loadtxt()
arr = np.loadtxt("test1.txt", delimiter=",", dtype=float)
rows = arr.shape[0]
cols = arr.shape[1]
arr = arr.flatten()

print_matrix(arr, rows,cols )