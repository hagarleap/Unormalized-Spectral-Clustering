import numpy as np 
import pandas as pd



def print_matrix(mat, rows, cols):
    
#    print(f"printing matrix, {rows}, {cols}")
    for i in range(rows):
        for j in range(cols):
            x = "%.4f" % round(mat[i,j], 4)
            if j==cols-1:
                print(j)
                print(x)
            else: 
                print(i)   
                print(f"{x}," , end="")
                

# using loadtxt()
arr = np.loadtxt("test1.txt", delimiter=",", dtype=float)
print("hi")
print_matrix(arr, arr.shape[0], arr.shape[1])