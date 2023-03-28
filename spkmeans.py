import getopt
import sys
import pandas as pd
import numpy as np
import math
import spkmeans_capi

#kmeans_pp has been updated since HW2! Here, instead of recieving csvs, we recieve a flat list that must be
# reshaped with numpy, then turned into U. The rest of the function works just as before afterwards.
def kmeans_pp(flat_mat, k, vectors_amt): 
    eps=0
    iter=300

    ##turn flat mat into U(n x k) mat##
    vectors = np.reshape(np.array(flat_mat),(vectors_amt+1,vectors_amt))
    vectors[:,vectors[0:1, :].argsort()[0]]
    vectors = vectors[1:, 0:k]        #now vectors == U!!! delete first row (that has eigenvalues) and keep first k column vectors
    keys = np.arange(0,vectors_amt)
    vector_len = k
    
    
    
    ##Start algorithm##
    
    ##step 1##
    centroid_keys = []
    centroids = [] ## python array of numpy arrays
    np.random.seed(0)
    curr_index = np.random.choice(vectors_amt)
    centroid_keys.append(curr_index)
    centroids.append(vectors[curr_index].tolist())

    
    for i in range(1,k):
        ##step 2##
        distances = np.array([min_distance(vector, centroids) for vector in vectors])
        
        ##step 3
        probabilities = np.divide(distances, distances.sum())
        
        ##step 2 but with prob function##
        curr_index = np.random.choice(vectors_amt, p=probabilities) 
        centroid_keys.append(curr_index)
        centroids.append(vectors[curr_index].tolist())
        
    
    for i in range(len(centroid_keys)-1):
        print(f"{int(centroid_keys[i])}," , end="")
    print(f"{int(centroid_keys[-1])}") 
   
    #convert vectors to 2Darray
    centroids = np.array(centroids)
    centroids = centroids.flatten()
    centroids = centroids.tolist()
    

    new_centroids = spkmeans_capi.cKmeans(k, iter, vector_len, vectors_amt, eps, flat_mat, centroids)

    if(new_centroids == None):
        print("An Error Has Occurred")
        exit()
   

    #print the centroids
    for i in range(len(new_centroids)):
        x = "%.4f" % round(new_centroids[i], 4)
        if(i%(vector_len)==vector_len-1):
            print(f"{x}")
        else:    
            print(f"{x}," , end="")
     
   

def euclidian_distance(vec1, vec2):
    sum = 0
    for i in range(np.size(vec1)):
        sum += (vec1[i]-vec2[i])**2
    return sum**(1/2)

def min_distance(vector, centroids):
    min_dis = float('inf')
    for centroid in centroids:
        curr_distance=euclidian_distance(vector,centroid)
        if(curr_distance<min_dis):
           min_dis=curr_distance
    return  min_dis    

def get_heuristic(flat_j_mat, vectors_amt):
    eigenvalues = np.sort(np.array(flat_j_mat[:vectors_amt])) #first row is e"e
    max_gap = 0
    max_k = 0
     
    for i in range((vectors_amt//2)):
        delta = abs(eigenvalues[i]-eigenvalues[i+1])
        if delta > max_gap:
            max_k = i+1
            max_gap = delta
    return max_k

def print_matrix(mat, rows, cols):

    for i in range(rows): 
        for j in range(cols):
            x = "%.4f" % round(mat[i*cols + j], 4)
            if j==cols-1:
                print(x)
            else: 
                print(f"{x}," , end="")
    

def main_func(goal, file_name, k_is_given, k=0):
    
    vectors = np.loadtxt(file_name, delimiter=",", dtype=float)
    
    #get dim2 and n
    vectors_amt = np.shape(vectors)[0]
    vector_len = np.shape(vectors)[1]
    #process to send to API
    vectors = vectors.flatten()
   
    vectors = vectors.tolist()
    

    if(goal == 'wam'):
        wam = spkmeans_capi.wam(vectors, vectors_amt, vector_len)
        if(wam == None):
            print("An Error Has Occurred")
            exit()
               
        print_matrix(wam, vectors_amt, vectors_amt)
        
    elif(goal == 'ddg'):
        
        ddg = spkmeans_capi.ddg(vectors, vectors_amt, vector_len)
        if(ddg == None):
            print("An Error Has Occurred")
            exit()
            
        print_matrix(ddg, vectors_amt, vectors_amt)
        
    elif(goal == 'gl'):
            
        gl = spkmeans_capi.gl(vectors, vectors_amt, vector_len)
        if(gl == None):
            print("An Error Has Occurred")
            exit()
            
        print_matrix(gl, vectors_amt, vectors_amt)
         
    elif(goal == 'jacobi'):
        jacobi = spkmeans_capi.jacobi(vectors, vectors_amt)
        if(jacobi == None):
            print("An Error Has Occurred")
            exit()
        jacobi_rows =  vectors_amt+1  
        print_matrix(jacobi, jacobi_rows, vectors_amt)
         
    elif(goal == 'spk'): 
        gl = spkmeans_capi.gl(vectors, vectors_amt, vector_len)
        jacobi = spkmeans_capi.jacobi(gl, vectors_amt)
        if(jacobi == None):
            print("An Error Has Occurred")
            exit()
            
        if(not k_is_given):
            k = get_heuristic(jacobi, vectors_amt)
            
        #send full matrix, turn it into U in kmeans
        kmeans_pp(jacobi, k, vectors_amt)
        
    else: #if goal is illegal
        print("An Error Has Occurred")
        exit()
            
         
            
argv = sys.argv[1:]
if len(argv)==3:
    k = int(argv[0])
    goal = argv[1]
    file_name = argv[2]
    main_func(goal, file_name, True,  k)

elif len(argv)==2:
    goal = argv[0]
    file_name = argv[1]
    main_func(goal, file_name, False)
    
else:
    print("An Error Has Occurred")
    exit()
    