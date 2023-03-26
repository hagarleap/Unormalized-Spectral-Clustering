import getopt
import sys
import pandas as pd
import numpy as np
import math
import spkmeans_capi

def kmeans_pp(df, k): #work on thisgit
    eps=0
    iter=300

    
    vectors_amt = df_1.shape[0]
  
    ##check that the amount of clusters is legal##
    if(k>vectors_amt):
        print("Invalid number of clusters!")
        exit()

    
    ##convert dataframe into numpy array, without keys!##
    keys = df_vectors[0].to_numpy()
    vectors = df_vectors.drop(0, axis=1).to_numpy()
    vector_len = np.size(vectors[0])
    
    ##Start algorithm##
    
    ##step 1##
    centroid_keys = []
    centroids = [] ## python array of numpy arrays
    np.random.seed(0)
    curr_index = np.random.choice(np.size(keys))
    centroid_keys.append(keys[curr_index])
    centroids.append(vectors[curr_index].tolist())

    
    for i in range(1,k):
        ##step 2##
        distances = np.array([min_distance(vector, centroids) for vector in vectors])
        
        ##step 3
        probabilities = np.divide(distances, distances.sum())
        
        ##step 2 but with prob function##
        curr_index = np.random.choice(np.size(keys), p=probabilities) 
        centroid_keys.append(keys[curr_index])
        centroids.append(vectors[curr_index].tolist())
        
    
    for i in range(len(centroid_keys)-1):
        print(f"{int(centroid_keys[i])}," , end="")
    print(f"{int(centroid_keys[-1])}") 
   
    #convert vectors to 2Darray
    vectors = vectors.flatten()
    vectors = vectors.tolist()
    centroids = np.array(centroids)
    centroids = centroids.flatten()
    centroids = centroids.tolist()
    

    new_centroids = kmeans_capi.cKmeans(k, iter, vector_len, vectors_amt, eps, vectors, centroids)

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


def main_func(goal, file_name, k_is_given, k=0):
           
    file = open(file_name, 'r')
    df = pd.read_csv(file_name,header=None)
    file.close()
    
    vectors = df.to_numpy()
    
    #get dim2 and n
    vectors_amt = np.shape(vectors)[0]
    vector_len = np.shape(vectors)[1]
    
    #process to send to API
    vectors = vectors.flatten()
    vectors = vectors.tolist()
    
    if(goal == 'wam'):
        wam = spkmeans_capi.wam(vectors, vectors_amt, vector_len)
        print_matrix(wam, vectors_amt, vector_len)
        
    elif(goal == 'ddg'):
        wam = spkmeans_capi.wam(vectors, vectors_amt, vector_len)
        ddg = spkmeans_capi.ddg(wam, vectors_amt)
        print_matrix(ddg, vectors_amt, vector_len)
        
    elif(goal == 'gl'):
        wam = spkmeans_capi.wam(vectors, vectors_amt, vector_len)
        ddg = spkmeans_capi.ddg(wam, vectors_amt)
        gl = spkmeans_capi.gl(wam, ddg, vectors_amt)
        print_matrix(gl, vectors_amt, vector_len)
         
    elif(goal == 'jacobi'):
        jacobi = spkmeans_capi.jacobi(vectors, vectors_amt)
        print_matrix(jacobi, vectors_amt, vector_len)
         
    elif(goal == 'spk'):
        wam = spkmeans_capi.wam(vectors, vectors_amt, vector_len)
        ddg = spkmeans_capi.ddg(wam, vectors_amt)
        gl = spkmeans_capi.gl(wam, ddg, vectors_amt)
        jacobi = spkmeans_capi.jacobi(gl, vectors_amt)
        if(not k_is_given):
            k = spkmeans_capi.get_heuristic(jacobi, vectors_amt)
        U = get_first_k_eigenvectors(jacobi)
        kmeans_pp(U, k)
    else:
            
 
# double** wam_func(double** data_matrix, int n, int dim2);
# double** ddg_func( double** weight_mat, int n);
# double** gl_func(double** weight_mat, double** diag_degree_mat, int n);
# double** jacobi_func(double** A, int n)          
            
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
    