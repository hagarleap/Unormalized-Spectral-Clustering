#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

//returns matrix of zeros
double** matrix_maker(int dim1, int dim2){
    double **mat;
    int i, j;
    mat = malloc(dim1*sizeof(double*));
    if(mat==NULL){
        return NULL;
    }
    for(i=0; i<dim1; i++){
        mat[i] = calloc(dim2, sizeof(double));
        if(mat==NULL){
            return NULL;
        }        
    }
}

//Recives 2 vectors. Returns the euclidian distance between them
double Squared_Euclidean_distance(double* vec1, double* vec2, int vector_len){
    double sum=0;
    int j;
    for(j=0; j<vector_len; j++)
      { 
        sum += pow((vec1[j]-vec2[j]), 2);
      }
    return sum;
}

//free matrix data!
void free_matrix(double** mat, int dim1){
    for(int i=0; i<dim1; i++){
        free(mat[i]);
    }
    free(mat);
}

// Calculate and output the Weighted Adjacency Matrix
double** wam(double** data_matrix, int n){
    double** weight_mat;
    int i, j;

    //create matrix
    weight_mat = matrix_maker(n,n);
    if(weight_mat==NULL){
        free_matrix(weight_mat, n, n);
        free_matrix(data_matrix, n, n);
        return NULL;
    }
    // calculate weights
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i!=j){
                weight_mat[i][j] = exp(-(Squared_Euclidean_distance(weight_mat[i], weight_mat[j], n))/2);
            }
    
        }
    }
    return weight_mat;
}

// Calculate and output the Diagonal Degree Matrix
 double** ddg(double** data_matrix, double** weight_mat, int n){
    double **diag_degree_mat;
    int i, j;

    //create matrix
    diag_degree_mat = matrix_maker(n,n);
    if(diag_degree_mat==NULL){
        free_matrix(diag_degree_mat, n, n);
        free_matrix(data_matrix, n, n);
        free_matrix(weight_mat, n, n);
        return NULL;
    }


    //create diagonal mat
    for(i=0; i<n; i++){
        for (j=0; j<n; j++){
            diag_degree_mat[i][i] = diag_degree_mat[i][i]+ weight_mat[i][j];
        }

    }
    return diag_degree_mat;
 }


//Calculate and output the Graph Laplacian
double** gl(double** weight_mat, double** diag_degree_mat, int n){
    double **laplacian_mat;
    int i, j;

    //create output matrix and 
    laplacian_mat = matrix_maker(n,n);
    if(laplacian_mat==NULL){
        free_matrix(diag_degree_mat, n, n);
        free_matrix(weight_mat, n, n);
        return NULL;
    }

    for(i=0; i<n; i++){
        for (j=0; j<n; j++){
            laplacian_mat[i][j]= diag_degree_mat[i][j]-weight_mat[i][j];
        }
    }
    return laplacian_mat;
}

double** jacobi(double** weight_mat, double** diag_degree_mat, int n){
    
}


double** wrapper(double** data_matrix, int n){

    //ddg side
    double **diag_degree_mat, **weight_mat;
    weight_mat = matrix_maker(n,n);
    if(weight_mat==NULL){
        free_matrix(data_matrix, n, n);
        return NULL;
    }
    diag_degree_mat = ddg_func(data_matrix, weight_mat, n);
    free_matrix(weight_mat, n, n);
    return diag_degree_mat;
 }

}



