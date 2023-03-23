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
double** wam_func(double** data_matrix, int n){
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
 double** ddg_func(double** data_matrix, double** weight_mat, int n){
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
double** gl_func(double** weight_mat, double** diag_degree_mat, int n){
    double **laplacian_mat;
    int i, j;

    //create output matrix 
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

/* get symmetric matrix A and return The i,j  for Aij is the off-diagonal element with the largest absolute value.*/
int* pivot_func(double** A, int n){
    double index[2];
    int i ,max = 0;
    for(i=0; i<n; i++){
         for (j=i+1; j<n; j++){
            if(abs(A[i][j])>max){
                max = A[i][j];
                index[0] = i;
                index[1] = j;
            }
        }
    }
    return index;
}

/* get matrix A and pivot indexes and return the A'*/
double** calculate_A_tag_matrix(double** A, int n, double c, double s){
    double *A_tag_coli, *A_tag_colj;
    

    A_tag_coli = calloc(n, sizeof(double));
        if(A_tag_coli==NULL){
            return NULL;
        }   

    A_tag_colj= calloc(n, sizeof(double));
        if(A_tag_colj==NULL){
            free(A_tag_coli);
            return NULL;
        }    

   
    for(r=0; r<n; r++){
        if(r!=index_i && r!=index_j){
            A_tag_coli[r] = c*(A[r][index_i]) - s*(A[r][index_j]);
            A_tag_colj[r] = c*(A[r][index_j]) + s*(A[r][index_i]);
        } 
             
    }
    A_tag_coli[index_i] = c*c*(A[index_i][index_i]) + s*s*(A[index_j][index_j]) - 2*s*c*(A[index_i][index_j]);
    A_tag_colj[index_j] = s*s*(A[index_i][index_i]) + c*c*(A[index_j][index_j]) + 2*s*c*(A[index_i][index_j]);
    A_tag_coli[index_j] = 0;
    A_tag_colj[index_i] = 0;

    
    for(r=0; r<n; r++){
        //replace row and col index_i
        A[index_i][r] = A_tag_coli[r];
        A[r][index_i] = A_tag_coli[r];
        //replace row and col index_j
        A[index_j][r] = A_tag_colj[r];
        A[r][index_j] = A_tag_colj[r];
    }

   free(A_tag_coli);
   free(A_tag_colj);

    return A;
}


calculate_new_V_matrix
off_func

 /* get symmetric matrix A and return the e eigenvalues and eigenvectors of a real A
 when A = laplacian_mat 
*/
double** jacobi_func(double** A, int n){
    int pivot[2];
    int r,sign_theta,index_i,index_j, cnt_num_rotation=0,i,j;
    double theta, s ,c ,t ,convergence, off_A;
    double **V_matrix , **res_matrix;
    
    

    //create V_matrix 
    V_matrix = matrix_maker(n,n);
    if(V_matrix==NULL){
        return NULL;
    }

    //create V_matrix 
    res_matrix = matrix_maker(n+1,n);
    if(res_matrix==NULL){
        free_matrix(V_matrix);
        return NULL;
    }

    
    //We will be using EPSILON = 1.0 × 10−5 OR maximum number of rotations = 100 
    while((cnt_num_rotation < 100) && (convergence>EPSILON)){
        off_A = off_func(A);
        pivot = pivot_func(A,n);
        index_i = pivot[0];
        index_j = pivot[1];
        theta= (A[index_j][index_j]-A[index_i][index_i])/(2*A[index_i][index_j]);
        //check this calc
        sign_theta = (theta >= 0) - (theta < 0);
        t = sign_theta / (abs(theta)+pow((theta*theta+1), 0.5));
        c = 1 / pow((t*t+1), 0.5);
        s = t*c;
        A = calculate_A_tag_matrix(A , n, c, s);
        V_matrix = calculate_new_V_matrix(V_matrix, n, c, s );
        convergence = off_A -off_func(A);
        cnt_num_rotation++;

    }
    
    //first row will be the eigenvalues, and than each col is the eigenvectors 
    for(i=0; i<n; i++){
        res_matrix[0][i] = A[i][i];
    }

    for(i=1; i<n+1; i++){
        for (j=0; j<n; j++){
            res_matrix[i][j] =V_matrix[i-1][j];
        }
    }

    free_matrix(V_matrix);
    return res_matrix;
}


double** wrapper(double** data_matrix, int n){

    //ddg side
   /* double **diag_degree_mat, **weight_mat;
    weight_mat = matrix_maker(n,n);
    if(weight_mat==NULL){
        free_matrix(data_matrix, n, n);
        return NULL;
    }
    diag_degree_mat = ddg_func(data_matrix, weight_mat, n);
    free_matrix(weight_mat, n, n);
    return diag_degree_mat;*/
 }





int main(int argc, char *argv[])
{
}