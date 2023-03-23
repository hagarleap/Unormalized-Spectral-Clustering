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
    if(mat==NULL){
        return;
    }
    for(int i=0; i<dim1; i++){
        free(mat[i]);
    }
    free(mat);
}

// Calculate and output the Weighted Adjacency Matrix
double** wam_func(double** data_matrix, int n, int dim2){
    double** weight_mat;
    int i, j;
    if(data_matrix==NULL){
        return NULL;
    }
    //create matrix
    weight_mat = matrix_maker(n,n);
    if(weight_mat==NULL){
        return NULL;
    }
    // calculate weights
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i!=j){
                weight_mat[i][j] = exp(-(Squared_Euclidean_distance(data_matrix[i], data_matrix[j], dim2))/2);
            }
    
        }
    }
    return weight_mat;
}

// Calculate and output the Diagonal Degree Matrix
 double** ddg_func( double** weight_mat, int n){
    double **diag_degree_mat;
    int i, j;

    if(weight_mat==NULL){
        return NULL;
    }

    //create matrix
    diag_degree_mat = matrix_maker(n,n);
    if(diag_degree_mat==NULL){
        free_matrix(weight_mat, n);
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

    if(weight_mat==NULL){
        return NULL;
    }
    if(diag_degree_mat==NULL){
        return NULL;
    }
    //create output matrix 
    laplacian_mat = matrix_maker(n,n);
    if(laplacian_mat==NULL){
        free_matrix(diag_degree_mat, n);
        free_matrix(weight_mat, n);
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
    int i, j,max = 0;
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
double** calculate_A_tag_matrix(double** A, int n, double c, double s, int index_i, int index_j){
    double *A_tag_coli, *A_tag_colj;
    int r;
    

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

//returns V*P^i
double** calculate_new_V_matrix(double** V_matrix, int n, double c, double s, int index_i, int index_j){
    // double **rotation_mat;
    // int i;

    // rotation_mat = matrix_maker(n,n);
    // if(rotation_mat==NULL){
    //     return NULL;
    // }

    // for(i=0;i<n;i++){
    //         rotation_mat[i][i] = 1;
    // }
    // rotation_mat[index_i][index_i] = c;
    // rotation_mat[index_j][index_j] = c;
    // rotation_mat[index_i][index_j] = s;
    // rotation_mat[index_j][index_i] = -s;

    double Vri, Vrj;
    int r;
    for(r=0;r<n;r++){
        Vri=V_matrix[r][index_i];
        Vrj=V_matrix[r][index_j];

        V_matrix[r][index_i] = Vri*c - s*Vrj;
        V_matrix[r][index_j] = Vrj*c + s*Vri;
    }
    return V_matrix;
}

double off_func(double **A, int n){
    int i, j;
    double res;
    for(i=0; i<n; i++){
         for (j=0; j<n; j++){
            if(i!=j){
               res = res + pow(A[i][j], 2); 
            } 
         }
    }
    return res;
}

 /*recieve symmetric matrix A and return the e eigenvalues and eigenvectors of a real A
 when A = laplacian_mat 
*/
double** jacobi_func(double** A, int n){
    int *pivot;
    int r,sign_theta,index_i,index_j, cnt_num_rotation=0,i,j;
    double theta, s ,c ,t ,convergence, off_A, EPSILON=0.00001;
    double **V_matrix , **res_matrix;
    
    if(A==NULL){
        return NULL;
    }
    
    //create V_matrix 
    V_matrix = matrix_maker(n,n);
    if(V_matrix==NULL){
        return NULL;
    }
    //initialize V matrix as unit matrix
    for(i=0;i<n;i++){
            V_matrix[i][i] = 1;
    }
    
    //create res_matrix 
    res_matrix = matrix_maker(n+1,n);
    if(res_matrix==NULL){
        free_matrix(V_matrix, n);
        return NULL;
    }

    
    //We will be using EPSILON = 1.0 × 10−5 OR maximum number of rotations = 100 
    while((cnt_num_rotation < 100) && (convergence>EPSILON)){
        off_A = off_func(A, n);
        pivot = pivot_func(A,n);
        index_i = pivot[0];
        index_j = pivot[1];
        theta= (A[index_j][index_j]-A[index_i][index_i])/(2*A[index_i][index_j]);
        //check this calc
        sign_theta = (theta >= 0) - (theta < 0);
        t = sign_theta / (abs(theta)+pow((theta*theta+1), 0.5));
        c = 1 / pow((t*t+1), 0.5);
        s = t*c;
        A = calculate_A_tag_matrix(A , n, c, s, index_i, index_j);
        if(A==NULL){
            return NULL;
        } 
        V_matrix = calculate_new_V_matrix(V_matrix, n, c, s, index_i, index_j);
        if(V_matrix==NULL){
            return NULL;
         }
        convergence = off_A - off_func(A, n);
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

    free_matrix(V_matrix, n);
    return res_matrix;
}

//get the goal and return the matrix
 double** wrapper(double** data_matrix, int n,int dim2,char * goal){
    double **wam_mat, **ddg_mat, **gl_mat, **jacobi_mat;
    int is_jacobi=0;
    if (!strcmp(goal,"wam"))
    {
        wam_mat = wam_func( data_matrix,n, dim2);
        wam_mat= non_neg_zero(wam_mat,n,is_jacobi);
        return wam_mat;
    }
    else if (!strcmp(goal,"ddg"))
    {
        wam_mat = wam_func(data_matrix,n, dim2);
        ddg_mat = ddg_func(wam_mat, n);
        free_matrix(wam_mat,n);
        ddg_mat= non_neg_zero(ddg_mat,n,is_jacobi);
        return ddg_mat;
    }
    else if (!strcmp(goal,"gl"))
    {
        wam_mat = wam_func(data_matrix,n, dim2);
        ddg_mat = ddg_func(wam_mat, n);
        gl_mat = gl_func(wam_mat, ddg_mat, n);
        free_matrix(wam_mat,n);
        free_matrix(ddg_mat,n);
        gl_mat= non_neg_zero(gl_mat,n,is_jacobi);
        return gl_mat;
    }
    else if (!strcmp(goal,"jacobi"))
    {
        wam_mat = wam_func( data_matrix, n, dim2);
        ddg_mat = ddg_func(   wam_mat, n);
        gl_mat = gl_func(wam_mat, ddg_mat, n);
        jacobi_mat = jacobi_func(gl_mat, n);
        free_matrix(wam_mat,n);
        free_matrix(ddg_mat,n);
        free_matrix(gl_mat,n);
        is_jacobi=1;
        jacobi_mat= non_neg_zero(jacobi_mat,n,is_jacobi);
        return jacobi_mat;
    }

    
 }


void print_matrix(double **goal_mat, int n,int is_jacobi){
    int i,j;
    if (is_jacobi==0){
        for(i=0; i<n; i++){
            for (j=0; j<n-1; j++){
                printf("%.4f,",  goal_mat[i][j]);  
            }
        printf("%.4f\n", goal_mat[i][j]);
        }
    }
    else{
       { for(i=0; i<n+1; i++)
            for (j=0; j<n-1; j++){
                printf("%.4f,",  goal_mat[i][j]);  
            }
        printf("%.4f\n", goal_mat[i][j]); 
        }
    }
}


double** non_neg_zero(double** mat,int n,int is_jacobi){
    int i,j;
    if(mat==NULL){
        return NULL;
    }

    if (is_jacobi==0){
        for(i=0; i<n; i++){
            for (j=0; j<n; j++){
                if(mat[i][j]==0){
                   mat[i][j]= 0;
                }  
            }
        }
    }
    else{
       { for(i=0; i<n+1; i++)
            if(mat[i][j]==0){
                   mat[i][j]= 0;
            }  
        }
    }
    return mat;
}



int main(int argc, char *argv[])
{
    double **goal_mat, **data_matrix, num;
    char *goal ,*file_name, ch, c;
    int is_jacobi, i=0 ,j=0 ,n=0 ,dim2=1;
    FILE *pf;

    goal = argv[1];
    file_name =  argv[2];

    //building data_matrix
    
    
    pf = fopen (file_name, "r");
    if (pf == NULL){
        printf("An Error Has Occurred");
        return 0;
    }

    /* figure out n */
    while((ch=fgetc(pf))!=EOF) 
    {
        if (ch=='\n')
        {
            n++;
        }
    }
    rewind(pf);
    
    /* figure out dim2 */
    while((ch=fgetc(pf))!='\n') 
    {
        if (ch==',')
        {
            dim2++;
        }
    }
    rewind(pf);


    data_matrix = matrix_maker(n,dim2);
    if(data_matrix==NULL){
        printf("An Error Has Occurred");
        return 0;
    }

    while (scanf("%lf%c", &num, &c) == 2)
    {        
        if (c == '\n')
        {
            data_matrix[i][j]=num;
            j=0;
            i++;
        }
        else{
            j++;
            data_matrix[i][j]=num;
        }
    }
    fclose (pf); 
    

    //print output
    goal_mat = wrapper( data_matrix ,n , dim2, goal);
    if(goal_mat==NULL){
        printf("An Error Has Occurred");
        return 0;
    }

    if (!strcmp(goal,"jacobi")){
        is_jacobi = 1;
        print_matrix(goal_mat,n, is_jacobi);
    }
    else{
        is_jacobi = 0;
        print_matrix(goal_mat,n, is_jacobi);
    }

    free_matrix(goal_mat,n);
    return 0;

}