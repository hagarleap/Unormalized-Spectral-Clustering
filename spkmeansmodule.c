# define PY_SSIZE_T_CLEAN
# include <Python.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "spkmeans.h"


/*Data structures*/
struct cord_node
{
    double value;
    struct cord_node *next;
};
struct vector_node
{
    struct vector_node *next;
    struct cord_node *cords;
};
struct dict_node
{
    struct dict_node *next;
    struct cord_node *centroid;
    struct cord_node *sum;
    int avg_divisor;
};

/*Receives sum field (cord_node type) from dict_centroid, and some other cord_node type.
Adds them and saves the result in the input sum field. */
void vector_addition(struct cord_node* closest_cluster_sum_vec,struct cord_node* vector,int vector_len){
    double cord_value_vector;
    int j;
    for(j=0; j<vector_len; j++)
      {
        cord_value_vector = vector->value;
        closest_cluster_sum_vec->value = closest_cluster_sum_vec->value + cord_value_vector;
        vector = vector->next; 
        closest_cluster_sum_vec = closest_cluster_sum_vec->next; 
      }

}

/*Recives 2 vectors. Returns the euclidian distance between them*/
double euclidian_distance(struct cord_node* vec1, struct cord_node* vec2, int vector_len){
    double sum=0;
    int j;
    double cord_value_vec1;
    double cord_value_vec2;
    for(j=0; j<vector_len; j++)
      {
        cord_value_vec1 = vec1->value;
        cord_value_vec2 = vec2->value;
        vec1 = vec1->next; 
        vec2 = vec2->next; 
        sum += pow((cord_value_vec1-cord_value_vec2), 2);
      }
    sum = pow(sum,0.5);
    return sum;
}

/*Creates the deltas linked list, makes it as long as amount of clusters there are.*/
struct cord_node* init_deltas(int K){
      struct cord_node *head_deltas_node, *curr_deltas_node, *prev_deltas_node;
      int i;
      head_deltas_node = malloc(sizeof(struct cord_node));
      if(head_deltas_node==NULL){return NULL;}
      
      curr_deltas_node = head_deltas_node;
      curr_deltas_node->next = NULL;

      prev_deltas_node = curr_deltas_node;
      for(i=0; i<K; i++)
      {
         curr_deltas_node->value = 1;
         curr_deltas_node->next = malloc(sizeof(struct cord_node));
         if((curr_deltas_node->next)==NULL){return NULL;}
         prev_deltas_node = curr_deltas_node;
         curr_deltas_node = curr_deltas_node->next; 
      }
      free(curr_deltas_node);
      prev_deltas_node->next=NULL; 
      return head_deltas_node;
}

/*Initializes the vector for the sum field. Filled with zeros, vector is as long as the input vectors.*/
struct cord_node* ZERO_vector(int vector_len){
      struct cord_node *head_zcord_node, *curr_zcord_node, *prev_zcord_node;
      int l;

      head_zcord_node = malloc(sizeof(struct cord_node));
      if(head_zcord_node==NULL){return NULL;}
      curr_zcord_node = head_zcord_node;
      curr_zcord_node->next = NULL;

      prev_zcord_node = curr_zcord_node;

      for(l=0; l<vector_len; l++)
      {
         curr_zcord_node->value = 0;
         curr_zcord_node->next = malloc(sizeof(struct cord_node));
         if((curr_zcord_node->next)==NULL){return NULL;}
         prev_zcord_node = curr_zcord_node;
         curr_zcord_node = curr_zcord_node->next; 
      }
      free(curr_zcord_node);
      prev_zcord_node->next=NULL; 
      return head_zcord_node;
      
}

/*Calculates new centroid using the avg_divisor and sum vector, saves it in the sum vector. Then calculates the 
delta difference between the new and old centroid. Finally it replaces the centroid field with the sum field,
and replaces the sum field with zeros. It also resets the avg_divisor field, and returns an integer 1 or 0.
0 means that all of the values in delta are lesser than epsilon. Otherwise, at least one is larger.*/
int update_centroid(struct dict_node *head_dict_centroid, struct cord_node *deltas,  int vector_len, double eps){
    int i =0;
    int max_delta_bigger_than_epsilon=0;
    struct cord_node *curr_sum_node;
    struct cord_node *curr_centroid_node;
    
    while(head_dict_centroid!=NULL){

        curr_sum_node = head_dict_centroid->sum;

        for(i=0;i<vector_len;i++){
            curr_sum_node->value = curr_sum_node->value /head_dict_centroid->avg_divisor;
            curr_sum_node = curr_sum_node->next;
        }
        head_dict_centroid->avg_divisor=0;
        deltas->value = euclidian_distance(head_dict_centroid->centroid, head_dict_centroid->sum ,vector_len);
        if((deltas->value) > eps){
            /*as long as ONE delta is bigger than epsilon, we want the while loop to keep going */
            max_delta_bigger_than_epsilon=1;
        }

        curr_centroid_node = head_dict_centroid->centroid;
        curr_sum_node = head_dict_centroid->sum;
        for (i=0; i<vector_len; i++){
            curr_centroid_node->value=curr_sum_node->value;
            curr_sum_node->value=0;
            curr_centroid_node= curr_centroid_node->next;
            curr_sum_node= curr_sum_node->next;
        }
        head_dict_centroid = head_dict_centroid->next;
    }
    return max_delta_bigger_than_epsilon;
}

/*Frees memory for an input cord node*/
void delete_cord_node(struct cord_node* cord_node){
    if (cord_node!=NULL){
        struct cord_node* next_cord;
        next_cord = NULL;
        if (cord_node != NULL)
        {
            next_cord= cord_node->next;
        }
        while(next_cord != NULL) {
            free(cord_node);
            cord_node = next_cord;
            next_cord = cord_node->next;
        }
        if (cord_node != NULL) { 
            free(cord_node);
        }  
    } 
}

/*Frees memory for an input vector node*/
void delete_vector_list( struct vector_node *vectors_list)
{
    if(vectors_list!=NULL){
        struct vector_node *curr_vector, *next_vector;
        curr_vector = vectors_list;
        next_vector = curr_vector->next;

        while (next_vector != NULL )
        {
            if (curr_vector->cords!=NULL){
                delete_cord_node(curr_vector->cords);
            }
            free(curr_vector);
            curr_vector = next_vector;
            next_vector = curr_vector->next;
        
        }
        if (curr_vector->cords!=NULL){
        delete_cord_node(curr_vector->cords);
        }
        free(curr_vector);
       }
}

/*Frees memory for an input dict node node*/
void delete_dict_list(struct dict_node *head_dict_centroid){
    if(head_dict_centroid!=NULL){
        struct dict_node *curr_dict, *next_dict;
        curr_dict = head_dict_centroid;
        next_dict = curr_dict->next;

        while (next_dict != NULL )
        {   
            if (curr_dict->centroid!=NULL){
            delete_cord_node(curr_dict->centroid);
            }
            if (curr_dict->sum!=NULL){
            delete_cord_node(curr_dict->sum);
            }
            free(curr_dict);
            curr_dict = next_dict;
            next_dict = curr_dict->next;
        
        }
        if (curr_dict->centroid!=NULL){
        delete_cord_node(curr_dict->centroid);
        }
        if (curr_dict->sum!=NULL){
        delete_cord_node(curr_dict->sum);
        }
        free(curr_dict);
    }
}

PyObject* kmeans(int K, int iter, int vector_len, int vectors_amt, double eps, PyObject *vectors, PyObject *centroids)
{
    int b, i,max_delta_bigger_than_epsilon=1, iter_count = 0, centroid_len_count=0,vectors_len_count=0,index=0;
    double  argmin, dist;

    PyObject* python_val;
    
    struct vector_node *vectors_list, *head_vec, *curr_vec, *prev_vec;
    struct cord_node *result_cord, *head_cord, *curr_cord, *head_cord2, *curr_cord2, *deltas, *delta_head_for_UC;
    struct dict_node *centroid_list_dict, *centroid_head_for_UC, *result, *head_dict_centroid;
    struct dict_node *curr_dict_centroid, *prev_dict_centroid, *closest_cluster;

    /*Simultaneously building from the text file the list of input vectors and the Dictionary
     that holds the centroids */
    head_dict_centroid = malloc(sizeof(struct dict_node));
    if(head_dict_centroid==NULL){
            return NULL;
         }
    curr_dict_centroid = head_dict_centroid;
    curr_dict_centroid->next = NULL;

    head_cord = malloc(sizeof(struct cord_node));
    if(head_cord==NULL){
            return NULL;

         }
    curr_cord = head_cord;
    curr_cord->next = NULL;

    head_vec = malloc(sizeof(struct vector_node));
    if(head_vec==NULL){
           return NULL;
         }
    curr_vec = head_vec;
    curr_vec->next = NULL;

    head_cord2 = malloc(sizeof(struct cord_node));
    if(head_cord2==NULL){
           return NULL;
         }
    curr_cord2 = head_cord2;
    curr_cord2->next = NULL;

    prev_dict_centroid = curr_dict_centroid;
    prev_vec = curr_vec;


    //creating vectors list as linked list
    for (i=0;i<(vectors_amt*vector_len);i++)
        {        
            if (vectors_len_count == (vector_len-1))
            {
                curr_cord->value =PyFloat_AsDouble(PyList_GetItem(vectors,i));
                curr_vec->cords = head_cord;
                curr_vec->next = malloc(sizeof(struct vector_node));
                if((curr_vec->next)==NULL){ return NULL;}
                prev_vec = curr_vec;
                curr_vec = curr_vec->next;
                curr_vec->next = NULL;
                curr_vec->cords = NULL;

                head_cord = malloc(sizeof(struct cord_node));
                if(head_cord==NULL){return NULL;}
                curr_cord = head_cord;
                curr_cord->next = NULL;
                vectors_len_count=0;
                continue;
                
            }

            curr_cord->value = PyFloat_AsDouble(PyList_GetItem(vectors,i));
            curr_cord->next = malloc(sizeof(struct cord_node));
            if((curr_cord->next)==NULL){
                return NULL;;
            }
            curr_cord = curr_cord->next;
            curr_cord->next = NULL;
            curr_cord->value = 0.0;
            vectors_len_count++;

        }


    //creating dict struct 
    for (i=0;i<(K*vector_len);i++)
        {        
            if (centroid_len_count == (vector_len-1))
            {

            curr_cord2->value = PyFloat_AsDouble(PyList_GetItem(centroids,i));
            curr_dict_centroid->centroid =head_cord2;
            curr_dict_centroid->sum = ZERO_vector(vector_len); 
            curr_dict_centroid->avg_divisor =0;
            curr_dict_centroid->next = malloc(sizeof(struct dict_node));
            if((curr_dict_centroid->next)==NULL){return NULL;}
            prev_dict_centroid = curr_dict_centroid;
            curr_dict_centroid = curr_dict_centroid->next;
            curr_dict_centroid->next = NULL;
            curr_dict_centroid->centroid = NULL;
            curr_dict_centroid->sum = NULL;

      
            head_cord2 = malloc(sizeof(struct cord_node));
            if(head_cord2==NULL){return NULL;}
            curr_cord2 = head_cord2;
            curr_cord2->next = NULL;   
            
            centroid_len_count=0;
            continue;
            
        }

        curr_cord2->value =PyFloat_AsDouble(PyList_GetItem(centroids,i));
        curr_cord2->next = malloc(sizeof(struct cord_node));
        if((curr_cord2->next)==NULL){return NULL;}
        curr_cord2 = curr_cord2->next;
        curr_cord2->next = NULL;
        curr_cord2->value = 0.0;
        centroid_len_count++;

    }


    free(curr_cord);
    free(curr_cord2);
    delete_vector_list(curr_vec);
    delete_dict_list(curr_dict_centroid);  

    prev_vec->next = NULL;
    prev_dict_centroid->next = NULL;   
    

    /*new head- for the lists of vector */
    vectors_list = head_vec;
    /*new head- for the list of centroid dict */
    centroid_list_dict = head_dict_centroid;
    deltas = init_deltas(K);

    /*Comparing and updating the centroids dictionary until we reach the maximum amount of iterations, or
    all deltas are less than epsilon.*/
     while ((iter_count<iter) && (max_delta_bigger_than_epsilon==1)){
        
        while(vectors_list != NULL){
            argmin = euclidian_distance(vectors_list->cords, centroid_list_dict->centroid, vector_len);
            closest_cluster = head_dict_centroid;
            centroid_list_dict= centroid_list_dict->next;
            while(centroid_list_dict != NULL){
                dist= euclidian_distance(vectors_list->cords, centroid_list_dict->centroid, vector_len);
                if (dist<argmin) {
                    argmin = dist;
                    /*need to verify it does not move */
                    closest_cluster = centroid_list_dict;
                }
                centroid_list_dict= centroid_list_dict->next;
            }
            closest_cluster->avg_divisor=  closest_cluster->avg_divisor + 1;
            vector_addition(closest_cluster->sum,vectors_list->cords, vector_len);
            vectors_list= vectors_list->next;
            centroid_list_dict = head_dict_centroid;
        }
        iter_count+=1;
        centroid_head_for_UC =head_dict_centroid;
        delta_head_for_UC = deltas;
        max_delta_bigger_than_epsilon = update_centroid(centroid_head_for_UC, delta_head_for_UC, vector_len, eps);
        vectors_list= head_vec;
     } 

    /*Get the copy of the head of the dict for final print */
    result = head_dict_centroid;
    
    python_val = PyList_New(K*vector_len); 
    while(result!=NULL){
        result_cord = result->centroid;
        for(b=0; b<vector_len; b++){
            PyList_SetItem(python_val, index, PyFloat_FromDouble(result_cord->value));
            result_cord = result_cord->next;
            index++;
        } 
        result = result->next;   
    }



    /*Free the following: head_vec, deltas, head_dict_centroid*/
    if(head_vec!=NULL){
        delete_vector_list(head_vec);
    }
    if(head_dict_centroid!=NULL){
        delete_dict_list(head_dict_centroid);
    }
    delete_cord_node(deltas);

    return python_val;
}


/////////////////////MODULE SECTION//////////////////////////

int counter = 0;
// wrapper function for kmeans
static PyObject* cKmeans(PyObject *self, PyObject *args) {
    int K, iter, vector_len, vectors_amt;
    double eps;
    PyObject *vectors, *centroids;

    if (!PyArg_ParseTuple(args, "iiiidOO", &K, &iter, &vector_len, &vectors_amt, &eps, &vectors, &centroids)) {
        return NULL;
    }

    return kmeans(K, iter, vector_len, vectors_amt, eps, vectors, centroids);
}

//wrapper function for C wrapper function
static PyObject* cWrapper(PyObject *self, PyObject *args) {}

// module's function table
static PyMethodDef kmeansMethods[] = {
    {"cKmeans",                   /* the Python method name that will be used */
      (PyCFunction) cKmeans, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters accepted for this function */
      PyDoc_STR("Returns centroids using a modified version of the k-means clustering program we created in HW 1, which uses new arguments: K, iter, vector_len, vectors_amt, eps, vectors, centroids.")}, /*  The docstring for the function */
    {NULL, NULL, 0, NULL} 
};

// modules definition
static struct PyModuleDef kmeans_Module = {
    PyModuleDef_HEAD_INIT,
    "kmeans_capi",     // name of module exposed to Python
    "Python wrapper for kmeans C extension library.", // module documentation
    -1,
    kmeansMethods
};

PyMODINIT_FUNC PyInit_kmeans_capi(void) {
    PyObject *m;
    m = PyModule_Create(&kmeans_Module);
    if (!m) {
        return NULL;
    }
    return m;
}