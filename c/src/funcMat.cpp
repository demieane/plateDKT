#ifndef MAIN_FEM_HEADER_FILE
    #define MAIN_FEM_HEADER_FILE

    #include "../include/mainFem.h"


#endif

//
template<class T>
void allocate1Darray(int rows, T **arrIn);
//
template<class T>
void allocate2Darray(int rows, int cols, T ***arrIn);



// from funcBLAS.c
/*=========================================================================================*/
/* Definition of the functions BELOW */
/*=========================================================================================*/
// Allocate 1-D array based on double pointer type
template<class T>
void allocate1Darray(int rows, T **arrIn){

    int i;//, j;

    // allocate local 2D array and pass the pointer to 
    // the pointer to double pointer, otherwise the main
    // can't access the memory.
    T *arrTemp = (T*)malloc(rows * sizeof(T));
    if(arrIn == NULL){
        printf("Memory allocation failed. allocate2Darray()");
        return;
    }

    //printf("\nInside allocate1Darray()..\n");
    // Note that arr[i][j] is same as *(*(arr+i)+j)
    for (i = 0; i < rows; i++){
        arrTemp[i] = 0.0;
        //printf("\n%f,",arrTemp[i]);
    }

    *arrIn = arrTemp; // this will do?
}


// Allocate 2-D array based on double pointer type
template<class T>
void allocate2Darray(int rows, int cols, T ***arrIn){

    int i, j;

    // allocate local 2D array and pass the pointer to 
    // the pointer to double pointer, otherwise the main
    // can't access the memory.
    T **arrTemp = (T**)malloc(rows * sizeof(T*));
    if(arrIn == NULL){
        printf("Memory allocation failed. allocate2Darray()");
        return;
    }

    for (i = 0; i < rows; i++){
        arrTemp[i] = (T*)malloc(cols * sizeof(T));
    }

    //printf("\nInside allocate2Darray()..\n");
    // Note that arr[i][j] is same as *(*(arr+i)+j)
    for (i = 0; i < rows; i++){
        for (j = 0; j < cols; j++){
            arrTemp[i][j] = 0.0;
            //printf("%f,",arrTemp[i][j]);
        } 
        //printf("\n");
    }

    *arrIn = arrTemp; // this will do?
}
