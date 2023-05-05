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
//
template<class T>
void squareMatInverse2(int rows, int cols, T **arrIn, T **arrOut);



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


#include <stdio.h>  /*This form is used for system header files.*/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include <cblas.h> // use -lblas 
#include <lapack.h> // use -llapack


//***********************************************************************************
// Interface for lapack functions
// sgetrf_() & sgetri_(): combination for inverse
//***********************************************************************************
template<class T>
void squareMatInverse2(int rows, int cols, T **arrIn, T **arrOut){
    //printf("\n------------------------------------");
    //printf("\n         Testing LAPACK           \n");

    // assuming that the given arrIn is a 2D array. 
    // transform it to 1D- array for lapack functions
    T *AA;
    AA = (T*)malloc((rows*cols) *sizeof(T));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            AA[i * cols + j] = arrIn[i][j];
            //printf("Local: %f, In: %f ", AA[i * cols + j],arrIn[i][j]);
        } 
        //printf("\n");
    }
    int *IPIV = (int*)malloc((rows) *sizeof(int));
    int INFO;
    //!     Factorize A
    #if PRECISION_MODE_FEM == 2 // SINGLE PRECISION
        sgetrf_(&rows,&cols,AA,&rows,IPIV,&INFO); //A = PLU
    #endif
    #if PRECISION_MODE_FEM == 1
        dgetrf_(&rows,&cols,AA,&rows,IPIV,&INFO); //A = PLU
    #endif

    if (INFO==0){
    //    !       Compute inverse of A
        int LWORK = 64*cols;
        T *WORK = (T*)malloc((LWORK) *sizeof(T));
        #if PRECISION_MODE_FEM == 2 // SINGLE PRECISION
            sgetri_(&rows, AA, &cols, IPIV, WORK, &LWORK, &INFO);
        #endif
        #if PRECISION_MODE_FEM == 1
            dgetri_(&rows, AA, &cols, IPIV, WORK, &LWORK, &INFO);
        #endif

        free(WORK);
    }
    
    //printf("\n INVERSE MATRIX ------------------\n\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++){
            //printf("%f ", AA[i * cols + j]);
            arrOut[i][j] = AA[i * cols + j];
        }
        //printf("\n");
    }
    //printf("\n----------------------------------\n\n");

    //for (int i = 0; i < rows; i++) {
    //    printf("%d ", IPIV[i]);
    //}
    //printf("\n%d,----------------------------------\n\n",INFO);

    free(IPIV);
    free(AA);
    //printf("         Inverse 2D OK..          \n");
    //printf("------------------------------------\n");
}