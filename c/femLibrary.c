#ifndef MAIN_FEM_HEADER_FILE
#define MAIN_FEM_HEADER_FILE

#include "mainFem.h"

#endif

#include <stdio.h>  /*This form is used for system header files.*/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cblas.h> // use -lblas 
#include <lapack.h> // use -llapack


//***********************************************************************************
// Interface for lapack functions
// sgetrf_() & sgetri_(): combination for inverse
//***********************************************************************************

void squareMatInverse2(int rows, int cols, float **arrIn, float **arrOut){
    printf("\n------------------------------------");
    printf("\n         Testing LAPACK           \n");

    // assuming that the given arrIn is a 2D array. 
    // transform it to 1D- array for lapack functions
    float *AA;
    AA = (float*)malloc((rows*cols) *sizeof(float));

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
    sgetrf_(&rows,&cols,AA,&rows,IPIV,&INFO); //A = PLU

    if (INFO==0){
    //    !       Compute inverse of A
        int LWORK = 64*cols;
        float *WORK = (float*)malloc((LWORK) *sizeof(float));
        sgetri_(&rows, AA, &cols, IPIV, WORK, &LWORK, &INFO);
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

    free(AA);
    printf("         Inverse 2D OK..          \n");
    printf("------------------------------------\n");
}

// Allocate 2-D array based on double pointer type
void allocate2Darray(int rows, int cols, float ***arrIn){

    int i, j;

    // allocate local 2D array and pass the pointer to 
    // the pointer to double pointer, otherwise the main
    // can't access the memory.
    float **arrTemp = (float**)malloc(rows * sizeof(float*));
    if(arrIn == NULL){
        printf("Memory allocation failed. allocate2Darray()");
        return;
    }

    for (i = 0; i < rows; i++){
        arrTemp[i] = (float*)malloc(cols * sizeof(float));
    }

    printf("\nInside allocate2Darray()..\n");
    // Note that arr[i][j] is same as *(*(arr+i)+j)
    for (i = 0; i < rows; i++){
        for (j = 0; j < cols; j++){
            arrTemp[i][j] = 99.0;
            printf("%f,",arrTemp[i][j]);
        } 
        printf("\n");
    }

    *arrIn = arrTemp; // this will do?
}