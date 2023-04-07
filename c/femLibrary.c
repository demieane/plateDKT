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
