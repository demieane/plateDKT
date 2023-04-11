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

void matMatMultiplication2(int rowsA, int colsA, int colsB, float **arrA, float **arrB, float **arrOut){
    //   SGEMM - perform one of the matrix-matrix operations   
    //   C := alpha*op( A )*op( B ) + beta*C,
    //   M number of rows (A)
    //   N number of columns for (B)

    float alpha = 1.0, beta = 0.0;
    int rowsB = colsA;
    //
    int rowsC = rowsA;
    int colsC = colsB;

    printf("\n------------------------------------");
    printf("\nMatMat Mult: [M x P] [P x N]    \n");
    // assuming that the given arrA, arrB are 2D arrays. 
    // transform it to 1D- array for lapack functions
    float *AA;
    AA = (float*)malloc((rowsA*colsA) *sizeof(float));

    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsA; j++){
            AA[i * colsA + j] = arrA[i][j];
            //printf("Local: %f, In: %f ", AA[i * cols + j],arrIn[i][j]);
        } 
        //printf("\n");
    }

    float *BB;
    BB = (float*)malloc((rowsB*colsB) *sizeof(float));
    for (int i = 0; i < rowsB; i++) {
        for (int j = 0; j < colsB; j++){
            BB[i * colsB + j] = arrB[i][j];
            //printf("Local: %f, In: %f ", AA[i * cols + j],arrIn[i][j]);
        } 
        //printf("\n");
    }

    float *CC;
    CC = (float*)malloc((rowsC*colsC) *sizeof(float));

    int LDA = colsA; // increment in the array (due to row major order)
    int LDB = colsB;
    int LDC = colsC;
    /*
    The leading dimension for a two-dimensional array is an 
    increment that is used to find the starting 
    point for the matrix elements. (length of the leading dimension - ROW)
    */
    cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, rowsC, colsC, rowsA,
     alpha, AA, LDA, BB, LDB, beta, CC, LDC);
    //cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 10, 9, 10, alpha,
    //        AA, 10, BB, 9, beta, CC, 9);
     
    
    //printf("\n INVERSE MATRIX ------------------\n\n");
    for (int i = 0; i < rowsC; i++) {
        for (int j = 0; j < colsC; j++){
            arrOut[i][j] = CC[i * colsC+ j];
            //printf("%f, ",arrOut[i][j]);
        }
        //printf("\n");
    }

    free(AA);
    free(BB);
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
            arrTemp[i][j] = 0.0;
            //printf("%f,",arrTemp[i][j]);
        } 
        //printf("\n");
    }

    *arrIn = arrTemp; // this will do?
}

// Allocate 1-D array based on double pointer type
void allocate1Darray(int rows, float **arrIn){

    int i, j;

    // allocate local 2D array and pass the pointer to 
    // the pointer to double pointer, otherwise the main
    // can't access the memory.
    float *arrTemp = (float*)malloc(rows * sizeof(float));
    if(arrIn == NULL){
        printf("Memory allocation failed. allocate2Darray()");
        return;
    }

    printf("\nInside allocate1Darray()..\n");
    // Note that arr[i][j] is same as *(*(arr+i)+j)
    for (i = 0; i < rows; i++){
        arrTemp[i] = 0.0;
        //printf("\n%f,",arrTemp[i]);
    }

    *arrIn = arrTemp; // this will do?
}