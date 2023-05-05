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
//
template<class T>
void matMatMultiplication2(int optionCalc, int rowsA, int colsA, int colsB, T alpha, T beta, T **arrA, T **arrB, T **arrOut);
//
template<class T>
void matSum2(T alpha, T beta, int rows, int cols, T **arrA, T **arrB, T **arrOut);
//
template<class T>
void deallocate2Darray(int rows, T ***arrIn);
//
template<class T>
void linearSystemSolve(int rowsA, int colsA, T **arrA, T **arrB, T **Usol);
    //

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
void deallocate2Darray(int rows, T **arrIn){

    for (int i = 0; i < rows; i++){
        free(arrIn[i]);
    }
    free(arrIn);
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


template<class T>
void matMatMultiplication2(int optionCalc, int rowsA, int colsA, int colsB, T alpha, T beta, T **arrA, T **arrB, T **arrOut){
    //   SGEMM - perform one of the matrix-matrix operations   
    //   C := alpha*op( A )*op( B ) + beta*C,
    //   M number of rows (A)
    //   N number of columns for (B)

    /*
    SUBROUTINE SGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
    */

    // optionCalc == 1: no transpose
    // optionCalc == 2: A transpose

    //float alpha = 1.0, beta = 0.0;
    int rowsB, rowsC, colsC;
    if (optionCalc == 1){
        rowsB = colsA;
        rowsC = rowsA;
        colsC = colsB;
        //printf("option (1)");
    }
    if (optionCalc == 2){
        rowsB = rowsA;
        rowsC = colsA; // transA
        colsC = colsB;
        //printf("option (2)");
    }
    //printf("%d, %d, %d, %d, %d, %d,", rowsA, colsA, rowsB, colsB, rowsC, colsC);

    //printf("\n------------------------------------");
    //printf("\n  MatMat Mult: [M x P] [P x N]      ");
    //printf("\n------------------------------------\n");
    // assuming that the given arrA, arrB are 2D arrays. 
    // transform it to 1D- array for lapack functions
    T *AA;
    AA = (T*)malloc((rowsA*colsA) *sizeof(T));

    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsA; j++){
            AA[i * colsA + j] = arrA[i][j];
            //printf("Local: %f, In: %f ", AA[i * colsA + j],arrA[i][j]);
        } 
        //printf("\n");
    }
    //printf("\n\n");
    //printf("Allocated A.. OK!\n");

    T *BB;
    BB = (T*)malloc((rowsB*colsB) *sizeof(T));
    for (int i = 0; i < rowsB; i++) {
        for (int j = 0; j < colsB; j++){
            BB[i * colsB + j] = arrB[i][j];
            //BB[i * colsB + j] = arrB[i][j]/1000;
            //printf("Local: %f, In: %f ", BB[i * colsB + j],arrB[i][j]);
        } 
        //printf("\n");
    }
    //printf("Allocated B.. OK!\n");

    T *CC;
    CC = (T*)malloc((rowsC*colsC) *sizeof(T));
    //printf("Allocated C.. OK!\n");

    for (int i = 0; i < rowsC; i++) {
        for (int j = 0; j < colsC; j++){
            arrOut[i][j] = CC[i * colsC+ j];
            //printf("%f, ",arrOut[i][j]);
        }
        //printf("\n");
    }

    if (optionCalc == 1){
        //#if DEBUG_ON   
        //printf("\n MAT MUL option (1)\n\n");
        //#endif
        int LDA = colsA; // increment in the array (due to row major order)
        int LDB = colsB;
        int LDC = colsC;

        //printf("%d, %d, %d,", LDA, LDB, LDC);

        int M,N,K;
        M = rowsA; //rows of op(A)
        N = colsB; //cols of op(B)
        K = rowsB; //rows of op(B)
        #if PRECISION_MODE_FEM == 2
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M, N, K,
            alpha, &(AA[0]), LDA, &(BB[0]), LDB, beta, &(CC[0]), LDC);
            //cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 10, 9, 10, alpha,
            //        AA, 10, BB, 9, beta, CC, 9);
        #endif
        #if PRECISION_MODE_FEM == 1
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M, N, K,
            alpha, &(AA[0]), LDA, &(BB[0]), LDB, beta, &(CC[0]), LDC);
            //cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 10, 9, 10, alpha,
            //        AA, 10, BB, 9, beta, CC, 9);
        #endif

    }
    if (optionCalc == 2){
        //#if DEBUG_ON   
        //printf("\n MAT MUL option (2) \n\n");
        //#endif
        int LDA = colsA; // increment in the array (due to row major order)
        int LDB = colsB;
        int LDC = colsC;
        /*
        The leading dimension for a two-dimensional array is an 
        increment that is used to find the starting 
        point for the matrix elements. (length of the leading dimension - ROW)
        */
        int M, N, K;
        M = colsA;
        N = colsB;
        K = rowsA;
        #if PRECISION_MODE_FEM == 2
            cblas_sgemm(CblasRowMajor,CblasTrans,CblasNoTrans, M, N, K,
            alpha, &(AA[0]), LDA, &(BB[0]), LDB, beta, &(CC[0]), LDC);
            //cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 10, 9, 10, alpha,
            //        AA, 10, BB, 9, beta, CC, 9);
        #endif
        #if PRECISION_MODE_FEM == 1
            cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans, M, N, K,
            alpha, &(AA[0]), LDA, &(BB[0]), LDB, beta, &(CC[0]), LDC);
            //cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 10, 9, 10, alpha,
            //        AA, 10, BB, 9, beta, CC, 9);
        #endif

    }

   
    //printf("Result in C.. OK!\n");
    for (int i = 0; i < rowsC; i++) {
        for (int j = 0; j < colsC; j++){
            arrOut[i][j] = CC[i * colsC+ j];
            //printf("%f, ",arrOut[i][j]);
        }
        //printf("\n");
    }


    free(AA);
    free(BB);
    free(CC);
//#if DEBUG_ON   
//    printf(" EXITING matMatMultiplication2...\n");
//#endif
}

template<class T>
void matSum2(T alpha, T beta, int rows, int cols, T **arrA, T **arrB, T **arrOut){

    for (int i=0;i<rows;i++){
        for (int j=0;j<cols;j++){
            arrOut[i][j]= arrOut[i][j] + alpha*arrA[i][j]  + beta*arrB[i][j];
        }
    }

/*
    for (int i=0;i<rows;i++){
        for (int j=0;j<cols;j++){
            printf("%f, ",arrOut[i][j]);
        }
        printf("\n");
    }
*/

}

template<class T>
void linearSystemSolve(int rowsA, int colsA, T **arrA, T **arrB, T **Usol){

    // assuming that the given arrA, arrB are 2D arrays. 
    // transform it to 1D- array for lapack functions
    T *AA;
    AA = (T*)malloc((rowsA*colsA) *sizeof(T));

    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsA; j++){
            AA[i * colsA + j] = arrA[i][j];
            //printf("Local: %f, In: %f ", AA[i * colsA + j],arrA[i][j]);
        } 
        //printf("\n");
    }
    //printf("\n\n");
    //printf("Allocated A.. OK!\n");

    T *BB;
    BB = (T*)malloc((rowsA*1) *sizeof(T));
    for (int i = 0; i < rowsA; i++) {
        BB[i] = arrB[i][0];
    }
    //printf("Allocated B.. OK!\n");

    int nrhs = 1; //number of rhs vectors B 
    int LDA = colsA;
    int LDB = rowsA;
    int info;

    int *IPIV = (int*)malloc((rowsA) *sizeof(int));

    // https://www.intel.com/content/www/us/en/docs/onemkl/code-samples-lapack/2022-1/sgesv-example-c.html
    #if PRECISION_MODE_FEM == 2
        sgesv_(&rowsA, &nrhs, AA, &LDA, IPIV , BB, &LDB, &info); //INTEL DOCS
    #endif
    #if PRECISION_MODE_FEM == 1
        dgesv_(&rowsA, &nrhs, AA, &LDA, IPIV , BB, &LDB, &info); //INTEL DOCS
    #endif

    if (info >= 0){
        printf("\n    Solution successfull");
    } 

    for (int i = 0; i < rowsA; i++) {
        Usol[i][0]=BB[i];
    }
    printf("\n    Transfered solution from B to Usol.. OK!\n");

    free(IPIV);
    free(BB);
    free(AA);


}