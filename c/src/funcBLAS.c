#ifndef MAIN_FEM_HEADER_FILE
#define MAIN_FEM_HEADER_FILE

#include "../include/mainFem.h"

#endif


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

void matMatMultiplication2(int optionCalc, int rowsA, int colsA, int colsB, float alpha, float beta, float **arrA, float **arrB, float **arrOut){
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
    float *AA;
    AA = (float*)malloc((rowsA*colsA) *sizeof(float));

    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsA; j++){
            AA[i * colsA + j] = arrA[i][j];
            //printf("Local: %f, In: %f ", AA[i * colsA + j],arrA[i][j]);
        } 
        //printf("\n");
    }
    //printf("\n\n");
    //printf("Allocated A.. OK!\n");

    float *BB;
    BB = (float*)malloc((rowsB*colsB) *sizeof(float));
    for (int i = 0; i < rowsB; i++) {
        for (int j = 0; j < colsB; j++){
            BB[i * colsB + j] = arrB[i][j];
            //BB[i * colsB + j] = arrB[i][j]/1000;
            //printf("Local: %f, In: %f ", BB[i * colsB + j],arrB[i][j]);
        } 
        //printf("\n");
    }
    //printf("Allocated B.. OK!\n");

    float *CC;
    CC = (float*)malloc((rowsC*colsC) *sizeof(float));
    //printf("Allocated C.. OK!\n");

    for (int i = 0; i < rowsC; i++) {
        for (int j = 0; j < colsC; j++){
            arrOut[i][j] = CC[i * colsC+ j];
            //printf("%f, ",arrOut[i][j]);
        }
        //printf("\n");
    }

    if (optionCalc == 1){
        #if DEBUG_ON   
        printf("\n MAT MUL option (1)\n\n");
        #endif
        int LDA = colsA; // increment in the array (due to row major order)
        int LDB = colsB;
        int LDC = colsC;

        //printf("%d, %d, %d,", LDA, LDB, LDC);

        int M,N,K;
        M = rowsA; //rows of op(A)
        N = colsB; //cols of op(B)
        K = rowsB; //rows of op(B)
        cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, M, N, K,
        alpha, &(AA[0]), LDA, &(BB[0]), LDB, beta, &(CC[0]), LDC);
        //cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 10, 9, 10, alpha,
        //        AA, 10, BB, 9, beta, CC, 9);
    }
    if (optionCalc == 2){
        #if DEBUG_ON   
        printf("\n MAT MUL option (2) \n\n");
        #endif
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
        cblas_sgemm(CblasRowMajor,CblasTrans,CblasNoTrans, M, N, K,
        alpha, &(AA[0]), LDA, &(BB[0]), LDB, beta, &(CC[0]), LDC);
        //cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, 10, 9, 10, alpha,
        //        AA, 10, BB, 9, beta, CC, 9);
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
#if DEBUG_ON   
    printf(" EXITING matMatMultiplication2...\n");
#endif
}

void matSum2(float alpha, float beta, int rows, int cols, float **arrA, float **arrB, float **arrOut){

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

// Allocate 1-D array based on double pointer type
void allocate1Darray(int rows, float **arrIn){

    int i;//, j;

    // allocate local 2D array and pass the pointer to 
    // the pointer to double pointer, otherwise the main
    // can't access the memory.
    float *arrTemp = (float*)malloc(rows * sizeof(float));
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

void linearSystemSolve(int rowsA, int colsA, float **arrA, float **arrB, float **Usol){

    // assuming that the given arrA, arrB are 2D arrays. 
    // transform it to 1D- array for lapack functions
    float *AA;
    AA = (float*)malloc((rowsA*colsA) *sizeof(float));

    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsA; j++){
            AA[i * colsA + j] = arrA[i][j];
            //printf("Local: %f, In: %f ", AA[i * colsA + j],arrA[i][j]);
        } 
        //printf("\n");
    }
    //printf("\n\n");
    printf("Allocated A.. OK!\n");

    float *BB;
    BB = (float*)malloc((rowsA*1) *sizeof(float));
    for (int i = 0; i < rowsA; i++) {
        BB[i] = arrB[i][0];
    }
    printf("Allocated B.. OK!\n");

    int nrhs = 1; //number of rhs vectors B 
    int LDA = colsA;
    int LDB = rowsA;
    int info;

    int *IPIV = (int*)malloc((rowsA) *sizeof(int));

    // https://www.intel.com/content/www/us/en/docs/onemkl/code-samples-lapack/2022-1/sgesv-example-c.html
    sgesv_(&rowsA, &nrhs, AA, &LDA, IPIV , BB, &LDB, &info); //INTEL DOCS

    if (info >= 0){
        printf("Solution successfull");
    } 

    for (int i = 0; i < rowsA; i++) {
        Usol[i][0]=BB[i];
    }
    printf("Transfered solution from B to Usol.. OK!\n");


}

void modalAnalysis_sggev(int N, float **arrA, float **arrB, float *eigVals){

    /* 
    INFO 
    https://netlib.org/lapack/lug/node35.html
    */

    // EXAMPLE FOR TEST
    /*
    N = 6;
       
    float A[6][6] = {{50.0, -60.0, 50.0, -27.0, 6.0, 6.0},
                    {38.0, -28.0, 27.0, -17, 5.0, 5.0},
                    {27.0, -17.0, 27.0, -17, 5.0, 5.0},
                    {27.0, -28.0, 38.0, -17, 5.0 ,5.0},
                    {27.0, -28.0, 27.0, -17, 16.0, 5.0},
                    {27.0, -28.0, 27.0, -17, 5.0, 16.0}};
                    
    float B[6][6] = {{16.0, 5.0, 5.0, 5.0, -6.0, 5.0},
                    {5.0, 16.0, 5.0, 5.0, -6.0, 5.0},
                    {5.0, 5.0, 16.0, 5.0, -6.0, 5.0},
                    {5.0, 5.0, 5.0, 16.0, -6.0, 5.0},
                    {5.0, 5.0, 5.0, 5.0, -6.0, 16.0},
                    {6.0, 6.0, 6.0, 6.0, -5.0, 6.0,}};
    */

    /* result for the above test case
        cc =

        0.5000 - 0.8660i
        0.5000 + 0.8660i
        0.5000 - 0.8660i
        0.5000 + 0.8660i
            Inf + 0.0000i
            Inf + 0.0000i
    */
    // square matrix N x N 
    float *AA;
    AA = (float*)malloc((N*N) *sizeof(float));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++){
            AA[i * N + j] = arrA[i][j];
            //AA[i * N + j] = A[i][j];
            //printf("Local: %f, In: %f ", AA[i * N + j],arrA[i][j]);
        } 
        //printf("\n");
    }
    //printf("\n\n");
    printf("\nAllocated A.. OK!\n");

    float *BB;
    BB = (float*)malloc((N*N) *sizeof(float));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++){
            BB[i * N + j] = arrB[i][j];
            //BB[i * N + j] = B[i][j];;
            //printf("Local: %f, In: %f ", BB[i * N + j],arrB[i][j]);
        } 
        //printf("\n");
    }
    printf("\nAllocated B.. OK!\n");

    char JOBVL = 'N';
    char JOBVR = 'N'; //the right generalized eigenvectors are computed.
    int LDA = N;
    int LDB = N;

    float *ALPHAR, *ALPHAI, *BETA;
    allocate1Darray(N, &ALPHAR);
    allocate1Darray(N, &ALPHAI);
    allocate1Darray(N, &BETA);

    int LDVL = 1;
    int LDVR = 1;
    float *VL, *VR, *WORK; // BUG
    int sizeVL = LDVL*N;
    int sizeVR = LDVR*N;
    allocate1Darray(sizeVL, &VL);
    allocate1Darray(sizeVR, &VR);

    int sizeWork = 8*N;
    allocate1Darray(sizeWork, &WORK);

    int INFO = 0;
    size_t dummy1;
    size_t dummy2;

    /* 
    ?GGEV: a simple driver that computes all the generalized eigenvalues of (A, B), 
    and optionally the left or right eigenvectors (or both);
    */
    sggev_(&JOBVL, &JOBVR, &N, AA, &LDA, BB, &LDB, ALPHAR,
        ALPHAI, BETA, VL, &LDVL, VR, &LDVR, WORK, &sizeWork, &INFO, dummy1, dummy2);

    printf("INFO = %d, \n", INFO);

    /* PROCESSING TO FIND THE SMALLEST EIGENVALUES */

    float *res = (float*)malloc((N) *sizeof(float));
    double kernel, r1, r2;
    for (int i = 0;i<N;i++){
        r1 = ALPHAR[i]/ BETA[i];
        r2 = ALPHAI[i]/ BETA[i];
        kernel = (double) pow(r1,2)+pow(r2,2);
        res[i] = sqrt(kernel); ///BETA[i]; //real eigenvalues
        //printf("ALPHAR[i]/BETA[i]=%f, ALPHAI[i]/BETA[i]=%f \n", ALPHAR[i]/BETA[i], ALPHAI[i]/BETA[i]);
    }
    char ID = 'I';
    slasrt_(&ID,&N,res,&INFO,dummy1);
    for (int i = 0;i<5;i++){
        printf("Res[i]=%f \n", res[i]);
    }

    /*
    xGGEVX: an expert driver that can additionally balance the matrix pair to improve 
    the conditioning of the eigenvalues and eigenvectors, and compute condition numbers 
    for the eigenvalues and/or left and right eigenvectors (or both). 
    */
   /* TO DO: DEBUG */
    char JOBVSL = 'N', JOBVSR = 'N';
    LDVR = N;
    char sortEigs = 'N';  //un-ordered
    LAPACK_S_SELECT3 selctG; // seg fault
    char sense = 'N';
    int sdimOut;
    float *rcondeOut = (float*)malloc((2) *sizeof(float));
    float *rcondvOut = (float*)malloc((2) *sizeof(float));

    int sizeIWork = 1;
    int *IWORK = (int*)malloc((sizeIWork) *sizeof(int));
    int *BWORK = (int*)malloc((N) *sizeof(int));
    size_t dummy3, dummy4;

    sggesx_(&JOBVSL, &JOBVSR, &sortEigs, selctG, &sense, &N, AA, &LDA, BB, &LDB, &sdimOut,
     ALPHAR, ALPHAI, BETA, VL, &LDVL, VR, &LDVR, rcondeOut,rcondvOut,
     WORK, &sizeWork, IWORK, &sizeIWork, BWORK, &INFO, dummy1, dummy2, dummy3, dummy4);

    printf("\n\n");
    for (int i = 0;i<N;i++){
        r1 = ALPHAR[i]/ BETA[i];
        r2 = ALPHAI[i]/ BETA[i];
        kernel = (double) pow(r1,2)+pow(r2,2);
        res[i] = sqrt(kernel); ///BETA[i]; //real eigenvalues
        //printf("ALPHAR[i]/BETA[i]=%f, ALPHAI[i]/BETA[i]=%f \n", ALPHAR[i]/BETA[i], ALPHAI[i]/BETA[i]);
    }

    slasrt_(&ID,&N,res,&INFO,dummy1);
    for (int i = 0;i<5;i++){
        printf("Resx[i]=%f \n", res[i]);
    }




}
