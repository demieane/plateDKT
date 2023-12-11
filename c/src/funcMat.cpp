#ifndef MAIN_FEM_HEADER_FILE
    #define MAIN_FEM_HEADER_FILE

    //#include "../include/mainFem.h"


#endif

#include <stdio.h>
#include <stdlib.h>
#include <cstddef> // for NULL
#include <math.h> 

//
template<class T>
T mypow(T base, T power);
//
template<class T>
T mysqrt(T base);
//
template<class T>
T myabs(T var);
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
void deallocate2Darray(int rows, T **arrIn); // BUG
//
template<class T>
void linearSystemSolve(int rowsA, int colsA, T **arrA, T **arrB, T **Usol);
//
template<class T>
void myeigs(int N, T **arrA, T **arrB, int n_eigs, T *eigVals);
// TIME INTEGRATION WITH CRANK-NICOLSON
template<class T>
void timeIntegrationCN(int d, T dt, T theta, int rowsColsG, T **G, T **Mglob_aug, T **Kglob_aug, T **C, T **u_t); 
// TIME INTEGRATION WITH NEWMARK
template<class T>
void timeIntegrationNewmark(); 
//
template<class T>
void writeMatrixInBinary(int rowsA, int colsA, T **arrA);
//
template<class T>
void writeMatrixInTextFile(int rowsA, int colsA, T **arrA);
//
// from funcBLAS.c
/*=========================================================================================*/
/* Definition of the functions BELOW */
/*=========================================================================================*/




template<class T>
T mypow(T base, T power){

    #if PRECISION_MODE_FEM == 1 //double
        return pow(base,power);
    #endif
    #if PRECISION_MODE_FEM == 2 //single
        return powf(base,power);
    #endif
}

template<class T>
T myabs(T var){

    #if PRECISION_MODE_FEM == 1 //double
        return fabs(var);
    #endif
    #if PRECISION_MODE_FEM == 2 //single
        return abs(var);
    #endif
}

template<class T>
T mysqrt(T base){

    #if PRECISION_MODE_FEM == 1 //double
        return sqrt(base);
    #endif
    #if PRECISION_MODE_FEM == 2 //single
        return sqrtf(base);
    #endif
}

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
        arrTemp[i] = (T)0;
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
    BB = (T*)malloc((rowsA) *sizeof(T));
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
        // This works OK for the single precision case.
        sgesv_(&rowsA, &nrhs, AA, &LDA, IPIV , BB, &LDB, &info); //INTEL DOCS
    #endif
    #if PRECISION_MODE_FEM == 1

    #if DEBUG_ON
        size_t dummyVar=1;
        const char NORM = 'I';
        T *workNorm;
        allocate1Darray<T>(rowsA,&workNorm);
        T ANORM;
        ANORM=dlange_(&NORM, &rowsA, &colsA, AA, &LDA, workNorm, dummyVar);
        printf("\n    ANORM=%f, ", ANORM);
        T RCOND;
        T *WORK;
        int *IWORK;
        allocate1Darray<T>(4*rowsA,&WORK);
        allocate1Darray<int>(rowsA,&IWORK);
        dgecon_(&NORM,&rowsA,AA,&LDA,&ANORM,&RCOND,WORK,IWORK,&info, dummyVar);
        printf("RCOND=%10.8f\n",RCOND);

        free(workNorm);
        free(WORK);
        free(IWORK);

        printf("\n    rowsA=%d,colsA=%d\n",rowsA,colsA);
        /*checking sparse of matrix*/
        int count = 0;
        T epsTEST = mypow<T>(10.0, -5);
        for (int i = 0; i < rowsA; i++){
            for (int j = 0; j < colsA; j++){
                if( myabs<T>(arrA[i][j])< epsTEST){
                    count = count + 1;
                }     
            }
        }
        int criterion = (rowsA * colsA)/2;
        if (count > ((rowsA * colsA)/2))
            printf("    arrA: Matrix is a sparse matrix, count=%d, criterion=%d\n",count, criterion);
        else
            printf("    arrA: Matrix is not sparse matrix, count=%d, criterion=%d\n",count, criterion);
    
        /*
        The reciprocal condition number is a scale-invariant measure 
        of how close a given matrix is to the set of singular matrices. 
        If C is near 0, the matrix is nearly singular and badly conditioned.
        If C is near 1.0, the matrix is well conditioned.
        */
    #endif
        //dgesv_(&rowsA, &nrhs, &(AA[0]), &LDA, IPIV , &(BB[0]), &LDB, &info); //INTEL DOCS
        dgesv_(&rowsA, &nrhs, AA, &LDA, IPIV , BB, &LDB, &info); //INTEL DOCS

        //size_t dummyVar = 0;
        //const char TRANS = 'N';
        //dgetrs_(&TRANS,&rowsA,&nrhs,AA,&LDA,IPIV,BB,&LDB,&info,dummyVar); //INTEL DOCS
    #endif


    if (info == 0){
        //printf("\n    Solution successfull\n");
 
        for (int i = 0; i < rowsA; i++) {
            Usol[i][0]=BB[i];
            //printf("%10.4f,%10.4f,\n",BB[i]/pow(10.0,8.0), Usol[i][0]/pow(10.0,8.0));
        }

/*
        for (int i=0;i<10;i++){
            printf("    %f,", Usol[i][0]);
        }
*/        
    }
    else{
        printf("\n    Solution failed in sgesv_: info=%d", info);
    } 


    free(IPIV);
    free(BB);
    free(AA);
    


}

template<class T>
void myeigs(int N, T **arrA, T **arrB, int n_eigs, T *eigVals){

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
    T *AA, *BB;
    AA = (T*)malloc((N*N) *sizeof(T));
    BB = (T*)malloc((N*N) *sizeof(T));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++){
            AA[i * N + j] = arrA[i][j];
            BB[i * N + j] = arrB[i][j];
        } 
    }

    char JOBVL = 'N';
    char JOBVR = 'N'; //the right generalized eigenvectors are computed.
    int LDA = N;
    int LDB = N;

    T *ALPHAR, *ALPHAI, *BETA;
    allocate1Darray<T>(N, &ALPHAR);
    allocate1Darray<T>(N, &ALPHAI);
    allocate1Darray<T>(N, &BETA);

    int LDVL = 1;
    int LDVR = 1;
    T **VL, **VR, *WORK; // BUG
    //int sizeVL = LDVL*N;
    //int sizeVR = LDVR*N;
    allocate2Darray<T>(LDVL, N, &VL);
    allocate2Darray<T>(LDVR, N, &VR);

    int sizeWork = 8*N;
    allocate1Darray<T>(sizeWork, &WORK);

    int INFO = 0;
    size_t dummy1 = 0;
    size_t dummy2 = 0;

    /* 
    ?GGEV: a simple driver that computes all the generalized eigenvalues of (A, B), 
    and optionally the left or right eigenvectors (or both);
    */
    #if PRECISION_MODE_FEM == 1
        dggev_(&JOBVL, &JOBVR, &N, AA, &LDA, BB, &LDB, ALPHAR,
            ALPHAI, BETA, VL[0], &LDVL, VR[0], &LDVR, WORK, &sizeWork, &INFO, dummy1, dummy2);
    #endif
    #if PRECISION_MODE_FEM == 2
        sggev_(&JOBVL, &JOBVR, &N, AA, &LDA, BB, &LDB, ALPHAR,
            ALPHAI, BETA, VL[0], &LDVL, VR[0], &LDVR, WORK, &sizeWork, &INFO, dummy1, dummy2);
    #endif

    if (INFO == 0) 
    {
        printf(", INFO: %d. \n", INFO);
        /* PROCESSING TO FIND THE SMALLEST EIGENVALUES */

        T *res = (T*)malloc((N) *sizeof(T));
        T kernel, r1, r2;
        for (int i = 0;i<N;i++){
            r1 = ALPHAR[i]/ BETA[i];
            r2 = ALPHAI[i]/ BETA[i];
            kernel = mysqrt<T>(mypow<T>(r1,2.0)+mypow<T>(r2,2.0));
            res[i] = mysqrt<T>(kernel)/2.0/M_PI;
        }
        char ID = 'I';
        #if PRECISION_MODE_FEM == 1
            dlasrt_(&ID,&N,res,&INFO,dummy1);
        #endif
        #if PRECISION_MODE_FEM == 2
            slasrt_(&ID,&N,res,&INFO,dummy1);
        #endif
        for (int i = 0;i<n_eigs;i++){
            eigVals[i]=res[i];
            printf("        Eigvals[%d]=%f \n",i, eigVals[i]);
        }
        

        free(res);
    }
    else 
    {
        printf("\n INFO: %d, \n", INFO);
    }
    
    
    
    free(AA);
    free(BB);
    free(WORK);
    free(ALPHAR);
    free(ALPHAI);
    free(BETA);
    
    deallocate2Darray<T>(LDVL, VL);
    deallocate2Darray<T>(LDVR, VR);

}

template<class T>
void timeIntegrationCN(int d, T dt, T theta, int rowsColsG, T **G, T **Mglob_aug, T **Kglob_aug, T **C, T **u_t){

// MAKE MATRICES FOR THE SECOND ORDER SYSTEM OF EQS. //
/*
    sizeM=size(Mglob,1);

    A = [Mglob, sparse(sizeM,sizeM); sparse(sizeM,sizeM), speye(sizeM,sizeM)];
    B = -[C, Kglob; -speye(sizeM,sizeM), sparse(sizeM,sizeM)];

    %lamda (1: implicit Euler, 1/2: Crank - Nicolson)

    AA =  A - theta*dt*B;
    BB =  A + (1 - theta)*dt*B;

    MAT = [a, b; 
           c, d]
*/

    int sz1 = rowsColsG/2;
    int sz2 = rowsColsG;

    T **A, **B, **AA, **BB;
    allocate2Darray<T>(sz2,sz2,&A);
    allocate2Darray<T>(sz2,sz2,&B);
    allocate2Darray<T>(sz2,sz2,&AA);
    allocate2Darray<T>(sz2,sz2,&BB);

    
    for (int i = 0;i<sz1;i++){
        for (int j = 0;j<sz1;j++){
            // a: part of matrix
            A[i][j] = Mglob_aug[i][j];
            B[i][j] = -C[i][j];
            // b: part of matrix
            B[i][j+sz1] = -Kglob_aug[i][j];
            // Ieye
            if (i == j){
                // c: part of matrix
               B[i+sz1][j] = 1.0;
                // d: part of matrix
               A[i+sz1][j+sz1] = 1.0;
            }
        }
    }

    //writeMatrixInBinary<T>(sz2, sz2, A);
    //writeMatrixInBinary<T>(sz2, sz2, B);
    //writeMatrixInBinary<T>(sz2, sz2, AA);
    //writeMatrixInBinary<T>(sz2, sz2, BB);

    //exit(55);

    //AA =  A - theta*dt*B;
    matSum2<T>(1.0, (-theta*dt), sz2, sz2, A, B, AA);
    //BB =  A + (1 - theta)*dt*B;
    matSum2<T>(1.0, ( (1.0-theta)*dt ), sz2, sz2, A, B, BB);

    //writeMatrixInBinary<T>(sz2, sz2, A);
    //writeMatrixInBinary<T>(sz2, sz2, B);
    //writeMatrixInBinary<T>(sz2, sz2, AA);
    //writeMatrixInBinary<T>(sz2, sz2, BB);

    //==================================================================
    //Q = (1 - theta)*dt*G(:,d-1) + (theta)*dt*G(:,d);
    T **Q;
    allocate2Darray<T>(sz2,1,&Q);

    for (int i = 0; i<sz2; i++){
        Q[i][0] = (1.0-theta)*dt*G[i][d-1] + theta*dt*G[i][d];
    }

    printf("\nQ[%d]=",d);
    for (int i = 0;i<10;i++){
        printf("%10.5f, ", Q[i][0]/pow(10.0,-6.0));
    }
    printf("\nQ[%d]=",d);
    for (int i = sz1;i<sz1+10;i++){
        printf("%10.5f, ", Q[i][0]/pow(10.0,-6.0));
    }

    //==================================================================

    printf("\nu_t[%d]=",d);
    for (int i = 0;i<10;i++){
        printf("%10.6f, ", u_t[i][d-1]);
    }

    T **rhs;//, **rhsDEBUG;
    allocate2Darray<T>(sz2,1,&rhs);
    //allocate2Darray<T>(sz2,1,&rhsDEBUG);

    for (int ii = 0; ii<sz2; ii++){
        for (int jj = 0; jj<sz2; jj++){      
            //rhs[ii][0] = rhs[ii][0] + BB[jj][ii]*u_t[jj][d-1];//matrix multiplication
            rhs[ii][0] = rhs[ii][0] + BB[ii][jj]*u_t[jj][d-1];//matrix multiplication
        }
        rhs[ii][0] = rhs[ii][0] + Q[ii][0];//matrix addition
    }

    printf("\nQ[%d]=",d);
    for (int i = 0;i<10;i++){
        printf("%10.5f, ", rhs[i][0]/pow(10.0,-6.0));
    }

/*
    int ii = 0;
    //for (int ii = 0; ii<sz2; ii++){
        for (int jj = 0; jj<sz2; jj++){      
            //rhs[ii][0] = rhs[ii][0] + BB[jj][ii]*u_t[jj][d-1];//matrix multiplication
            rhsDEBUG[ii][0] = rhsDEBUG[ii][0] + BB[jj][ii];//matrix multiplication
        }
        //rhs[ii][0] = rhs[ii][0] + Q[ii][0];//matrix addition
    //}

    printf("\nrhs = %f\n", rhsDEBUG[ii][0]);
*/
    //exit(55);
/*
    printf("\nutemp2[%d]=",d);
    for (int i = 0;i<10;i++){    
        printf("%10.5f, ", rhs[i][0]/pow(10.0,-6.0));
    }
*/    
    //==================================================================
    T **Usol;
    allocate2Darray(sz2,1,&Usol);

/*
    FILE *file;
    file = fopen("test_lin_solve.bin", "rb"); // r for read, b for binary
    //printf("file: %d", file);
    int rows_dummy1, cols_dummy2;
    fread(&(rows_dummy1), sizeof(T) , 1, file);
    fread(&(cols_dummy2), sizeof(T) , 1, file);
    
    for (int i=0;i<rows_dummy1;i++){
        for (int j=0;j<cols_dummy2;j++){
            fread(&(AA[i][j]), sizeof(T), 1, file);
        } 
    }
    printf("AA: %d,%d", rows_dummy1, cols_dummy2);
    //
    int rows_dummy3, cols_dummy4;
    fread(&(rows_dummy3), sizeof(T) , 1, file);
    fread(&(cols_dummy4), sizeof(T) , 1, file);
    for (int i=0;i<rows_dummy3;i++){
        for (int j=0;j<cols_dummy4;j++){
            fread(&(rhs[i][j]), sizeof(T), 1, file);
        } 
    }
    fclose(file);
*/
    writeMatrixInTextFile<T>(sz2, sz2, AA);
    writeMatrixInTextFile<T>(sz2, 1, rhs);



    linearSystemSolve<T>(sz2, sz2, AA, rhs, &(Usol[0]));

    //if (d == 2){
    //writeMatrixInBinary<T>(sz2, sz2, AA);
    //writeMatrixInBinary<T>(sz2, 1, rhs);

    //}

    for (int i = 0; i<sz2; i++){
        u_t[i][d]=Usol[i][0];
    }

    printf("\nu[%d]=",d);
    for (int i = 0;i<10;i++){
        printf("    %10.6f, ",Usol[i][0]);
    }

    printf("\nu[%d]=",d);
    for (int i = sz1-4;i<sz1+10;i++){
        printf("    %10.15f, ",Usol[i][0]);
    }

    //==================================================================

    deallocate2Darray<T>(sz2,A);
    deallocate2Darray<T>(sz2,B);
    deallocate2Darray<T>(sz2,AA);
    deallocate2Darray<T>(sz2,BB);
    //
    deallocate2Darray<T>(rowsColsG,Q);
    deallocate2Darray<T>(rowsColsG,rhs);
    deallocate2Darray<T>(rowsColsG,Usol);

    exit(55);
}


template<class T>
void writeMatrixInBinary(int rowsA, int colsA, T **arrA){

    //If the file exists it just overwrites it (delete DEBUG file
    //prior to each execution)

    FILE *fileOut;
	fileOut = fopen("../c/OUTDATA_FEM_DEBUG.bin", "ab"); // w for write, b for binary

    fwrite(&rowsA, sizeof(int), 1, fileOut);
    fwrite(&colsA, sizeof(int), 1, fileOut);

    for (int i = 0; i < rowsA; i++){
        for (int j = 0; j < colsA; j++){
            fwrite(&(arrA[i][j]), sizeof(T), 1, fileOut);
        }
    }

    fclose(fileOut);
    printf("\n    Exiting writeMatrixInBinary...");

}

template<class T>
void writeMatrixInTextFile(int rowsA, int colsA, T **arrA){

    //If the file exists it just overwrites it (delete DEBUG file
    //prior to each execution)

    FILE *fileOut;
	fileOut = fopen("../c/OUTDATA_FEM_DEBUG.txt", "ab"); // w for write, b for binary

    fprintf(fileOut, "%d\n", rowsA);
    fprintf(fileOut, "%d\n", colsA);

    for (int i = 0; i < rowsA; i++){
        for (int j = 0; j < colsA; j++){
            fprintf(fileOut, "%.15f\n", arrA[i][j]);
        }
    }

    fclose(fileOut);
    printf("\n    Exiting writeMatrixInTextFile...");

}