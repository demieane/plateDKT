#ifndef MAIN_FEM_HEADER_FILE
#define MAIN_FEM_HEADER_FILE

#include "mainFem.h"

#endif

#include <math.h>

void assignRowArrayMatrixG_DST(int rowID, float xsi, float eta, float **array){

    float rowMat[10] = {1.0, xsi, eta, xsi*eta, pow(xsi,2), pow(eta,2),
     pow(xsi,2)*eta, pow(eta,2)*xsi, pow(xsi,3), pow(eta,3)};

    //for (int i=0;i<10;i++){
    //    printf("%f, ", rowMat[i]);
    //}

    //printf("\nInside assignRowArrayMatrixG_DST()... \n\n");
    int i = rowID -1;
    for (int j=0;j<10;j++){
        array[i][j]=rowMat[j];
        //printf("%d, %f, ",j, array[I][j]);
    }
    //printf("\n\n");
}

void assignRowArrayMatrixG_DKT(int rowID, float xsi, float eta, float **array){

    float rowMat[6] = {1.0, xsi, eta, xsi*eta, pow(xsi,2), pow(eta,2)};

    //for (int i=0;i<6;i++){
    //    printf("%f, ", rowMat[i]);
    //}

    //printf("\nInside assignRowArrayMatrixG_DKT()... \n\n");
    int i = rowID -1;
    for (int j=0;j<6;j++){
        array[i][j]=rowMat[j];
    //    printf("%d, %f, ",j, array[I][j]);
    }
    //printf("\n\n");
}