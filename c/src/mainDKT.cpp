/*
 COMMENT: Make sure that the INDATA_FEM.bin is of the same precision type 
 as the preprocessor directives used in the code. Else, it gives instant
 segmentation fault. 

 COMPLILE using:
 >> gcc -Wall -g3 -fsanitize=address mainDKT.cpp -lm

*/
#include "../include/mainDKT.h"

#include "../src/funcFEM.cpp" // Karperaki functions for DKT fem
#include "../src/funcFSI.cpp" // Functions used for the coupling of bem - fem 
#include "../src/funcMat.cpp" // Functions used to facilitate martix, vector operations in c


/*=========================================================================================*/
/* MAIN PROGRAM BELOW */
/*=========================================================================================*/
int main(int argc, char **argv){

    printf("    RUNNING IN MODE: %d (1. DOUBLE, 2. SINGLE)", PRECISION_MODE_FEM);

    clock_t tstart, tend;
    tstart = clock();

    /* read input parameters from a file */
    /* Boundary conditions nodes, dofs */
    struct InDataRecFem<mytype> inDataFem;
    CuFEMNum2DReadInData<mytype>(&inDataFem);

    /* Create or load from matlab IEN, ID, LM */
    struct triangleDKT<mytype> wingMeshFem;  
    ConnectivityFEM_IEN_ID_LM<mytype>(&inDataFem, &wingMeshFem); // BUG FOUND IN PREVIOUS VERSIONS in IEN_3

    /* Gauss integration function - read about it */
    mytype xw[GaussIntegrPoints][3]; // {xg,yg,wg}
    TriGaussPoints<mytype>(xw);

    /* Distributed properties */
    mytype *distrLoad = NULL, *distrThick = NULL;
    allocate1Darray<mytype>(wingMeshFem.Nelem,&distrLoad);
    allocate1Darray<mytype>(wingMeshFem.Nelem,&distrThick);

    if (inDataFem.LL == 3){
        /* DISTRIBUTED LOAD & THICKNESS CASE */
        int nd = inDataFem.sizexcp;
        int ni = wingMeshFem.Nelem; //size(xm)

        mytype p1 = 10.55;
        mytype p2 = 10.55;
        mytype *pparam1, *pparam2;
        pparam1 = &p1;  
        pparam2 = &p2; 

        shepard_interp_2d<mytype>(nd, inDataFem.xcp, inDataFem.ycp, inDataFem.fcp, 
        pparam1, ni, wingMeshFem.xm, wingMeshFem.ym, distrLoad);

        shepard_interp_2d<mytype>(nd, inDataFem.xcp, inDataFem.ycp, inDataFem.tcp, 
        pparam2, ni, wingMeshFem.xm, wingMeshFem.ym, distrThick);

        printf("distrThick[i], distrLoad[i]\n");
        for (int i = 0; i<10; i++){
            printf("%f, %f\n", distrThick[i], distrLoad[i]);
        }
    }

    mytype BeSt[3][3] = {0};
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            printf("%f,",BeSt[i][j]);
        } 
        printf("\n");
    }
    if (inDataFem.LL==2 || inDataFem.LL==1){
        BendingStiffness<mytype>(inDataFem.E, inDataFem.v, inDataFem.h, BeSt);
    }

    
    /* DKT */
    TrigElCoefsDKT<mytype>(&inDataFem, &wingMeshFem);
    //
    LNShapeFunDST<mytype>(GaussIntegrPoints, xw, &wingMeshFem);
    //
    LNShapeFunMassDST<mytype>(GaussIntegrPoints, xw, &wingMeshFem);
    //
    matrixG(&wingMeshFem);
    




    
    tend = clock();

    double cpu_time_used = ((double) (tend-tstart))/ CLOCKS_PER_SEC;
    printf("\n----\n Elapsed time [s]: %f\n----\n", cpu_time_used);
    printf("In Matlab the same operations using vectorization take 2.0383 sec.\n");


    // DE-ALLOCATE MEMORY TO RESOLVE MEMORY LEAKS
    free(distrLoad);
    free(distrThick);
    freeInDataRecFem(&inDataFem);
    freetriangleDKT(GaussIntegrPoints,&wingMeshFem);


    return 0;

}


