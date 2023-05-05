/*
 COMMENT: Make sure that the INDATA_FEM.bin is of the same precision type 
 as the preprocessor directives used in the code. Else, it gives instant
 segmentation fault. 

 COMPLILE using:
 >> gcc -Wall -g3 -fsanitize=address mainDKT.cpp -lm -lblas -llapack

*/
#include "../include/mainDKT.h"

#include "../src/funcFEM.cpp" // Karperaki functions for DKT fem
#include "../src/funcFSI.cpp" // Functions used for the coupling of bem - fem 

#ifndef FUNCMAT
    #define FUNCMAT

    #include "../src/funcMat.cpp" // Functions used to facilitate martix, vector operations in c
#endif

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
    LNShapeFunDST<mytype>(GaussIntegrPoints, xw, &wingMeshFem);
    LNShapeFunMassDST<mytype>(GaussIntegrPoints, xw, &wingMeshFem);
    matrixG<mytype>(&wingMeshFem);
    //
    //int rows, cols;
    squareMatInverse2<mytype>(10, 10, wingMeshFem.GGDST, wingMeshFem.GGin);
    squareMatInverse2<mytype>(6, 6, wingMeshFem.GGDKT, wingMeshFem.GGin2);

    printf("\n GGin \n");
    for (int i=0;i<10;i++){
        for (int j=0;j<10;j++){
            printf("%f,",wingMeshFem.GGin[i][j]);
        }
        printf("\n");
    }
    printf("\n GGin2 \n");
    for (int i=0;i<6;i++){
        for (int j=0;j<6;j++){
            printf("%f,",wingMeshFem.GGin2[i][j]);
        }
        printf("\n");
    }
    //************************************************************************************
    //  DKT PLATE SOLVER: LOCAL MATRIX (mloc, kloc, floc)
    //************************************************************************************
    if (inDataFem.LL == 1){
         printf("\n    P_load=%f in [Pa]", inDataFem.P_load);
    }
    /*
    initialize structure that contains all the arrays
    needed for the final global matrix assembly
    */
    struct femArraysDKT<mytype> elemFemArr;

    allocate2Darray<mytype>(wingMeshFem.GEN, 1, &(elemFemArr.Fglob)); //[GEN x 1]
    allocate2Darray<mytype>(10, 9, &(elemFemArr.Hm)); 
    allocate2Darray<mytype>(10, 9, &(elemFemArr.HW)); 
    allocate2Darray<mytype>(9, 9, &(elemFemArr.kloc));
    allocate2Darray<mytype>(9, 9, &(elemFemArr.mloc));
    //==========================================
    allocate2Darray<mytype>(9, 1, &(elemFemArr.floc)); // TODO : floc, floc1 uniform/distributed load
    allocate2Darray<mytype>(10, 1, &(elemFemArr.floc1)); // point load
    //==========================================
    allocate2Darray<mytype>(6, 9, &(elemFemArr.Hxx));
    allocate2Darray<mytype>(6, 9, &(elemFemArr.Hyy));
    //
    allocate2Darray<mytype>(1, 9, &(elemFemArr.Hx)); // [1 x 9]
    allocate2Darray<mytype>(1, 9, &(elemFemArr.Hy)); // [1 x 9]
    allocate1Darray<mytype>(9, &(elemFemArr.Hx_xsi)); // [1 x 9]
    allocate1Darray<mytype>(9, &(elemFemArr.Hx_eta)); // [1 x 9]
    allocate1Darray<mytype>(9, &(elemFemArr.Hy_xsi)); // [1 x 9]
    allocate1Darray<mytype>(9, &(elemFemArr.Hy_eta)); // [1 x 9]
    //
    allocate2Darray<mytype>(1, 10, &(elemFemArr.LW)); //[1 x 10]
    allocate1Darray<mytype>(6, &(elemFemArr.L)); 
    //
    allocate2Darray<mytype>(3, 9, &(elemFemArr.Bb));
    allocate2Darray<mytype>(81, wingMeshFem.Nelem, &(elemFemArr.Mg));
    allocate2Darray<mytype>(81, wingMeshFem.Nelem, &(elemFemArr.Kg));

    //for each triangle in the mesh
    for (int kk = 0;kk<wingMeshFem.Nelem;kk++){
        massHmDKT<mytype>(kk, &wingMeshFem, &elemFemArr); // Hm, HW
        rotationMass2<mytype>(kk, &wingMeshFem, &elemFemArr); // Hxx, Hyy

    }




    
    tend = clock();

    double cpu_time_used = ((double) (tend-tstart))/ CLOCKS_PER_SEC;
    printf("\n----\n Elapsed time [s]: %f\n----\n", cpu_time_used);
    printf("In Matlab the same operations using vectorization take 2.0383 sec.\n");


    // DE-ALLOCATE MEMORY TO RESOLVE MEMORY LEAKS
    free(distrLoad);
    free(distrThick);
    freeInDataRecFem(&inDataFem);
    freetriangleDKT(GaussIntegrPoints,&wingMeshFem);

    freefemArraysDKT(&wingMeshFem, &elemFemArr);


    return 0;

}


