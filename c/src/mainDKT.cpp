/*
 COMMENT: Make sure that the INDATA_FEM.bin is of the same precision type 
 as the preprocessor directives used in the code. Else, it gives instant
 segmentation fault. 

 COMPLILE using:
 >> gcc -g mainDKT.cpp -lm

*/
#include "../include/mainDKT.h"

#include "../src/funcFEM.cpp"
#include "../src/funcFSI.cpp"
#include "../src/funcMat.cpp"


/*=========================================================================================*/
/* MAIN PROGRAM BELOW */
/*=========================================================================================*/
int main(int argc, char **argv){

    clock_t tstart, tend;
    tstart = clock();

    /* read input parameters from a file */
    /* Boundary conditions nodes, dofs */
    struct InDataRecFem<mytype> inDataFem;
    CuFEMNum2DReadInData(&inDataFem );
    /* Create or load from matlab IEN, ID, LM */
    struct triangleDKT<mytype> wingMeshFem;  
    ConnectivityFEM_IEN_ID_LM(&inDataFem, &wingMeshFem); // BUG FOUND IN PREVIOUS VERSIONS in IEN_3

    /* Gauss integration function - read about it */
    mytype xw[GaussIntegrPoints][3]; // {xg,yg,wg}
    TriGaussPoints(xw);

    /* Distributed properties */
    mytype *distrLoad, *distrThick;
    allocate1Darray(wingMeshFem.Nelem,&distrLoad);
    allocate1Darray(wingMeshFem.Nelem,&distrThick);

    if (inDataFem.LL == 3){
        /* DISTRIBUTED LOAD & THICKNESS CASE */
        int nd = inDataFem.sizexcp;
        int ni = wingMeshFem.Nelem; //size(xm)

        mytype p1 = 10.55;
        mytype p2 = 10.55;
        mytype *pparam1, *pparam2;
        pparam1 = &p1;  
        pparam2 = &p2; 

        shepard_interp_2d(nd, inDataFem.xcp, inDataFem.ycp, inDataFem.fcp, 
        pparam1, ni, wingMeshFem.xm, wingMeshFem.ym, distrLoad);

        shepard_interp_2d(nd, inDataFem.xcp, inDataFem.ycp, inDataFem.tcp, 
        pparam2, ni, wingMeshFem.xm, wingMeshFem.ym, distrThick);

        printf("distrThick[i], distrLoad[i]\n");
        for (int i = 0; i<10; i++){
            printf("%f, %f\n", distrThick[i], distrLoad[i]);
        }
    }
    float **BeSt;
    allocate2Darray(3, 3, &BeSt);
    if (inDataFem.LL==2 || inDataFem.LL==1){
        BendingStiffness(inDataFem.E, inDataFem.v, inDataFem.h, BeSt);
    }


    



    
    tend = clock();

    double cpu_time_used = ((double) (tend-tstart))/ CLOCKS_PER_SEC;
    printf("\n----\n Elapsed time [s]: %f\n----\n", cpu_time_used);
    printf("In Matlab the same operations using vectorization take 2.0383 sec.\n");


    freeInDataRecFem(&inDataFem);
    int Ng = GaussIntegrPoints;
    //freetriangleDKT(Ng,&wingMeshFem);

    return 0;

}


