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

    mytype **BeSt = NULL;
    allocate2Darray<mytype>(3, 3, &(BeSt)); //[GEN x 1]
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
        //------------------------------------------------------------->>
        // for each gauss point
        for (int ii = 0; ii<GaussIntegrPoints; ii++){
            ShapeFunDKT2<mytype>(ii, kk, &wingMeshFem, &elemFemArr);
            pseudoMassDKT<mytype>(ii, kk, &wingMeshFem, &elemFemArr); // not exactly used (only LW)

            // matrix addition needed 
            // kb=kb+Area(kk)*xw(ii,3)*(Bb'*BeSt2(:,:,kk)*Bb);
            if (inDataFem.LL == 3.0){
                BendingStiffness<mytype>(inDataFem.E, inDataFem.v, distrThick[kk], BeSt);
            }
            mytype **kb;
            allocate2Darray<mytype>(9, 3, &kb);
            matMatMultiplication2<mytype>(2, 3, 9, 3, 1.0, 0.0, elemFemArr.Bb, BeSt, kb);

            mytype **kb1;
            allocate2Darray<mytype>(9,9,&kb1);
            mytype var1 = wingMeshFem.area[kk] * xw[ii][2];
            matMatMultiplication2<mytype>(1, 9, 3, 9, var1, 0.0, kb, elemFemArr.Bb, kb1);

            matSum2<mytype>(1.0, 0.0, 9, 9, kb1, kb1, elemFemArr.kloc); // kloc = kloc + kb1

//==========================DE - ALLOCATE ---->
            for (int i=0;i<9;i++){
                free(kb[i]);
            }
            free(kb);
            //
            for (int i=0;i<9;i++){
                free(kb1[i]);
            }
            free(kb1);   
//==========================DE - ALLOCATE ---->

            if (inDataFem.LL == 1){
                // floc1=floc1+Area(kk)*xw(ii,3)*(LW');
                for (int i=0;i<10;i++){
                    elemFemArr.floc1[i][0]= elemFemArr.floc1[i][0]+ wingMeshFem.area[kk]*xw[ii][2]*elemFemArr.LW[0][i]; 
                }
            }

                // mlocTEMP = (HW'*(LW'*LW)*HW)+
                // txxBEM(kk)^2/12*(Hx'*Hx)+
                // txxBEM(kk)^2/12*(Hy'*Hy)
            
            mytype **term1, **term2, **term3, **term4, **term5;
            allocate2Darray<mytype>(9,9,&term1);
            allocate2Darray<mytype>(9,9,&term2);
            allocate2Darray<mytype>(10,10,&term3);
            allocate2Darray<mytype>(9,10,&term4);
            allocate2Darray<mytype>(9,9,&term5);
  
            mytype var0;
            var0 = pow(inDataFem.h,2)/12.0;
            if (inDataFem.LL==3){
                var0 = pow(distrThick[kk],2)/12.0;
            }
            
            matMatMultiplication2<mytype>(2, 1, 9, 9, var0, 0.0, elemFemArr.Hx, elemFemArr.Hx, term1); //Hx'*Hx
            //
            matMatMultiplication2<mytype>(2, 1, 9, 9, var0, 0.0, elemFemArr.Hy, elemFemArr.Hy, term2); //Hy'*Hy
            //
            matMatMultiplication2<mytype>(2, 1, 10, 10, 1.0, 0.0, elemFemArr.LW, elemFemArr.LW, term3); //LW'*LW
            //
            matMatMultiplication2<mytype>(2, 10, 9, 10, 1.0, 0.0, elemFemArr.HW, term3, term4); //HW'*(LW'*LW) -> [10 x 9] [10 x 10]
            //
            matMatMultiplication2<mytype>(1, 9, 10, 9, 1.0, 0.0, term4, elemFemArr.HW, term5); //(HW'*(LW'*LW)*HW) -> [9 x 10] [10 x 9]

            //  mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*((HW'*(LW'*LW)*HW)+txxBEM(kk)^2/12*(Hx'*Hx)+txxBEM(kk)^2/12*(Hy'*Hy));
            mytype **sum1, **sum2;
            allocate2Darray<mytype>(9,9,&sum1);
            allocate2Darray<mytype>(9,9,&sum2);
            float varmloc = inDataFem.mass*inDataFem.h*wingMeshFem.area[kk]*xw[ii][2];
            if (inDataFem.LL==3){
                varmloc = inDataFem.mass*distrThick[kk]*wingMeshFem.area[kk]*xw[ii][2];
            }
            matSum2<mytype>(varmloc, varmloc, 9, 9, term1, term2, sum1); // C = C + a * A + b * B
            matSum2<mytype>(1.0, varmloc, 9, 9, sum1, term5, sum2);
            matSum2<mytype>(1.0, 0.0, 9, 9, sum2, sum2, elemFemArr.mloc);

//==========================DE - ALLOCATE ---->
            for (int i=0;i<9;i++){
                free(term1[i]);
                free(term2[i]);
                free(term4[i]);
                free(term5[i]);
                free(sum1[i]);
                free(sum2[i]);
            }
            free(term1);
            free(term2);
            free(term4);
            free(term5);
            free(sum1);
            free(sum2);
            //
            for (int i=0;i<10;i++){
                free(term3[i]);
            }
            free(term3);
//==========================DE - ALLOCATE ---->
  
            // mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*(mlocTEMP);

        }
        //------------------------------------------------------------->> for each gauss point
        // lumped mass approach for the uniform load
        float lumpedMass[9] = {1, 0, 0, 1, 0, 0, 1, 0, 0};
        if (inDataFem.LL == 2){
            // version - 1
            for (int i=0;i<9;i++){
                for (int j=0;j<1;j++){
                    elemFemArr.floc[i][j] = elemFemArr.floc[i][j] + wingMeshFem.area[kk]*inDataFem.P_load/3.0*lumpedMass[i];
                }
            }
            // version - 2 (EQUIVALENT) floc=P*HW'*floc1;
            //matMatMultiplication2(2, 10, 9, inDataFem.P_load, 1.0, 0.0, elemFemArr.HW, elemFemArr.floc1, elemFemArr.floc); //HW'*floc1
        }
        if (inDataFem.LL == 3){
            // version - 1
            //printf(" LLL = 3, %f,",distrLoad[kk]);
            for (int i=0;i<9;i++){
                for (int j=0;j<1;j++){
                    elemFemArr.floc[i][j] = elemFemArr.floc[i][j] + wingMeshFem.area[kk]*distrLoad[kk]/3.0*lumpedMass[i];
                }
            }
        }

        //Mg(:,kk)=[mloc(:,1);mloc(:,2);mloc(:,3);mloc(:,4);mloc(:,5);mloc(:,6);mloc(:,7);mloc(:,8);mloc(:,9)];
        //Kg(:,kk)=[kloc(:,1);kloc(:,2);kloc(:,3);kloc(:,4);kloc(:,5);kloc(:,6);kloc(:,7);kloc(:,8);kloc(:,9)];
        int cntMg = 0;
        for (int i=0;i<9;i++){
            for (int j=0;j<9;j++){
                elemFemArr.Mg[cntMg][kk] = elemFemArr.mloc[j][i];
                elemFemArr.Kg[cntMg][kk] = elemFemArr.kloc[j][i];
                cntMg++;
                // i++ will increment the value of i, but return the original value that i held before being incremented.
            }
        }

        if (inDataFem.LL == 2 || inDataFem.LL==3){
            for (int q=0;q<9;q++){
                int cntFglob = wingMeshFem.LM[q][kk]-1;
                elemFemArr.Fglob[cntFglob][0] = elemFemArr.Fglob[cntFglob][0] + elemFemArr.floc[q][0];
                //Fglob(LM(q,kk))=Fglob(LM(q,kk))+floc(q);
            }
        }

        // re-initialize mloc, kloc, floc, floc1
        for (int i = 0;i<9;i++){
            elemFemArr.floc[i][0] = 0;
            for (int j = 0;j<9;j++){
                elemFemArr.mloc[i][j] = 0;
                elemFemArr.kloc[i][j] = 0;
            }
        }
        for (int i = 0;i<10;i++){
            //printf("%f, ",elemFemArr.floc1[i][0] );
            elemFemArr.floc1[i][0] = 0;
        }


    }

    printf("\n    Calculated Mg(:,kk), Kg(:,kk), Fglob(kk)");

    //************************************************************************************
    //  DKT PLATE SOLVER: GLOBAL MATRIX ASSEMBLY (Mglob, Kglob, Fglob)
    //************************************************************************************


    
    tend = clock();

    double cpu_time_used = ((double) (tend-tstart))/ CLOCKS_PER_SEC;
    printf("\n----\n Elapsed time [s]: %f\n----\n", cpu_time_used);
    printf("In Matlab the same operations using vectorization take 2.0383 sec.\n");


    // DE-ALLOCATE MEMORY TO RESOLVE MEMORY LEAKS
    free(distrLoad);
    free(distrThick);
    //
    for (int i = 0;i<3;i++){
        free(BeSt[i]);
    }
    free(BeSt);
    //
    freeInDataRecFem(&inDataFem);
    freetriangleDKT(GaussIntegrPoints,&wingMeshFem);

    freefemArraysDKT(&wingMeshFem, &elemFemArr);


    return 0;

}


