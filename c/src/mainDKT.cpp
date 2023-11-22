/*
 COMMENT: Make sure that the INDATA_FEM.bin is of the same precision type 
 as the preprocessor directives used in the code. Else, it gives instant
 segmentation fault. 

 COMPLILE using:
 >> gcc -Wall -g3 -fsanitize=address ./src/mainDKT.cpp ./src/funcFEM.cpp ./src/funcMat.cpp -o mainDKT_CPP -lm -lblas -llapack

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

    printf("\n----\n RUNNING IN MODE: %d (1. DOUBLE, 2. SINGLE) \n----\n", PRECISION_MODE_FEM);

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

        #if PRECISION_MODE_FEM == 2 
            mytype p1 = 2.5;
            mytype p2 = 2.5;
        #endif
        #if PRECISION_MODE_FEM == 1 
            //mytype p1 = 20.55;
            //mytype p2 = 20.55;
            mytype p1 = 20.55;
            mytype p2 = 20.55;
        #endif

        mytype *pparam1, *pparam2;
        pparam1 = &p1;  
        pparam2 = &p2; 

        shepard_interp_2d<mytype>(nd, inDataFem.xcp, inDataFem.ycp, inDataFem.fcp, 
        pparam1, ni, wingMeshFem.xm, wingMeshFem.ym, distrLoad);

        shepard_interp_2d<mytype>(nd, inDataFem.xcp, inDataFem.ycp, inDataFem.tcp, 
        pparam2, ni, wingMeshFem.xm, wingMeshFem.ym, distrThick);
    }
    mytype **BeSt;
    allocate2Darray<mytype>(3, 3, &(BeSt)); //[GEN x 1]

    if (inDataFem.LL==2 || inDataFem.LL==1 ){
        BendingStiffness<mytype>(inDataFem.E, inDataFem.v, inDataFem.h, BeSt);
    }

/*
    printf("\n");
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            printf("%f,",BeSt[i][j]);
        } 
        printf("\n");
    }
*/

    /* DKT */
    TrigElCoefsDKT<mytype>(&inDataFem, &wingMeshFem);
    LNShapeFunDST<mytype>(GaussIntegrPoints, xw, &wingMeshFem);
    LNShapeFunMassDST<mytype>(GaussIntegrPoints, xw, &wingMeshFem);
    matrixG<mytype>(&wingMeshFem);
    //
    //int rows, cols;
    squareMatInverse2<mytype>(10, 10, wingMeshFem.GGDST, wingMeshFem.GGin);
    squareMatInverse2<mytype>(6, 6, wingMeshFem.GGDKT, wingMeshFem.GGin2);

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
    //
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

    printf("\n\n");
    //for each triangle in the mesh
    for (int kk = 0;kk<wingMeshFem.Nelem;kk++){   
    //for (int kk = 0;kk<1;kk++){       
        massHmDKT<mytype>(kk, &wingMeshFem, &elemFemArr); // Hm, HW
        rotationMass2<mytype>(kk, &wingMeshFem, &elemFemArr); // Hxx, Hyy
        //------------------------------------------------------------->>
        // for each gauss point
        //for (int ii = 0; ii<GaussIntegrPoints; ii++){
        for (int ii = 0; ii<GaussIntegrPoints; ii++){
        //for (int ii = 0; ii<1; ii++){    
            ShapeFunDKT2<mytype>(ii, kk, &wingMeshFem, &elemFemArr);
            pseudoMassDKT<mytype>(ii, kk, &wingMeshFem, &elemFemArr); // not exactly used (only LW)

            // matrix addition needed 
            // kb=kb+Area(kk)*xw(ii,3)*(Bb'*BeSt2(:,:,kk)*Bb);
            if (inDataFem.LL == 3.0){
                BendingStiffness<mytype>(inDataFem.E, inDataFem.v, distrThick[kk], BeSt);
            }

/*
            for (int i=0;i<3;i++){
                for (int j=0;j<3;j++){
                    printf("%f, ", BeSt[i][j]);
                }
                printf("\n");
            }
*/
            
            mytype **kb;
            allocate2Darray<mytype>(9, 3, &kb); //Bb'*BeSt2(:,:,kk)
            matMatMultiplication2<mytype>(2, 3, 9, 3, 1.0, 0.0, elemFemArr.Bb, BeSt, kb);

            mytype **kb1;
            allocate2Darray<mytype>(9,9,&kb1);
            mytype var1 = wingMeshFem.area[kk] * xw[ii][2];
            matMatMultiplication2<mytype>(1, 9, 3, 9, var1, 0.0, kb, elemFemArr.Bb, kb1);

/*
            for (int i=0;i<9;i++){
                for (int j=0;j<9;j++){
                    printf("%10.4f, ",kb1[i][j]/pow(10.0,5.0));
                }
                printf("\n");
            }
            printf("\n---->");

*/
            matSum2<mytype>(1.0, 0.0, 9, 9, kb1, kb1, elemFemArr.kloc); // kloc = kloc + kb1

/*
            for (int i=0;i<9;i++){
                for (int j=0;j<9;j++){
                    printf("%10.4f, ",elemFemArr.kloc[i][j]/pow(10.0,6.0));
                }
                printf("\n");
            }
            printf("\n---->");
*/            

//==========================DE - ALLOCATE ---->
            deallocate2Darray<mytype>(9,kb);
            deallocate2Darray<mytype>(9,kb1);
            //printf("It worked!");
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
                var0 = mypow<mytype>(distrThick[kk],2)/12.0;
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
            mytype varmloc = inDataFem.mass*inDataFem.h*wingMeshFem.area[kk]*xw[ii][2];
            if (inDataFem.LL==3){
                varmloc = inDataFem.mass*distrThick[kk]*wingMeshFem.area[kk]*xw[ii][2];
            }
            matSum2<mytype>(varmloc, varmloc, 9, 9, term1, term2, sum1); // C = C + a * A + b * B
            matSum2<mytype>(1.0, varmloc, 9, 9, sum1, term5, sum2);
            matSum2<mytype>(1.0, 0.0, 9, 9, sum2, sum2, elemFemArr.mloc);

//==========================DE - ALLOCATE ---->
            deallocate2Darray<mytype>(9,term1);
            deallocate2Darray<mytype>(9,term2);
            deallocate2Darray<mytype>(10,term3);
            deallocate2Darray<mytype>(9,term4);
            deallocate2Darray<mytype>(9,term5);
            deallocate2Darray<mytype>(9,sum1);
            deallocate2Darray<mytype>(9,sum2);
//==========================DE - ALLOCATE ---->
  
            // mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*(mlocTEMP);

        }

/*
        for (int i=0;i<9;i++){
            for (int j=0;j<9;j++){
                printf("%10.4f, ",elemFemArr.kloc[i][j]/pow(10.0,6.0));
            }
            printf("\n");
        }
        printf("\n---->");
*/        

        //------------------------------------------------------------->> for each gauss point
        // lumped mass approach for the uniform load
        mytype lumpedMass[9] = {1, 0, 0, 1, 0, 0, 1, 0, 0};
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

/*
        for (int i = 0;i<9;i++){
            for (int j = 0;j<1;j++){
                printf("%10.10f,",elemFemArr.floc[i][0]);
            }
            printf("\n");
        }
        printf("\n");printf("\n");
*/

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
/*
    for (int i = 0;i<10;i++){
        for (int j = 0;j<10;j++){
            printf("%10.4f,",elemFemArr.Kg[i][j]/pow(10.0,6.0));
        }
        printf("\n");
    }

    printf("--->");
    for (int i = 0;i<10;i++){
        printf("%10.4f,\n",elemFemArr.Fglob[i][0]);
    }
*/

    printf("    Calculated Mg(:,kk), Kg(:,kk), Fglob(kk)");




/*
    printf("    Fglob(kk)");
    for (int i=0;i<10;i++){
        printf("\n%f, ", elemFemArr.Fglob[i][0]);
    }
*/    

    //************************************************************************************
    //  DKT PLATE SOLVER: GLOBAL MATRIX ASSEMBLY (Mglob, Kglob, Fglob)
    //************************************************************************************
    int **iii, **rr, **iii_col, **rr_col;
    allocate2Darray<int>(9,9,&iii);
    allocate2Darray<int>(9,9,&rr);
    allocate2Darray<int>(81,1,&iii_col);
    allocate2Darray<int>(81,1,&rr_col);

    for (int i=0;i<9;i++){
        for (int j=0;j<9;j++){
            iii[i][j] = i;//from 0-8 (instead of 1-9 in matlab)
            rr[i][j] = j;
        }
    }

    int cnt = 0;
    for (int i=0;i<9;i++){
        for (int j=0;j<9;j++){
            iii_col[cnt][0]=iii[j][i];
            rr_col[cnt][0]=rr[j][i];
            cnt++;
        }
    }

    int **Ig, **Jg;
    allocate2Darray<int>(81,wingMeshFem.Nelem,&Ig);//LM(iii(:),:);
    allocate2Darray<int>(81,wingMeshFem.Nelem,&Jg);//LM(rr(:),:);

    int indexIg, indexJg;
    for (int i=0;i<81;i++){
        for (int j=0;j<wingMeshFem.Nelem;j++){
            indexIg = iii_col[i][0];
            indexJg = rr_col[i][0];
            Ig[i][j] = wingMeshFem.LM[indexIg][j]-1;
            Jg[i][j] = wingMeshFem.LM[indexJg][j]-1;
        }
    }   

    mytype **Kglob, **Mglob;
    allocate2Darray<mytype>(wingMeshFem.GEN,wingMeshFem.GEN,&Kglob);
    allocate2Darray<mytype>(wingMeshFem.GEN,wingMeshFem.GEN,&Mglob);

    int indexKi, indexKj;
    for (int i=0;i<81;i++){
        for (int j=0;j<wingMeshFem.Nelem;j++){
            indexKi = Ig[i][j];
            indexKj = Jg[i][j];
            Kglob[indexKi][indexKj] = Kglob[indexKi][indexKj] + elemFemArr.Kg[i][j];
            Mglob[indexKi][indexKj] = Mglob[indexKi][indexKj] + elemFemArr.Mg[i][j];
        }
    }
/*
    printf("\n----\n");
    printf("Mglob[i][j]\n");
    for (int i = 0;i<5;i++){
        for (int j = 0;j<5;j++){
            printf("%f, ",Mglob[i][j]);
        }
        printf("\n");
    }
    printf("----\n");

    printf("\n----\n");
    printf("Kglob[i][j]\n");
    for (int i = 0;i<5;i++){
        for (int j = 0;j<5;j++){
            printf("%f, ",Kglob[i][j]);
        }
        printf("\n");
    }
    printf("----\n");
    exit(55);
*/

    printf("\n    Kglob, Mglob OK... DKT PLATE SOLVER: AUGMENTED GLOBAL MATRIX (for BCs)\n");

    //************************************************************************************
    //  DKT PLATE SOLVER: AUGMENTED GLOBAL MATRIX (for BCs)
    //************************************************************************************
    mytype **kkk, **mmm;
    allocate2Darray<mytype>(inDataFem.sizeBdofs,wingMeshFem.GEN,&kkk);
    allocate2Darray<mytype>(inDataFem.sizeBdofs,wingMeshFem.GEN,&mmm);

    int index_kkk;
    for (int j=0;j<inDataFem.sizeBdofs;j++){
        index_kkk = inDataFem.Bdofs[j]-1;
        kkk[j][index_kkk] = 1.0;
        //kkk(j,:)=[zeros(1,Bdofs(j)-1) 1 zeros(1,Dofs-Bdofs(j))]; %matlab code sample
    }

    //Kglob=[Kglob kkk'; kkk zeros(length(BBnodes))];  %matlab code sample
    //Mglob=[Mglob mmm'; mmm zeros(length(BBnodes))];  %matlab code sample
    mytype **Kglob_aug, **Mglob_aug; // augmented
    int sizeKMglob_aug = wingMeshFem.GEN+inDataFem.sizeBdofs;
    allocate2Darray<mytype>(sizeKMglob_aug,sizeKMglob_aug,&Kglob_aug);
    allocate2Darray<mytype>(sizeKMglob_aug,sizeKMglob_aug,&Mglob_aug);

    for (int i=0;i<wingMeshFem.GEN;i++){
        for (int j=0;j<wingMeshFem.GEN;j++){
            Kglob_aug[i][j] = Kglob[i][j];
            Mglob_aug[i][j] = Mglob[i][j]; //OK
        }
    }

    int indexKaug;
    for (int i=0;i<inDataFem.sizeBdofs;i++){
        for (int j=0;j<wingMeshFem.GEN;j++){
            indexKaug = wingMeshFem.GEN + i;
            Kglob_aug[indexKaug][j] = kkk[i][j];
            Kglob_aug[j][indexKaug] = kkk[i][j];
        }
    }

/*
    printf("\n----\n");
    printf("Mglob_aug[i][j]\n");
    for (int i = 0;i<5;i++){
        for (int j = 0;j<5;j++){
            printf("%f, ",Mglob_aug[i][j]);
        }
        printf("\n");
    }
    printf("----\n");

    printf("\n----\n");
    printf("Kglob_aug[i][j]\n");
    for (int i = 0;i<5;i++){
        for (int j = 0;j<5;j++){
            printf("%f, ",Kglob_aug[i][j]);
        }
        printf("\n");
    }
    printf("----\n");
*/

    //************************************************************************************
    //  DKT PLATE SOLVER: SOLUTION OPTIONS (1. EIGEN, 2. STATIC, 3. DYNAMIC)
    //************************************************************************************
    printf("    Starting linear system solution (Kglob_aug, Mglob_aug are dense matrices!)\n");

    int indexPointLoad;
    if (inDataFem.LL == 1){
        indexPointLoad = wingMeshFem.ID[0][inDataFem.P_node-1]-1;
        elemFemArr.Fglob[indexPointLoad][0] = inDataFem.P_load;

        //printf("index point load %d", indexPointLoad);
    }

    mytype **Usol;
    allocate2Darray<mytype>(sizeKMglob_aug,1,&Usol);

    mytype **Fglob_aug;
    allocate2Darray<mytype>(sizeKMglob_aug,1,&Fglob_aug);
    for (int i=0;i<wingMeshFem.GEN;i++){
        Fglob_aug[i][0]=elemFemArr.Fglob[i][0];
    }

    if (inDataFem.LL==2){
        printf("    UNIFORM LOAD: P=%f [Pa]\n",inDataFem.P_load);
    }

    // solve linear system of eqs. using LAPACK sgels_ function
    //linearSystemSolve<mytype>(sizeKMglob_aug, sizeKMglob_aug, Kglob_aug, Fglob_aug, Usol);
    //printf("    SOLVED LINEAR SYSTEM\n");

    //************************************************************************************
    //  DKT PLATE SOLVER: OUTPUT BINARY FILE for Matlab Post-Processor
    //************************************************************************************
    //int optionSelect = 0;

    //CuFEMNum2DWriteDataInBinary<mytype>(sizeKMglob_aug, 1, Usol, wingMeshFem.GEN);

    //CuFEMNum2DWriteMatrix<mytype>(sizeKMglob_aug, sizeKMglob_aug, Kglob_aug, Mglob_aug, Fglob_aug);

    //----
    // TODO: MODAL ANALYSIS
    #if (MODAL_ANALYSIS == 1)
    //  TO DO : FIX - I GET INFINITE EIGENVALUES
        printf("\n    Performing MODAL ANALYSIS");
        // Generalized Nonsymmetric Eigenvalue Problems TO DO: It doesnt work as matlab
        
        //http://matlab.izmiran.ru/help/techdoc/ref/eig.html
        // Real nonsymmetric A, real general B: sggev() from LAPACK
        mytype *eigVals;
        int n_eigs = 6;
        allocate1Darray<mytype>(n_eigs,&eigVals);
        printf(" (w/2/pi)");
        myeigs<mytype>(sizeKMglob_aug, Kglob_aug, Mglob_aug, n_eigs, eigVals);
    #endif
    //----

    exit(55);

    #if (DYNAMIC_ANALYSIS == 1)
        // select method for time integration
        printf("\n    Starting DYNAMIC ANALYSIS using Crank-Nicolson.\n");

        if (inDataFem.LL == 3){
            int d = 0; // index for time integration
            mytype t = 0.0;
            mytype dt = inDataFem.dt;
            mytype Tp = 2*M_PI/inDataFem.omega3; // period of motion
            int NtimeSteps = ceil((2*Tp)/dt)+1;
            //int NtimeSteps = ceil((inDataFem.Nper*Tp)/dt)+1;
            printf("    dt=%10.4f, Tp=%10.4f, NtimeSteps=%d \n", dt, Tp, NtimeSteps);

            //==========DAMPING MATRIX============
            mytype **Cdamp; // Rayleigh damping matrix
            allocate2Darray<mytype>(sizeKMglob_aug,sizeKMglob_aug,&Cdamp);

            mytype a=0, b=0;
            RayleighDampingCoefs<mytype>(&a, &b); // TO DO (based on eigenfrequencies)
            matSum2<mytype>(a, b, sizeKMglob_aug, sizeKMglob_aug, Mglob_aug, Kglob_aug, Cdamp); //C = a*Mglob+b*Kglob;

            printf("    Rayleigh damping coefs. a=%10.4f, b=%10.4f\n",a,b);

            printf("\nCdamp: \n");
            for (int i=0;i<10;i++){
                for (int j=0;j<10;j++){
                    printf("%10.4f,", Cdamp[i][j]/mypow<mytype>(10.0,3.0));
                }
                printf("\n");
            }
            //printf("\nCdamp(end,end)=%10.4f\n", Cdamp[sizeKMglob_aug-1][sizeKMglob_aug-1]/mypow<mytype>(10.0,3.0));
            //==========DAMPING MATRIX============

            //writeMatrixInBinary(sizeKMglob_aug, sizeKMglob_aug, Cdamp);

            // Create loading vector G(:,d) for the current time step
            mytype **G;
            int sz1 = sizeKMglob_aug;
            int sz2 = 2*sz1;
            allocate2Darray<mytype>(sz2, NtimeSteps, &G); //[G(:,d), G(:,d+1)]
            
            createRHS<mytype>(&inDataFem, &wingMeshFem, &elemFemArr,
                             distrLoad, G, d); //G(:,d)

            //writeMatrixInBinary(sz2, NtimeSteps, G);

            printf("\n    sizeKMglob_aug=%d, GEN=%d \n", sizeKMglob_aug, wingMeshFem.GEN);
            printf("    sz1 = %d, sz2 = %d\n",sz1,sz2);
            printf("\n    G(1:10,%d)=\n",d);
            for (int i=0;i<10;i++){
                printf("    %10.8f, ",G[i][d]);
            }

            //==========TIME INTEGRATION============
            mytype theta = 0.5; //Crank-Nicolson
            //
            mytype **u_t; // u(:,d)
            allocate2Darray(sz2,NtimeSteps,&u_t); //u=[qdot;q]
            //
            mytype **q; //displacement unknown vector (similar to Ustatic)
            mytype **qdot; //velocity
            mytype **qdot2, **qdot2_buffer; //acceleration
            allocate2Darray<mytype>(sz1,NtimeSteps,&q); 
            allocate2Darray<mytype>(sz1,NtimeSteps,&qdot); 
            allocate2Darray<mytype>(sz1,NtimeSteps,&qdot2); 
            allocate2Darray<mytype>(sz1,1,&qdot2_buffer); 

            mytype beta = 0.25, gamma = 0.5;
            /*====initialization=======================

            AA=Mglob + gamma*ddt*C + ddt^2*beta*Kglob;
            BB=Fglob_t - C*qdot(:,1)- Kglob*q(:,1);
            qdot2(:,1)=AA\BB;
            
            ===========================================*/
            mytype **AA,**BB;
            allocate2Darray<mytype>(sz1,sz1,&AA); 
            allocate2Darray<mytype>(sz1,1,&BB); 
            //printf("\nsz1=%d",sz1);

            for (int i=0; i<sz1;i++){
                for (int j=0; j<sz1; j++){
                    AA[i][j] = Mglob_aug[i][j] + gamma*dt*Cdamp[i][j] + mypow<mytype>(dt,2.0)*beta*Kglob_aug[i][j];
                    BB[i][0] = BB[i][0] -Cdamp[i][j]*qdot[j][0] - Kglob_aug[i][j]*q[j][0];
                }
                BB[i][0] = BB[i][0] + G[i][0];//maybe problematic but we will see
            }

            for (int i = 0;i<sz1;i++){
                qdot2_buffer[i][0] = 0; // re-initialize
            }
            linearSystemSolve<mytype>(sz1, sz1, AA, BB, qdot2_buffer);

            //printf("qdot2_buffer[295][0]=%10.15f, ",qdot2_buffer[295][0]);

            for (int i = 0;i<sz1;i++){
                qdot2[i][d]=qdot2_buffer[i][0];
            }

            mytype **pr_vel, **pr_disp;
            allocate2Darray<mytype>(sz1,1,&pr_vel);
            allocate2Darray<mytype>(sz1,1,&pr_disp);

            mytype suma1, suma2;
            for (int d = 1; d< NtimeSteps-1 ; d++){  
            //for (int d = 0; d< 100 ; d++){   
                printf("\n    d (time) = %d\n",d);   
                t = t + dt; // t in [sec]

                createRHS<mytype>(&inDataFem, &wingMeshFem, &elemFemArr,
                            distrLoad, G, d+1);//G(:,d+1)
/*
                printf("\n    G(1:10,%d)=\n",d);
                for (int i=0;i<10;i++){
                    printf("%10.8f, ",G[i][d+1]);
                }            
*/
                /* 
                // DOES NOT WORK!!!! ILL-CONDITIONED DENSE SYSTEM (dgesv_) fails
                timeIntegrationCN((d+1), dt, theta, sz2, G, Mglob_aug, Kglob_aug, Cdamp, u_t); // TIME INTEGRATION WITH CRANK-NICOLSON 
                */

                // TO DO: Use newmark since it is the same with crank-nicolson

                /*
                pr_vel = qdot(:,d)+(1-gamma)*ddt*qdot2(:,d);% + gamma*hhh*qdot2(:,d);
                pr_disp = q(:,d)+ddt*qdot(:,d)+ddt^2*(1/2-beta)*qdot2(:,d);%+hhh^2*beta*qdot2(:,d);

                AA = Mglob + gamma*ddt*C + ddt^2*beta*Kglob;
                BB = Fglob_t - C*pr_vel- Kglob*pr_disp;
                qdot2(:,d+1) = AA\BB; 
                */
                for (int i = 0; i <sz1;i++) {
                    pr_vel[i][0] = qdot[i][d] + (1.0-gamma)*dt*qdot2[i][d];
                    pr_disp[i][0] = q[i][d] + dt*qdot[i][d] + mypow<mytype>(dt, 2.0)*(0.5-beta)*qdot2[i][d];
                }
/*
                printf("\n    pr_vel=\n");
                for (int i=0;i<10;i++){
                    printf("%10.8f, ",pr_vel[i][0]);
                }

                printf("\n    pr_disp=\n");
                for (int i=0;i<10;i++){
                    printf("%10.8f, ",pr_disp[i][0]);
                }
*/
                // re - initialize BB 
                for (int i=0; i<sz1;i++){
                    BB[i][0] = 0.0;
                }

                //printf("%\n");
                for (int i=0; i<sz1;i++){
                    suma1 = 0;
                    suma2 = 0;
                    for (int j=0; j<sz1;j++){
                        //AA[i][j] = Mglob_aug[i][j] + gamma*dt*Cdamp[i][j] + mypow<mytype>(dt,2.0)*beta*Kglob_aug[i][j];
                        suma1 = suma1 -Cdamp[i][j]*pr_vel[j][0];
                        suma2 = suma2 -Kglob_aug[i][j]*pr_disp[j][0];  
                        //printf("%10.8f,%10.8f,%10.8f\n", -Cdamp[i][j], pr_vel[j][0]/pow(10.0,-6.0), suma1);
                    }
                    //printf("%f", suma1);
                    BB[i][0] = suma1 + suma2 + G[i][d+1];//maybe problematic but we will see
                }
/*
                printf("\n   BB=\n");
                // re - initialize BB 
                for (int i=0; i<10;i++){
                    printf("%10.15f, ", BB[i][0]);
                }
*/
                //printf("qdot2[295][0]=%10.15f, ",qdot2[295][0]);

                for (int i = 0;i<sz1;i++){
                    qdot2_buffer[i][0] = 0; // re-initialize
                }
                linearSystemSolve<mytype>(sz1, sz1, AA, BB, qdot2_buffer);

                for (int i = 0;i<sz1;i++){
                    qdot2[i][d+1]=qdot2_buffer[i][0];
                }
/*
                printf("\n qdot2=\n");
                for (int i = 0;i<10;i++){
                    printf("%f,",qdot2[i][d+1]);
                }
*/

                for (int i=0; i<sz1; i++){
                    q[i][d+1] = pr_disp[i][0] + mypow<mytype>(dt, 2.0)*beta*qdot2[i][d+1];
                    qdot[i][d+1] = pr_vel[i][0] + gamma*dt*qdot2[i][d+1];
                }
/*
                printf("\n q=\n");
                for (int i = 0;i<10;i++){
                    printf("%f,",q[i][d+1]/pow(10.0,-7));
                }

                printf("\n qdot=\n");
                for (int i = 0;i<10;i++){
                    printf("%f,",qdot[i][d+1]/pow(10.0,-4));
                }

*/              
                // return solution to the u_t vector
                for (int i=0; i<sz1; i++){
                    u_t[i][d+1] = qdot[i][d+1];
                    u_t[i+sz1][d+1] = q[i][d+1];
                }

/*
                printf("\n\nu=\n");
                for (int i = 0;i<10;i++){
                    //for (int j=0;j<1;j++){
                    printf("    %10.8f, ",u_t[i][d+1]);
                    //}
                    printf("\n");
                }
*/
                //exit(55);

                //timeIntegrationNewmark<mytype>(); // TIME INTEGRATION WITH CRANK-NICOLSON 
                //u(:,d+1) = timeIntegration(u, d+1, GEN, Mglob, Kglob, C, G, ddt, theta); %[w,bx,by,lambda]

            }

            CuFEMNum2DWriteDataInBinary<mytype>(sz2, NtimeSteps, u_t, wingMeshFem.GEN);

            printf("\n\nG=\n");
            for (int i = 0;i<10;i++){
                for (int j=0;j<10;j++){
                    printf("    %10.8f, ",G[i][j]);
                }
                printf("\n");
            }

            printf("\n\nu=\n");
            for (int i = 0;i<10;i++){
                for (int j=0;j<10;j++){
                    if (u_t[i][j]>0.0){
                        printf("    %10.9f, ",u_t[i][j]);
                    }
                    else{
                        printf("    %10.8f, ",u_t[i][j]);
                    }
                    
                }
                printf("\n");
            }

            printf("\n\nu(sz1-4:sz1+10,1:10)=\n");
            for (int i = sz1-4;i<sz1+10;i++){
                for (int j=0;j<10;j++){
                    if (u_t[i][j]>0.0){
                        printf("    %10.9f, ",u_t[i][j]);
                    }
                    else{
                        printf("    %10.8f, ",u_t[i][j]);
                    }
                }
                printf("\n");
            }

            printf("UP TP HERE..\n");

            deallocate2Darray<mytype>(sz1,Cdamp);
            deallocate2Darray<mytype>(sz1,AA); 
            deallocate2Darray<mytype>(sz1,BB);
            //
            deallocate2Darray<mytype>(sz2,u_t);  
            deallocate2Darray<mytype>(sz2,G); //[G(:,d), G(:,d+1)] 
            deallocate2Darray<mytype>(sz1,q); 
            deallocate2Darray<mytype>(sz1,qdot); 
            deallocate2Darray<mytype>(sz1,qdot2); 
            deallocate2Darray<mytype>(sz1,qdot2_buffer);
            deallocate2Darray<mytype>(sz1,pr_vel);
            deallocate2Darray<mytype>(sz1,pr_disp);
           

        }
        else{
            printf("    Attempting DYNAMIC ANALYSIS with uniform or point load. Please support time\n"
                    "dependent forcing data");
        }     

    #endif


    deallocate2Darray<int>(9,iii);
    deallocate2Darray<int>(9,rr);

    deallocate2Darray(inDataFem.sizeBdofs,kkk);
    deallocate2Darray(inDataFem.sizeBdofs,mmm);
    //
    deallocate2Darray<mytype>(sizeKMglob_aug,Kglob_aug);
    deallocate2Darray<mytype>(sizeKMglob_aug,Mglob_aug);

    deallocate2Darray(sizeKMglob_aug,Usol);
    deallocate2Darray<mytype>(sizeKMglob_aug,Fglob_aug);

    tend = clock();

    double cpu_time_used = ((double) (tend-tstart))/ CLOCKS_PER_SEC;
    printf("\n----\n Elapsed time [s]: %f\n----\n", cpu_time_used);
    printf("In Matlab the same operations using vectorization take 2.0383 sec.\n");
    
    // DE-ALLOCATE MEMORY TO RESOLVE MEMORY LEAKS
    free(distrLoad);
    free(distrThick);
    deallocate2Darray<mytype>(3,BeSt);
    //
    freeInDataRecFem(&inDataFem);
    freetriangleDKT(GaussIntegrPoints,&wingMeshFem);
    freefemArraysDKT<mytype>(&wingMeshFem, &elemFemArr);    
    deallocate2Darray<mytype>(wingMeshFem.GEN,Kglob);
    deallocate2Darray<mytype>(wingMeshFem.GEN,Mglob);
    //
    deallocate2Darray<int>(81,iii_col);
    deallocate2Darray<int>(81,rr_col);
    deallocate2Darray<int>(81,Ig);//LM(iii(:),:);
    deallocate2Darray<int>(81,Jg);//LM(rr(:),:);

    return 0;

}


