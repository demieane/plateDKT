/*This form is used for header files of your own program.
It searches for a file named 'file' in the directory containing the current file.
You can prepend directories to this list with the -I 
option while compiling your source code. */
#include "../include/mainFem.h"
#include <time.h>
/* 
Search file for TODO
Search file for COMMENT
*/

/* TODO: single or double precision? I have the project in C <float> & <double> */

/*=========================================================================================*/
/*  The main program follows. */
/*=========================================================================================*/
int main(int argc, char **argv){

    printf("\n----\n plateDKT solver in single precision");
    printf("\n----\n\n");
    /* Variable declarations */
    clock_t tstart, tend;
    int Ng = 3; // Gauss integration points
    struct InDataRecFem inDataFem;
    struct triangleDKT wingMeshFem;   

    tstart = clock();

    //************************************************************************************
    //  DKT PLATE SOLVER: PREPARATIONS
    //************************************************************************************
    /* read input parameters from a file */
    /* Boundary conditions nodes, dofs */
    CuFEMNum2DReadInData( &inDataFem );
    /* Create or load from matlab IEN, ID, LM */
    ConnectivityFEM_IEN_ID_LM( &inDataFem, &wingMeshFem ); // BUG FOUND IN PREVIOUS VERSIONS in IEN_3

    float *distrLoad, *distrThick;
    allocate1Darray(wingMeshFem.Nelem,&distrLoad);
    allocate1Darray(wingMeshFem.Nelem,&distrThick);

    if (inDataFem.LL == 3){
        printf(" SHEPARD'S INTERP \n");

        /* DISTRIBUTED LOAD & THICKNESS CASE */
        int nd = inDataFem.sizexcp;
        int ni = wingMeshFem.Nelem; //size(xm)

        float p1 = 10.55, p2 = 10.55;
        float *pparam1, *pparam2;
        pparam1 = &p1;  
        pparam2 = &p2; 

        shepard_interp_2d(nd, inDataFem.xcp, inDataFem.ycp, inDataFem.fcp, 
        pparam1, ni, wingMeshFem.xm, wingMeshFem.ym, distrLoad);

        shepard_interp_2d(nd, inDataFem.xcp, inDataFem.ycp, inDataFem.tcp, 
        pparam2, ni, wingMeshFem.xm, wingMeshFem.ym, distrThick);
    }
   
    /* Gauss integration function - read about it */
    float xw[Ng][3]; // {xg,yg,wg}
    TriGaussPoints(Ng, xw);

    float **BeSt;
    allocate2Darray(3, 3, &BeSt);
    if (inDataFem.LL==2 || inDataFem.LL==1){
        BendingStiffness(inDataFem.E, inDataFem.v, inDataFem.h, BeSt);
    }

    /* DKT */
    TrigElCoefsDKT(&inDataFem, &wingMeshFem);
    LNShapeFunDST(Ng, xw, &wingMeshFem);
    LNShapeFunMassDST(Ng, xw, &wingMeshFem);
    matrixG(&wingMeshFem);

    //int rows, cols;
    squareMatInverse2(10, 10, wingMeshFem.GGDST, wingMeshFem.GGin);
    squareMatInverse2(6, 6, wingMeshFem.GGDKT, wingMeshFem.GGin2);


#if DEBUG_ON
    printf("wingMeshFem.GGDST\n");
    for (int i=0;i<10;i++){
        for (int j=0;j<10;j++){
            printf("%f, ",wingMeshFem.GGDST[i][j]);
        }
        printf("\n\n");
    }

    printf("wingMeshFem.GGin\n");
    for (int i=0;i<10;i++){
        for (int j=0;j<10;j++){
            printf("%f, ",wingMeshFem.GGin[i][j]);
        }
        printf("\n");
    }

    printf("wingMeshFem.GGDKT\n");
    for (int i=0;i<6;i++){
        for (int j=0;j<6;j++){
            printf("%f, ",wingMeshFem.GGDKT[i][j]);
        }
        printf("\n\n");
    }

    printf("wingMeshFem.GGin2\n");
    for (int i=0;i<6;i++){
        for (int j=0;j<6;j++){
            printf("%f, ",wingMeshFem.GGin2[i][j]);
        }
        printf("\n");
    }
#endif


#if DEBUG_ON
    printf("from matrixG() \n");
    for (int i=0;i<10;i++){
        for (int j=0;j<10;j++){
            printf("%f, ", wingMeshFem.GGDST[i][j]);
        }
        printf("\n\n");
    }
    printf("from matrixG() \n");
    for (int i=0;i<6;i++){
        for (int j=0;j<6;j++){
            printf("%f, ", wingMeshFem.GGDKT[i][j]);
        }
        printf("\n\n");
    }
#endif

#if DEBUG_ON
    for (int i=0;i<3;i++){
        //for (int j=0;j<5;j++){
            int j = wingMeshFem.NN-1;
            //printf("wingMeshFem->IEN[%d][%d]= %d\n,   ", i,j, wingMeshFem->IEN[i][j]);
            printf("wingMeshFem.ID[%d][%d]= %d\n,   ", i,j, wingMeshFem.ID[i][j]);
        //}
    }

    for (int i=0;i<3;i++){
        //for (int j=0;j<5;j++){
            int j = wingMeshFem.Nelem-1;
            //printf("wingMeshFem->IEN[%d][%d]= %d\n,   ", i,j, wingMeshFem->IEN[i][j]);
            printf("wingMeshFem.IEN[%d][%d]= %d\n,   ", i,j, wingMeshFem.IEN[i][j]);
        //}
    }
#endif
    
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
    struct femArraysDKT elemFemArr;

    allocate2Darray(wingMeshFem.GEN, 1, &(elemFemArr.Fglob)); //[GEN x 1]
    allocate2Darray(10, 9, &(elemFemArr.Hm)); 
    allocate2Darray(10, 9, &(elemFemArr.HW)); 
    //
    allocate2Darray(9, 9, &(elemFemArr.kloc));
    allocate2Darray(9, 9, &(elemFemArr.mloc));
    //==========================================
    allocate2Darray(9, 1, &(elemFemArr.floc)); // TODO : floc, floc1 uniform/distributed load
    allocate2Darray(10, 1, &(elemFemArr.floc1)); // point load
    //==========================================
    //
    allocate2Darray(6, 9, &(elemFemArr.Hxx));
    allocate2Darray(6, 9, &(elemFemArr.Hyy));
    //
    allocate2Darray(1, 9, &(elemFemArr.Hx)); // [1 x 9]
    allocate2Darray(1, 9, &(elemFemArr.Hy)); // [1 x 9]

    /* NEW */
    allocate1Darray(9, &(elemFemArr.Hx_xsi)); // [1 x 9]
    allocate1Darray(9, &(elemFemArr.Hx_eta)); // [1 x 9]
    allocate1Darray(9, &(elemFemArr.Hy_xsi)); // [1 x 9]
    allocate1Darray(9, &(elemFemArr.Hy_eta)); // [1 x 9]

    allocate2Darray(1, 10, &(elemFemArr.LW)); //[1 x 10]
    allocate1Darray(6, &(elemFemArr.L)); 
    //
    allocate2Darray(3, 9, &(elemFemArr.Bb));
    allocate2Darray(81, wingMeshFem.Nelem, &(elemFemArr.Mg));
    allocate2Darray(81, wingMeshFem.Nelem, &(elemFemArr.Kg));


#if DEBUG_ON
    float **arrIn; 
    allocate2Darray(2, 2, &arrIn); // TODO
    printf("\n (arrIn==NULL) = %d\n", (arrIn==NULL));

    printf("arrIn: %f:",arrIn[0][0]);


    allocate2Darray(3, 9, &(elemFemArr.Bb)); // TODO
    printf("\n (arrIn==NULL) = %d\n", (elemFemArr.Bb==NULL));

    printf("elemFemArr.Bb: %f:",elemFemArr.Bb[0][0]);


    allocate1Darray(9, &(elemFemArr.Fglob)); // TODO
    printf("\n (elemFemArr.Fglob==NULL) = %d\n", (elemFemArr.Fglob==NULL));

    printf("elemFemArr.Fglob: %f:",elemFemArr.Fglob[0]);
#endif

    //for each triangle in the mesh
    for (int kk = 0;kk<wingMeshFem.Nelem;kk++){
    //for (int kk = 0;kk<1;kk++){   

        massHmDKT(kk, &wingMeshFem, &elemFemArr); // Hm, HW
        rotationMass2(kk, &wingMeshFem, &elemFemArr); // Hxx, Hyy

#if DEBUG_ON
        printf("\nPrinting HW...(in main)\n");
        for (int i=0;i<10;i++){
            for (int j=0;j<9;j++){
                printf("%f, ", elemFemArr.HW[i][j]);
            }
            printf("\n");
        }

        printf("\nPrinting GGin...(in main)\n");
        for (int i=0;i<10;i++){
            for (int j=0;j<9;j++){
                printf("%f, ", wingMeshFem.GGin[i][j]);
            }
            printf("\n");
        }
#endif
        //------------------------------------------------------------->>
        // for each gauss point
        //for (int ii = 0; ii<1; ii++){
        for (int ii = 0; ii<Ng; ii++){

            //printf(" ENTERING ShapeFunDKT2\n");
            ShapeFunDKT2(ii, kk, &wingMeshFem, &elemFemArr);
            //printf(" EXITING ShapeFunDKT2\n");
            pseudoMassDKT(ii, kk, &wingMeshFem, &elemFemArr); // not exactly used (only LW)

#if DEBUG_ON
            //The C compiler can glue adjacent string literals into one
            printf("\n------------------------------\n"
                   "  kloc calculations: area(kk)=%f, gauss weight xw(ii,3)=%f\n"
                   "------------------------------\n"
                   ,wingMeshFem.area[kk], xw[ii][2]);
            // TODO: Use the available CBLAS & LAPACKE routines 
#endif            
            // matrix addition needed 
            // kb=kb+Area(kk)*xw(ii,3)*(Bb'*BeSt2(:,:,kk)*Bb);
            if (inDataFem.LL == 3.0){
                BendingStiffness(inDataFem.E, inDataFem.v, distrThick[kk], BeSt);
            }
            float **kb;
            allocate2Darray(9, 3, &kb);
            matMatMultiplication2(2, 3, 9, 3, 1.0, 0.0, elemFemArr.Bb, BeSt, kb);

            float **kb1;
            allocate2Darray(9,9,&kb1);
            float var1 = wingMeshFem.area[kk] * xw[ii][2];
            matMatMultiplication2(1, 9, 3, 9, var1, 0.0, kb, elemFemArr.Bb, kb1);

            matSum2(1.0, 0.0, 9, 9, kb1, kb1, elemFemArr.kloc); // kloc = kloc + kb1

            free(kb1);
            free(kb);

#if DEBUG_ON
            printf("\n kloc (in main)\n");
            for (int i=0;i<9;i++){
                for (int j=0;j<9;j++){
                    printf("%f,",elemFemArr.kloc[i][j]);
                }
                printf("\n");
            }
#endif         

            if (inDataFem.LL == 1){
                // floc1=floc1+Area(kk)*xw(ii,3)*(LW');
                for (int i=0;i<10;i++){
                    elemFemArr.floc1[i][0]= elemFemArr.floc1[i][0]+ wingMeshFem.area[kk]*xw[ii][2]*elemFemArr.LW[0][i]; 
                }
            }

#if DEBUG_ON
            printf("\n------------------------------\n"
                   "  mloc calculations\n"
                   "------------------------------\n");
#endif 
                // mlocTEMP = (HW'*(LW'*LW)*HW)+
                // txxBEM(kk)^2/12*(Hx'*Hx)+
                // txxBEM(kk)^2/12*(Hy'*Hy)
            
            float **term1, **term2, **term3, **term4, **term5;
            allocate2Darray(9,9,&term1);
            allocate2Darray(9,9,&term2);
            allocate2Darray(10,10,&term3);
            allocate2Darray(9,10,&term4);
            allocate2Darray(9,9,&term5);

            //printf("h = %f,\n",inDataFem.h);    
            float var0;
            var0 = pow(inDataFem.h,2)/12.0;
            if (inDataFem.LL==3){
                var0 = pow(distrThick[kk],2)/12.0;
            }
            //printf("txx^2/12*1000 = %f\n",var0*1000);
            
            matMatMultiplication2(2, 1, 9, 9, var0, 0.0, elemFemArr.Hx, elemFemArr.Hx, term1); //Hx'*Hx
            //
            matMatMultiplication2(2, 1, 9, 9, var0, 0.0, elemFemArr.Hy, elemFemArr.Hy, term2); //Hy'*Hy
            //
            matMatMultiplication2(2, 1, 10, 10, 1.0, 0.0, elemFemArr.LW, elemFemArr.LW, term3); //LW'*LW
            //
            matMatMultiplication2(2, 10, 9, 10, 1.0, 0.0, elemFemArr.HW, term3, term4); //HW'*(LW'*LW) -> [10 x 9] [10 x 10]
            //
            matMatMultiplication2(1, 9, 10, 9, 1.0, 0.0, term4, elemFemArr.HW, term5); //(HW'*(LW'*LW)*HW) -> [9 x 10] [10 x 9]

            //  mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*((HW'*(LW'*LW)*HW)+txxBEM(kk)^2/12*(Hx'*Hx)+txxBEM(kk)^2/12*(Hy'*Hy));
            float **sum1, **sum2;
            allocate2Darray(9,9,&sum1);
            allocate2Darray(9,9,&sum2);
            float varmloc = inDataFem.mass*inDataFem.h*wingMeshFem.area[kk]*xw[ii][2];
            if (inDataFem.LL==3){
                varmloc = inDataFem.mass*distrThick[kk]*wingMeshFem.area[kk]*xw[ii][2];
            }
            matSum2(varmloc, varmloc, 9, 9, term1, term2, sum1); // C = C + a * A + b * B
            matSum2(1.0, varmloc, 9, 9, sum1, term5, sum2);
            matSum2(1.0, 0.0, 9, 9, sum2, sum2, elemFemArr.mloc);

            free(term1);
            free(term2);
            free(term3);
            free(term4);
            free(term5);
            free(sum1);
            free(sum2);
#if DEBUG_ON
            printf("\n mloc (in main)\n");
            for (int i=0;i<9;i++){
                for (int j=0;j<9;j++){
                    printf("%f,",elemFemArr.mloc[i][j]);
                }
                printf("\n");
            }
#endif            
            // mloc=mloc+m*txxBEM(kk)*Area(kk)*xw(ii,3)*(mlocTEMP);

        }
        //------------------------------------------------------------->> for each gauss point

#if DEBUG_ON
        printf("\n------------------------------\n"
                "  floc calculations\n"
                "------------------------------\n");
#endif
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

#if DEBUG_ON
        printf("\n");
        for (int i=0;i<9;i++){
            for (int j=0;j<1;j++){
                printf("%f,",elemFemArr.floc[i][j]); 
            }
            printf("\n");
        }

        printf("\n");
        for (int i=0;i<10;i++){
            for (int j=0;j<1;j++){
                printf("%f,",elemFemArr.floc1[i][j]); 
            }
            printf("\n");
        }

        printf("\n");
        for (int i=0;i<10;i++){
            for (int j=0;j<9;j++){
                printf("%f,",elemFemArr.HW[i][j]); 
            }
            printf("\n");
        }
#endif

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

#if DEBUG_ON
        printf("\nMg(:,k)\n");
        for (int j=0;j<81;j++){
            printf("%f,\n",elemFemArr.Mg[j][kk]); 
        }

        printf("\nKg(:,k)\n");
        for (int j=0;j<81;j++){
            printf("%f,\n",elemFemArr.Kg[j][kk]); 
        }
#endif

        if (inDataFem.LL == 2 || inDataFem.LL==3){
            for (int q=0;q<9;q++){
                int cntFglob = wingMeshFem.LM[q][kk]-1;
                elemFemArr.Fglob[cntFglob][0] = elemFemArr.Fglob[cntFglob][0] + elemFemArr.floc[q][0];
                //Fglob(LM(q,kk))=Fglob(LM(q,kk))+floc(q);
            }
        }

        //printf("\nFglob...\n");
        //for (int j=0;j<20;j++){
        //    printf("j=%d, %f,\n",j, elemFemArr.Fglob[j][0]); 
        //}

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

#if DEBUG_ON 
    printf("\nMg(:,k)\n");
    for (int j=0;j<81;j++){
        printf("%f,\n",elemFemArr.Mg[j][0]); 
    }

    printf("\nKg(:,k)\n");
    for (int j=0;j<81;j++){
        printf("%f,\n",elemFemArr.Kg[j][0]); 
    }
#endif

//exit(5);
    //************************************************************************************
    //  DKT PLATE SOLVER: GLOBAL MATRIX ASSEMBLY (Mglob, Kglob, Fglob)
    //************************************************************************************
    float **iii, **rr, **iii_col, **rr_col;
    allocate2Darray(9,9,&iii);
    allocate2Darray(9,9,&rr);
    allocate2Darray(81,1,&iii_col);
    allocate2Darray(81,1,&rr_col);

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

    free(iii);
    free(rr);

    float **Ig, **Jg;
    allocate2Darray(81,wingMeshFem.Nelem,&Ig);//LM(iii(:),:);
    allocate2Darray(81,wingMeshFem.Nelem,&Jg);//LM(rr(:),:);

    int indexIg, indexJg;
    for (int i=0;i<81;i++){
        for (int j=0;j<wingMeshFem.Nelem;j++){
            indexIg = iii_col[i][0];
            indexJg = rr_col[i][0];
            Ig[i][j] = wingMeshFem.LM[indexIg][j]-1;
            Jg[i][j] = wingMeshFem.LM[indexJg][j]-1;
        }
    }   

    float **Kglob, **Mglob;
    allocate2Darray(wingMeshFem.GEN,wingMeshFem.GEN,&Kglob);
    allocate2Darray(wingMeshFem.GEN,wingMeshFem.GEN,&Mglob);

    int indexKi, indexKj;
    for (int i=0;i<81;i++){
        for (int j=0;j<wingMeshFem.Nelem;j++){
            indexKi = Ig[i][j];
            indexKj = Jg[i][j];
            Kglob[indexKi][indexKj] = Kglob[indexKi][indexKj] + elemFemArr.Kg[i][j];
            Mglob[indexKi][indexKj] = Mglob[indexKi][indexKj] + elemFemArr.Mg[i][j];
        }
    }

#if DEBUG_ON
    //printf("Kglob = \n");
    printf("\n\nMglob = \n");
    for (int i=0;i<10;i++){
        for (int j=0;j<10;j++){
            //printf("%f, ",Kglob[i][j]);
            printf("%f, ",Mglob[i][j]);
        }
        printf("\n");
    }

    printf("\n\nKglob = \n");
    for (int i=0;i<10;i++){
        for (int j=0;j<10;j++){
            printf("%f, ",Kglob[i][j]);
        }
        printf("\n");
    }

    printf("iii_cols = \n");
    for (int i=0;i<81;i++){
        printf("%f, ",iii_col[i][0]);
    }

    printf("\n rr_cols = \n");
    for (int i=0;i<81;i++){
        printf("%d, %f, ",i, rr_col[i][0]);
    }
#endif    
    
//exit(10);
    printf("\n    Kglob, Mglob OK... DKT PLATE SOLVER: AUGMENTED GLOBAL MATRIX (for BCs)\n");
    //************************************************************************************
    //  DKT PLATE SOLVER: AUGMENTED GLOBAL MATRIX (for BCs)
    //************************************************************************************
/*
    for (int i=0; i<inDataFem.sizeBdofs;i++){
        printf("%d, ", inDataFem.Bdofs[i]);
    }
*/
    //printf("\n\ninDataFem.sizeBdofs = %d",inDataFem.sizeBdofs );
    //printf("\n\ninDataFem.sizeBdofs = %d",inDataFem.sizeBBnodes );

    
    float **kkk, **mmm;
    allocate2Darray(inDataFem.sizeBdofs,wingMeshFem.GEN,&kkk);
    allocate2Darray(inDataFem.sizeBdofs,wingMeshFem.GEN,&mmm);
    
    //printf("\n\ninDataFem.sizeBdofs = %d",inDataFem.sizeBdofs );
    
    int index_kkk;
    for (int j=0;j<inDataFem.sizeBdofs;j++){
        index_kkk = inDataFem.Bdofs[j]-1;
        kkk[j][index_kkk] = 1.0;
        //kkk(j,:)=[zeros(1,Bdofs(j)-1) 1 zeros(1,Dofs-Bdofs(j))]; %matlab code sample
    }
    
/*
    printf("\nkkk = \n");
    for (int i=0;i<2;i++){
        for (int j=0;j<10;j++){
            printf("%f, ",kkk[i][j]);
        }
        printf("\n");
    }
*/    

    //Kglob=[Kglob kkk'; kkk zeros(length(BBnodes))];  %matlab code sample
    //Mglob=[Mglob mmm'; mmm zeros(length(BBnodes))];  %matlab code sample
    //exit(5);
    float **Kglob_aug, **Mglob_aug; // augmented
    int sizeKMglob_aug = wingMeshFem.GEN+inDataFem.sizeBdofs;
    allocate2Darray(sizeKMglob_aug,sizeKMglob_aug,&Kglob_aug);
    allocate2Darray(sizeKMglob_aug,sizeKMglob_aug,&Mglob_aug);

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

#if DEBUG_ON
    printf("Kglob sample (kkk)\n");
    //for (int i=wingMeshFem.GEN;i<sizeKMglob_aug;i++){
    for (int i=1;i<10;i++){    
        for (int j=0;j<10;j++){
            printf("%f, ",Kglob_aug[i][j]);
        }
        printf("\n");
    }

    printf("Mglob sample (kkk')\n");
    //for (int i=wingMeshFem.GEN;i<sizeKMglob_aug;i++){
    for (int i=1;i<10;i++){    
        for (int j=0;j<10;j++){
            //printf("%f, ",Mglob_aug[j][i]);
            printf("%f, ",Mglob_aug[i][j]);
        }
        printf("\n");
    }
#endif

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

    float **Usol;
    allocate2Darray(sizeKMglob_aug,1, &Usol);

    float **Fglob_aug;
    allocate2Darray(sizeKMglob_aug,1,&Fglob_aug);
    for (int i=0;i<wingMeshFem.GEN;i++){
        Fglob_aug[i][0]=elemFemArr.Fglob[i][0];
    }
    
/*
    for (int i=0;i<sizeKMglob_aug;i++){
        printf("%f, ", Fglob_aug[i][0]);
    }
*/
    if (inDataFem.LL==2){
        printf("    UNIFORM LOAD: P=%f [Pa]\n",inDataFem.P_load);
    }

    // solve linear system of eqs. using LAPACK sgels_ function
    linearSystemSolve(sizeKMglob_aug, sizeKMglob_aug, Kglob_aug, Fglob_aug, Usol);

    //printf("    Usol...\n");
    //for (int i=0;i<20;i++){
    //    printf("%f, ", Usol[i][0]);
    //}

    //sgels_();
    free(Fglob_aug);
    //}

    //************************************************************************************
    //  DKT PLATE SOLVER: OUTPUT BINARY FILE for Matlab Post-Processor
    //************************************************************************************
    //int optionSelect = 0;
    CuFEMNum2DWriteDataInBinary(sizeKMglob_aug, 1, Usol, wingMeshFem.GEN);

    #if (MODAL_ANALYSIS == 1)
    //  TO DO : FIX - I GET INFINITE EIGENVALUES
        printf("\n\nPerforming MODAL ANALYSIS");
        // Generalized Nonsymmetric Eigenvalue Problems
        
        //http://matlab.izmiran.ru/help/techdoc/ref/eig.html
        // Real nonsymmetric A, real general B: sggev() from LAPACK
        float *eigVals;
        allocate1Darray(5,&eigVals);
        modalAnalysis_sggev(sizeKMglob_aug, Kglob_aug, Mglob_aug, eigVals);

        CuFEMNum2DWriteKglobMglobBCs(sizeKMglob_aug, sizeKMglob_aug, Kglob_aug, Mglob_aug);

    #endif


    //************************************************************************************
    //  DKT PLATE SOLVER: CLEAN UP AND EXIT
    //************************************************************************************

    /* TODO free pointers - after using the malloc() */
    free(distrLoad);
    free(distrThick);
    free(inDataFem.pp[0]);
    free(inDataFem.pp[1]);
    free(inDataFem.tt[0]);
    free(inDataFem.tt[1]);
    //free(inDataFem.tt);
    //free(inDataFem.ee);
    //free(wingMeshFem.ID);
    //free(wingMeshFem.IEN);
    free(Kglob);
    free(Mglob);
    free(Kglob_aug);
    free(Mglob_aug);
    free(Usol);


    tend = clock();

    double cpu_time_used = ((double) (tend-tstart))/ CLOCKS_PER_SEC;
    printf("\n----\n Elapsed time [s]: %f\n----\n", cpu_time_used);
    printf("In Matlab the same operations using vectorization take 2.0383 sec.\n");
    return 0;
}


