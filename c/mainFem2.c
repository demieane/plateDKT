/*This form is used for header files of your own program.
It searches for a file named 'file' in the directory containing the current file.
You can prepend directories to this list with the -I 
option while compiling your source code. */
#include "mainFem.h"
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

    clock_t tstart, tend;

    tstart = clock();
    /* Do the work. */

    struct InDataRecFem inDataFem;
    struct triangleDKT wingMeshFem;
    int Ng = 3; // Gauss integration points

    // Preparing to run a script for the purpose of scatter interpolation //
    char command[] = "python3 ";
    char scriptName[] = "dataExchange_BEM_FEM.py";
    // concatenates str1 and str2
    // the resultant string is stored in str1.
    // strcat(str1, str2);
    strcat(command, scriptName);

    //************************************************************************************
    //  DKT PLATE SOLVER: PREPARATIONS
    //************************************************************************************
    /* read input parameters from a file */
    /* Boundary conditions nodes, dofs */
    CuFEMNum2DReadInData( &inDataFem );

    /* TODO: Create output function to check whether the BCs are correct in a figure */

    /*if the structure is given as a reference*/
    //printf("Accessing data structure: %f\n", (&inDataFem)->cRoot);
    /*if the structure is given as a name*/ 
    //printf("Accessing data structure: %f\n", inDataFem.cRoot); 

    //printf("inDataFem.pp[0][256]: %f\n", inDataFem.pp[0][250]); 
    //printf("inDataFem.pp[1][256]: %f\n", inDataFem.pp[1][250]); 

    /* Create or load from matlab IEN, ID, LM */
    ConnectivityFEM_IEN_ID_LM( &inDataFem, &wingMeshFem );

    /* Gauss integration function - read about it */
    float xw[Ng][3]; // {xg,yg,wg}
    TriGaussPoints(Ng, xw);

#if DEBUG_ON
    for (int i=0;i<Ng;i++){
        for (int j=0;j<3;j++){
            printf("MAIN: xw [%d]:%f,",j,xw[i][j]);
        }
        printf("\n");
    }
#endif

    /* Bending stiffness for each triangle - using python??? or constant thickness?? */
    system(command);

    printf("BeSt");
    //float BeSt[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    float **BeSt;
    allocate2Darray(3, 3, &BeSt);
    BendingStiffness(inDataFem.E, inDataFem.v, inDataFem.h, BeSt);
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            printf("%f, ", BeSt[i][j]);
        }
        printf("\n");
    }

    /* DKT */
    TrigElCoefsDKT(&inDataFem, &wingMeshFem);


#if DEBUG_ON
    printf("l23 1: %f, Nelem: %f\n", wingMeshFem.l23[0],wingMeshFem.l23[wingMeshFem.Nelem-1]);
    printf("l31 1: %f, Nelem: %f\n", wingMeshFem.l31[0],wingMeshFem.l31[wingMeshFem.Nelem-1]);
    printf("l12 1: %f, Nelem: %f\n", wingMeshFem.l12[0],wingMeshFem.l12[wingMeshFem.Nelem-1]);
    //
    printf("y23 1: %f, Nelem: %f\n", wingMeshFem.y23[0],wingMeshFem.y23[wingMeshFem.Nelem-1]);
    printf("y31 1: %f, Nelem: %f\n", wingMeshFem.y31[0],wingMeshFem.y31[wingMeshFem.Nelem-1]);
    printf("y12 1: %f, Nelem: %f\n", wingMeshFem.y12[0],wingMeshFem.y12[wingMeshFem.Nelem-1]);
    //
    printf("x23 1: %f, Nelem: %f\n", wingMeshFem.x23[0],wingMeshFem.x23[wingMeshFem.Nelem-1]);
    printf("x31 1: %f, Nelem: %f\n", wingMeshFem.x31[0],wingMeshFem.x31[wingMeshFem.Nelem-1]);
    printf("x12 1: %f, Nelem: %f\n", wingMeshFem.x12[0],wingMeshFem.x12[wingMeshFem.Nelem-1]);
    //
    printf("C4 1: %f, Nelem: %f\n", wingMeshFem.C4[0],wingMeshFem.C4[wingMeshFem.Nelem-1]);
    printf("C5 1: %f, Nelem: %f\n", wingMeshFem.C5[0],wingMeshFem.C5[wingMeshFem.Nelem-1]);
    printf("C6 1: %f, Nelem: %f\n", wingMeshFem.C6[0],wingMeshFem.C6[wingMeshFem.Nelem-1]);
    //
    printf("c4 1: %f, Nelem: %f\n", wingMeshFem.c4[0],wingMeshFem.c4[wingMeshFem.Nelem-1]);
    printf("c5 1: %f, Nelem: %f\n", wingMeshFem.c5[0],wingMeshFem.c5[wingMeshFem.Nelem-1]);
    printf("c6 1: %f, Nelem: %f\n", wingMeshFem.c6[0],wingMeshFem.c6[wingMeshFem.Nelem-1]);
#endif

    LNShapeFunDST(Ng, xw, &wingMeshFem);

    LNShapeFunMassDST(Ng, xw, &wingMeshFem);

    matrixG(&wingMeshFem);


    printf("Testing BLAS\n\n");

    // EXAMPLE FROM cblas.h
    float a[2] = {-1, -1};
    float summ;
    summ = cblas_sasum(2, a, 1);
    printf("cblas_sasum: %f",summ);

    //int rows, cols;
    squareMatInverse2(10, 10, wingMeshFem.GGDST, wingMeshFem.GGin);
    squareMatInverse2(6, 6, wingMeshFem.GGDKT, wingMeshFem.GGin2);


#if DEBUG_ON
    for (int i=0;i<10;i++){
        for (int j=0;j<10;j++){
            printf("%f, ",wingMeshFem.GGDST[i][j]);
        }
        printf("\n\n");
    }

    for (int i=0;i<10;i++){
        for (int j=0;j<10;j++){
            printf("%f, ",wingMeshFem.GGin[i][j]);
        }
        printf("\n");
    }

    for (int i=0;i<6;i++){
        for (int j=0;j<6;j++){
            printf("%f, ",wingMeshFem.GGDKT[i][j]);
        }
        printf("\n\n");
    }

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
    printf("%f =>  P_load in [Pa]", inDataFem.P_load);

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
    //for (int kk = 0;kk<2;kk++){   
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
            // floc1=floc1+Area(kk)*xw(ii,3)*(LW');
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
            float var0 = pow(inDataFem.h,2)/12.0;
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
            for (int i=0;i<9;i++){
                for (int j=0;j<1;j++){
                    elemFemArr.floc[i][j] = elemFemArr.floc[i][j] + wingMeshFem.area[kk]*inDataFem.P_load/3.0*lumpedMass[i];
                }
            }
        }
#if DEBUG_ON
        for (int i=0;i<9;i++){
            for (int j=0;j<1;j++){
                printf("%f,",elemFemArr.floc[i][j]); 
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

        for (int q=0;q<9;q++){
            int cntFglob = wingMeshFem.LM[q][kk]-1;
            elemFemArr.Fglob[cntFglob][0] = elemFemArr.Fglob[cntFglob][0] + elemFemArr.floc[q][0];
            //Fglob(LM(q,kk))=Fglob(LM(q,kk))+floc(q);
        }

        //printf("\nFglob...\n");
        //for (int j=0;j<20;j++){
        //    printf("j=%d, %f,\n",j, elemFemArr.Fglob[j][0]); 
        //}

        // re-initialize mloc, kloc, floc
        for (int i = 0;i<9;i++){
            elemFemArr.floc[i][0] = 0;
            for (int j = 0;j<9;j++){
                elemFemArr.mloc[i][j] = 0;
                elemFemArr.kloc[i][j] = 0;
            }
        }
   
    }
    printf("\nCalculated Mg(:,kk), Kg(:,kk), Fglob(kk)");

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

    //************************************************************************************
    //  DKT PLATE SOLVER: GLOBAL MATRIX ASSEMBLY (Mglob, Kglob, Fglob)
    //************************************************************************************

    


    //************************************************************************************
    //  DKT PLATE SOLVER: AUGMENTED GLOBAL MATRIX (for BCs)
    //************************************************************************************


    //************************************************************************************
    //  DKT PLATE SOLVER: SOLUTION OPTIONS (1. EIGEN, 2. STATIC, 3. DYNAMIC)
    //************************************************************************************

    //************************************************************************************
    //  DKT PLATE SOLVER: OUTPUT BINARY FILE for Matlab Post-Processor
    //************************************************************************************

    //************************************************************************************
    //  DKT PLATE SOLVER: CLEAN UP AND EXIT
    //************************************************************************************

    /* TODO free pointers - after using the malloc() */
    free(inDataFem.pp[0]);
    free(inDataFem.pp[1]);
    free(inDataFem.tt[0]);
    free(inDataFem.tt[1]);
    //free(inDataFem.tt);
    //free(inDataFem.ee);
    //free(wingMeshFem.ID);
    //free(wingMeshFem.IEN);

    tend = clock();

    double cpu_time_used = ((double) (tend-tstart))/ CLOCKS_PER_SEC;
    printf("\n----\n Elapsed time [s]: %f\n----\n", cpu_time_used);
    printf("In Matlab the same operations using vectorization take 0.2493 sec.\n");
    return 0;
}


