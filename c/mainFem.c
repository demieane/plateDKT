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
    float BeSt[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    BendingStiffness(inDataFem.E, inDataFem.v, inDataFem.h, BeSt);
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            printf("%f\n", BeSt[i][j]);
        }
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
    allocate2Darray(10, 1, &(elemFemArr.floc));
    //
    allocate2Darray(6, 9, &(elemFemArr.Hxx));
    allocate2Darray(6, 9, &(elemFemArr.Hyy));
    //
    allocate1Darray(9, &(elemFemArr.Hx)); // [1 x 9]

    /* NEW */
    allocate1Darray(9, &(elemFemArr.Hx_xsi)); // [1 x 9]
    allocate1Darray(9, &(elemFemArr.Hx_eta)); // [1 x 9]
    allocate1Darray(9, &(elemFemArr.Hy_xsi)); // [1 x 9]
    allocate1Darray(9, &(elemFemArr.Hy_eta)); // [1 x 9]

    allocate1Darray(9, &(elemFemArr.Hy)); 
    allocate1Darray(10, &(elemFemArr.LW)); 
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
    //for (int kk = 0;kk<wingMeshFem.Nelem;kk++){
    for (int kk = 0;kk<1;kk++){   

        massHmDKT(kk, &wingMeshFem, &elemFemArr); // Hm, HW
        rotationMass2(kk, &wingMeshFem, &elemFemArr); // Hxx, Hyy

//#if DEBUG_ON
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
//#endif
        // for each gauss point
        for (int ii = 0; ii<1; ii++){
        //for (int ii = 0; ii<Ng; ii++){
            ShapeFunDKT2(ii, kk, &wingMeshFem, &elemFemArr);
            pseudoMassDKT(ii, kk, &wingMeshFem, &elemFemArr); // not exactly used (only LW)

            // TODO: Use the available CBLAS & LAPACKE routines 
            
            //matrix addition needed 
            // kb=kb+Area(kk)*xw(ii,3)*(Bb'*BeSt2(:,:,kk)*Bb);

            // kb_temp = Bb'*BeSt2(:,:,kk)
            //matMatMultiplication2(9, 3, 3, float **arrA, float **arrB, float **arrOut)

        }


    }

    

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


