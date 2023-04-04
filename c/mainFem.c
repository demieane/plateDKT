#include"mainFem.h"
/* 
Search file for TODO
Search file for COMMENT
*/

/*=========================================================================================*/
/* The main program follows. */
/*=========================================================================================*/
int main(int argc, char **argv){

    struct InDataRecFem inDataFem;
    /* read input parameters from a file */
    CuFEMNum2DReadInData( &inDataFem );
    /*if the structure is given as a reference*/
    //printf("Accessing data structure: %f\n", (&inDataFem)->cRoot);
    /*if the structure is given as a name*/ 
    //printf("Accessing data structure: %f\n", inDataFem.cRoot); 

    //printf("inDataFem.pp[0][256]: %f\n", inDataFem.pp[0][250]); 
    //printf("inDataFem.pp[1][256]: %f\n", inDataFem.pp[1][250]); 

    /* Create or load from matlab IEN, ID, LM */
    struct triangleDKT wingMeshFem;
    ConnectivityFEM_IEN_ID_LM( &inDataFem, &wingMeshFem );

#if DEBUG
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
    /* Boundary conditions nodes */

    /* Create output function to check whether the BCs are correct in a figure */

    /* Gauss integration function - read about it */

    /* Bending stiffness for each triangle - using python??? or constant thickness?? */


    /* TODO free pointers - after using the malloc() */
    free(inDataFem.pp[0]);
    free(inDataFem.pp[1]);
    free(inDataFem.tt[0]);
    free(inDataFem.tt[1]);
    //free(inDataFem.tt);
    //free(inDataFem.ee);
    //free(wingMeshFem.ID);
    //free(wingMeshFem.IEN);

    return 0;
}


