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
    struct triangleDKT wingMeshFem;
    int Ng = 6; // Gauss integration points

    // Preparing to run a script for the purpose of scatter interpolation //
    char command[] = "python3 ";
    char scriptName[] = "dataExchange_BEM_FEM.py";
    strcat(command, scriptName);

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

    /* Bending stiffness for each triangle - using python??? or constant thickness?? */
    system(command);

//#if DEBUG
    for (int i=0;i<Ng;i++){
        for (int j=0;j<3;j++){
            printf("xw [%d]:%f,",j,xw[i][j]);
        }
        printf("\n");
    }
//#endif
    

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


