#include <stdio.h>COMMENT
#include <stdlib.h>

/* 
Search file for TODO
Search file for COMMENT
*/
/* Declarations for data structures and functions */

struct InDataRecFem{
    float cRoot;
    float span;
    float U;
    float densfluid;
    /*material properties*/
    float mass, E, v, h; 
    /*boundary conditions: 1- SS, 2-CC */
    unsigned int CC; 
    /*type of load: 1- concentrated load, 2- uniform load, 3- distributed load*/
    unsigned int LL; 
    /* Delaunay triangulation - unstructured mesh 
            pp (mesh nodes) = [x1,x2,x3,x4...;
                               y1,y2,y3,y4...]; nodal coordinates
    
            ee (edges): 
            ee(1,k) is the index of the first point in mesh edge k
            ee(2,k) is the index of the second point in mesh edge k 
            ee(3,k) is the parameter value of the first point in edge k
            ee(4,k) is the parameter value at the second poiny of edge k
            ee(5,k) is the ID of the geometric edge containing the mesh edge
            ee(6,k),ee(7,k) subdomain numbers
    
            tt (triangles)
    */
    //float *ee; /*2-D array*/
    unsigned int pp_rows; // dummy
    unsigned int pp_cols;
    float *pp[2]; /*2-D array*/
    unsigned int tt_rows; // dummy
    unsigned int tt_cols;
    int *tt[4]; /*2-D array*/
    unsigned int ee_rows; // dummy
    unsigned int ee_cols;
    float *ee[7]; /*2-D array*/
};

struct triangleDKT{
    /*Information about each triangle in the mesh*/
    /*
    BeSt2
    thickness
    forcing 
    */
   int Nelem;
   int NN;
   int GEN;
   int *ID[3];  // [3,NN]
   int *IEN[3]; // [3,Nelem]
   int *LM[9];  // [9,Nelem]
};
//
void CuFEMNum2DReadInData(struct InDataRecFem *inDataFem );
//
void ConnectivityFEM_IEN_ID_LM(struct InDataRecFem *inDataFem, struct triangleDKT *wingMeshFem );

/* The main program follows. */

int main(int argc, char **argv){

    struct InDataRecFem inDataFem;
    /* read input parameters from a file */
    CuFEMNum2DReadInData( &inDataFem );
    /*if the structure is given as a reference*/
    printf("Accessing data structure: %f\n", (&inDataFem)->cRoot);
    /*if the structure is given as a name*/ 
    printf("Accessing data structure: %f\n", inDataFem.cRoot); 

    printf("inDataFem.pp[0][256]: %f\n", inDataFem.pp[0][250]); 
    printf("inDataFem.pp[1][256]: %f\n", inDataFem.pp[1][250]); 

    /* Create or load from matlab IEN, ID, LM */
    struct triangleDKT wingMeshFem;
    ConnectivityFEM_IEN_ID_LM( &inDataFem, &wingMeshFem );

    for (int i=0;i<3;i++)
    {
        for (int j=0;j<6;j++)
        {
            printf("wingMeshFem->IEN[%d][%d]= %d\n,   ", i,j, wingMeshFem.IEN[i][j]);
        }
    }
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

void CuFEMNum2DReadInData(struct InDataRecFem *inDataFem ){
    FILE *file;
	file = fopen("../c/INDATA_FEM.bin", "rb");
    fread(&(inDataFem->cRoot), sizeof(float) , 1, file);
    fread(&(inDataFem->span), sizeof(float) , 1, file);
    fread(&(inDataFem->U), sizeof(float) , 1, file);
    fread(&(inDataFem->densfluid), sizeof(float) , 1, file);
    fread(&(inDataFem->mass), sizeof(float) , 1, file);
    fread(&(inDataFem->E), sizeof(float) , 1, file);
    fread(&(inDataFem->v), sizeof(float) , 1, file);
    fread(&(inDataFem->h), sizeof(float) , 1, file);
    fread(&(inDataFem->CC), sizeof(int) , 1, file);
    fread(&(inDataFem->LL), sizeof(int) , 1, file);
    /* TODO: Create a separate function for reading 2-D matrices*/
    //---------------------------------------------------------------------------->> pp
    fread(&(inDataFem->pp_rows), sizeof(int) , 1, file);
    fread(&(inDataFem->pp_cols), sizeof(int) , 1, file);
    //
    for (int i=0;i<inDataFem->pp_rows;i++){
        inDataFem->pp[i] = (float*)malloc(inDataFem->pp_cols *sizeof(float));
    }
    // Note that arr[i][j] is same as *(*(arr+i)+j)
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < inDataFem->pp_cols; j++){
            fread(&(inDataFem->pp[i][j]), sizeof(float), 1, file);
        }
    }
    //---------------------------------------------------------------------------->> tt 
    fread(&(inDataFem->tt_rows), sizeof(int) , 1, file);
    fread(&(inDataFem->tt_cols), sizeof(int) , 1, file);
    //
    for (int i=0;i<inDataFem->tt_rows;i++){
        inDataFem->tt[i] = (int*)malloc(inDataFem->tt_cols *sizeof(int));
    }
    //
    for (int i = 0; i < inDataFem->tt_rows; i++){
        for (int j = 0; j < inDataFem->tt_cols; j++){
            fread(&(inDataFem->tt[i][j]), sizeof(int), 1, file);
        }
    }
    //---------------------------------------------------------------------------->> ee
    fread(&(inDataFem->ee_rows), sizeof(int) , 1, file);
    fread(&(inDataFem->ee_cols), sizeof(int) , 1, file);
    //
    for (int i=0;i<inDataFem->ee_rows;i++){
        inDataFem->ee[i] = (float*)malloc(inDataFem->ee_cols *sizeof(float));
    }
    //
    for (int i = 0; i < inDataFem->ee_rows; i++){
        for (int j = 0; j < inDataFem->ee_cols; j++){
            fread(&(inDataFem->ee[i][j]), sizeof(float), 1, file);
        }
    }
    //---------------------------------------------------------------------------->>
    fclose(file);

    /* Printing data */
    printf("cRoot = %f\n",inDataFem->cRoot );
	printf("span  = %f\n",inDataFem->span );
	printf("U = %f\n",inDataFem->U );
	printf("densfluid = %f\n",inDataFem->densfluid );
	printf("mass = %f\n",inDataFem->mass );
	printf("E = %f\n",inDataFem->E );
	printf("v = %f\n",inDataFem->v );
    printf("h = %f\n",inDataFem->h );
    printf("CC = %u\n",inDataFem->CC );
    printf("LL = %u\n",inDataFem->LL );
    printf("pp: rows = %u, cols=%u \n",inDataFem->pp_rows, inDataFem->pp_cols );

    for (int i = 0; i < inDataFem->pp_rows; i++){
        for (int j = 0; j < 5; j++){
            printf("inDataFem->pp[%d][%d]= %f,   ", i,j,inDataFem->pp[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    printf("tt: rows = %u, cols=%u \n",inDataFem->tt_rows, inDataFem->tt_cols );

    for (int i = 0; i < inDataFem->tt_rows; i++){
        for (int j = 0; j < 5; j++){
            printf("inDataFem->tt[%d][%d]= %d,   ", i,j,inDataFem->tt[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
    printf("ee: rows = %u, cols=%u \n",inDataFem->ee_rows, inDataFem->ee_cols );

    for (int i = 0; i < inDataFem->ee_rows; i++){
        for (int j = 0; j < 5; j++){
            printf("inDataFem->ee[%d][%d]= %f,   ", i,j,inDataFem->ee[i][j]);
        }
        printf("\n");
    }
}

void ConnectivityFEM_IEN_ID_LM(struct InDataRecFem *inDataFem, struct triangleDKT *wingMeshFem ){
    /* Use the information on Delaunay triangulation from matlab to
    generate connectivity arrays ID, IEN, LM*/

    /*
     COMMENT: I keep the values of ID, IEN, LM the same as matlab, but for the 
    final assembly the number MUST be minus 1. In C matrices start with 0 index!!!!!
    */

    /* IEN = tt(1:3,:) by selecting a triangle of the unstructured mesh
     - pick a column - it gives you the index numbers of the triangle nodes on the 
      global numbering */
    for (int i=0;i<3;i++){
        wingMeshFem->IEN[i] = (int*)malloc(inDataFem->tt_cols *sizeof(int));
    }
    //
    //for (int i = 0; i < inDataFem->tt_rows-1; i++)
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < inDataFem->tt_cols; j++)
        {
            wingMeshFem->IEN[i][j]=inDataFem->tt[i][j];
        }
    }
    wingMeshFem->Nelem=inDataFem->tt_cols;
    printf("\n Nelem=%d\n",wingMeshFem->Nelem);

    /* ID For each node (column), gives you the dof number */
    wingMeshFem->NN=inDataFem->pp_cols;
    printf("\n NN=%d\n",wingMeshFem->NN);

    for (int i=0;i<3;i++){
        wingMeshFem->ID[i] = (int*)malloc(inDataFem->pp_cols *sizeof(int));
    }
    /*
    ID(1,:)=1:3:3*NN-2;
    ID(2,:)=2:3:3*NN-1;
    ID(3,:)=3:3:3*NN;
    */
    int cnt=0;
    for (int j = 1; j < 3*wingMeshFem->NN-2+1; j=j+3){
        wingMeshFem->ID[0][cnt]= j;
        cnt = cnt + 1;
    }
    cnt=0;
    for (int j = 2; j < 3*wingMeshFem->NN-1+1; j=j+3){
        wingMeshFem->ID[1][cnt]= j;
        cnt = cnt + 1;
    }
    cnt=0;
    for (int j = 3; j < 3*wingMeshFem->NN+1; j=j+3){
        wingMeshFem->ID[2][cnt]= j;
        cnt = cnt + 1;
    }
    wingMeshFem->GEN=wingMeshFem->ID[2][wingMeshFem->NN-1]; /*or GEN=max(max(LM)); */
    printf("\n GEN=%d\n",wingMeshFem->GEN);

    for (int i = 0; i < 3; i++){
        for (int j = 250; j < 253; j++){
            printf("ID[%d][%d]= %d,   ", i,j,wingMeshFem->ID[i][j]);
        }
        printf("\n");
    }
    printf("ID[%d][%d]= %d,   \n", 0,wingMeshFem->NN -1 ,wingMeshFem->ID[0][wingMeshFem->NN -1]);

    /*
    LM=zeros(9,Nelem);
    for k=1:Nelem
       for i=1:3
           for j=1:3
               P=(3)*(j-1)+i; %the 9 dofs per triangle
               LM(P,k)=ID(i,IEN(j,k));
           end
       end
    end
    */
    for (int i=0;i<9;i++){
        wingMeshFem->LM[i] = (int*)malloc(wingMeshFem->Nelem *sizeof(int));
    }
    int Pindex;
    int i,j,k; // for each triangle
    for (k = 0; k< wingMeshFem->Nelem; k++){
        for (i=0; i<3; i++){
            for (j=0; j<3; j++){
                Pindex=(3)*(j-1)+i; //the 9 dofs per triangle
                wingMeshFem->LM[Pindex][k]=wingMeshFem->ID[i][wingMeshFem->IEN[j][k]];
            }
        }
    }

}

