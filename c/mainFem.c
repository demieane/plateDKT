#include <stdio.h>
#include <stdlib.h>

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
    unsigned int pp_rows;
    unsigned int pp_cols;
    float *pp[2]; /*2-D array*/
    //float *tt; /*2-D array*/
};

struct triangleDKT{
    /*Information about each triangle in the mesh*/
    /*
    BeSt2
    thickness
    forcing 
    */
};

void CuFEMNum2DReadInData(struct InDataRecFem *inDataFem );

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

    /* Boundary conditions nodes */

    /* Gauss integration function */

    /* Bending stiffness for each triangle */

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
    //
    fread(&(inDataFem->pp_rows), sizeof(int) , 1, file);
    fread(&(inDataFem->pp_cols), sizeof(int) , 1, file);

    //inDataFem->pp = (float*)malloc(inDataFem->pp_rows * sizeof(float));
    for (int i=0;i<2;i++){
        inDataFem->pp[i] = (float*)malloc(inDataFem->pp_cols *sizeof(float));
    }

    // Note that arr[i][j] is same as *(*(arr+i)+j)
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < inDataFem->pp_cols; j++){
            fread(&(inDataFem->pp[i][j]), sizeof(float), 1, file);
        }
    }

    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 10; j++){
            printf("inDataFem->pp[%d][%d]= %f,   ", i,j,inDataFem->pp[i][j]);
        }
        printf("\n");
    }

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

}
