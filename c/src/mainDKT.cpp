/*
 COMMENT: Make sure that the INDATA_FEM.bin is of the same precision type 
 as the preprocessor directives used in the code. Else, it gives instant
 segmentation fault. 

*/

#include<stdio.h>
#include<stdlib.h> //malloc
#include<time.h>

#ifndef DEBUG_ON
    #define DEBUG_ON 0 /*allow printf for debugging purposes*/
#endif

#ifndef PRECISION_MODE_FEM
    #define PRECISION_MODE_FEM 2 /* 1. DOUBLE, 2. SINGLE */
    #if PRECISION_MODE_FEM == 1
        typedef double mytype;
    #endif
    #if PRECISION_MODE_FEM == 2
        typedef float mytype;
    #endif
#endif



#define Ng 3 /* Gauss Integration Points */

/*=========================================================================================*/
/* Declarations for data structures and functions */
/*=========================================================================================*/
template<class T>
struct InDataRecFem{
    /* UNIT SYSTEM SI */
    T cRoot;
    T span;
    T U;
    T densfluid;
    /*material properties (mass: kg/m3, Young's modulus: Pa, Poisson's ratio, (average)thickness: m)*/
    T mass, E, v, h; 
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
    unsigned int pp_rows; // dummy
    unsigned int pp_cols;
    T *pp[2]; /*2-D array*/
    unsigned int tt_rows; // dummy
    unsigned int tt_cols;
    int *tt[4]; /*2-D array*/
    unsigned int ee_rows; // dummy
    unsigned int ee_cols;
    T *ee[7]; /*2-D array*/
    unsigned int sizeBBnodes;
    int *BBnodes; /* id of nodes affected by boundary conditions*/
    unsigned int sizeBdofs;
    int *Bdofs; /* global numbering of dofs affected by boundary conditions*/

    T P_load; // Pa: positive values point towards the positive Z-axis (reverse for ANSYS) 
    // ONLY FOR CONCENTRATED LOAD <--BELOW
    T P_xy[2];
    int P_node; 
    // ONLY FOR DISTRIBUTED PROPERTIES LOAD/THICKNESS <--BELOW
    int sizexcp;
    T *xcp, *ycp, *fcp, *tcp;
};

template<class T>
struct triangleDKT{
    /*Information about each triangle in the mesh*/
    /*
    BeSt2
    thickness
    forcing 
    xm, ym : center of triangles
    */
    int Nelem; // number of triangles
    int NN; // number of nodes
    int GEN; // number of dofs (system of eqs. before BCs)
    int *ID[3];  // [3,NN]
    int *IEN[3]; // [3,Nelem]
    int *LM[9];  // [9,Nelem]
    //
    T *xm; // x barycentric coordinate [Nelem]
    T *ym; // y barycentric coordinate [Nelem]

    /* output of TrigElCoefsDKT() */ 
    T *l23, *l31, *l12; // [1 x Nelem]
    T *y12, *y31, *y23;
    T *x12, *x31, *x23;
    T *area;
    T *a4, *a5, *a6, *b4, *b5, *b6, *c4, *c5, *c6;
    T *d4, *d5, *d6, *e4, *e5, *e6;
    T *C4, *C5, *C6, *S4, *S5, *S6;
    /*
    The C language is case-sensitive. This means that all language keywords,
    identifiers, function names, and other variables
    must be entered with consistent letter capitalization. 
    */
    /* output of LNShapeFunDST() */
    T **SF, **DxsiSF, **DetaSF; //[Ng x 6]
    T *D2xsiSF, *D2xsietaSF, *D2etaSF; // [1 x 6]

    /* output of LNShapeFunMassDST() */
    T **SFm, **DxsiSFm, **DetaSFm; //[Ng x 3]

    T **GGDST, **GGDKT; //[10 x 10]  
    T **GGin, **GGin2; //inverse of above
};

/*=========================================================================================*/
/* Functions prototypes*/
/*=========================================================================================*/
template<class T>
void CuFEMNum2DReadInData(struct InDataRecFem<T> *inDataFem );
//
template<class T>
void ConnectivityFEM_IEN_ID_LM(struct InDataRecFem<T> *inDataFem, struct triangleDKT<T> *wingMeshFem );


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

    



    
    tend = clock();

    double cpu_time_used = ((double) (tend-tstart))/ CLOCKS_PER_SEC;
    printf("\n----\n Elapsed time [s]: %f\n----\n", cpu_time_used);
    printf("In Matlab the same operations using vectorization take 2.0383 sec.\n");
    return 0;

}





/*=========================================================================================*/
/* Definition of the functions BELOW */
/*=========================================================================================*/
template<class T>
void CuFEMNum2DReadInData(struct InDataRecFem<T> *inDataFem ){
    printf("\n    Entering CuFEMNum2DReadInData().\n");
    FILE *file;
	file = fopen("../INDATA_FEM.bin", "rb"); // r for read, b for binary
    fread(&(inDataFem->cRoot), sizeof(T) , 1, file);
    fread(&(inDataFem->span), sizeof(T) , 1, file);
    fread(&(inDataFem->U), sizeof(T) , 1, file);
    fread(&(inDataFem->densfluid), sizeof(T) , 1, file);
    fread(&(inDataFem->mass), sizeof(T) , 1, file);
    fread(&(inDataFem->E), sizeof(T) , 1, file);
    fread(&(inDataFem->v), sizeof(T) , 1, file);
    fread(&(inDataFem->h), sizeof(T) , 1, file);
    fread(&(inDataFem->CC), sizeof(int) , 1, file);
    fread(&(inDataFem->LL), sizeof(int) , 1, file);
    /* TODO: Create a separate function for reading 2-D matrices*/
    //---------------------------------------------------------------------------->> pp
    fread(&(inDataFem->pp_rows), sizeof(int) , 1, file);
    fread(&(inDataFem->pp_cols), sizeof(int) , 1, file);
    //
    for (int i=0;i<inDataFem->pp_rows;i++){
        inDataFem->pp[i] = (T*)malloc(inDataFem->pp_cols *sizeof(T));
    }
    // Note that arr[i][j] is same as *(*(arr+i)+j)
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < inDataFem->pp_cols; j++){
            fread(&(inDataFem->pp[i][j]), sizeof(T), 1, file);
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
        inDataFem->ee[i] = (T*)malloc(inDataFem->ee_cols *sizeof(T));
    }
    //
    for (int i = 0; i < inDataFem->ee_rows; i++){
        for (int j = 0; j < inDataFem->ee_cols; j++){
            fread(&(inDataFem->ee[i][j]), sizeof(T), 1, file);
        }
    }
    //---------------------------------------------------------------------------->>
    fread(&(inDataFem->sizeBBnodes), sizeof(int) , 1, file);
    //printf("sizeBBnodes = %d\n",inDataFem->sizeBBnodes );
    inDataFem->BBnodes = (int*)malloc(inDataFem->sizeBBnodes *sizeof(int));
    for (int i=0;i<inDataFem->sizeBBnodes;i++){
        fread(&(inDataFem->BBnodes[i]), sizeof(int), 1, file);
        //printf("i=%d,BBnodes[i]=%d\n", i,inDataFem->BBnodes[i]);
    }
    
    fread(&(inDataFem->sizeBdofs), sizeof(int) , 1, file);
    inDataFem->Bdofs = (int*)malloc(inDataFem->sizeBdofs *sizeof(int));
    for (int i=0;i<inDataFem->sizeBdofs;i++){
        fread(&(inDataFem->Bdofs[i]), sizeof(int), 1, file);
        //printf("i=%d,Bdofs[i]=%d\n", i,inDataFem->Bdofs[i]);
    }
    //printf("BBnodes = %d, Bdofs=%d\n",inDataFem->sizeBBnodes, inDataFem->sizeBdofs );
    if (inDataFem->LL == 2){
        printf("    UNIFORM LOAD CASE LL == 2\n");
        /* UNIFORM LOAD CASE */
        fread(&(inDataFem->P_load), sizeof(T) , 1, file);
    }
    if (inDataFem->LL == 1){
        printf("    POINT LOAD CASE LL == 1\n");
        /* POINT LOAD CASE */
        fread(&(inDataFem->P_load), sizeof(T) , 1, file);
        fread(&(inDataFem->P_xy[0]), sizeof(T) , 1, file);
        fread(&(inDataFem->P_xy[1]), sizeof(T) , 1, file);
        printf("    Px=%f, Py=%f\n\n",inDataFem->P_xy[0],inDataFem->P_xy[1]);
        fread(&(inDataFem->P_node), sizeof(int) , 1, file);
        printf("    P_NODE=%d \n\n",inDataFem->P_node);
    }
    if (inDataFem->LL == 3){
        printf("    DISTRIBUTED LOAD CASE LL == 3\n");
        /* DISTRIBUTED LOAD CASE */ /* TODO: use pre-processor directives instead*/
        fread(&(inDataFem->sizexcp), sizeof(int), 1, file);
        printf("sizexcp = %d\n", inDataFem->sizexcp);
        inDataFem->xcp = (T*)malloc(inDataFem->sizexcp *sizeof(T));
        inDataFem->ycp = (T*)malloc(inDataFem->sizexcp *sizeof(T));
        inDataFem->fcp = (T*)malloc(inDataFem->sizexcp *sizeof(T));
        inDataFem->tcp = (T*)malloc(inDataFem->sizexcp *sizeof(T));
        for (int i = 0;i<inDataFem->sizexcp;i++){
            fread(&(inDataFem->xcp[i]),sizeof(T),1,file);
        }
        for (int i = 0;i<inDataFem->sizexcp;i++){
            fread(&(inDataFem->ycp[i]),sizeof(T),1,file);
        }
        for (int i = 0;i<inDataFem->sizexcp;i++){
            fread(&(inDataFem->fcp[i]),sizeof(T),1,file);
        }
        for (int i = 0;i<inDataFem->sizexcp;i++){
            fread(&(inDataFem->tcp[i]),sizeof(T),1,file);
        }
    }
    fclose(file);

#if DEBUG_ON
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
#endif
        printf("    Exiting CuFEMNum2DReadInData().  OK.\n");
}

template<class T>
void ConnectivityFEM_IEN_ID_LM(struct InDataRecFem<T> *inDataFem, struct triangleDKT<T> *wingMeshFem ){
    printf("\n    Exiting ConnectivityFEM_IEN_ID_LM().  \n");
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
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < inDataFem->tt_cols; j++){
            wingMeshFem->IEN[i][j]=inDataFem->tt[i][j];
        }
    }
    wingMeshFem->Nelem=inDataFem->tt_cols;
    /* ID For each node (column), gives you the dof number */
    wingMeshFem->NN=inDataFem->pp_cols;

    /*
    ID(1,:)=1:3:3*NN-2;
    ID(2,:)=2:3:3*NN-1;
    ID(3,:)=3:3:3*NN;
    */
    for (int i=0;i<3;i++){
        wingMeshFem->ID[i] = (int*)malloc(inDataFem->pp_cols *sizeof(int));
    }
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
    int aa,bb;
    for (k = 0; k< wingMeshFem->Nelem; k++){
        for (i=0; i<3; i++){
            for (j=0; j<3; j++){
                Pindex=(3*j+i); //the 9 dofs per triangle
                aa = i;
                bb = wingMeshFem->IEN[j][k]-1; //SINCE IT COMES FROM MATLAB NUMBERING
                wingMeshFem->LM[Pindex][k]= wingMeshFem->ID[aa][bb];
            }
        }
    }

#if DEBUG_ON_ON
    for (int i=0;i<3;i++){
    //for (int j=0;j<5;j++){
        int j = wingMeshFem->Nelem-1;
        //printf("wingMeshFem->IEN[%d][%d]= %d\n,   ", i,j, wingMeshFem->IEN[i][j]);
        printf("wingMeshFem->IEN[%d][%d]= %d\n,   ", i,j, wingMeshFem->IEN[i][j]);
    //}
    }

    printf("\n NN=%d\n",wingMeshFem->NN);
    printf("\n Nelem=%d\n",wingMeshFem->Nelem);
    printf("\n GEN=%d\n",wingMeshFem->GEN);

    for (int i = 0; i < 3; i++){
        for (int j = 250; j < 253; j++){
            printf("ID[%d][%d]= %d,   ", i,j,wingMeshFem->ID[i][j]);
        }
        printf("\n");
    }
    printf("ID[%d][%d]= %d,   \n", 0,wingMeshFem->NN -1 ,wingMeshFem->ID[0][wingMeshFem->NN -1]);

    for (i=0; i<9; i++){
        int j =wingMeshFem->Nelem-1;//for (j=0; j<3; j++){
            printf("LM[%d][%d]= %d\n\n", i,j,wingMeshFem->LM[i][j]);
        //}
    }
#endif

    printf("    Calculating barycentric coordinates.\n");
    // allocate memory
    wingMeshFem->xm = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->ym = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    
    int IEN_1, IEN_2, IEN_3;    
    //i = wingMeshFem->Nelem-1;
    for (int i=0;i<wingMeshFem->Nelem;i++){ 
        IEN_1 = wingMeshFem->IEN[0][i] - 1;
        IEN_2 = wingMeshFem->IEN[1][i] - 1;
        IEN_3 = wingMeshFem->IEN[2][i] - 1;
        //printf("%d, %d, %d, \n",IEN_1, IEN_2, IEN_3);
        wingMeshFem->xm[i] = (1.0/3.0)*(inDataFem->pp[0][IEN_1] + inDataFem->pp[0][IEN_2] + inDataFem->pp[0][IEN_3]);
        wingMeshFem->ym[i] = (1.0/3.0)*(inDataFem->pp[1][IEN_1] + inDataFem->pp[1][IEN_2] + inDataFem->pp[1][IEN_3]);
        //printf("xm(%d)=%f, ym(%d)=%f\n",i, wingMeshFem->xm[i],i, wingMeshFem->ym[i]);
    }
    
    printf("    Exiting ConnectivityFEM_IEN_ID_LM. OK.\n");
}