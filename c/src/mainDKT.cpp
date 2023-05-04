/*
 COMMENT: Make sure that the INDATA_FEM.bin is of the same precision type 
 as the preprocessor directives used in the code. Else, it gives instant
 segmentation fault. 

 COMPLILE using:
 >> gcc -g mainDKT.cpp -lm

*/

#include<stdio.h>
#include<stdlib.h> //malloc
#include<time.h>
#include<math.h>

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

#define GaussIntegrPoints 3 /* Gauss Integration Points */

/*=========================================================================================*/
/* Declarations for data structures and functions */
/*=========================================================================================*/
template<class T>
struct InDataRecFem{
    /* UNIT SYSTEM SI */
    int modeFem;
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
void freeInDataRecFem(struct InDataRecFem<T> *inDataFem);
//
template<class T>
void ConnectivityFEM_IEN_ID_LM(struct InDataRecFem<T> *inDataFem, struct triangleDKT<T> *wingMeshFem );
//
template<class T>
void freetriangleDKT(int Ng, struct triangleDKT<T> *wingMeshFem);
//
template<class T>
void TriGaussPoints(T xw[GaussIntegrPoints][3]);
//
template<class T>
void BendingStiffness(float E, float v, float tx, T **BeSt);
//
template<class T>
void allocate1Darray(int rows, T **arrIn);
//
template<class T>
void allocate2Darray(int rows, int cols, T ***arrIn);
//
template<class T>
void shepard_interp_2d(int nd, T *xd, T *yd, T *zd,
    T *p, int ni, T *xi, T *yi, T *zi);

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

    /* Gauss integration function - read about it */
    mytype xw[GaussIntegrPoints][3]; // {xg,yg,wg}
    TriGaussPoints(xw);

    /* Distributed properties */
    mytype *distrLoad, *distrThick;
    allocate1Darray(wingMeshFem.Nelem,&distrLoad);
    allocate1Darray(wingMeshFem.Nelem,&distrThick);

    if (inDataFem.LL == 3){
        /* DISTRIBUTED LOAD & THICKNESS CASE */
        int nd = inDataFem.sizexcp;
        int ni = wingMeshFem.Nelem; //size(xm)

        mytype p1 = 10.55;
        mytype p2 = 10.55;
        mytype *pparam1, *pparam2;
        pparam1 = &p1;  
        pparam2 = &p2; 

        shepard_interp_2d(nd, inDataFem.xcp, inDataFem.ycp, inDataFem.fcp, 
        pparam1, ni, wingMeshFem.xm, wingMeshFem.ym, distrLoad);

        shepard_interp_2d(nd, inDataFem.xcp, inDataFem.ycp, inDataFem.tcp, 
        pparam2, ni, wingMeshFem.xm, wingMeshFem.ym, distrThick);

        printf("distrThick[i], distrLoad[i]\n");
        for (int i = 0; i<10; i++){
            printf("%f, %f\n", distrThick[i], distrLoad[i]);
        }
    }
    float **BeSt;
    allocate2Darray(3, 3, &BeSt);
    if (inDataFem.LL==2 || inDataFem.LL==1){
        BendingStiffness(inDataFem.E, inDataFem.v, inDataFem.h, BeSt);
    }


    



    
    tend = clock();

    double cpu_time_used = ((double) (tend-tstart))/ CLOCKS_PER_SEC;
    printf("\n----\n Elapsed time [s]: %f\n----\n", cpu_time_used);
    printf("In Matlab the same operations using vectorization take 2.0383 sec.\n");


    freeInDataRecFem(&inDataFem);
    int Ng = GaussIntegrPoints;
    //freetriangleDKT(Ng,&wingMeshFem);

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
    fread(&(inDataFem->modeFem), sizeof(int) , 1, file);
    if ((inDataFem->modeFem != PRECISION_MODE_FEM)){
        printf("    Compile code with correct precision mode. Enjoy the seg fault :) \n");   
    }
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
        printf("    sizexcp = %d\n", inDataFem->sizexcp);
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
void freeInDataRecFem(struct InDataRecFem<T> *inDataFem){
/*
    float *pp[2]; 
    int *tt[4]; 
    float *ee[7]; 
    int *BBnodes; 
    int *Bdofs; 
    float *xcp, *ycp, *fcp, *tcp;
*/
    for (int i=0;i<inDataFem->pp_rows;i++){
        free(inDataFem->pp[i]);
    }
    for (int i=0;i<inDataFem->tt_rows;i++){
        free(inDataFem->tt[i]);
    }
    for (int i=0;i<inDataFem->ee_rows;i++){
        free(inDataFem->ee[i]);
    }
    if (inDataFem->LL == 3){
        free(inDataFem->xcp);
        free(inDataFem->ycp);
        free(inDataFem->fcp);
        free(inDataFem->tcp);
    }
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

template<class T>
void freetriangleDKT(int Ng, struct triangleDKT<T> *wingMeshFem){
/*
struct triangleDKT  
{
    int Nelem; // number of triangles
    int NN; // number of nodes
    int GEN; // number of dofs (system of eqs. before BCs)
    int *ID[3];  // [3,NN]
    int *IEN[3]; // [3,Nelem]
    int *LM[9];  // [9,Nelem]
    float *xm; // x barycentric coordinate [Nelem]
    float *ym; // y barycentric coordinate [Nelem]
    float *l23, *l31, *l12; // [1 x Nelem]
    float *y12, *y31, *y23;
    float *x12, *x31, *x23;
    float *area;
    float *a4, *a5, *a6, *b4, *b5, *b6, *c4, *c5, *c6;
    float *d4, *d5, *d6, *e4, *e5, *e6;
    float *C4, *C5, *C6, *S4, *S5, *S6;
    float **SF, **DxsiSF, **DetaSF; //[Ng x 6]
    float *D2xsiSF, *D2xsietaSF, *D2etaSF; // [1 x 6]
    float **SFm, **DxsiSFm, **DetaSFm; //[Ng x 3]
    float **GGDST, **GGDKT; //[10 x 10]  
    float **GGin, **GGin2; //inverse of above
}
*/
    for (int i=0;i<3;i++){
        free(wingMeshFem->ID[i]);
        free(wingMeshFem->IEN[i]);
    }
    for (int i=0;i<9;i++){
        free(wingMeshFem->LM[i]);
    }
    free(wingMeshFem->xm);
    free(wingMeshFem->ym);
    //    
    free(wingMeshFem->l12);
    free(wingMeshFem->l23);
    free(wingMeshFem->l31);
    //
    free(wingMeshFem->x12);
    free(wingMeshFem->x23);
    free(wingMeshFem->x31);
    //
    free(wingMeshFem->y12);
    free(wingMeshFem->y23);
    free(wingMeshFem->y31);
    //
    free(wingMeshFem->area);
    //
    free(wingMeshFem->a4);
    free(wingMeshFem->a5);
    free(wingMeshFem->a6);
    //
    free(wingMeshFem->b4);
    free(wingMeshFem->b5);
    free(wingMeshFem->b6);
    //
    free(wingMeshFem->c4);
    free(wingMeshFem->c5);
    free(wingMeshFem->c6);
    //
    free(wingMeshFem->e4);
    free(wingMeshFem->e5);
    free(wingMeshFem->e6);
    //
    free(wingMeshFem->d4);
    free(wingMeshFem->d5);
    free(wingMeshFem->d6);
    //
    free(wingMeshFem->C4);
    free(wingMeshFem->C5);
    free(wingMeshFem->C6);
    //
    free(wingMeshFem->S4);
    free(wingMeshFem->S5);
    free(wingMeshFem->S6);
    //
    for (int i=0;i<Ng;i++){
        free(wingMeshFem->SF[i]);
        free(wingMeshFem->DxsiSF[i]);
        free(wingMeshFem->DetaSF[i]);
        free(wingMeshFem->SFm[i]);
        free(wingMeshFem->DxsiSFm[i]);
        free(wingMeshFem->DetaSFm[i]);
    }
    free(wingMeshFem->SF);
    free(wingMeshFem->DxsiSF);
    free(wingMeshFem->DetaSF);
    free(wingMeshFem->SFm);
    free(wingMeshFem->DxsiSFm);
    free(wingMeshFem->DetaSFm);
    // 
    free(wingMeshFem->D2xsiSF);
    free(wingMeshFem->D2xsietaSF);
    free(wingMeshFem->D2etaSF);

    for (int i = 0;i<10;i++){
        free(wingMeshFem->GGDST[i]);//[10 x 10]    
        free(wingMeshFem->GGin[i]);
    }
    for (int i = 0;i<6;i++){
        free(wingMeshFem->GGDKT[i]);//[6 X 6] 
        free(wingMeshFem->GGin2[i]); //inverse of above
    }
    free(wingMeshFem->GGDST);
    free(wingMeshFem->GGDKT); //[10 x 10]  
    free(wingMeshFem->GGin);
    free(wingMeshFem->GGin2); //inverse of above
}



template<class T>
void TriGaussPoints(T xw[GaussIntegrPoints][3]){
    /* TODO : Add more options for the gauss integration */
    int Ng = GaussIntegrPoints;
    int Mcol=Ng;
    int Ncol=3;

    if (Ng==1){
        T xw_temp[1][3] = {{0.33333333333333, 0.33333333333333, 1.00000000000000}};

        for (int i=0;i<Mcol;i++){
            for (int j=0;j<Ncol;j++){
                //printf("xw_temp [%d]:%f,",j,xw_temp[i][j]);
                xw[i][j]=xw_temp[i][j];
                //printf("xw [%d]:%f,",j,xw[i][j]);
            }
        }
    }
    
    if (Ng==3){
        T xw_temp[3][3]={{0.16666666666667, 0.16666666666667, 0.33333333333333},
                            {0.16666666666667, 0.66666666666667, 0.33333333333333},
                            {0.66666666666667, 0.16666666666667, 0.33333333333333}};

        for (int i=0;i<Mcol;i++){
            for (int j=0;j<Ncol;j++){
            //printf("xw_temp [%d]:%f,",j,xw_temp[i][j]);
            xw[i][j]=xw_temp[i][j];
            //printf("xw [%d]:%f,",j,xw[i][j]);
            }
        }
    }

    if (Ng==4){
        T xw_temp[4][3]={ {0.33333333333333, 0.33333333333333, -0.56250000000000},
                                {0.20000000000000, 0.20000000000000, 0.52083333333333},
                                {0.20000000000000, 0.60000000000000, 0.52083333333333},
                                {0.60000000000000, 0.20000000000000, 0.52083333333333}};


        for (int i=0;i<Mcol;i++){
            for (int j=0;j<Ncol;j++){
                //printf("xw_temp [%d]:%f,",j,xw_temp[i][j]);
                xw[i][j]=xw_temp[i][j];
                //printf("xw [%d]:%f,",j,xw[i][j]);
            }
        }
    }

    if (Ng==6){
        T xw_temp[6][3]={{0.44594849091597, 0.44594849091597, 0.22338158967801},
                                {0.44594849091597, 0.10810301816807, 0.22338158967801},
                                {0.10810301816807, 0.44594849091597, 0.22338158967801},
                                {0.09157621350977, 0.09157621350977, 0.10995174365532},
                                {0.09157621350977, 0.81684757298046, 0.10995174365532},
                                {0.81684757298046, 0.09157621350977, 0.10995174365532}};


        for (int i=0;i<Mcol;i++){
            for (int j=0;j<Ncol;j++){
                //printf("xw_temp [%d]:%f,",j,xw_temp[i][j]);
                xw[i][j]=xw_temp[i][j];
                //printf("xw [%d]:%f,",j,xw[i][j]);
            }
        }
    }
}

template<class T>
void BendingStiffness(float E, float v, float tx, T **BeSt){
  
    T la = (E*pow(tx,3.0))/(12*(1-pow(v,2)));

    BeSt[0][0] = la;
    BeSt[0][1] = v*la;
    BeSt[1][0] = v*la;
    BeSt[1][1] = la;
    BeSt[2][2] = (1.0-v)*la/2.0; 

}


// from funcBLAS.c

// Allocate 1-D array based on double pointer type
template<class T>
void allocate1Darray(int rows, T **arrIn){

    int i;//, j;

    // allocate local 2D array and pass the pointer to 
    // the pointer to double pointer, otherwise the main
    // can't access the memory.
    T *arrTemp = (T*)malloc(rows * sizeof(T));
    if(arrIn == NULL){
        printf("Memory allocation failed. allocate2Darray()");
        return;
    }

    //printf("\nInside allocate1Darray()..\n");
    // Note that arr[i][j] is same as *(*(arr+i)+j)
    for (i = 0; i < rows; i++){
        arrTemp[i] = 0.0;
        //printf("\n%f,",arrTemp[i]);
    }

    *arrIn = arrTemp; // this will do?
}


// Allocate 2-D array based on double pointer type
template<class T>
void allocate2Darray(int rows, int cols, T ***arrIn){

    int i, j;

    // allocate local 2D array and pass the pointer to 
    // the pointer to double pointer, otherwise the main
    // can't access the memory.
    T **arrTemp = (T**)malloc(rows * sizeof(T*));
    if(arrIn == NULL){
        printf("Memory allocation failed. allocate2Darray()");
        return;
    }

    for (i = 0; i < rows; i++){
        arrTemp[i] = (T*)malloc(cols * sizeof(T));
    }

    //printf("\nInside allocate2Darray()..\n");
    // Note that arr[i][j] is same as *(*(arr+i)+j)
    for (i = 0; i < rows; i++){
        for (j = 0; j < cols; j++){
            arrTemp[i][j] = 0.0;
            //printf("%f,",arrTemp[i][j]);
        } 
        //printf("\n");
    }

    *arrIn = arrTemp; // this will do?
}


// FSI 

template<class T>
void shepard_interp_2d(int nd, T *xd, T *yd, T *zd,
    T *p, int ni, T *xi, T *yi, T *zi){

    printf("\n    Entering SHEPARD INTERP, p=%f...\n",*p);
    
    //printf("nd=%d, ni=%d\n", nd, ni);
    //printf("p=%f\n",*p);

    int z;
    T suma, s;
    T dotproc;

    T *w = ( T * ) malloc ( nd * sizeof ( T ) );

    for (int i=0;i<ni;i++){
        if (abs(*p) < 0.01){
            for ( int j = 0; j < nd; j++ ){
                w[j] = 1.0 / ( double ) ( nd );
            }
            printf("here...\n");
        }
        else{
            //w = zeros ( nd, 1 );
            //for ( int k = 0; k < nd; k++ ){
            //    w[k] = 0.0;
            //}

            z = -1;
            for ( int j = 0; j < nd; j++ ){
                w[j] = sqrt ( pow ( (xi[i] - xd[j]), 2 )
                            + pow ( (yi[i] - yd[j]), 2 ) );
                //printf("w[%d]=%f\n",j,w[j]);
                //printf("%f,%f,%f,%f,%f\n",xi[i],xd[j],yi[i],yd[j],w[j]);
                if ( w[j] == 0.0 ){
                    z = j;
                    break;
                }
            }
            if ( z != -1 ){
                for ( int j = 0; j < nd; j++ ){
                    w[j] = 0.0;
                }
                w[z] = 1.0;
            }
            else{
                for (int j = 0; j < nd; j++ ){
                    w[j] = 1.0 / pow ( w[j], *p );
                }
                suma = 0.0;
                for (int k=0;k<nd;k++){
                    suma = suma + w[k];
                }
                s = suma;
                //s = r8vec_sum ( nd, w );
                for ( int j = 0; j < nd; j++ ){
                    w[j] = w[j] / s;
                }
            }
        }
        dotproc = 0.0;
        for (int k=0;k<nd;k++){
            dotproc = dotproc + w[k]*zd[k];
        }
        zi[i] = dotproc;
        //zi[i] = r8vec_dot_product ( nd, w, zd );
    }

    printf("    Exiting SHEPARD INTERP. OK.\n");
}