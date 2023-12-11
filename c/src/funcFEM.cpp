#ifndef MAIN_FEM_HEADER_FILE
    #define MAIN_FEM_HEADER_FILE

    //#include "../include/mainFem.h"

    #ifndef GaussIntegrPoints
        #define GaussIntegrPoints 3 /* Gauss Integration Points */
    #endif

    #ifndef PRECISION_MODE_FEM
        #define PRECISION_MODE_FEM 1
    #endif

#endif

#ifndef FUNCMAT
    #define FUNCMAT

    #include "../src/funcMat.cpp" // Functions used to facilitate martix, vector operations in c
#endif

#include<stdio.h>
#include<math.h>



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
    int CC; 
    /*type of load: 1- concentrated load, 2- uniform load, 3- distributed load*/
    int LL; 
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
    int pp_rows; // dummy
    int pp_cols;
    T *pp[2] = {0}; /*2-D array*/
    int tt_rows; // dummy
    int tt_cols;
    int *tt[4] = {0}; /*2-D array*/
    int ee_rows; // dummy
    int ee_cols;
    T *ee[7] = {0}; /*2-D array*/
    int sizeBBnodes;
    int *BBnodes = NULL; /* id of nodes affected by boundary conditions*/
    int sizeBdofs;
    int *Bdofs = NULL; /* global numbering of dofs affected by boundary conditions*/

    T P_load; // Pa: positive values point towards the positive Z-axis (reverse for ANSYS) 
    // ONLY FOR CONCENTRATED LOAD <--BELOW
    T P_xy[2] = {0};
    int P_node; 
    // ONLY FOR DISTRIBUTED PROPERTIES LOAD/THICKNESS <--BELOW
    int sizexcp;
    T *xcp = NULL, *ycp = NULL, *fcp = NULL, *tcp = NULL;
    T omega3;
    T dt;
    T Nper;
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
    int *ID[3] = {0};  // [3,NN]
    int *IEN[3]= {0}; // [3,Nelem]
    int *LM[9]= {0};  // [9,Nelem]
    //
    T *xm = NULL; // x barycentric coordinate [Nelem]
    T *ym = NULL; // y barycentric coordinate [Nelem]

    /* output of TrigElCoefsDKT() */ 
    T *l23 = NULL, *l31 = NULL, *l12 = NULL; // [1 x Nelem]
    T *y12 = NULL, *y31 = NULL, *y23 = NULL;
    T *x12 = NULL, *x31 = NULL, *x23 = NULL;
    T *area = NULL;
    T *a4 = NULL, *a5 = NULL, *a6 = NULL,
      *b4 = NULL, *b5 = NULL, *b6 = NULL,
      *c4 = NULL, *c5 = NULL, *c6 = NULL;
    T *d4 = NULL, *d5 = NULL, *d6 = NULL,
      *e4 = NULL, *e5 = NULL, *e6 = NULL;
    T *C4 = NULL, *C5 = NULL, *C6 = NULL,
      *S4 = NULL, *S5 = NULL, *S6 = NULL;
    /*
    The C language is case-sensitive. This means that all language keywords,
    identifiers, function names, and other variables
    must be entered with consistent letter capitalization. 
    */
    /* output of LNShapeFunDST() */
    T **SF = NULL, **DxsiSF = NULL, **DetaSF = NULL; //[Ng x 6]
    T *D2xsiSF = NULL, *D2xsietaSF = NULL, *D2etaSF = NULL;// [1 x 6]
    /* output of LNShapeFunMassDST() */
    T **SFm = NULL;
    T **DxsiSFm = NULL;
    T **DetaSFm = NULL; //[Ng x 3]

    T **GGDST = NULL, **GGDKT = NULL; //[10 x 10]  
    T **GGin = NULL, **GGin2 = NULL; //inverse of above
};

template<class T>
struct femArraysDKT{
    /* global intermediate matrices Mg, Kg, Fglob */
    T **Fglob; //[GEN x 1] global
    T **Hm, **HW; //[10 x 9] overwrite massHmDKT()
    T **kloc, **mloc;//[9 x 9]
    T **floc; // [9 x 1]
    T **floc1; // [10 x 1]
    T **Hxx, **Hyy; // [6 x 9] overwrite rotationMass2()
    T **Hx, **Hy; // [1 x 9] from ShapeFunDKT2()
    /*
    NEW
    */
    T *Hx_xsi, *Hx_eta, *Hy_xsi, *Hy_eta; // [1 x 9] from ShapeFunDKT2()

    T **Bb; // [3 x 9] from ShapeFunDKT2()
    T **LW; // [1 x 10] from pseudoMassDKT()
    T *L; // [1 x 6] from pseudoMassDKT()
    T **Mg, **Kg; // [81 x Nelem] pre-assembly matrices
};


/*=========================================================================================*/
/* Function prototypes*/
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
void BendingStiffness(T E, T v, T tx, T **BeSt);
//
template<class T>
void TrigElCoefsDKT(struct InDataRecFem<T> *inDataFem, struct triangleDKT<T> *wingMeshFem);
//
// Calculation of shape functions and their derivatives at gauss points on
// the parent element
template<class T>
void LNShapeFunDST(int Ng, T xw[GaussIntegrPoints][3], struct triangleDKT<T> *wingMeshFem);
//
template<class T>
void LNShapeFunMassDST(int Ng, T xw[GaussIntegrPoints][3], struct triangleDKT<T> *wingMeshFem);
//
// ax, ay, bx,by are constant for constant h (independednt of î,ç)
template<class T>
void matrixG(struct triangleDKT<T> *wingMeshFem);
//
// utilities
template<class T>
void assignRowArrayMatrixG_DST(int rowID, T xsi, T eta, T **array);
//
template<class T>
void assignRowArrayMatrixG_DKT(int rowID, T xsi, T eta, T **array);
//
template<class T>
void freefemArraysDKT(struct triangleDKT<T> *wingMeshFem, struct femArraysDKT<T> *elemFemArr);  
//
template<class T>
void massHmDKT(int kk, struct triangleDKT<T> *wingMeshFem, struct femArraysDKT<T> *elemFemArr);
//
template<class T>
void rotationMass2(int kk, struct triangleDKT<T> *wingMeshFem, struct femArraysDKT<T> *elemFemArr);
//
template<class T>
void ShapeFunDKT2(int ii, int kk, struct triangleDKT<T> *wingMeshFem, struct femArraysDKT<T> *elemFemArr);
//
template<class T>
void pseudoMassDKT(int ii, int kk, struct triangleDKT<T> *wingMeshFem, struct femArraysDKT<T> *elemFemArr);
//
template<class T>
void CuFEMNum2DWriteDataInBinary(int rows, int cols, T **Usol, int GEN);
//
template<class T>
void CuFEMNum2DWriteMatrix(int rows, int cols, T **K, T **M, T **F);
//
//-----------------added in 10/05/2023
template<class T>
void RayleighDampingCoefs(T *a, T *b); // TO DO
//
template<class T>
void createRHS(struct InDataRecFem<T> *inDataFem, 
                struct triangleDKT<T> *wingMeshFem,
                struct femArraysDKT<T> *elemFemArr,
                 T *distrLoad, T **G, int d);

/*=========================================================================================*/
/* Definition of the functions BELOW */
/*=========================================================================================*/
template<class T>
void CuFEMNum2DReadInData(struct InDataRecFem<T> *inDataFem ){
    printf("\n    Entering CuFEMNum2DReadInData().\n");
    FILE *file;
    #if PRECISION_MODE_FEM == 1
	    file = fopen("../c/INDATA_FEM_double.bin", "rb"); // r for read, b for binary
    #endif
    #if PRECISION_MODE_FEM == 2
	    file = fopen("../c/INDATA_FEM_single.bin", "rb"); // r for read, b for binary
    #endif
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
    #if (DYNAMIC_ANALYSIS == 1)
        /* DYNAMIC ANALYSIS */
        fread(&(inDataFem->omega3), sizeof(T) , 1, file);
        fread(&(inDataFem->dt), sizeof(T) , 1, file);
        fread(&(inDataFem->Nper), sizeof(T) , 1, file);
    #endif
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
    free(inDataFem->BBnodes);
    free(inDataFem->Bdofs);
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
        wingMeshFem->IEN[i] = (int*)malloc( (inDataFem->tt_cols) *sizeof(int));
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
void BendingStiffness(T E, T v, T tx, T **BeSt){
  
    T la = (E*mypow<T>(tx,3.0))/(12*(1-mypow<T>(v,2)));

    BeSt[0][0] = la;
    BeSt[0][1] = v*la;
    BeSt[1][0] = v*la;
    BeSt[1][1] = la;
    BeSt[2][2] = (1.0-v)*la/2.0; 

}

template<class T>
void TrigElCoefsDKT(struct InDataRecFem<T> *inDataFem, struct triangleDKT<T> *wingMeshFem){
    /* memory allocation */
    wingMeshFem->l23 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->l31 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->l12 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    //
    wingMeshFem->y12 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->y31 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->y23 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    //
    wingMeshFem->x12 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->x31 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->x23 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    //
    wingMeshFem->area = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    //
    wingMeshFem->a4 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->a5 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->a6 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    //
    wingMeshFem->b4 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->b5 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->b6 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    //
    wingMeshFem->c4 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->c5 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->c6 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    //
    wingMeshFem->d4 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->d5 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->d6 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    //
    wingMeshFem->e4 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->e5 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->e6 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    //
    wingMeshFem->C4 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->C5 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->C6 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    //
    wingMeshFem->S4 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->S5 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));
    wingMeshFem->S6 = (T*)malloc(wingMeshFem->Nelem *sizeof(T));

    int IEN_1, IEN_2, IEN_3;    
    //i = wingMeshFem->Nelem-1;
    for (int i=0;i<wingMeshFem->Nelem;i++){ 
        IEN_1 = wingMeshFem->IEN[0][i] - 1;
        IEN_2 = wingMeshFem->IEN[1][i] - 1;
        IEN_3 = wingMeshFem->IEN[2][i] - 1;
        //printf("%d, %d, %d, \n",IEN_1, IEN_2, IEN_3);

                //y23=y(IEN(2,:))-y(IEN(3,:));  % 2-3
                //y31=y(IEN(3,:))-y(IEN(1,:));  % 3-1
                //y12=y(IEN(1,:))-y(IEN(2,:));  % 1-2

        wingMeshFem->y23[i] = (inDataFem->pp[1][IEN_2] - inDataFem->pp[1][IEN_3]);
        wingMeshFem->y31[i] = (inDataFem->pp[1][IEN_3] - inDataFem->pp[1][IEN_1]);
        wingMeshFem->y12[i] = (inDataFem->pp[1][IEN_1] - inDataFem->pp[1][IEN_2]);

                //x23=x(IEN(2,:))-x(IEN(3,:));% 2-3
                //x31=x(IEN(3,:))-x(IEN(1,:));% 3-1
                //x12=x(IEN(1,:))-x(IEN(2,:));% 1-2

        wingMeshFem->x23[i] = (inDataFem->pp[0][IEN_2] - inDataFem->pp[0][IEN_3]);
        wingMeshFem->x31[i] = (inDataFem->pp[0][IEN_3] - inDataFem->pp[0][IEN_1]);
        wingMeshFem->x12[i] = (inDataFem->pp[0][IEN_1] - inDataFem->pp[0][IEN_2]);

        // Length of triagle sides with corresponding (s) node
        wingMeshFem->l23[i]=mysqrt<T>(mypow<T>(wingMeshFem->x23[i],2.0)+mypow<T>(wingMeshFem->y23[i],2.0));
        wingMeshFem->l31[i]=mysqrt<T>(mypow<T>(wingMeshFem->x31[i],2.0)+mypow<T>(wingMeshFem->y31[i],2.0));
        wingMeshFem->l12[i]=mysqrt<T>(mypow<T>(wingMeshFem->x12[i],2.0)+mypow<T>(wingMeshFem->y12[i],2.0));

        wingMeshFem->a4[i]=-wingMeshFem->x23[i]/mypow<T>(wingMeshFem->l23[i],2);
        wingMeshFem->a5[i]=-wingMeshFem->x31[i]/mypow<T>(wingMeshFem->l31[i],2);
        wingMeshFem->a6[i]=-wingMeshFem->x12[i]/mypow<T>(wingMeshFem->l12[i],2);

        wingMeshFem->b4[i]=(3.0/4.0)*(wingMeshFem->x23[i])*(wingMeshFem->y23[i])/mypow<T>(wingMeshFem->l23[i],2.0);
        wingMeshFem->b5[i]=(3.0/4.0)*(wingMeshFem->x31[i])*(wingMeshFem->y31[i])/mypow<T>(wingMeshFem->l31[i],2.0);
        wingMeshFem->b6[i]=(3.0/4.0)*(wingMeshFem->x12[i])*(wingMeshFem->y12[i])/mypow<T>(wingMeshFem->l12[i],2.0);

        wingMeshFem->c4[i]=((1.0/4.0)*mypow<T>(wingMeshFem->x23[i],2.0)-(1.0/2.0)*mypow<T>(wingMeshFem->y23[i],2.0))/mypow<T>(wingMeshFem->l23[i],2.0);
        wingMeshFem->c5[i]=((1.0/4.0)*mypow<T>(wingMeshFem->x31[i],2.0)-(1.0/2.0)*mypow<T>(wingMeshFem->y31[i],2.0))/mypow<T>(wingMeshFem->l31[i],2.0);
        wingMeshFem->c6[i]=((1.0/4.0)*mypow<T>(wingMeshFem->x12[i],2.0)-(1.0/2.0)*mypow<T>(wingMeshFem->y12[i],2.0))/mypow<T>(wingMeshFem->l12[i],2.0);

        wingMeshFem->d4[i]=-wingMeshFem->y23[i]/mypow<T>(wingMeshFem->l23[i],2);
        wingMeshFem->d5[i]=-wingMeshFem->y31[i]/mypow<T>(wingMeshFem->l31[i],2);
        wingMeshFem->d6[i]=-wingMeshFem->y12[i]/mypow<T>(wingMeshFem->l12[i],2);

        wingMeshFem->e4[i]=((1.0/4.0)*mypow<T>(wingMeshFem->y23[i],2.0)-(1.0/2.0)*mypow<T>(wingMeshFem->x23[i],2.0))/mypow<T>(wingMeshFem->l23[i],2.0);
        wingMeshFem->e5[i]=((1.0/4.0)*mypow<T>(wingMeshFem->y31[i],2.0)-(1.0/2.0)*mypow<T>(wingMeshFem->x31[i],2.0))/mypow<T>(wingMeshFem->l31[i],2.0);
        wingMeshFem->e6[i]=((1.0/4.0)*mypow<T>(wingMeshFem->y12[i],2.0)-(1.0/2.0)*mypow<T>(wingMeshFem->x12[i],2.0))/mypow<T>(wingMeshFem->l12[i],2.0);

        wingMeshFem->area[i]=(1.0/2.0)*(wingMeshFem->x31[i]*wingMeshFem->y12[i]-wingMeshFem->x12[i]*wingMeshFem->y31[i]);

        wingMeshFem->C4[i]=-wingMeshFem->y23[i]/wingMeshFem->l23[i];
        wingMeshFem->C5[i]=-wingMeshFem->y31[i]/wingMeshFem->l31[i];
        wingMeshFem->C6[i]=-wingMeshFem->y12[i]/wingMeshFem->l12[i];

        wingMeshFem->S4[i]=wingMeshFem->x23[i]/wingMeshFem->l23[i];
        wingMeshFem->S5[i]=wingMeshFem->x31[i]/wingMeshFem->l31[i];
        wingMeshFem->S6[i]=wingMeshFem->x12[i]/wingMeshFem->l12[i];

    }
#if DEBUG_ON
    printf("    Exiting TrigElCoefsDKT...\n\n");
#endif
}

template<class T>
void LNShapeFunDST(int Ng, T xw[GaussIntegrPoints][3], struct triangleDKT<T> *wingMeshFem){
    // Ng, xw are given data based on which we will fill up some matrices
    //int Ng = GaussIntegrPoints;

#if DEBUG_ON    
    printf("    Ng: %d\n", Ng);

    for (int i=0;i<Ng;i++){
        for (int j=0;j<3;j++){
            printf("    xw [%d]:%f,",j,xw[i][j]);
        }
        printf("\n");
    }
#endif

    // Initialize matrices that are [Ng x 6]   
    wingMeshFem->SF = (T**)malloc(Ng *sizeof(T*)); // pointer array with Ng rows
    wingMeshFem->DxsiSF = (T**)malloc(Ng *sizeof(T*)); // pointer array with Ng rows
    wingMeshFem->DetaSF = (T**)malloc(Ng *sizeof(T*)); // pointer array with Ng rows
    for (int i=0;i<Ng;i++){
        wingMeshFem->SF[i] = (T*)malloc(6 *sizeof(T));
        wingMeshFem->DxsiSF[i] = (T*)malloc(6 *sizeof(T));
        wingMeshFem->DetaSF[i] = (T*)malloc(6 *sizeof(T));
    }

    T xg, yg;
    int i;
    for (i=0; i<Ng; i++){
        
        xg = xw[i][0];
        yg = xw[i][1];

        wingMeshFem->SF[i][0]=(1.0 - xg-yg)*(1.0-2.0*xg-2.0*yg); //checked
        wingMeshFem->SF[i][1]=xg*(2.0*xg-1.0);//checked
        wingMeshFem->SF[i][2]=yg*(2.0*yg-1.0);  
        wingMeshFem->SF[i][3]=4.0*xg*yg;
        wingMeshFem->SF[i][4]=4.0*yg*(1.0-xg-yg);
        wingMeshFem->SF[i][5]=4.0*xg*(1.0-xg-yg);
        
        // 1st derivative --> î
        wingMeshFem->DxsiSF[i][0]=-3.0+4.0*(xg+yg);//checked
        wingMeshFem->DxsiSF[i][1]=4.0*xg-1.0;//checked
        wingMeshFem->DxsiSF[i][2]=0.0;//checked
        wingMeshFem->DxsiSF[i][3]=4.0*yg;//checked                  
        wingMeshFem->DxsiSF[i][4]=-4.0*yg;//checked
        wingMeshFem->DxsiSF[i][5]=4.0-8.0*xg-4.0*yg;
        
        //1st derivative ç
        wingMeshFem->DetaSF[i][0]=-3.0+4.0*(xg+yg);
        wingMeshFem->DetaSF[i][1]=0.0;
        wingMeshFem->DetaSF[i][2]=4.0*yg-1.0;
        wingMeshFem->DetaSF[i][3]=4.0*xg;        
        wingMeshFem->DetaSF[i][4]=4.0-4.0*xg-8.0*yg;
        wingMeshFem->DetaSF[i][5]=-4.0*xg; 
    }


#if DEBUG_ON
    for (int i=0;i<Ng;i++){
        for (int j=0;j<6;j++){
            //printf("    SF[%d][%d]=%f, ", i,j,wingMeshFem->SF[i][j] );
            //printf("    DxsiSF[%d][%d]=%f, ", i,j,wingMeshFem->DxsiSF[i][j] );
            printf("    DetaSF[%d][%d]=%f, ", i,j,wingMeshFem->DetaSF[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");

#endif

    //float *SF, *DxsiSF, *DetaSF; //[Ng x 6]
    //    T *D2xsiSF = NULL, *D2xsietaSF = NULL, *D2etaSF = NULL;// [1 x 6]
    wingMeshFem->D2xsiSF = (T*)malloc(6 *sizeof(T)); // [1 x 6]
    wingMeshFem->D2xsietaSF = (T*)malloc(6 *sizeof(T)); // [1 x 6]
    wingMeshFem->D2etaSF = (T*)malloc(6 *sizeof(T)); // [1 x 6]

    // 2nd derivative--> î    %checked
    wingMeshFem->D2xsiSF[0]=4.0;
    wingMeshFem->D2xsiSF[1]=4.0;
    wingMeshFem->D2xsiSF[2]=0.0;
    wingMeshFem->D2xsiSF[3]=0.0;
    wingMeshFem->D2xsiSF[4]=0.0;
    wingMeshFem->D2xsiSF[5]=-8.0;

    // mixed derivative--> î,ç symmetric!%checked
    wingMeshFem->D2xsietaSF[0]=4.0;
    wingMeshFem->D2xsietaSF[1]=0.0;
    wingMeshFem->D2xsietaSF[2]=0.0;
    wingMeshFem->D2xsietaSF[3]=4.0;
    wingMeshFem->D2xsietaSF[4]=-4.0;
    wingMeshFem->D2xsietaSF[5]=-4.0;
    
    // 2nd derivative ç    %checked
    wingMeshFem->D2etaSF[0]=4.0;
    wingMeshFem->D2etaSF[1]=0.0;
    wingMeshFem->D2etaSF[2]=4.0;
    wingMeshFem->D2etaSF[3]=0.0;
    wingMeshFem->D2etaSF[4]=-8.0;
    wingMeshFem->D2etaSF[5]=0.0;

#if DEBUG_ON
    for (int j=0;j<6;j++){
            //printf("D2xsiSF[%d]=%f, ", j,wingMeshFem->D2xsiSF[j]);
            //printf("D2xsietaSF[%d]=%f, ", j,wingMeshFem->D2xsietaSF[j]);
            printf("    D2etaSF[%d]=%f, ", j,wingMeshFem->D2etaSF[j]);
    }
#endif
#if DEBUG_ON
    printf("\n    Exiting LNShapeFunDST...\n");
#endif
}

//
template<class T>
void LNShapeFunMassDST(int Ng, T xw[GaussIntegrPoints][3], struct triangleDKT<T> *wingMeshFem){
    // Ng, xw are given data based on which we will fill up some matrices


#if DEBUG_ON    
    printf("Ng: %d\n", Ng);

    for (int i=0;i<Ng;i++){
        for (int j=0;j<3;j++){
            printf("xw [%d]:%f,",j,xw[i][j]);
        }
        printf("\n");
    }
#endif

    // Initialize matrices that are Ng x 3
    wingMeshFem->SFm = (T**)malloc(Ng *sizeof(T*)); // pointer array with Ng rows
    wingMeshFem->DxsiSFm = (T**)malloc(Ng *sizeof(T*)); // pointer array with Ng rows
    wingMeshFem->DetaSFm = (T**)malloc(Ng *sizeof(T*)); // pointer array with Ng rows
    for (int i=0;i<Ng;i++){
        wingMeshFem->SFm[i] = (T*)malloc(3 *sizeof(T));
        wingMeshFem->DxsiSFm[i] = (T*)malloc(3 *sizeof(T));
        wingMeshFem->DetaSFm[i] = (T*)malloc(3 *sizeof(T));
    }

    T xg, yg;
    int i;
    for (i=0; i<Ng; i++){
        
        xg = xw[i][0];
        yg = xw[i][1];

        wingMeshFem->SFm[i][0]=(1.0-xg-yg); //checked
        wingMeshFem->SFm[i][1]=xg;//checked
        wingMeshFem->SFm[i][2]=yg;

        // 1st derivative --> î
        wingMeshFem->DxsiSFm[i][0]=-1.0; //checked
        wingMeshFem->DxsiSFm[i][1]=1.0; //checked
        wingMeshFem->DxsiSFm[i][2]=0.0; //checked

        //1st derivative ç
        wingMeshFem->DetaSFm[i][0]=-1.0;
        wingMeshFem->DetaSFm[i][1]=0.0;
        wingMeshFem->DetaSFm[i][2]=1.0;

    }

#if DEBUG_ON
    for (int i=0;i<Ng;i++){
        for (int j=0;j<3;j++){
            //printf("SFm[%d][%d]=%f, ", i,j,wingMeshFem->SFm[i][j] );
            //printf("DxsiSFm[%d][%d]=%f, ", i,j,wingMeshFem->DxsiSFm[i][j] );
            printf("DetaSFm[%d][%d]=%f, ", i,j,wingMeshFem->DetaSFm[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");

#endif
#if DEBUG_ON
    printf("EXITING LNShapeFunMassDST...\n\n");
#endif

}

template<class T>
void matrixG(struct triangleDKT<T> *wingMeshFem){

    int r = 10, c = 10, i, j;
    T xsi, eta, rowID;
 
    wingMeshFem->GGDST = (T**)malloc(r * sizeof(T*));
    wingMeshFem->GGin = (T**)malloc(r * sizeof(T*));//inverse matrix
    for (i = 0; i < r; i++){
        wingMeshFem->GGDST[i] = (T*)malloc(c * sizeof(T));
        wingMeshFem->GGin[i] = (T*)malloc(c * sizeof(T));
    }

    // Note that arr[i][j] is same as *(*(arr+i)+j)
    for (i = 0; i < r; i++){
        for (j = 0; j < c; j++){
                wingMeshFem->GGDST[i][j] = 0.0;
                wingMeshFem->GGin[i][j] = 0.0;
                //printf("%f, ", wingMeshFem->GGDST[i][j]);
        }
        //printf("\n\n");
    }

    r = 6; c = 6;
    wingMeshFem->GGDKT = (T**)malloc(r * sizeof(T*));
    wingMeshFem->GGin2 = (T**)malloc(r * sizeof(T*));
    for (i = 0; i < r; i++){
        wingMeshFem->GGDKT[i] = (T*)malloc(c * sizeof(T));
        wingMeshFem->GGin2[i] = (T*)malloc(c * sizeof(T));
    }

    // Note that arr[i][j] is same as *(*(arr+i)+j)
    for (i = 0; i < r; i++){
        for (j = 0; j < c; j++){
                wingMeshFem->GGDKT[i][j] = 0.0;
                wingMeshFem->GGin2[i][j] = 0.0;
                //printf("%f, ", wingMeshFem->GGDKT[i][j]);
        }
        //printf("\n\n");
    }

    assignRowArrayMatrixG_DST<T>(rowID = 1.0, xsi=0.0, eta = 0.0, wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST<T>(rowID = 2.0, xsi=1.0, eta = 0.0, wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST<T>(rowID = 3.0, xsi=0.0, eta = 1.0, wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST<T>(rowID = 4.0, xsi=(1.0/3.0), eta = (1.0/3.0), wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST<T>(rowID = 5.0, xsi=(2.0/3.0), eta = (1.0/3.0), wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST<T>(rowID = 6.0, xsi=(1.0/3.0), eta = (2.0/3.0), wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST<T>(rowID = 7.0, xsi=0.0, eta = (2.0/3.0), wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST<T>(rowID = 8.0, xsi=0.0, eta = (1.0/3.0), wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST<T>(rowID = 9.0, xsi=(1.0/3.0), eta = 0.0, wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST<T>(rowID = 10.0, xsi=(2.0/3.0), eta = 0.0, wingMeshFem->GGDST);
    //
    assignRowArrayMatrixG_DKT<T>(rowID = 1.0, xsi=0.0, eta = 0.0, wingMeshFem->GGDKT);
    assignRowArrayMatrixG_DKT<T>(rowID = 2.0, xsi=1.0, eta = 0.0, wingMeshFem->GGDKT);
    assignRowArrayMatrixG_DKT<T>(rowID = 3.0, xsi=0.0, eta = 1.0, wingMeshFem->GGDKT);
    assignRowArrayMatrixG_DKT<T>(rowID = 4.0, xsi=(1.0/2.0), eta = (1.0/2.0), wingMeshFem->GGDKT);
    assignRowArrayMatrixG_DKT<T>(rowID = 5.0, xsi=0.0, eta = (1.0/2.0), wingMeshFem->GGDKT);
    assignRowArrayMatrixG_DKT<T>(rowID = 6.0, xsi=(1.0/2.0), eta = 0.0, wingMeshFem->GGDKT);
}


/* 
COMMENT: The following functions are not dependent on my new datatypes! Therefore,
they can be included without the structure declaration in a separate .c file! That is why in the BEM code
the DATASTRUCTURES are defined in multiple destinations within preprocessor directives.
*/
template<class T>
void assignRowArrayMatrixG_DST(int rowID, T xsi, T eta, T **array){

    T rowMat[10] = {1.0, xsi, eta, xsi*eta, mypow<T>(xsi,2.0), mypow<T>(eta,2.0),
     mypow<T>(xsi,2)*eta, mypow<T>(eta,2)*xsi, mypow<T>(xsi,3.0), mypow<T>(eta,3.0)};

    //for (int i=0;i<10;i++){
    //    printf("%f, ", rowMat[i]);
    //}

    //printf("\nInside assignRowArrayMatrixG_DST()... \n\n");
    int i = rowID -1;
    for (int j=0;j<10;j++){
        array[i][j]=rowMat[j];
        //printf("%d, %f, ",j, array[I][j]);
    }
    //printf("\n\n");
}

template<class T>
void assignRowArrayMatrixG_DKT(int rowID, T xsi, T eta, T **array){

    T rowMat[6] = {1.0, xsi, eta, xsi*eta, mypow<T>(xsi,2.0), mypow<T>(eta,2.0)};

    //for (int i=0;i<6;i++){
    //    printf("%f, ", rowMat[i]);
    //}

    //printf("\nInside assignRowArrayMatrixG_DKT()... \n\n");
    int i = rowID -1;
    for (int j=0;j<6;j++){
        array[i][j]=rowMat[j];
    //    printf("%d, %f, ",j, array[I][j]);
    }
    //printf("\n\n");
}

template<class T>
void freefemArraysDKT(struct triangleDKT<T> *wingMeshFem, struct femArraysDKT<T> *elemFemArr){

/*
struct femArraysDKT
{
    float **Fglob; //[GEN x 1] global
    float **Hm, **HW; //[10 x 9] overwrite massHmDKT()
    float **kloc, **mloc;//[9 x 9]
    float **floc; // [9 x 1]
    float **floc1; // [10 x 1]
    float **Hxx, **Hyy; // [6 x 9] overwrite rotationMass2()
    float **Hx, **Hy; // [1 x 9] from ShapeFunDKT2()
    float *Hx_xsi, *Hx_eta, *Hy_xsi, *Hy_eta; // [1 x 9] from ShapeFunDKT2()
    float **Bb; // [3 x 9] from ShapeFunDKT2()
    float **LW; // [1 x 10] from pseudoMassDKT()
    float *L; // [1 x 6] from pseudoMassDKT()
    float **Mg, **Kg; // [81 x Nelem] pre-assembly matrices
}
*/
// How to free double pointer.
for (int i=0;i<wingMeshFem->GEN;i++){
    free(elemFemArr->Fglob[i]);
}
free(elemFemArr->Fglob);
//
for (int i=0;i<10;i++){
    //printf("elemFemArr->Hm[%d]=%f\n",i,elemFemArr->Hm[i][0]);
    free(elemFemArr->Hm[i]);
    free(elemFemArr->HW[i]);
}
free(elemFemArr->Hm);
free(elemFemArr->HW);
//
for (int i=0;i<9;i++){
    free(elemFemArr->kloc[i]);
    free(elemFemArr->mloc[i]);
}
free(elemFemArr->kloc);
free(elemFemArr->mloc);
//
for (int i=0;i<9;i++){
    free(elemFemArr->floc[i]);
}
free(elemFemArr->floc);
//
for (int i=0;i<10;i++){
    free(elemFemArr->floc1[i]);
}
free(elemFemArr->floc1);
//
for (int i=0;i<6;i++){
    free(elemFemArr->Hxx[i]);
    free(elemFemArr->Hyy[i]);
}
free(elemFemArr->Hxx);
free(elemFemArr->Hyy);
//
for (int i=0;i<1;i++){
    free(elemFemArr->Hx[i]);
    free(elemFemArr->Hy[i]);
    free(elemFemArr->LW[i]);
}
free(elemFemArr->Hx);
free(elemFemArr->Hy);
free(elemFemArr->LW);// [1 x 10] from pseudoMassDKT()
//
free(elemFemArr->Hx_xsi);
free(elemFemArr->Hx_eta);
free(elemFemArr->Hy_xsi);
free(elemFemArr->Hy_eta);
for (int i=0;i<3;i++){
    free(elemFemArr->Bb[i]);
}
free(elemFemArr->Bb); // [3 x 9] from ShapeFunDKT2()
free(elemFemArr->L); // [1 x 6] from pseudoMassDKT()
//

for (int i=0;i<81;i++){
    free(elemFemArr->Mg[i]);// [81 x Nelem] pre-assembly matrices
    free(elemFemArr->Kg[i]);
}
free(elemFemArr->Mg);
free(elemFemArr->Kg);


}

template<class T>
void massHmDKT(int kk, struct triangleDKT<T> *wingMeshFem, struct femArraysDKT<T> *elemFemArr){

    elemFemArr->Hm[0][0] = 1.0;// Hm(1,:)
    elemFemArr->Hm[1][3] = 1.0;// Hm(2,:)
    elemFemArr->Hm[2][6] = 1.0;// Hm(3,:)

    elemFemArr->Hm[3][0] = (1.0/3.0);// Hm(4,:)
    elemFemArr->Hm[3][3] = (1.0/3.0);
    elemFemArr->Hm[3][6] = (1.0/3.0);

    //printf("\nPrinting Hm...\n");
    //for (int i=0;i<10;i++){
    //    for (int j=0;j<9;j++){
    //        printf("%f, ", elemFemArr->Hm[i][j]);
    //    }
    //    printf("\n");
    //}

    T Hm_5[] = {(T)0.0, (T)0.0, (T)0.0,
     (T)(20.0/27.0), (T)4.0*(wingMeshFem->l23[kk])/(T)27.0*(wingMeshFem->S4[kk]), (T)-4.0*(wingMeshFem->l23[kk])/(T)27.0*(wingMeshFem->C4[kk]),
     (T)(7.0/27.0), (T)-2.0*(wingMeshFem->l23[kk])/(T)27.0*(wingMeshFem->S4[kk]), (T)2.0*(wingMeshFem->l23[kk])/(T)27.0*(wingMeshFem->C4[kk])};

    T Hm_6[] = {(T)0.0, (T)0.0, (T)0.0,
     (T)(7.0/27.0), (T)2.0*(wingMeshFem->l23[kk])/(T)27.0*(wingMeshFem->S4[kk]), (T)-2.0*(wingMeshFem->l23[kk])/(T)27.0*(wingMeshFem->C4[kk]),
     (T)(20.0/27.0), (T)-4.0*(wingMeshFem->l23[kk])/(T)27.0*(wingMeshFem->S4[kk]), (T)4.0*(wingMeshFem->l23[kk])/(T)27.0*(wingMeshFem->C4[kk])};

    T Hm_7[] = {(T)(7.0/27.0), (T)-2.0*wingMeshFem->l31[kk]/(T)27.0*wingMeshFem->S5[kk], (T)2.0*wingMeshFem->l31[kk]/(T)27.0*wingMeshFem->C5[kk],
     (T)0.0, (T)0.0, (T)0.0,
     (T)(20.0/27.0), (T)4.0*wingMeshFem->l31[kk]/(T)27.0*wingMeshFem->S5[kk], (T)-4.0*wingMeshFem->l31[kk]/(T)27.0*wingMeshFem->C5[kk]};

    T Hm_8[] = {(T)(20.0/27.0), -(T)4.0*wingMeshFem->l31[kk]/(T)27.0*wingMeshFem->S5[kk], (T)4.0*wingMeshFem->l31[kk]/(T)27.0*wingMeshFem->C5[kk],
     (T)0.0, (T)0.0, (T)0.0,
     (7.0/27.0), (T)2.0*wingMeshFem->l31[kk]/(T)27.0*wingMeshFem->S5[kk], -(T)2.0*wingMeshFem->l31[kk]/(T)27.0*wingMeshFem->C5[kk]};

    T Hm_9[] = {(T)(20.0/27.0), (T)4.0*wingMeshFem->l12[kk]/(T)27.0*wingMeshFem->S6[kk], -(T)4.0*wingMeshFem->l12[kk]/(T)27.0*wingMeshFem->C6[kk],
    (T)(7.0/27.0), -(T)2.0*wingMeshFem->l12[kk]/(T)27.0*wingMeshFem->S6[kk], (T)2.0*wingMeshFem->l12[kk]/(T)27.0*wingMeshFem->C6[kk],
      (T)0.0, (T)0.0, (T)0.0};

    T Hm_10[] = {(T)(7.0/27.0), (T)2.0*wingMeshFem->l12[kk]/(T)27.0*wingMeshFem->S6[kk], -(T)2.0*wingMeshFem->l12[kk]/(T)27.0*wingMeshFem->C6[kk],
     (T)(20.0/27.0), -(T)4.0*wingMeshFem->l12[kk]/(T)27.0*wingMeshFem->S6[kk], (T)4.0*wingMeshFem->l12[kk]/(T)27.0*wingMeshFem->C6[kk],
      (T)0.0, (T)0.0, (T)0.0};

    for (int j=0; j<9; j++){
        elemFemArr->Hm[4][j]=Hm_5[j];
        elemFemArr->Hm[5][j]=Hm_6[j];
        elemFemArr->Hm[6][j]=Hm_7[j];
        elemFemArr->Hm[7][j]=Hm_8[j];
        elemFemArr->Hm[8][j]=Hm_9[j];
        elemFemArr->Hm[9][j]=Hm_10[j];
    }

    //printf("\nPrinting Hm...\n");
    //for (int i=0;i<10;i++){
    //    for (int j=0;j<9;j++){
    //        printf("%f, ", elemFemArr->Hm[i][j]);
    //    }
    //    printf("\n");
    //}

    //int myVariable = -5;
    /*
    with op( A ) an m by k matrix, op( B )  a  k by n matrix
    
    HW=GGin*Hm;  % [10 x 10] x[10 x 9] 
    */
    int M = 10, N = 10, K = 9;
    int optionCalc = 1; //without transpose
    T alpha = 1.0, beta = 0.0;
    matMatMultiplication2<T>(optionCalc, M, N, K, alpha, beta, wingMeshFem->GGin, elemFemArr->Hm, elemFemArr->HW);

//#if DEBUG_ON
//    printf("\n EXITING massHmDKT...");
//#endif
}

template<class T>
void rotationMass2(int kk, struct triangleDKT<T> *wingMeshFem, struct femArraysDKT<T> *elemFemArr){

// Hxx, Hyy [6 x 9]
T A4 = wingMeshFem->a4[kk];
T A6 = wingMeshFem->a6[kk];
T A5 = wingMeshFem->a5[kk];
//
T B4 = wingMeshFem->b4[kk];
T B5 = wingMeshFem->b5[kk];
T B6 = wingMeshFem->b6[kk];
//
T C4 = wingMeshFem->c4[kk];
T C5 = wingMeshFem->c5[kk];
T C6 = wingMeshFem->c6[kk];
//
T D4 = wingMeshFem->d4[kk];
T D5 = wingMeshFem->d5[kk];
T D6 = wingMeshFem->d6[kk];
//
T E4 = wingMeshFem->e4[kk];
T E5 = wingMeshFem->e5[kk];
T E6 = wingMeshFem->e6[kk];
//
T Hxx1[]={(T)0.0, (T)6.0*A6, -(T)6.0*A5, (T)6.0*A5-(T)6.0*A6, -(T)6.0*A6, (T)6.0*A5}; //columns of Hxx
T Hxx2[]={(T)0.0, (T)4.0*B6, (T)4.0*B5, -(T)4.0*B5-(T)4.0*B6, -(T)4.0*B6, -(T)4*B5};  
T Hxx3[]={(T)1.0, -(T)4.0*C6-(T)3, -(T)3.0-(T)4.0*C5, (T)4.0+(T)4.0*C5+(T)4.0*C6, (T)4.0*C6+(T)2.0, (T)4.0*C5+(T)2.0};
T Hxx4[]={(T)0.0, -(T)6.0*A6, (T)0.0, (T)6.0*A4+(T)6.0*A6, (T)6.0*A6, (T)0.0};
T Hxx5[]={(T)0.0, (T)4.0*B6, (T)0.0, (T)4.0*B4-(T)4.0*B6, -(T)4.0*B6, (T)0.0};
T Hxx6[]={(T)0.0, -(T)1.0-(T)4.0*C6, (T)0.0, (T)4.0*C6-(T)4.0*C4, (T)2.0+(T)4.0*C6, (T)0.0};
T Hxx7[]={(T)0.0, (T)0.0, (T)6.0*A5, -(T)6.0*A5-(T)6.0*A4, (T)0.0, -(T)6.0*A5};
T Hxx8[]={(T)0.0, (T)0.0, (T)4.0*B5, (T)4.0*B4-(T)4.0*B5, (T)0.0, -(T)4.0*B5};
T Hxx9[]={(T)0.0, (T)0.0, -(T)4.0*C5-(T)1.0, (T)4.0*C5-(T)4.0*C4, (T)0.0, (T)4.0*C5+(T)2.0};

for (int i=0; i<6; i++){
    elemFemArr->Hxx[i][0]=Hxx1[i];
    elemFemArr->Hxx[i][1]=Hxx2[i];
    elemFemArr->Hxx[i][2]=Hxx3[i];
    elemFemArr->Hxx[i][3]=Hxx4[i];
    elemFemArr->Hxx[i][4]=Hxx5[i];
    elemFemArr->Hxx[i][5]=Hxx6[i];
    elemFemArr->Hxx[i][6]=Hxx7[i];
    elemFemArr->Hxx[i][7]=Hxx8[i];
    elemFemArr->Hxx[i][8]=Hxx9[i];
}
/*
printf("\n");
for (int i=0; i<6; i++){
    for (int j=0; j<9; j++){
        printf("%f, ", elemFemArr->Hxx[i][j]);
    }
    printf("\n");
}
*/

// now we do the same for Hyy
T Hyy1[]={(T)0.0, (T)6.0*D6, (T)-6.0*D5, (T)6.0*D5-(T)6.0*D6, (T)-6.0*D6, (T)6.0*D5}; //column
T Hyy2[]={(T)-1.0, (T)4,0*E6+(T)3.0, (T)4.0*E5+(T)3.0, -(T)4.0*E5-(T)4.0*E6-(T)4.0, -(T)4.0*E6-(T)2.0, (T)-4.0*E5-(T)2.0};
T Hyy3[]={(T)0.0, (T)-4.0*B6, (T)-4.0*B5, (T)4.0*B5+(T)4.0*B6, (T)4.0*B6, (T)4.0*B5};
T Hyy4[]={(T)0.0, (T)-6.0*D6, (T)0.0, (T)6.0*D4+(T)6.0*D6, (T)6.0*D6, (T)0.0};
T Hyy5[]={(T)0.0, (T)1.0+(T)4.0*E6, (T)0.0, (T)4.0*E4-(T)4.0*E6, -(T)4.0*E6-(T)2.0, (T)0.0};
T Hyy6[]={(T)0.0, (T)-4.0*B6, (T)0.0, (T)-4.0*B4+(T)4.0*B6, (T)4.0*B6, (T)0.0};
T Hyy7[]={(T)0.0, (T)0.0, (T)6.0*D5, (T)-6.0*D5-(T)6.0*D4, 0.0, -(T)6.0*D5};
T Hyy8[]={(T)0.0, (T)0.0, (T)1.0+(T)4.0*E5, (T)4.0*E4-(T)4.0*E5, (T)0.0, -(T)4.0*E5-(T)2.0};
T Hyy9[]={(T)0.0, (T)0.0, (T)-4.0*B5, (T)-4.0*B4+(T)4.0*B5, (T)0.0, (T)4.0*B5};

for (int i=0; i<6; i++){
    elemFemArr->Hyy[i][0]=Hyy1[i];
    elemFemArr->Hyy[i][1]=Hyy2[i];
    elemFemArr->Hyy[i][2]=Hyy3[i];
    elemFemArr->Hyy[i][3]=Hyy4[i];
    elemFemArr->Hyy[i][4]=Hyy5[i];
    elemFemArr->Hyy[i][5]=Hyy6[i];
    elemFemArr->Hyy[i][6]=Hyy7[i];
    elemFemArr->Hyy[i][7]=Hyy8[i];
    elemFemArr->Hyy[i][8]=Hyy9[i];
}
/*
printf("\n");
for (int i=0; i<6; i++){
    for (int j=0; j<9; j++){
        printf("%f, ", elemFemArr->Hyy[i][j]);
    }
    printf("\n");
}
*/
}

template<class T>
void ShapeFunDKT2(int ii, int kk, struct triangleDKT<T> *wingMeshFem, struct femArraysDKT<T> *elemFemArr){

/* output of LNShapeFunDST() */
//float **SF, **DxsiSF, **DetaSF; //[Ng x 6]
//float *D2xsiSF, *D2xsietaSF, *D2etaSF; // [1 x 6]

/* output of LNShapeFunMassDST() */
//float **SFm, **DxsiSFm, **DetaSFm; //[Ng x 3]

T S4=wingMeshFem->S4[kk];
T S5=wingMeshFem->S5[kk];
T S6=wingMeshFem->S6[kk];
T C4=wingMeshFem->C4[kk];
T C5=wingMeshFem->C5[kk];
T C6=wingMeshFem->C6[kk];
//
T l31=wingMeshFem->l31[kk];
T l12=wingMeshFem->l12[kk];
T l23=wingMeshFem->l23[kk];

//size Hx [1 x 9]
elemFemArr->Hx[0][0]=S5*3.0/(2.0*l31)*(wingMeshFem->SF[ii][4])-S6*3.0/(2.0*l12)*(wingMeshFem->SF[ii][5]);
elemFemArr->Hx[0][1]=-(wingMeshFem->SF[ii][4])*(3.0/4.0*S5*C5)-(wingMeshFem->SF[ii][5])*(3.0/4.0*S6*C6);
elemFemArr->Hx[0][2]=(wingMeshFem->SF[ii][0])+(0.5*mypow<T>(C5,2)-0.25*mypow<T>(S5,2))*(wingMeshFem->SF[ii][4])+(0.5*mypow<T>(C6,2)-0.25*mypow<T>(S6,2))*(wingMeshFem->SF[ii][5]);
//
elemFemArr->Hx[0][3]=S6*3.0/(2.0*l12)*(wingMeshFem->SF[ii][5])-S4*3.0/(2.0*l23)*(wingMeshFem->SF[ii][3]);
elemFemArr->Hx[0][4]=-(wingMeshFem->SF[ii][3])*(3.0/4.0*S4*C4)-(wingMeshFem->SF[ii][5])*(3.0/4.0*S6*C6);
elemFemArr->Hx[0][5]=(wingMeshFem->SF[ii][1])+(0.5*mypow<T>(C4,2)-0.25*mypow<T>(S4,2))*(wingMeshFem->SF[ii][3])+(0.5*mypow<T>(C6,2)-0.25*mypow<T>(S6,2))*(wingMeshFem->SF[ii][5]);
//
elemFemArr->Hx[0][6]=S4*3.0/(2.0*l23)*(wingMeshFem->SF[ii][3])-S5*3.0/(2.0*l31)*(wingMeshFem->SF[ii][4]);
elemFemArr->Hx[0][7]=-(wingMeshFem->SF[ii][3])*(3.0/4.0*S4*C4)-(wingMeshFem->SF[ii][4])*(3.0/4.0*S5*C5);
elemFemArr->Hx[0][8]=(wingMeshFem->SF[ii][2])+(0.5*mypow<T>(C4,2)-0.25*mypow<T>(S4,2))*(wingMeshFem->SF[ii][3])+(0.5*mypow<T>(C5,2)-0.25*mypow<T>(S5,2))*(wingMeshFem->SF[ii][4]);

/*
#if DEBUG_ON
printf("ShapeFunDKT2 ...\n");
printf("\nelemFemArr->Hx...\n");
for (int i=0;i<9;i++){
    printf("%f, ",elemFemArr->Hx[0][i]);
}
#endif
*/
//size Hx_xsi [1 x 9]
elemFemArr->Hx_xsi[0]=S5*3.0/(2.0*l31)*(wingMeshFem->DxsiSF[ii][4])-S6*3.0/(2.0*l12)*(wingMeshFem->DxsiSF[ii][5]);
elemFemArr->Hx_xsi[1]=-(wingMeshFem->DxsiSF[ii][4])*(3.0/4.0*S5*C5)-(wingMeshFem->DxsiSF[ii][5])*(3.0/4.0*S6*C6);
elemFemArr->Hx_xsi[2]=(wingMeshFem->DxsiSF[ii][0])+(0.5*mypow<T>(C5,2)-0.25*mypow<T>(S5,2))*(wingMeshFem->DxsiSF[ii][4])+(0.5*mypow<T>(C6,2)-0.25*mypow<T>(S6,2))*(wingMeshFem->DxsiSF[ii][5]);
//
elemFemArr->Hx_xsi[3]=S6*3.0/(2.0*l12)*(wingMeshFem->DxsiSF[ii][5])-S4*3.0/(2.0*l23)*(wingMeshFem->DxsiSF[ii][3]);
elemFemArr->Hx_xsi[4]=-(wingMeshFem->DxsiSF[ii][3])*(3.0/4.0*S4*C4)-(wingMeshFem->DxsiSF[ii][5])*(3.0/4.0*S6*C6);
elemFemArr->Hx_xsi[5]=(wingMeshFem->DxsiSF[ii][1])+(0.5*mypow<T>(C4,2)-0.25*mypow<T>(S4,2))*(wingMeshFem->DxsiSF[ii][3])+(0.5*mypow<T>(C6,2)-0.25*mypow<T>(S6,2))*(wingMeshFem->DxsiSF[ii][5]);
//
elemFemArr->Hx_xsi[6]=S4*3.0/(2.0*l23)*(wingMeshFem->DxsiSF[ii][3])-S5*3.0/(2.0*l31)*(wingMeshFem->DxsiSF[ii][4]);
elemFemArr->Hx_xsi[7]=-(wingMeshFem->DxsiSF[ii][3])*(3.0/4.0*S4*C4)-(wingMeshFem->DxsiSF[ii][4])*(3.0/4.0*S5*C5);
elemFemArr->Hx_xsi[8]=(wingMeshFem->DxsiSF[ii][2])+(0.5*mypow<T>(C4,2)-0.25*mypow<T>(S4,2))*(wingMeshFem->DxsiSF[ii][3])+(0.5*mypow<T>(C5,2)-0.25*mypow<T>(S5,2))*(wingMeshFem->DxsiSF[ii][4]);

/*
#if DEBUG_ON
printf("\n\nelemFemArr->Hx_xsi...\n");
for (int i=0;i<9;i++){
    printf("%f, ",elemFemArr->Hx_xsi[i]);
}
#endif
*/
elemFemArr->Hx_eta[0]=S5*3.0/(2.0*l31)*(wingMeshFem->DetaSF[ii][4])-S6*3.0/(2.0*l12)*(wingMeshFem->DetaSF[ii][5]);
elemFemArr->Hx_eta[1]=-(wingMeshFem->DetaSF[ii][4])*(3.0/4.0*S5*C5)-(wingMeshFem->DetaSF[ii][5])*(3.0/4.0*S6*C6);
elemFemArr->Hx_eta[2]=(wingMeshFem->DetaSF[ii][0])+(0.5*mypow<T>(C5,2)-0.25*mypow<T>(S5,2))*(wingMeshFem->DetaSF[ii][4])+(0.5*mypow<T>(C6,2)-0.25*mypow<T>(S6,2))*(wingMeshFem->DetaSF[ii][5]);
//
elemFemArr->Hx_eta[3]=S6*3.0/(2.0*l12)*(wingMeshFem->DetaSF[ii][5])-S4*3.0/(2.0*l23)*(wingMeshFem->DetaSF[ii][3]);
elemFemArr->Hx_eta[4]=-(wingMeshFem->DetaSF[ii][3])*(3.0/4.0*S4*C4)-(wingMeshFem->DetaSF[ii][5])*(3.0/4.0*S6*C6);
elemFemArr->Hx_eta[5]=(wingMeshFem->DetaSF[ii][1])+(0.5*mypow<T>(C4,2)-0.25*mypow<T>(S4,2))*(wingMeshFem->DetaSF[ii][3])+(0.5*mypow<T>(C6,2)-0.25*mypow<T>(S6,2))*(wingMeshFem->DetaSF[ii][5]);
//
elemFemArr->Hx_eta[6]=S4*3.0/(2.0*l23)*(wingMeshFem->DetaSF[ii][3])-S5*3.0/(2.0*l31)*(wingMeshFem->DetaSF[ii][4]);
elemFemArr->Hx_eta[7]=-(wingMeshFem->DetaSF[ii][3])*(3.0/4.0*S4*C4)-(wingMeshFem->DetaSF[ii][4])*(3.0/4.0*S5*C5);
elemFemArr->Hx_eta[8]=(wingMeshFem->DetaSF[ii][2])+(0.5*mypow<T>(C4,2)-0.25*mypow<T>(S4,2))*(wingMeshFem->DetaSF[ii][3])+(0.5*mypow<T>(C5,2)-0.25*mypow<T>(S5,2))*(wingMeshFem->DetaSF[ii][4]);

/*
#if DEBUG_ON
printf("\n\nelemFemArr->Hx_eta...\n");
for (int i=0;i<9;i++){
    printf("%f, ",elemFemArr->Hx_eta[i]);
}
#endif
*/
elemFemArr->Hy[0][0]=-C5*3.0/(2.0*l31)*(wingMeshFem->SF[ii][4])+C6*3.0/(2.0*l12)*(wingMeshFem->SF[ii][5]);
elemFemArr->Hy[0][1]=-(wingMeshFem->SF[ii][0])-(0.5*mypow<T>(S5,2)-0.25*mypow<T>(C5,2))*(wingMeshFem->SF[ii][4])-(0.5*mypow<T>(S6,2)-0.25*mypow<T>(C6,2))*(wingMeshFem->SF[ii][5]);
elemFemArr->Hy[0][2]=(wingMeshFem->SF[ii][4])*(3.0/4.0*S5*C5)+(wingMeshFem->SF[ii][5])*(3.0/4.0*S6*C6);
//
elemFemArr->Hy[0][3]=-C6*3/(2*l12)*(wingMeshFem->SF[ii][5])+C4*3.0/(2.0*l23)*(wingMeshFem->SF[ii][3]);
elemFemArr->Hy[0][4]=-(wingMeshFem->SF[ii][1])-(0.5*mypow<T>(S4,2)-0.25*mypow<T>(C4,2))*(wingMeshFem->SF[ii][3])-(0.5*mypow<T>(S6,2)-0.25*mypow<T>(C6,2))*(wingMeshFem->SF[ii][5]);
elemFemArr->Hy[0][5]=(wingMeshFem->SF[ii][3])*(3.0/4.0*S4*C4)+(wingMeshFem->SF[ii][5])*(3.0/4.0*S6*C6);
//
//-C4*3/(2*l23(k))*SF(ii,4)+C5*3/(2*l31(k))*SF(ii,5);
elemFemArr->Hy[0][6]=-C4*3.0/(2.0*l23)*(wingMeshFem->SF[ii][3])+C5*3.0/(2.0*l31)*(wingMeshFem->SF[ii][4]);
elemFemArr->Hy[0][7]=-(wingMeshFem->SF[ii][2])-(0.5*mypow<T>(S4,2)-0.25*mypow<T>(C4,2))*(wingMeshFem->SF[ii][3])-(0.5*mypow<T>(S5,2)-0.25*mypow<T>(C5,2))*(wingMeshFem->SF[ii][4]);
elemFemArr->Hy[0][8]=(wingMeshFem->SF[ii][3])*(3.0/4.0*S4*C4)+(wingMeshFem->SF[ii][4])*(3.0/4.0*S5*C5);

/*
#if DEBUG_ON
printf("\n\nelemFemArr->Hy...\n");
for (int i=0;i<9;i++){
    printf("%f, ",elemFemArr->Hy[0][i]);
}
#endif
*/
elemFemArr->Hy_xsi[0]=-C5*3.0/(2.0*l31)*(wingMeshFem->DxsiSF[ii][4])+C6*3.0/(2.0*l12)*(wingMeshFem->DxsiSF[ii][5]);
elemFemArr->Hy_xsi[1]=-(wingMeshFem->DxsiSF[ii][0])-(0.5*mypow<T>(S5,2)-0.25*mypow<T>(C5,2))*(wingMeshFem->DxsiSF[ii][4])-(0.5*mypow<T>(S6,2)-0.25*mypow<T>(C6,2))*(wingMeshFem->DxsiSF[ii][5]);
elemFemArr->Hy_xsi[2]=(wingMeshFem->DxsiSF[ii][4])*(3.0/4.0*S5*C5)+(wingMeshFem->DxsiSF[ii][5])*(3.0/4.0*S6*C6);
//
elemFemArr->Hy_xsi[3]=-C6*3.0/(2.0*l12)*(wingMeshFem->DxsiSF[ii][5])+C4*3.0/(2.0*l23)*(wingMeshFem->DxsiSF[ii][3]);
elemFemArr->Hy_xsi[4]=-(wingMeshFem->DxsiSF[ii][1])-(0.5*mypow<T>(S4,2)-0.25*mypow<T>(C4,2))*(wingMeshFem->DxsiSF[ii][3])-(0.5*mypow<T>(S6,2)-0.25*mypow<T>(C6,2))*(wingMeshFem->DxsiSF[ii][5]);
elemFemArr->Hy_xsi[5]=(wingMeshFem->DxsiSF[ii][3])*(3.0/4.0*S4*C4)+(wingMeshFem->DxsiSF[ii][5])*(3.0/4.0*S6*C6);
//
elemFemArr->Hy_xsi[6]=-C4*3.0/(2.0*l23)*(wingMeshFem->DxsiSF[ii][3])+C5*3.0/(2.0*l31)*(wingMeshFem->DxsiSF[ii][4]);
elemFemArr->Hy_xsi[7]=-(wingMeshFem->DxsiSF[ii][2])-(0.5*mypow<T>(S4,2)-0.25*mypow<T>(C4,2))*(wingMeshFem->DxsiSF[ii][3])-(0.5*mypow<T>(S5,2)-0.25*mypow<T>(C5,2))*(wingMeshFem->DxsiSF[ii][4]);
elemFemArr->Hy_xsi[8]=(wingMeshFem->DxsiSF[ii][3])*(3.0/4.0*S4*C4)+(wingMeshFem->DxsiSF[ii][4])*(3.0/4.0*S5*C5);

/*
#if DEBUG_ON
printf("\n\nelemFemArr->Hy_xsi...\n");
for (int i=0;i<9;i++){
    printf("%f, ",elemFemArr->Hy_xsi[i]);
}
#endif
*/
elemFemArr->Hy_eta[0]=-C5*3.0/(2.0*l31)*(wingMeshFem->DetaSF[ii][4])+C6*3.0/(2.0*l12)*(wingMeshFem->DetaSF[ii][5]);
elemFemArr->Hy_eta[1]=-(wingMeshFem->DetaSF[ii][0])-(0.5*mypow<T>(S5,2)-0.25*mypow<T>(C5,2))*(wingMeshFem->DetaSF[ii][4])-(0.5*mypow<T>(S6,2)-0.25*mypow<T>(C6,2))*(wingMeshFem->DetaSF[ii][5]);
elemFemArr->Hy_eta[2]=(wingMeshFem->DetaSF[ii][4])*(3.0/4.0*S5*C5)+(wingMeshFem->DetaSF[ii][5])*(3.0/4.0*S6*C6);
//
elemFemArr->Hy_eta[3]=-C6*3.0/(2.0*l12)*(wingMeshFem->DetaSF[ii][5])+C4*3.0/(2.0*l23)*(wingMeshFem->DetaSF[ii][3]);
elemFemArr->Hy_eta[4]=-(wingMeshFem->DetaSF[ii][1])-(0.5*mypow<T>(S4,2)-0.25*mypow<T>(C4,2))*(wingMeshFem->DetaSF[ii][3])-(0.5*mypow<T>(S6,2)-0.25*mypow<T>(C6,2))*(wingMeshFem->DetaSF[ii][5]);
elemFemArr->Hy_eta[5]=(wingMeshFem->DetaSF[ii][3])*(3.0/4.0*S4*C4)+(wingMeshFem->DetaSF[ii][5])*(3.0/4.0*S6*C6);
//
elemFemArr->Hy_eta[6]=-C4*3.0/(2.0*l23)*(wingMeshFem->DetaSF[ii][3])+C5*3.0/(2.0*l31)*(wingMeshFem->DetaSF[ii][4]);
elemFemArr->Hy_eta[7]=-(wingMeshFem->DetaSF[ii][2])-(0.5*mypow<T>(S4,2)-0.25*mypow<T>(C4,2))*(wingMeshFem->DetaSF[ii][3])-(0.5*mypow<T>(S5,2)-0.25*mypow<T>(C5,2))*(wingMeshFem->DetaSF[ii][4]);
elemFemArr->Hy_eta[8]=(wingMeshFem->DetaSF[ii][3])*(3.0/4.0*S4*C4)+(wingMeshFem->DetaSF[ii][4])*(3.0/4.0*S5*C5);

/*    
#if DEBUG_ON
printf("\n\nelemFemArr->Hy_eta...\n");
for (int i=0;i<9;i++){
    printf("%f, ",elemFemArr->Hy_eta[i]);
}
#endif
*/
T y31=wingMeshFem->y31[kk];
T y12=wingMeshFem->y12[kk];
T x31=wingMeshFem->x31[kk];
T x12=wingMeshFem->x12[kk];

for (int i = 0;i<9;i++){
    elemFemArr->Bb[0][i]=1.0/(2.0*wingMeshFem->area[kk])*(y31*(elemFemArr->Hx_xsi[i])+y12*(elemFemArr->Hx_eta[i]));
    elemFemArr->Bb[1][i]=1.0/(2.0*wingMeshFem->area[kk])*(-x31*(elemFemArr->Hy_xsi[i])-x12*(elemFemArr->Hy_eta[i]));
    elemFemArr->Bb[2][i]=1.0/(2.0*wingMeshFem->area[kk])*(-x31*(elemFemArr->Hx_xsi[i])-x12*(elemFemArr->Hx_eta[i])
    +y31*(elemFemArr->Hy_xsi[i])+y12*(elemFemArr->Hy_eta[i]));
}

/*
#if DEBUG_ON
printf("\n\nelemFemArr->Bb...\n");
for (int i=0;i<3;i++){
    for (int j=0;j<9;j++){
        printf("%f, ",elemFemArr->Bb[i][j]);
    }
    printf("\n");
}
#endif
*/
}

template<class T>
void pseudoMassDKT(int ii, int kk, struct triangleDKT<T> *wingMeshFem, struct femArraysDKT<T> *elemFemArr){

T xsi = wingMeshFem->SFm[ii][1];
T eta = wingMeshFem->SFm[ii][2];

T LWtemp[]={(T)1.0, xsi, eta, xsi*eta, mypow<T>(xsi,2.0), mypow<T>(eta,2.0), 
            mypow<T>(xsi,2.0)*eta, mypow<T>(eta,2.0)*xsi, mypow<T>(xsi,3.0), mypow<T>(eta,3.0)}; //[1 x 10]
//T Ltemp[]={(T)1.0, xsi, eta, xsi*eta, mypow<T>(xsi,2.0), mypow<T>(eta,2.0)}; // [1 x 6]

for (int i = 0; i<10; i++){
    elemFemArr->LW[0][i] = LWtemp[i];
}

//printf("\n\nelemFemArr->LW...\n");
//for (int i=0;i<10;i++){
//    printf("%f, ",elemFemArr->LW[0][i]);
//}

}

template<class T>
void CuFEMNum2DWriteDataInBinary(int rows, int cols, T **Usol, int GEN){

    printf("\n    Entering CuFEMNum2DWriteDataInBinary...  (STATIC SOLUTION):");

    FILE *fileOut;
    #if PRECISION_MODE_FEM == 1
	    fileOut = fopen("../c/OUTDATA_FEM_double.bin", "wb"); // w for write, b for binary
    #endif
    #if PRECISION_MODE_FEM == 2
	    fileOut = fopen("../c/OUTDATA_FEM_single.bin", "wb"); // w for write, b for binary
    #endif

    fwrite(&GEN, sizeof(int), 1, fileOut);
    fwrite(&rows, sizeof(int), 1, fileOut);
    fwrite(&cols, sizeof(int), 1, fileOut);

    printf("    rows = %d, cols=%d\n",rows,cols);

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            fwrite(&(Usol[i][j]), sizeof(T), 1, fileOut);
        }
    }

    fclose(fileOut);
    printf("\n    EXITING CuFEMNum2DWriteDataInBinary...");
}

//============================== FEM FUNCTIONS ==================================================//
template<class T>
void CuFEMNum2DWriteMatrix(int rows, int cols, T **K, T **M, T **F){

    FILE *fileOut;
	fileOut = fopen("../c/OUTDATA_FEM_Kglob_Mglob_BCs.bin", "wb"); // w for write, b for binary

    fwrite(&rows, sizeof(int), 1, fileOut);
    fwrite(&cols, sizeof(int), 1, fileOut);

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            fwrite(&(K[i][j]), sizeof(T), 1, fileOut);
        }
    }

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            fwrite(&(M[i][j]), sizeof(T), 1, fileOut);
        }
    }

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < 1; j++){
            fwrite(&(F[i][j]), sizeof(T), 1, fileOut);
        }
    }

    fclose(fileOut);
    printf("\n    Exiting CuFEMNum2DWriteKglobMglobBCs...");

}

//======================DYNAMIC ANALYSIS==============================/

template<class T>
void RayleighDampingCoefs(T *a, T *b){
    
    *a = 0.166328454612080;
    *b = 2.473011484105673*mypow<T>(10.0,-4.0);


}

template<class T>
void createRHS(struct InDataRecFem<T> *inDataFem, 
                struct triangleDKT<T> *wingMeshFem,
                struct femArraysDKT<T> *elemFemArr,
                 T *distrLoad, T **G, int d){

    // Make sure it is initialized with zeros.
    for (int i = 0;i<wingMeshFem->GEN;i++){
        elemFemArr->Fglob[i][0] = 0;
    }

    T lumpedMass[9] = {1, 0, 0, 1, 0, 0, 1, 0, 0};
    T q = 0.0;
    T w3 = inDataFem->omega3;
    //T t = d*inDataFem->dt;
    T t = (d+1)*inDataFem->dt; // to match matlab
    int cntFglob;

    if (inDataFem->LL == 2){
        q = (inDataFem->P_load)*sin(w3*t);
    }

    //printf("\n cos(w3*t)=%f \n",cos(w3*t));
    for (int kk = 0;kk<wingMeshFem->Nelem;kk++){

        if (inDataFem->LL == 3){
            q = distrLoad[kk]*sin(w3*t);
            /*
            if (kk == 0){
                printf("\n    d = %4d, t=%10.4f, sin(w3 t)=%10.4f,\n",d, t, sin(w3*t));
            }
            */
        }
        else{
            printf("Problem with load case! Inside createFglob\n");
            exit(55);
        }

        //overwrites
        for (int i=0;i<9;i++){
            elemFemArr->floc[i][0] = wingMeshFem->area[kk]*q/3.0*lumpedMass[i];
        }
        // version - 2 (EQUIVALENT) floc=P*HW'*floc1;
        //matMatMultiplication2(2, 10, 9, inDataFem.P_load, 1.0, 0.0, elemFemArr.HW, elemFemArr.floc1, elemFemArr.floc); //HW'*floc1

        for (int q=0;q<9;q++){
            cntFglob = wingMeshFem->LM[q][kk]-1;
            elemFemArr->Fglob[cntFglob][0] = elemFemArr->Fglob[cntFglob][0] + elemFemArr->floc[q][0];
            //Fglob(LM(q,kk))=Fglob(LM(q,kk))+floc(q);
        }


        //for (int i = 0;i<9;i++){
        //    elemFemArr->floc[i][0] = 0;
        //}
        /*
        for (int i = 0;i<10;i++){
            //printf("%f, ",elemFemArr.floc1[i][0] );
            elemFemArr->floc1[i][0] = 0;
        }
        */

    }

    for (int i=0;i<wingMeshFem->GEN;i++){
        G[i][d]=elemFemArr->Fglob[i][0];
    }
    
}
