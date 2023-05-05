#ifndef MAIN_FEM_HEADER_FILE
    #define MAIN_FEM_HEADER_FILE

    #include "../include/mainFem.h"

    #ifndef GaussIntegrPoints
        #define GaussIntegrPoints 3 /* Gauss Integration Points */
    #endif
    

#endif




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
    T **SF, **DxsiSF, **DetaSF; //[Ng x 6]
    T *D2xsiSF, *D2xsietaSF, *D2etaSF; // [1 x 6]
    /* output of LNShapeFunMassDST() */
    T **SFm = NULL;
    T **DxsiSFm = NULL;
    T **DetaSFm = NULL; //[Ng x 3]

    T **GGDST, **GGDKT; //[10 x 10]  
    T **GGin, **GGin2; //inverse of above
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
void LNShapeFunDST(T xw[GaussIntegrPoints][3], struct triangleDKT<T> *wingMeshFem);
//



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
void BendingStiffness(T E, T v, T tx, T **BeSt){
  
    T la = (E*pow(tx,3.0))/(12*(1-pow(v,2)));

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
        wingMeshFem->l23[i]=sqrt(pow(wingMeshFem->x23[i],2)+pow(wingMeshFem->y23[i],2));
        wingMeshFem->l31[i]=sqrt(pow(wingMeshFem->x31[i],2)+pow(wingMeshFem->y31[i],2));
        wingMeshFem->l12[i]=sqrt(pow(wingMeshFem->x12[i],2)+pow(wingMeshFem->y12[i],2));

        wingMeshFem->a4[i]=-wingMeshFem->x23[i]/pow(wingMeshFem->l23[i],2);
        wingMeshFem->a5[i]=-wingMeshFem->x31[i]/pow(wingMeshFem->l31[i],2);
        wingMeshFem->a6[i]=-wingMeshFem->x12[i]/pow(wingMeshFem->l12[i],2);

        wingMeshFem->b4[i]=(3.0/4.0)*(wingMeshFem->x23[i])*(wingMeshFem->y23[i])/pow(wingMeshFem->l23[i],2);
        wingMeshFem->b5[i]=(3.0/4.0)*(wingMeshFem->x31[i])*(wingMeshFem->y31[i])/pow(wingMeshFem->l31[i],2);
        wingMeshFem->b6[i]=(3.0/4.0)*(wingMeshFem->x12[i])*(wingMeshFem->y12[i])/pow(wingMeshFem->l12[i],2);

        wingMeshFem->c4[i]=((1.0/4.0)*pow(wingMeshFem->x23[i],2)-(1.0/2.0)*pow(wingMeshFem->y23[i],2))/pow(wingMeshFem->l23[i],2);
        wingMeshFem->c5[i]=((1.0/4.0)*pow(wingMeshFem->x31[i],2)-(1.0/2.0)*pow(wingMeshFem->y31[i],2))/pow(wingMeshFem->l31[i],2);
        wingMeshFem->c6[i]=((1.0/4.0)*pow(wingMeshFem->x12[i],2)-(1.0/2.0)*pow(wingMeshFem->y12[i],2))/pow(wingMeshFem->l12[i],2);

        wingMeshFem->d4[i]=-wingMeshFem->y23[i]/pow(wingMeshFem->l23[i],2);
        wingMeshFem->d5[i]=-wingMeshFem->y31[i]/pow(wingMeshFem->l31[i],2);
        wingMeshFem->d6[i]=-wingMeshFem->y12[i]/pow(wingMeshFem->l12[i],2);

        wingMeshFem->e4[i]=((1.0/4.0)*pow(wingMeshFem->y23[i],2)-(1.0/2.0)*pow(wingMeshFem->x23[i],2))/pow(wingMeshFem->l23[i],2);
        wingMeshFem->e5[i]=((1.0/4.0)*pow(wingMeshFem->y31[i],2)-(1.0/2.0)*pow(wingMeshFem->x31[i],2))/pow(wingMeshFem->l31[i],2);
        wingMeshFem->e6[i]=((1.0/4.0)*pow(wingMeshFem->y12[i],2)-(1.0/2.0)*pow(wingMeshFem->x12[i],2))/pow(wingMeshFem->l12[i],2);

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
void LNShapeFunDST(T xw[GaussIntegrPoints][3], struct triangleDKT<T> *wingMeshFem){
    // Ng, xw are given data based on which we will fill up some matrices
    int Ng = GaussIntegrPoints;

#if DEBUG_ON    
    printf("    Ng: %d\n", Ng);

    for (int i=0;i<Ng;i++){
        for (int j=0;j<3;j++){
            printf("    xw [%d]:%f,",j,xw[i][j]);
        }
        printf("\n");
    }
#endif

    // Initialize matrices that are Ng x 6
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
            //printf("SF[%d][%d]=%f, ", i,j,wingMeshFem->SF[i][j] );
            //printf("DxsiSF[%d][%d]=%f, ", i,j,wingMeshFem->DxsiSF[i][j] );
            printf("    DetaSF[%d][%d]=%f, ", i,j,wingMeshFem->DetaSF[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");

#endif

    //float *SF, *DxsiSF, *DetaSF; //[Ng x 6]
    //float *D2xsiSF, *D2xsietaSF, *D2etaSF; // [1 x 6]

    wingMeshFem->D2xsiSF = (T*)malloc(Ng *sizeof(T)); // [1 x 6]
    wingMeshFem->D2xsietaSF = (T*)malloc(Ng *sizeof(T)); // [1 x 6]
    wingMeshFem->D2etaSF = (T*)malloc(Ng *sizeof(T)); // [1 x 6]

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

