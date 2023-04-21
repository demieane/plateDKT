#ifndef MAIN_FEM_HEADER_FILE
#define MAIN_FEM_HEADER_FILE

/*Place the entire header file contents below within the #ifndef...#endif preprocessor directives */

#include <stdio.h>  /*This form is used for system header files.*/
#include <stdlib.h>
#include <string.h>
#include <math.h>

//#include <clapacke.h>
#include <cblas.h> // use -lblas 
#include <lapack.h> // use -llapack

/* suppress or not execution times (custom profiler) */
#define DEBUG_ON 0 /*allow printf for DEBUG_ONging purposes*/
#ifndef DEBUG_ON
    #define DEBUG_ON 1
#endif

/*=========================================================================================*/
/* Declarations for data structures and functions */
/*=========================================================================================*/
struct InDataRecFem{
    /* UNIT SYSTEM SI */
    float cRoot;
    float span;
    float U;
    float densfluid;
    /*material properties (mass: kg/m3, Young's modulus: Pa, Poisson's ratio, (average)thickness: m)*/
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
    unsigned int sizeBBnodes;
    int *BBnodes; /* id of nodes affected by boundary conditions*/
    unsigned int sizeBdofs;
    int *Bdofs; /* global numbering of dofs affected by boundary conditions*/

    float P_load; // Pa: positive values point towards the positive Z-axis (reverse for ANSYS) 
};

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
    float *xm; // x barycentric coordinate [Nelem]
    float *ym; // y barycentric coordinate [Nelem]

    /* output of TrigElCoefsDKT() */ 
    float *l23, *l31, *l12; // [1 x Nelem]
    float *y12, *y31, *y23;
    float *x12, *x31, *x23;
    float *area;
    float *a4, *a5, *a6, *b4, *b5, *b6, *c4, *c5, *c6;
    float *d4, *d5, *d6, *e4, *e5, *e6;
    float *C4, *C5, *C6, *S4, *S5, *S6;
    /*
    The C language is case-sensitive. This means that all language keywords,
    identifiers, function names, and other variables
    must be entered with consistent letter capitalization. 
    */
    /* output of LNShapeFunDST() */
    float **SF, **DxsiSF, **DetaSF; //[Ng x 6]
    float *D2xsiSF, *D2xsietaSF, *D2etaSF; // [1 x 6]

    /* output of LNShapeFunMassDST() */
    float **SFm, **DxsiSFm, **DetaSFm; //[Ng x 3]

    float **GGDST, **GGDKT; //[10 x 10]  
    float **GGin, **GGin2; //inverse of above
};

struct femArraysDKT{
    /* global intermediate matrices Mg, Kg, Fglob */
    float **Fglob; //[GEN x 1] global
    float **Hm, **HW; //[10 x 9] overwrite massHmDKT()
    float **kloc, **mloc, **floc; //[9 x 9]
    float **Hxx, **Hyy; // [6 x 9] overwrite rotationMass2()
    float *Hx, *Hy; // [1 x 9] from ShapeFunDKT2()
    /*
    NEW
    */
    float *Hx_xsi, *Hx_eta, *Hy_xsi, *Hy_eta; // [1 x 9] from ShapeFunDKT2()

    float **Bb; // [3 x 9] from ShapeFunDKT2()
    float *LW; // [1 x 10] from pseudoMassDKT()
    float *L; // [1 x 6] from pseudoMassDKT()
    float **Mg, **Kg; // [81 x Nelem] pre-assembly matrices
};

//
void CuFEMNum2DReadInData(struct InDataRecFem *inDataFem );
//
void ConnectivityFEM_IEN_ID_LM(struct InDataRecFem *inDataFem, struct triangleDKT *wingMeshFem );

void TriGaussPoints(int Ng, float xw[Ng][3]);

void BendingStiffness(float E, float v, float tx, float **BeSt);

//--------------------------- 05/04/2023 ADDED
void TrigElCoefsDKT(struct InDataRecFem *inDataFem, struct triangleDKT *wingMeshFem);

// Calculation of shape functions and their derivatives at gauss points on
// the parent element
void LNShapeFunDST(int Ng, float xw[Ng][3], struct triangleDKT *wingMeshFem);

void LNShapeFunMassDST(int Ng, float xw[Ng][3], struct triangleDKT *wingMeshFem);

// ax, ay, bx,by are constant for constant h (independednt of î,ç)
void matrixG(struct triangleDKT *wingMeshFem);

// utilities
void assignRowArrayMatrixG_DST(int rowID, float xsi, float eta, float **array);
//
void assignRowArrayMatrixG_DKT(int rowID, float xsi, float eta, float **array);

//---------------------------07/04/2023 ADDED
void squareMatInverse2(int rows, int cols, float **arrIn, float **arrOut);
//
void allocate2Darray(int rows, int cols, float ***arrIn);
//---------------------------11/04/2023 ADDED
void allocate1Darray(int rows, float **arrIn);
//
void massHmDKT(int kk, struct triangleDKT *wingMeshFem, struct femArraysDKT *elemFemArr);
//
void matMatMultiplication2(int optionCalc, int rowsA, int colsA, int colsB, float alpha, float beta, float **arrA, float **arrB, float **arrOut);
//
void rotationMass2(int kk, struct triangleDKT *wingMeshFem, struct femArraysDKT *elemFemArr);
//
void ShapeFunDKT2(int ii, int kk, struct triangleDKT *wingMeshFem, struct femArraysDKT *elemFemArr);
//
void pseudoMassDKT(int ii, int kk, struct triangleDKT *wingMeshFem, struct femArraysDKT *elemFemArr);
//---------------------------21/04/2023 ADDED
void matSum2(float alpha, float beta, int rows, int cols, float **arrA, float **arrB, float **arrOut);
//


/*=========================================================================================*/
/* Definition of the functions follows */
/*=========================================================================================*/
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
    fread(&(inDataFem->sizeBBnodes), sizeof(int) , 1, file);
    inDataFem->BBnodes = (int*)malloc(inDataFem->sizeBBnodes *sizeof(int));
    for (int i=0;i<inDataFem->sizeBBnodes;i++){
        fread(&(inDataFem->BBnodes[i]), sizeof(int), 1, file);
        //printf("i=%d,BBnodes[i]=%d\n", i,inDataFem->BBnodes[i]);
    }
    //
    fread(&(inDataFem->sizeBdofs), sizeof(int) , 1, file);
    inDataFem->Bdofs = (int*)malloc(inDataFem->sizeBdofs *sizeof(int));
    for (int i=0;i<inDataFem->sizeBdofs;i++){
        fread(&(inDataFem->Bdofs[i]), sizeof(int), 1, file);
        //printf("i=%d,Bdofs[i]=%d\n", i,inDataFem->Bdofs[i]);
    }
    //printf("BBnodes = %d, Bdofs=%d\n",inDataFem->sizeBBnodes, inDataFem->sizeBdofs );

    fread(&(inDataFem->P_load), sizeof(float) , 1, file);

    fclose(file);

#if DEBUG_ON_ON
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
    printf("EXITING CuFEMNum2DReadInData...\n\n");
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

    printf("Calculating barycentric coordinates.\n");
    // allocate memory
    wingMeshFem->xm = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->ym = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    
    int IEN_1, IEN_2, IEN_3;    
    //i = wingMeshFem->Nelem-1;
    for (int i=0;i<wingMeshFem->Nelem;i++){ 
        IEN_1 = wingMeshFem->IEN[0][i] - 1;
        IEN_2 = wingMeshFem->IEN[1][i] - 1;
        IEN_3 = wingMeshFem->IEN[2][i] - 1;
        //printf("%d, %d, %d, \n",IEN_1, IEN_2, IEN_3);
        wingMeshFem->xm[i] = (1.0/3.0)*(inDataFem->pp[0][IEN_1] + inDataFem->pp[0][IEN_1] + inDataFem->pp[0][IEN_2]);
        wingMeshFem->ym[i] = (1.0/3.0)*(inDataFem->pp[1][IEN_1] + inDataFem->pp[1][IEN_1] + inDataFem->pp[1][IEN_2]);
        //printf("xm(%d)=%f, ym(%d)=%f\n",i, wingMeshFem->xm[i],i, wingMeshFem->ym[i]);
    }
    
    printf("EXITING ConnectivityFEM_IEN_ID_LM...\n\n");
}


void TriGaussPoints(int Ng, float xw[Ng][3]){
    
    int Mcol=Ng;
    int Ncol=3;
    //float xw[Ng][3];

    if (Ng==1){
        float xw_temp[1][3] = {0.33333333333333, 0.33333333333333, 1.00000000000000};

        for (int i=0;i<Mcol;i++){
            for (int j=0;j<Ncol;j++){
                //printf("xw_temp [%d]:%f,",j,xw_temp[i][j]);
                xw[i][j]=xw_temp[i][j];
                //printf("xw [%d]:%f,",j,xw[i][j]);
            }
        }
    }
    
    if (Ng==3){
        float xw_temp[3][3]={{0.16666666666667, 0.16666666666667, 0.33333333333333},
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
        float xw_temp[4][3]={ {0.33333333333333, 0.33333333333333, -0.56250000000000},
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
        float xw_temp[6][3]={{0.44594849091597, 0.44594849091597, 0.22338158967801},
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

    /* TODO : Add more options for the gauss integration */
    


void BendingStiffness(float E, float v, float tx, float **BeSt){
  
    float la = (E*pow(tx,3.0))/(12*(1-pow(v,2))); // TODO work with matrices [1,Nelem]
    //printf("---> %f, %f, %f, %f\n\n",E,E*pow(tx,3),la, (1.0-v)*la/2.0);
    /*
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            BeSt[i][j] = 0;
            //printf("%f\n", BeSt[i][j]);
        }
    }
    */

    BeSt[0][0] = la;
    BeSt[0][1] = v*la;
    BeSt[1][0] = v*la;
    BeSt[1][1] = la;
    BeSt[2][2] = (1.0-v)*la/2.0; 
    /*
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            printf("Inside %f\n", BeSt[i][j]);
        }
    }    
    */
}

void TrigElCoefsDKT(struct InDataRecFem *inDataFem, struct triangleDKT *wingMeshFem){
    /* memory allocation */
    wingMeshFem->l23 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->l31 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->l12 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    //
    wingMeshFem->y12 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->y31 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->y23 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    //
    wingMeshFem->x12 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->x31 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->x23 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    //
    wingMeshFem->area = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    //
    wingMeshFem->a4 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->a5 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->a6 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    //
    wingMeshFem->b4 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->b5 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->b6 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    //
    wingMeshFem->c4 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->c5 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->c6 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    //
    wingMeshFem->d4 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->d5 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->d6 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    //
    wingMeshFem->e4 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->e5 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->e6 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    //
    wingMeshFem->C4 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->C5 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->C6 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    //
    wingMeshFem->S4 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->S5 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));
    wingMeshFem->S6 = (float*)malloc(wingMeshFem->Nelem *sizeof(float));


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

    printf("EXITING TrigElCoefsDKT...\n\n");

}

void LNShapeFunDST(int Ng, float xw[Ng][3], struct triangleDKT *wingMeshFem){
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

    // Initialize matrices that are Ng x 6
    wingMeshFem->SF = (float**)malloc(Ng *sizeof(float*)); // pointer array with Ng rows
    wingMeshFem->DxsiSF = (float**)malloc(Ng *sizeof(float*)); // pointer array with Ng rows
    wingMeshFem->DetaSF = (float**)malloc(Ng *sizeof(float*)); // pointer array with Ng rows
    for (int i=0;i<Ng;i++){
        wingMeshFem->SF[i] = (float*)malloc(6 *sizeof(float));
        wingMeshFem->DxsiSF[i] = (float*)malloc(6 *sizeof(float));
        wingMeshFem->DetaSF[i] = (float*)malloc(6 *sizeof(float));
    }

    float xg, yg;
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
            printf("DetaSF[%d][%d]=%f, ", i,j,wingMeshFem->DetaSF[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");

#endif
    //float *SF, *DxsiSF, *DetaSF; //[Ng x 6]
    //float *D2xsiSF, *D2xsietaSF, *D2etaSF; // [1 x 6]

    wingMeshFem->D2xsiSF = (float*)malloc(Ng *sizeof(float)); // [1 x 6]
    wingMeshFem->D2xsietaSF = (float*)malloc(Ng *sizeof(float)); // [1 x 6]
    wingMeshFem->D2etaSF = (float*)malloc(Ng *sizeof(float)); // [1 x 6]

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
            printf("D2etaSF[%d]=%f, ", j,wingMeshFem->D2etaSF[j]);
    }
#endif

    printf("EXITING LNShapeFunDST...\n\n");

}

void LNShapeFunMassDST(int Ng, float xw[Ng][3], struct triangleDKT *wingMeshFem){
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
    wingMeshFem->SFm = (float**)malloc(Ng *sizeof(float*)); // pointer array with Ng rows
    wingMeshFem->DxsiSFm = (float**)malloc(Ng *sizeof(float*)); // pointer array with Ng rows
    wingMeshFem->DetaSFm = (float**)malloc(Ng *sizeof(float*)); // pointer array with Ng rows
    for (int i=0;i<Ng;i++){
        wingMeshFem->SFm[i] = (float*)malloc(3 *sizeof(float));
        wingMeshFem->DxsiSFm[i] = (float*)malloc(3 *sizeof(float));
        wingMeshFem->DetaSFm[i] = (float*)malloc(3 *sizeof(float));
    }

    float xg, yg;
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

    printf("EXITING LNShapeFunMassDST...\n\n");


}

void matrixG(struct triangleDKT *wingMeshFem){

    int r = 10, c = 10, i, j;
    float xsi, eta, rowID;
 
    wingMeshFem->GGDST = (float**)malloc(r * sizeof(float*));
    wingMeshFem->GGin = (float**)malloc(r * sizeof(float*));//inverse matrix
    for (i = 0; i < r; i++){
        wingMeshFem->GGDST[i] = (float*)malloc(c * sizeof(float));
        wingMeshFem->GGin[i] = (float*)malloc(c * sizeof(float));
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
    wingMeshFem->GGDKT = (float**)malloc(r * sizeof(float*));
    wingMeshFem->GGin2 = (float**)malloc(r * sizeof(float*));
    for (i = 0; i < r; i++){
        wingMeshFem->GGDKT[i] = (float*)malloc(c * sizeof(float));
        wingMeshFem->GGin2[i] = (float*)malloc(c * sizeof(float));
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

    assignRowArrayMatrixG_DST(rowID = 1.0, xsi=0.0, eta = 0.0, wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST(rowID = 2.0, xsi=1.0, eta = 0.0, wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST(rowID = 3.0, xsi=0.0, eta = 1.0, wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST(rowID = 4.0, xsi=(1.0/3.0), eta = (1.0/3.0), wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST(rowID = 5.0, xsi=(2.0/3.0), eta = (1.0/3.0), wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST(rowID = 6.0, xsi=(1.0/3.0), eta = (2.0/3.0), wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST(rowID = 7.0, xsi=0.0, eta = (2.0/3.0), wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST(rowID = 8.0, xsi=0.0, eta = (1.0/3.0), wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST(rowID = 9.0, xsi=(1.0/3.0), eta = 0.0, wingMeshFem->GGDST);
    assignRowArrayMatrixG_DST(rowID = 10.0, xsi=(2.0/3.0), eta = 0.0, wingMeshFem->GGDST);
    //
    assignRowArrayMatrixG_DKT(rowID = 1.0, xsi=0.0, eta = 0.0, wingMeshFem->GGDKT);
    assignRowArrayMatrixG_DKT(rowID = 2.0, xsi=1.0, eta = 0.0, wingMeshFem->GGDKT);
    assignRowArrayMatrixG_DKT(rowID = 3.0, xsi=0.0, eta = 1.0, wingMeshFem->GGDKT);
    assignRowArrayMatrixG_DKT(rowID = 4.0, xsi=(1.0/2.0), eta = (1.0/2.0), wingMeshFem->GGDKT);
    assignRowArrayMatrixG_DKT(rowID = 5.0, xsi=0.0, eta = (1.0/2.0), wingMeshFem->GGDKT);
    assignRowArrayMatrixG_DKT(rowID = 6.0, xsi=(1.0/2.0), eta = 0.0, wingMeshFem->GGDKT);
}

/* 
COMMENT: The following functions are not dependent on my new datatypes! Therefore,
they can be included without the structure declaration in a separate .c file! That is why in the BEM code
the DATASTRUCTURES are defined in multiple destinations within preprocessor directives.
*/
void assignRowArrayMatrixG_DST(int rowID, float xsi, float eta, float **array){

    float rowMat[10] = {1.0, xsi, eta, xsi*eta, pow(xsi,2), pow(eta,2),
     pow(xsi,2)*eta, pow(eta,2)*xsi, pow(xsi,3), pow(eta,3)};

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

void assignRowArrayMatrixG_DKT(int rowID, float xsi, float eta, float **array){

    float rowMat[6] = {1.0, xsi, eta, xsi*eta, pow(xsi,2), pow(eta,2)};

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


void massHmDKT(int kk, struct triangleDKT *wingMeshFem, struct femArraysDKT *elemFemArr){

    elemFemArr->Hm[0][0] = 1.0;// Hm(1,:)
    elemFemArr->Hm[1][3] = 1.0;// Hm(2,:)
    elemFemArr->Hm[2][6] = 1.0;// Hm(3,:)

    elemFemArr->Hm[3][0] = (1.0/3.0);// Hm(4,:)
    elemFemArr->Hm[3][3] = (1.0/3.0);
    elemFemArr->Hm[3][6] = (1.0/3.0);

    printf("\nPrinting Hm...\n");
    for (int i=0;i<10;i++){
        for (int j=0;j<9;j++){
            printf("%f, ", elemFemArr->Hm[i][j]);
        }
        printf("\n");
    }

    float Hm_5[] = {0.0, 0.0, 0.0,
     (20.0/27.0), 4.0*(wingMeshFem->l23[kk])/27.0*(wingMeshFem->S4[kk]), -4.0*(wingMeshFem->l23[kk])/27.0*(wingMeshFem->C4[kk]),
     (7.0/27.0), -2.0*(wingMeshFem->l23[kk])/27.0*(wingMeshFem->S4[kk]), 2.0*(wingMeshFem->l23[kk])/27.0*(wingMeshFem->C4[kk])};

    float Hm_6[] = {0.0, 0.0, 0.0,
     (7.0/27.0), 2.0*(wingMeshFem->l23[kk])/27.0*(wingMeshFem->S4[kk]), -2.0*(wingMeshFem->l23[kk])/27.0*(wingMeshFem->C4[kk]),
     (20.0/27.0), -4.0*(wingMeshFem->l23[kk])/27.0*(wingMeshFem->S4[kk]), 4.0*(wingMeshFem->l23[kk])/27.0*(wingMeshFem->C4[kk])};

    float Hm_7[] = {(7.0/27.0), -2.0*wingMeshFem->l31[kk]/27.0*wingMeshFem->S5[kk], 2.0*wingMeshFem->l31[kk]/27.0*wingMeshFem->C5[kk],
     0.0, 0.0, 0.0,
     (20.0/27.0), 4.0*wingMeshFem->l31[kk]/27.0*wingMeshFem->S5[kk], -4.0*wingMeshFem->l31[kk]/27.0*wingMeshFem->C5[kk]};

    float Hm_8[] = {(20.0/27.0), -4.0*wingMeshFem->l31[kk]/27.0*wingMeshFem->S5[kk], 4.0*wingMeshFem->l31[kk]/27.0*wingMeshFem->C5[kk],
     0.0, 0.0, 0.0,
     (7.0/27.0), 2.0*wingMeshFem->l31[kk]/27.0*wingMeshFem->S5[kk], -2.0*wingMeshFem->l31[kk]/27.0*wingMeshFem->C5[kk]};

    float Hm_9[] = {(20.0/27.0), 4.0*wingMeshFem->l12[kk]/27.0*wingMeshFem->S6[kk], -4.0*wingMeshFem->l12[kk]/27.0*wingMeshFem->C6[kk],
     (7.0/27.0), -2.0*wingMeshFem->l12[kk]/27.0*wingMeshFem->S6[kk], 2.0*wingMeshFem->l12[kk]/27.0*wingMeshFem->C6[kk],
      0.0, 0.0, 0.0};

    float Hm_10[] = {(7.0/27.0), 2.0*wingMeshFem->l12[kk]/27.0*wingMeshFem->S6[kk], -2.0*wingMeshFem->l12[kk]/27.0*wingMeshFem->C6[kk],
     (20.0/27.0), -4.0*wingMeshFem->l12[kk]/27.0*wingMeshFem->S6[kk], 4.0*wingMeshFem->l12[kk]/27.0*wingMeshFem->C6[kk],
      0.0, 0.0, 0.0};

    for (int j=0; j<9; j++){
        elemFemArr->Hm[4][j]=Hm_5[j];
        elemFemArr->Hm[5][j]=Hm_6[j];
        elemFemArr->Hm[6][j]=Hm_7[j];
        elemFemArr->Hm[7][j]=Hm_8[j];
        elemFemArr->Hm[8][j]=Hm_9[j];
        elemFemArr->Hm[9][j]=Hm_10[j];
    }

    printf("\nPrinting Hm...\n");
    for (int i=0;i<10;i++){
        for (int j=0;j<9;j++){
            printf("%f, ", elemFemArr->Hm[i][j]);
        }
        printf("\n");
    }

    int myVariable = -5;
    /*
    with op( A ) an m by k matrix, op( B )  a  k by n matrix
    
    HW=GGin*Hm;  % [10 x 10] x[10 x 9] 
    */
    int M = 10, N = 10, K = 9;
    int optionCalc = 1; //without transpose
    float alpha = 1.0, beta = 0.0;
    matMatMultiplication2(optionCalc, M, N, K, alpha, beta, wingMeshFem->GGin, elemFemArr->Hm, elemFemArr->HW);

#if DEBUG_ON
    printf("\n EXITING massHmDKT...");
#endif
}

void rotationMass2(int kk, struct triangleDKT *wingMeshFem, struct femArraysDKT *elemFemArr){

// Hxx, Hyy [6 x 9]
float A4 = wingMeshFem->a4[kk];
float A6 = wingMeshFem->a6[kk];
float A5 = wingMeshFem->a5[kk];
//
float B4 = wingMeshFem->b4[kk];
float B5 = wingMeshFem->b5[kk];
float B6 = wingMeshFem->b6[kk];
//
float C4 = wingMeshFem->c4[kk];
float C5 = wingMeshFem->c5[kk];
float C6 = wingMeshFem->c6[kk];
//
float D4 = wingMeshFem->d4[kk];
float D5 = wingMeshFem->d5[kk];
float D6 = wingMeshFem->d6[kk];
//
float E4 = wingMeshFem->e4[kk];
float E5 = wingMeshFem->e5[kk];
float E6 = wingMeshFem->e6[kk];
//
float Hxx1[]={0.0, 6.0*A6, -6.0*A5, 6.0*A5-6.0*A6, -6.0*A6, 6.0*A5}; //columns of Hxx
float Hxx2[]={0.0, 4.0*B6, 4.0*B5, -4.0*B5-4.0*B6, -4.0*B6, -4*B5};  
float Hxx3[]={1.0, -4.0*C6-3, -3.0-4.0*C5, 4.0+4.0*C5+4.0*C6, 4.0*C6+2.0, 4.0*C5+2.0};
float Hxx4[]={0.0, -6.0*A6, 0.0, 6.0*A4+6.0*A6, 6.0*A6, 0.0};
float Hxx5[]={0.0, 4.0*B6, 0.0, 4.0*B4-4.0*B6, -4.0*B6, 0.0};
float Hxx6[]={0.0, -1.0-4.0*C6, 0.0, 4.0*C6-4.0*C4, 2.0+4.0*C6, 0.0};
float Hxx7[]={0.0, 0.0, 6.0*A5, -6.0*A5-6.0*A4, 0.0, -6.0*A5};
float Hxx8[]={0.0, 0.0, 4.0*B5, 4.0*B4-4.0*B5, 0.0, -4.0*B5};
float Hxx9[]={0.0, 0.0, -4.0*C5-1.0, 4.0*C5-4.0*C4, 0.0, 4.0*C5+2.0};

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
float Hyy1[]={0.0, 6.0*D6, -6.0*D5, 6.0*D5-6.0*D6, -6.0*D6, 6.0*D5}; //column
float Hyy2[]={-1.0, 4,0*E6+3.0, 4.0*E5+3.0, -4.0*E5-4.0*E6-4.0, -4.0*E6-2.0, -4.0*E5-2.0};
float Hyy3[]={0.0, -4.0*B6, -4.0*B5, 4.0*B5+4.0*B6, 4.0*B6, 4.0*B5};
float Hyy4[]={0.0, -6.0*D6, 0.0, 6.0*D4+6.0*D6, 6.0*D6, 0.0};
float Hyy5[]={0.0, 1.0+4.0*E6, 0.0, 4.0*E4-4.0*E6, -4.0*E6-2.0, 0.0};
float Hyy6[]={0.0, -4.0*B6, 0.0, -4.0*B4+4.0*B6, 4.0*B6, 0.0};
float Hyy7[]={0.0, 0.0, 6.0*D5, -6.0*D5-6.0*D4, 0.0, -6.0*D5};
float Hyy8[]={0.0, 0.0, 1.0+4.0*E5, 4.0*E4-4.0*E5, 0.0, -4.0*E5-2.0};
float Hyy9[]={0.0, 0.0, -4.0*B5, -4.0*B4+4.0*B5, 0.0, 4.0*B5};

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

void ShapeFunDKT2(int ii, int kk, struct triangleDKT *wingMeshFem, struct femArraysDKT *elemFemArr){

/* output of LNShapeFunDST() */
//float **SF, **DxsiSF, **DetaSF; //[Ng x 6]
//float *D2xsiSF, *D2xsietaSF, *D2etaSF; // [1 x 6]

/* output of LNShapeFunMassDST() */
//float **SFm, **DxsiSFm, **DetaSFm; //[Ng x 3]

float S4=wingMeshFem->S4[kk];
float S5=wingMeshFem->S5[kk];
float S6=wingMeshFem->S6[kk];
float C4=wingMeshFem->C4[kk];
float C5=wingMeshFem->C5[kk];
float C6=wingMeshFem->C6[kk];
//
float l31=wingMeshFem->l31[kk];
float l12=wingMeshFem->l12[kk];
float l23=wingMeshFem->l23[kk];

//size Hx [1 x 9]
elemFemArr->Hx[0]=S5*3.0/(2.0*l31)*(wingMeshFem->SF[ii][4])-S6*3.0/(2.0*l12)*(wingMeshFem->SF[ii][5]);
elemFemArr->Hx[1]=-(wingMeshFem->SF[ii][4])*(3.0/4.0*S5*C5)-(wingMeshFem->SF[ii][5])*(3.0/4.0*S6*C6);
elemFemArr->Hx[2]=(wingMeshFem->SF[ii][0])+(0.5*pow(C5,2)-0.25*pow(S5,2))*(wingMeshFem->SF[ii][4])+(0.5*pow(C6,2)-0.25*pow(S6,2))*(wingMeshFem->SF[ii][5]);
//
elemFemArr->Hx[3]=S6*3.0/(2.0*l12)*(wingMeshFem->SF[ii][5])-S4*3.0/(2.0*l23)*(wingMeshFem->SF[ii][3]);
elemFemArr->Hx[4]=-(wingMeshFem->SF[ii][3])*(3.0/4.0*S4*C4)-(wingMeshFem->SF[ii][5])*(3.0/4.0*S6*C6);
elemFemArr->Hx[5]=(wingMeshFem->SF[ii][1])+(0.5*pow(C4,2)-0.25*pow(S4,2))*(wingMeshFem->SF[ii][3])+(0.5*pow(C6,2)-0.25*pow(S6,2))*(wingMeshFem->SF[ii][5]);
//
elemFemArr->Hx[6]=S4*3.0/(2.0*l23)*(wingMeshFem->SF[ii][3])-S5*3.0/(2.0*l31)*(wingMeshFem->SF[ii][4]);
elemFemArr->Hx[7]=-(wingMeshFem->SF[ii][3])*(3.0/4.0*S4*C4)-(wingMeshFem->SF[ii][4])*(3.0/4.0*S5*C5);
elemFemArr->Hx[8]=(wingMeshFem->SF[ii][2])+(0.5*pow(C4,2)-0.25*pow(S4,2))*(wingMeshFem->SF[ii][3])+(0.5*pow(C5,2)-0.25*pow(S5,2))*(wingMeshFem->SF[ii][4]);

#if DEBUG_ON
printf("ShapeFunDKT2 ...\n");
printf("\nelemFemArr->Hx...\n");
for (int i=0;i<9;i++){
    printf("%f, ",elemFemArr->Hx[i]);
}
#endif

//size Hx_xsi [1 x 9]
elemFemArr->Hx_xsi[0]=S5*3.0/(2.0*l31)*(wingMeshFem->DxsiSF[ii][4])-S6*3.0/(2.0*l12)*(wingMeshFem->DxsiSF[ii][5]);
elemFemArr->Hx_xsi[1]=-(wingMeshFem->DxsiSF[ii][4])*(3.0/4.0*S5*C5)-(wingMeshFem->DxsiSF[ii][5])*(3.0/4.0*S6*C6);
elemFemArr->Hx_xsi[2]=(wingMeshFem->DxsiSF[ii][0])+(0.5*pow(C5,2)-0.25*pow(S5,2))*(wingMeshFem->DxsiSF[ii][4])+(0.5*pow(C6,2)-0.25*pow(S6,2))*(wingMeshFem->DxsiSF[ii][5]);
//
elemFemArr->Hx_xsi[3]=S6*3.0/(2.0*l12)*(wingMeshFem->DxsiSF[ii][5])-S4*3.0/(2.0*l23)*(wingMeshFem->DxsiSF[ii][3]);
elemFemArr->Hx_xsi[4]=-(wingMeshFem->DxsiSF[ii][3])*(3.0/4.0*S4*C4)-(wingMeshFem->DxsiSF[ii][5])*(3.0/4.0*S6*C6);
elemFemArr->Hx_xsi[5]=(wingMeshFem->DxsiSF[ii][1])+(0.5*pow(C4,2)-0.25*pow(S4,2))*(wingMeshFem->DxsiSF[ii][3])+(0.5*pow(C6,2)-0.25*pow(S6,2))*(wingMeshFem->DxsiSF[ii][5]);
//
elemFemArr->Hx_xsi[6]=S4*3.0/(2.0*l23)*(wingMeshFem->DxsiSF[ii][3])-S5*3.0/(2.0*l31)*(wingMeshFem->DxsiSF[ii][4]);
elemFemArr->Hx_xsi[7]=-(wingMeshFem->DxsiSF[ii][3])*(3.0/4.0*S4*C4)-(wingMeshFem->DxsiSF[ii][4])*(3.0/4.0*S5*C5);
elemFemArr->Hx_xsi[8]=(wingMeshFem->DxsiSF[ii][2])+(0.5*pow(C4,2)-0.25*pow(S4,2))*(wingMeshFem->DxsiSF[ii][3])+(0.5*pow(C5,2)-0.25*pow(S5,2))*(wingMeshFem->DxsiSF[ii][4]);

#if DEBUG_ON
printf("\n\nelemFemArr->Hx_xsi...\n");
for (int i=0;i<9;i++){
    printf("%f, ",elemFemArr->Hx_xsi[i]);
}
#endif

elemFemArr->Hx_eta[0]=S5*3.0/(2.0*l31)*(wingMeshFem->DetaSF[ii][4])-S6*3.0/(2.0*l12)*(wingMeshFem->DetaSF[ii][5]);
elemFemArr->Hx_eta[1]=-(wingMeshFem->DetaSF[ii][4])*(3.0/4.0*S5*C5)-(wingMeshFem->DetaSF[ii][5])*(3.0/4.0*S6*C6);
elemFemArr->Hx_eta[2]=(wingMeshFem->DetaSF[ii][0])+(0.5*pow(C5,2)-0.25*pow(S5,2))*(wingMeshFem->DetaSF[ii][4])+(0.5*pow(C6,2)-0.25*pow(S6,2))*(wingMeshFem->DetaSF[ii][5]);
//
elemFemArr->Hx_eta[3]=S6*3.0/(2.0*l12)*(wingMeshFem->DetaSF[ii][5])-S4*3.0/(2.0*l23)*(wingMeshFem->DetaSF[ii][3]);
elemFemArr->Hx_eta[4]=-(wingMeshFem->DetaSF[ii][3])*(3.0/4.0*S4*C4)-(wingMeshFem->DetaSF[ii][5])*(3.0/4.0*S6*C6);
elemFemArr->Hx_eta[5]=(wingMeshFem->DetaSF[ii][1])+(0.5*pow(C4,2)-0.25*pow(S4,2))*(wingMeshFem->DetaSF[ii][3])+(0.5*pow(C6,2)-0.25*pow(S6,2))*(wingMeshFem->DetaSF[ii][5]);
//
elemFemArr->Hx_eta[6]=S4*3.0/(2.0*l23)*(wingMeshFem->DetaSF[ii][3])-S5*3.0/(2.0*l31)*(wingMeshFem->DetaSF[ii][4]);
elemFemArr->Hx_eta[7]=-(wingMeshFem->DetaSF[ii][3])*(3.0/4.0*S4*C4)-(wingMeshFem->DetaSF[ii][4])*(3.0/4.0*S5*C5);
elemFemArr->Hx_eta[8]=(wingMeshFem->DetaSF[ii][2])+(0.5*pow(C4,2)-0.25*pow(S4,2))*(wingMeshFem->DetaSF[ii][3])+(0.5*pow(C5,2)-0.25*pow(S5,2))*(wingMeshFem->DetaSF[ii][4]);

#if DEBUG_ON
printf("\n\nelemFemArr->Hx_eta...\n");
for (int i=0;i<9;i++){
    printf("%f, ",elemFemArr->Hx_eta[i]);
}
#endif

elemFemArr->Hy[0]=-C5*3.0/(2.0*l31)*(wingMeshFem->SF[ii][4])+C6*3.0/(2.0*l12)*(wingMeshFem->SF[ii][5]);
elemFemArr->Hy[1]=-(wingMeshFem->SF[ii][0])-(0.5*pow(S5,2)-0.25*pow(C5,2))*(wingMeshFem->SF[ii][4])-(0.5*pow(S6,2)-0.25*pow(C6,2))*(wingMeshFem->SF[ii][5]);
elemFemArr->Hy[2]=(wingMeshFem->SF[ii][4])*(3.0/4.0*S5*C5)+(wingMeshFem->SF[ii][5])*(3.0/4.0*S6*C6);
//
elemFemArr->Hy[3]=-C6*3/(2*l12)*(wingMeshFem->SF[ii][5])+C4*3.0/(2.0*l23)*(wingMeshFem->SF[ii][3]);
elemFemArr->Hy[4]=-(wingMeshFem->SF[ii][1])-(0.5*pow(S4,2)-0.25*pow(C4,2))*(wingMeshFem->SF[ii][3])-(0.5*pow(S6,2)-0.25*pow(C6,2))*(wingMeshFem->SF[ii][5]);
elemFemArr->Hy[5]=(wingMeshFem->SF[ii][3])*(3.0/4.0*S4*C4)+(wingMeshFem->SF[ii][5])*(3.0/4.0*S6*C6);
//
elemFemArr->Hy[6]=-C4*3.0/(2.0*l23)*(wingMeshFem->SF[ii][3])+C5*3.0/(2.0*l31)*(wingMeshFem->SF[ii][5]);
elemFemArr->Hy[7]=-(wingMeshFem->SF[ii][2])-(0.5*pow(S4,2)-0.25*pow(C4,2))*(wingMeshFem->SF[ii][3])-(0.5*pow(S5,2)-0.25*pow(C5,2))*(wingMeshFem->SF[ii][4]);
elemFemArr->Hy[8]=(wingMeshFem->SF[ii][3])*(3.0/4.0*S4*C4)+(wingMeshFem->SF[ii][4])*(3.0/4.0*S5*C5);

#if DEBUG_ON
printf("\n\nelemFemArr->Hy...\n");
for (int i=0;i<9;i++){
    printf("%f, ",elemFemArr->Hy[i]);
}
#endif

elemFemArr->Hy_xsi[0]=-C5*3.0/(2.0*l31)*(wingMeshFem->DxsiSF[ii][4])+C6*3.0/(2.0*l12)*(wingMeshFem->DxsiSF[ii][5]);
elemFemArr->Hy_xsi[1]=-(wingMeshFem->DxsiSF[ii][0])-(0.5*pow(S5,2)-0.25*pow(C5,2))*(wingMeshFem->DxsiSF[ii][4])-(0.5*pow(S6,2)-0.25*pow(C6,2))*(wingMeshFem->DxsiSF[ii][5]);
elemFemArr->Hy_xsi[2]=(wingMeshFem->DxsiSF[ii][4])*(3.0/4.0*S5*C5)+(wingMeshFem->DxsiSF[ii][5])*(3.0/4.0*S6*C6);
//
elemFemArr->Hy_xsi[3]=-C6*3.0/(2.0*l12)*(wingMeshFem->DxsiSF[ii][5])+C4*3.0/(2.0*l23)*(wingMeshFem->DxsiSF[ii][3]);
elemFemArr->Hy_xsi[4]=-(wingMeshFem->DxsiSF[ii][1])-(0.5*pow(S4,2)-0.25*pow(C4,2))*(wingMeshFem->DxsiSF[ii][3])-(0.5*pow(S6,2)-0.25*pow(C6,2))*(wingMeshFem->DxsiSF[ii][5]);
elemFemArr->Hy_xsi[5]=(wingMeshFem->DxsiSF[ii][3])*(3.0/4.0*S4*C4)+(wingMeshFem->DxsiSF[ii][5])*(3.0/4.0*S6*C6);
//
elemFemArr->Hy_xsi[6]=-C4*3.0/(2.0*l23)*(wingMeshFem->DxsiSF[ii][3])+C5*3.0/(2.0*l31)*(wingMeshFem->DxsiSF[ii][4]);
elemFemArr->Hy_xsi[7]=-(wingMeshFem->DxsiSF[ii][2])-(0.5*pow(S4,2)-0.25*pow(C4,2))*(wingMeshFem->DxsiSF[ii][3])-(0.5*pow(S5,2)-0.25*pow(C5,2))*(wingMeshFem->DxsiSF[ii][4]);
elemFemArr->Hy_xsi[8]=(wingMeshFem->DxsiSF[ii][3])*(3.0/4.0*S4*C4)+(wingMeshFem->DxsiSF[ii][4])*(3.0/4.0*S5*C5);

#if DEBUG_ON
printf("\n\nelemFemArr->Hy_xsi...\n");
for (int i=0;i<9;i++){
    printf("%f, ",elemFemArr->Hy_xsi[i]);
}
#endif

elemFemArr->Hy_eta[0]=-C5*3.0/(2.0*l31)*(wingMeshFem->DetaSF[ii][4])+C6*3.0/(2.0*l12)*(wingMeshFem->DetaSF[ii][5]);
elemFemArr->Hy_eta[1]=-(wingMeshFem->DetaSF[ii][0])-(0.5*pow(S5,2)-0.25*pow(C5,2))*(wingMeshFem->DetaSF[ii][4])-(0.5*pow(S6,2)-0.25*pow(C6,2))*(wingMeshFem->DetaSF[ii][5]);
elemFemArr->Hy_eta[2]=(wingMeshFem->DetaSF[ii][4])*(3.0/4.0*S5*C5)+(wingMeshFem->DetaSF[ii][5])*(3.0/4.0*S6*C6);
//
elemFemArr->Hy_eta[3]=-C6*3.0/(2.0*l12)*(wingMeshFem->DetaSF[ii][5])+C4*3.0/(2.0*l23)*(wingMeshFem->DetaSF[ii][3]);
elemFemArr->Hy_eta[4]=-(wingMeshFem->DetaSF[ii][1])-(0.5*pow(S4,2)-0.25*pow(C4,2))*(wingMeshFem->DetaSF[ii][3])-(0.5*pow(S6,2)-0.25*pow(C6,2))*(wingMeshFem->DetaSF[ii][5]);
elemFemArr->Hy_eta[5]=(wingMeshFem->DetaSF[ii][3])*(3.0/4.0*S4*C4)+(wingMeshFem->DetaSF[ii][5])*(3.0/4.0*S6*C6);
//
elemFemArr->Hy_eta[6]=-C4*3.0/(2.0*l23)*(wingMeshFem->DetaSF[ii][3])+C5*3.0/(2.0*l31)*(wingMeshFem->DetaSF[ii][4]);
elemFemArr->Hy_eta[7]=-(wingMeshFem->DetaSF[ii][2])-(0.5*pow(S4,2)-0.25*pow(C4,2))*(wingMeshFem->DetaSF[ii][3])-(0.5*pow(S5,2)-0.25*pow(C5,2))*(wingMeshFem->DetaSF[ii][4]);
elemFemArr->Hy_eta[8]=(wingMeshFem->DetaSF[ii][3])*(3.0/4.0*S4*C4)+(wingMeshFem->DetaSF[ii][4])*(3.0/4.0*S5*C5);
    
#if DEBUG_ON
printf("\n\nelemFemArr->Hy_eta...\n");
for (int i=0;i<9;i++){
    printf("%f, ",elemFemArr->Hy_eta[i]);
}
#endif

float y31=wingMeshFem->y31[kk];
float y12=wingMeshFem->y12[kk];
float x31=wingMeshFem->x31[kk];
float x12=wingMeshFem->x12[kk];

for (int i = 0;i<9;i++){
    elemFemArr->Bb[0][i]=1.0/(2.0*wingMeshFem->area[kk])*(y31*(elemFemArr->Hx_xsi[i])+y12*(elemFemArr->Hx_eta[i]));
    elemFemArr->Bb[1][i]=1.0/(2.0*wingMeshFem->area[kk])*(-x31*(elemFemArr->Hy_xsi[i])-x12*(elemFemArr->Hy_eta[i]));
    elemFemArr->Bb[2][i]=1.0/(2.0*wingMeshFem->area[kk])*(-x31*(elemFemArr->Hx_xsi[i])-x12*(elemFemArr->Hx_eta[i])
    +y31*(elemFemArr->Hy_xsi[i])+y12*(elemFemArr->Hy_eta[i]));
}

#if DEBUG_ON
printf("\n\nelemFemArr->Bb...\n");
for (int i=0;i<3;i++){
    for (int j=0;j<9;j++){
        printf("%f, ",elemFemArr->Bb[i][j]);
    }
    printf("\n");
}
#endif
}


void pseudoMassDKT(int ii, int kk, struct triangleDKT *wingMeshFem, struct femArraysDKT *elemFemArr){

float xsi = wingMeshFem->SFm[ii][1];
float eta = wingMeshFem->SFm[ii][2];

float LWtemp[]={1.0, xsi, eta, xsi*eta, pow(xsi,2), pow(eta,2), pow(xsi,2)*eta, pow(eta,2)*xsi, pow(xsi,3), pow(eta,3)}; //[1 x 10]
float Ltemp[]={1.0, xsi, eta, xsi*eta, pow(xsi,2), pow(eta,2)}; // [1 x 6]

for (int i = 0; i<10; i++){
    elemFemArr->LW[i] = LWtemp[i];
}

printf("\n\nelemFemArr->LW...\n");
for (int i=0;i<10;i++){
    printf("%f, ",elemFemArr->LW[i]);
}

}


#endif