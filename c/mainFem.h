#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* suppress or not execution times (custom profiler) */
#define DEBUG_ON 0 /*allow printf for debugging purposes*/
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


};
//
void CuFEMNum2DReadInData(struct InDataRecFem *inDataFem );
//
void ConnectivityFEM_IEN_ID_LM(struct InDataRecFem *inDataFem, struct triangleDKT *wingMeshFem );

void TriGaussPoints(int Ng, float xw[Ng][3]);

void BendingStiffness(float E, float v, float tx, float BeSt[3][3]);

//--------------------------- 05/04/2023 ADDED
void TrigElCoefsDKT(struct InDataRecFem *inDataFem, struct triangleDKT *wingMeshFem);

// Calculation of shape functions and their derivatives at gauss points on
// the parent element
void LNShapeFunDST(int Ng, float xw[Ng][3], struct triangleDKT *wingMeshFem);

void LNShapeFunMassDST(int Ng, float xw[Ng][3], struct triangleDKT *wingMeshFem);

// ax, ay, bx,by are constant for constant h (independednt of î,ç)
void matrixG(struct triangleDKT *wingMeshFem);
//---------------------------

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

#if DEBUG_ON
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
    


void BendingStiffness(float E, float v, float tx, float BeSt[3][3]){
  
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

#if DEBUG    
    printf("Ng: %d\n", Ng);

    for (int i=0;i<Ng;i++){
        for (int j=0;j<3;j++){
            printf("xw [%d]:%f,",j,xw[i][j]);
        }
        printf("\n");
    }
#endif

    // Initialize matrices that are Ng x 6
    wingMeshFem->SF = (float**)malloc(Ng *sizeof(float)); // pointer array with Ng rows
    wingMeshFem->DxsiSF = (float**)malloc(Ng *sizeof(float)); // pointer array with Ng rows
    wingMeshFem->DetaSF = (float**)malloc(Ng *sizeof(float)); // pointer array with Ng rows
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

#if DEBUG
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

#if DEBUG
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


#if DEBUG    
    printf("Ng: %d\n", Ng);

    for (int i=0;i<Ng;i++){
        for (int j=0;j<3;j++){
            printf("xw [%d]:%f,",j,xw[i][j]);
        }
        printf("\n");
    }
#endif

    // Initialize matrices that are Ng x 3
    wingMeshFem->SFm = (float**)malloc(Ng *sizeof(float)); // pointer array with Ng rows
    wingMeshFem->DxsiSFm = (float**)malloc(Ng *sizeof(float)); // pointer array with Ng rows
    wingMeshFem->DetaSFm = (float**)malloc(Ng *sizeof(float)); // pointer array with Ng rows
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

#if DEBUG
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

    float **GGDST, **GGDKT; //[10 x 10]  
    int M=10,N=10;

    // Initialize matrices that are Ng x 3
    wingMeshFem->GGDST = (float**)malloc(M *sizeof(float)); // pointer array with M=10 rows
    wingMeshFem->GGDKT = (float**)malloc(M *sizeof(float)); // pointer array with M=10 rows
    for (int i=0;i<M;i++){
        wingMeshFem->SFm[i] = (float*)malloc(N *sizeof(float));
        wingMeshFem->DxsiSFm[i] = (float*)malloc(N *sizeof(float));
        wingMeshFem->DetaSFm[i] = (float*)malloc(N *sizeof(float));
    }

    float xsi=0.0;
    float eta=0.0;
    /*
    GGDST(1,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

    xsi=1; eta=0;
    GGDST(2,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

    xsi=0; eta=1;
    GGDST(3,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

    xsi=1/3; eta=1/3;
    GGDST(4,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

    xsi=2/3; eta=1/3;
    GGDST(5,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

    xsi=1/3; eta=2/3;
    GGDST(6,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

    xsi=0; eta=2/3;
    GGDST(7,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

    xsi=0; eta=1/3;
    GGDST(8,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

    xsi=1/3; eta=0;
    GGDST(9,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];

    xsi=2/3; eta=0;
    GGDST(10,:)=[1 xsi eta xsi*eta xsi^2 eta^2 xsi^2*eta eta^2*xsi xsi^3 eta^3];
    */

    /*
    xsi=0; eta=0;
    GGDKT(1,:)=[1 xsi eta xsi*eta xsi^2 eta^2];

    xsi=1; eta=0;
    GGDKT(2,:)=[1 xsi eta xsi*eta xsi^2 eta^2 ];

    xsi=0; eta=1;
    GGDKT(3,:)=[1 xsi eta xsi*eta xsi^2 eta^2 ];

    xsi=1/2; eta=1/2;
    GGDKT(4,:)=[1 xsi eta xsi*eta xsi^2 eta^2 ];

    xsi=0; eta=1/2;
    GGDKT(5,:)=[1 xsi eta xsi*eta xsi^2 eta^2];

    xsi=1/2; eta=0;
    GGDKT(6,:)=[1 xsi eta xsi*eta xsi^2 eta^2 ];
    */

}