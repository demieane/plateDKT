#ifndef MAIN_FEM_HEADER_FILE
    #define MAIN_FEM_HEADER_FILE

    #include "../include/mainFem.h"

    #ifndef PRECISION_MODE_FEM
        #define PRECISION_MODE_FEM 1 /* 1. DOUBLE, 2. SINGLE */
        #if PRECISION_MODE_FEM == 1
            typedef double mytype;
        #endif
    #if PRECISION_MODE_FEM == 2
            typedef float mytype;
        #endif
    #endif


#endif

#include <stdio.h>
#include <stdlib.h> //malloc
#include <math.h> 

#ifndef FUNCMAT
    #define FUNCMAT
    #include "../src/funcMat.cpp" // Functions used to facilitate martix, vector operations in c
#endif

//
template<class T>
void shepard_interp_2d(int nd, T *xd, T *yd, T *zd,
    T *p, int ni, T *xi, T *yi, T *zi);


/*=========================================================================================*/
/* Definition of the functions BELOW */
/*=========================================================================================*/

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
                //w[j] = 1.0 / ( double ) ( nd );
                w[j] = 1.0/(mytype) nd;
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
                w[j] = sqrt ( mypow<mytype> ( (xi[i] - xd[j]), 2.0 )
                            + mypow<mytype> ( (yi[i] - yd[j]), 2.0 ) );
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
                    w[j] = 1.0 / mypow<mytype>( w[j], *p );
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

    free(w);
    printf("    Exiting SHEPARD INTERP. OK.");
}
