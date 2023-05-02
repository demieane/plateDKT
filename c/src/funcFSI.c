#include <stdio.h>
#include <stdlib.h> //malloc
#include <math.h> 

void shepard_interp_2d(int nd, float *xd, float *yd, float *zd,
    float p, int ni, float *xi, float *yi, float *zi);

void shepard_interp_2d(int nd, float *xd, float *yd, float *zd,
    float p, int ni, float *xi, float *yi, float *zi){
/*
  Purpose:

    SHEPARD_INTERP_2D evaluates a 2D Shepard interpolant.

    Modified:

    02 October 2012

  Author:

    John Burkardt

  Reference:

    Donald Shepard,
    A two-dimensional interpolation function for irregularly spaced data,
    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
    ACM, pages 517-524, 1969.

  Parameters:

    Input, int ND, the number of data points.

    Input, float XD[ND], YD[ND], the data points.

    Input, float ZD[ND], the data values.

    Input, float P, the power.

    Input, int NI, the number of interpolation points.

    Input, float XI[NI], YI[NI], the interpolation points.

    Output, float ZI[NI], the interpolated values.

  Comments: 

    All matrices need to be allocated prior to the function call.

*/

    int i,j;
    int z;
    float s = 0.0;
    float dotproc;
    float *w;
    w = ( float * ) malloc ( nd * sizeof ( float ) ); //weights

    // for each point at which the interpolation is intended
    for (i = 0;i<ni;i++){
        if (p == 0.0){
            for (j=0; j<nd;j++){
                w[j] = 1.0/ (float) nd;
            }
        }
        else{
            z = -1.0; //flag
            for (j = 0;j<nd;j++){
                w[j] = sqrt ( pow ( xi[i] - xd[j], 2 ) + pow ( yi[i] - yd[j], 2 ) );
                if ( w[j] == 0.0 ){
                    z = j;
                    break;
                }
            }
            if ( z != -1 ){
                for ( j = 0; j < nd; j++ ){
                    w[j] = 0.0; // initialize
                }
                w[z] = 1.0;
            }
            else{  
                
                for ( j = 0; j < nd; j++ ){
                    w[j] = 1.0 / pow ( w[j], p );
                    s = s + w[j];
                }
                for ( j = 0; j < nd; j++ ){
                    w[j] = w[j] / s;
                }
                s  = 0.0;
            }
        }

        dotproc = 0.0;
        for ( int k = 0;k<nd;k++){
            dotproc = dotproc + w[k]*zd[k];
        }
        zi[i] = dotproc;
    }
    
    free ( w );


}

