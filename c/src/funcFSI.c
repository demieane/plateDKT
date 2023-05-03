#include <stdio.h>
#include <stdlib.h> //malloc
#include <math.h> 

void shepard_interp_2d(int nd, float *xd, float *yd, float *zd,
    float *p, int ni, float *xi, float *yi, float *zi);

/* single precision plays a role here as it changes the results a bit
compared to matlab */
void shepard_interp_2d(int nd, float *xd, float *yd, float *zd,
    float *p, int ni, float *xi, float *yi, float *zi){

    printf("\n\nSHEPARD INTERP...\n");
    
    printf("nd=%d, ni=%d\n", nd, ni);
    printf("p=%f\n",*p);

    int z;
    float suma, s;
    float dotproc;

    float *w = ( float * ) malloc ( nd * sizeof ( float ) );

    for (int i=0;i<1;i++){
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
            for ( int j = 0; j < 10; j++ ){
                w[j] = sqrt ( pow ( (xi[i] - xd[j]), 2 )
                            + pow ( (yi[i] - yd[j]), 2 ) );
                //printf("w[%d]=%f\n",j,w[j]);
                printf("%f,%f,%f,%f,%f\n",xi[i],xd[j],yi[i],yd[j],w[j]);
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

    printf("\n\nEXITING SHEPARD INTERP...\n");
}