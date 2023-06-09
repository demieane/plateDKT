#ifndef MAIN_FEM_HEADER_FILE
    #define MAIN_FEM_HEADER_FILE

    #include<stdio.h>
    #include<stdlib.h> //malloc
    #include<time.h>
    #include<math.h>

    #ifndef DEBUG_ON
        #define DEBUG_ON 0 /*allow printf for debugging purposes*/
    #endif

    #ifndef MODAL_ANALYSIS
        #define MODAL_ANALYSIS 1 /* Find eigenfrequencies */
    #endif

    #ifndef DYNAMIC_ANALYSIS
        #define DYNAMIC_ANALYSIS 1 /* Time-varying load */

        #ifndef TIME_MARCHING_METHOD
            #define TIME_MARCHING_METHOD 1 /* TODO 1. Newmark, 2. Crank-Nicolson */
        #endif

    #endif

    #ifndef PRECISION_MODE_FEM
        #define PRECISION_MODE_FEM 1 /* 1. DOUBLE, 2. SINGLE */
        #if PRECISION_MODE_FEM == 1
            typedef double mytype;
        #endif
        #if PRECISION_MODE_FEM == 2
            typedef float mytype;
        #endif
    #endif

    #ifndef GaussIntegrPoints
        #define GaussIntegrPoints 3 /* Gauss Integration Points */
    #endif


#endif

