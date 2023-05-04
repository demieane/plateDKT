/* The template feature is only available in cpp codes, not in c. */

#include<stdio.h>
#include<float.h>

#ifndef MODE
    #define MODE 1 /* 1. double, 2. single */
#endif

#if MODE == 1
    typedef double mytype;
#endif
#if MODE == 2
    typedef float mytype;
#endif

template<class T>
struct mydata
{
    T var1;
    T var2;
    /* data */
};

int main(){

    mytype a = 1.0/3.0, b=4.2;

    printf("res=%10.10f\n", a+b);
    mydata<double> demy;

    demy.var1 = 1.0;
    demy.var2 = 1.5/3.3333;

    printf("%10.10f\n", demy.var1+demy.var2);



    return 0;
}