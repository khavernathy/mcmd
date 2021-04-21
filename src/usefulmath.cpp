#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

double dddotprod( double * a, double * b ) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double iidotprod( int * a, int * b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double * crossprod( double * a, double * b) {
    static double output[3];

    output[0] = a[1]*b[2] - a[2]*b[1];
    output[1] = a[2]*b[0] - a[0]*b[2];
    output[2] = a[0]*b[1] - a[1]*b[0];

    return output;
}

// custom erf^-1(x)
// http://stackoverflow.com/questions/27229371/inverse-error-function-in-c
double erfInverse(double x) {

    double tt1, tt2, lnx, sgn;
    sgn = (x < 0) ? -1.0 : 1.0;

    x = (1 -x)*(1 + x);
    lnx = log(x);

    tt1 = 2/(M_PI*0.147) + 0.5 * lnx;
    tt2 = 1/(0.147) * lnx;

    return (sgn*sqrt(-tt1 + sqrt(tt1*tt1 - tt2)));
}

double factorial(double f) {
    if ( f==0 )
        return 1.0;
    return(f * factorial(f - 1));
}


