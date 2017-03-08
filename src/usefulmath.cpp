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


