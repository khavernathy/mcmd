#include <string>
#include <algorithm>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

double * rotatePoint(System &system, double x, double y, double z, int plane, double angle) {

    double finalx, finaly, finalz;
    // the function takes an ANGLE 0->360, trig functions take rads
    angle = angle*M_PI/180.0;
    if (plane == 0) {
        finalx = x;
        finaly = y*cos(angle) -z*sin(angle);
        finalz = y*sin(angle) + z*cos(angle);
    } else if (plane == 1) {
        finalx = x*cos(angle) + z*sin(angle);
        finaly = y;
        finalz = -x*sin(angle) + z*cos(angle);
    } else if (plane == 2) {
        finalx = x*cos(angle) - y*sin(angle);
        finaly = x*sin(angle) + y*cos(angle);
        finalz = z;
    }

    static double output[3];
    output[0] = finalx;
    output[1] = finaly;
    output[2] = finalz;
    return output;
}

double * rotatePointRadians(System &system, double x, double y, double z, int plane, double angle) {

    double finalx, finaly, finalz;
    if (plane == 0) {
        finalx = x;
        finaly = y*cos(angle) -z*sin(angle);
        finalz = y*sin(angle) + z*cos(angle);
    } else if (plane == 1) {
        finalx = x*cos(angle) + z*sin(angle);
        finaly = y;
        finalz = -x*sin(angle) + z*cos(angle);
    } else if (plane == 2) {
        finalx = x*cos(angle) - y*sin(angle);
        finaly = x*sin(angle) + y*cos(angle);
        finalz = z;
    }

    static double output[3];
    output[0] = finalx;
    output[1] = finaly;
    output[2] = finalz;
    return output;
}

