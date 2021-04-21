#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

using namespace std;

// LENNARD-JONES
// ==================================
// Lorentz-Berthelot
double lj_lb_eps(double e1, double e2) {
    if (e1 != e2) return sqrt(e1*e2);
    else return e1;
}
double lj_lb_sig(double s1, double s2) {
    if (s1 != s2) return 0.5*(s1 + s2);
    else return s1;
}

// W-H and H are not being used now but code is provided for future use.

// Waldman-Hagler
double lj_wh_sig(double s1, double s2) {
    double si6 = s1;
    si6 *= si6*si6;
    si6 *= si6;

    double sj6 = s2;
    sj6 *= sj6*sj6;
    sj6 *= sj6;

    return pow(0.5*(si6+sj6), 1./6.);
}
double lj_wh_eps(double e1, double e2, double s1, double s2) {
    double si3 = s1*s1*s1;
    double sj3 = s2*s2*s2;
    double si6 = si3*si3;
    double sj6 = sj3*sj3;

    return sqrt(e1*e2) * 2.0*si3*sj3/(si6+sj6);
}

// Halgren
double lj_h_sig(double s1, double s2) {
    double si2 = s1*s1;
    double sj2 = s2*s2;
    return (si2*s1 + sj2*s2)/(si2 + sj2);
}
double lj_h_eps(double e1, double e2) {
    return 4*e1*e2 / pow(sqrt(e1) + sqrt(e2), 2);
}



// TANG TOENNIES
// =======================================
double tt_sigma(double sigi, double sigj) {
    return 0.5*(sigi + sigj);
}
double tt_b(double bi, double bj) {
    return 2.0*bi*bj/(bi + bj);
}
double tt_c6(double c6i, double c6j) {
    return sqrt(c6i*c6j)*0.021958709/(3.166811429*0.000001); // Convert H*Bohr^6 to K*Angstrom^6, etc
}
double tt_c8(double c8i, double c8j) {
    return sqrt(c8i*c8j)*0.0061490647/(3.166811429*0.000001); // Convert H*Bohr^6 to K*Angstrom^6, etc
}
double tt_c10(double c6, double c8) {
    if (c6 != 0 && c8 != 0)
        return 49./40.*c8*c8/c6;
    else
        return 0.0;
}


