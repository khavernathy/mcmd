#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>


// PHAST2 REP/DISP ALGORITHM ISN'T GOOD YET
/*
double factorial(double f) {
	if ( f==0 )
		return 1.0;
	return(f * factorial(f - 1));
}

double tang_toennies(double n, double b, double r) { // for use in PHAST2 model (for He)
	double sum=0;
	for (double k=0; k<=n; k++) {
		sum += exp(-b*r)*pow((b*r),k)/(factorial(k));
	}
	return 1.0-sum;
}
*/
