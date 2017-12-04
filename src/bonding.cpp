/* Douglas M Franz
 * Space group, USF, 2017
 * 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

using namespace std;

/*
 * Essentially everything below comes from the parameters and
 * functional forms prescribed in UFF, via
 * J. Am. Chem. Soc. Vol. 114, No. 25, 1992
 * */

// get the total potential from bond stretches
// via the Morse potential
double stretch_energy(System &system) {
    return 0;
}

// get the total potential from angle bends
// via simple Fourier small cosine expansion
double angle_bend_energy(System &system) {
    return 0;
}

// get the total potential from torsions
// via simple Fourier small cosine expansion
double torsions_energy(System &system) {
    return 0;
}




