#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

using namespace std;

int pickRandomAtom(System &system) {
    return floor(getrand()*(double)system.molecules[0].atoms.size());
}

// Optimize the molecule via MM forcefield(s)
void optimize(System &system) {
    printf("hi\n");

    // iterate monte-carlo style pertubations until
    // the molecular energy is minimized
    int converged = 0;
    double error_tolerance = 0.0001;
    double current_error;
    int randatom;
    while (!converged) {
        randatom = pickRandomAtom(system);
        

        return; // for now.
        
        
        if (current_error < error_tolerance) converged=1;
    } // end while not converged
} // end optimize
