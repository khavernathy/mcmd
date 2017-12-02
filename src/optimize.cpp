#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

using namespace std;

int pickRandomAtom(System &system) {
    return floor(getrand()*(double)system.molecules[0].atoms.size());
}

void perturbAtom(System &system, int i, int j) {
    
}

// Optimize the molecule via MM forcefield(s)
void optimize(System &system) {
    printf("hi\n");

    int i;
    // print out all the bonds
    printf("Bond summary: \n============= \n");
    for (i=0; i<system.molecules[0].atoms.size(); i++) {
        printf("Atom %i bonds:\n", i);
        for (std::map<int,double>::iterator it=system.molecules[0].atoms[i].bonds.begin(); it!=system.molecules[0].atoms[i].bonds.end(); ++it)
            std::cout << "    " << system.molecules[0].atoms[i].name.c_str() << "-" << system.molecules[0].atoms[it->first].name.c_str() << " " << it->first << " => " << it->second << '\n';
    }

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
