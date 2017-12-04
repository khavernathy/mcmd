#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

using namespace std;

int pickRandomAtom(System &system) {
    return floor(getrand()*(double)system.molecules[0].atoms.size());
}

void perturbAtom(System &system, int i) {
    double rand[3] = {getrand()*2-1, getrand()*2-1, getrand()*2-1}; // 3 vals from -1 -> +1
    for (int n=0;n<3;n++)
        system.molecules[0].atoms[i].pos[n] += rand[n];
    return;
}

// Optimize the molecule (ID=0) via MM forcefield(s)
void optimize(System &system) {

    int i;
    // print out all the bonds
    printf("Dynamically-found Bonds Summary: \n============= \n");
    for (i=0; i<system.molecules[0].atoms.size(); i++) {
        printf("Atom %i (UFF: %s) bonds:\n", i, system.molecules[0].atoms[i].UFFlabel.c_str());
        for (std::map<int,double>::iterator it=system.molecules[0].atoms[i].bonds.begin(); it!=system.molecules[0].atoms[i].bonds.end(); ++it)
            std::cout << "    " << system.molecules[0].atoms[i].name.c_str() << "-" << system.molecules[0].atoms[it->first].name.c_str() << " " << it->first << " => " << it->second << '\n';
    }

    // iterate monte-carlo style pertubations until
    // the molecular energy is minimized
    int converged = 0;
    double error_tolerance = 0.00000001;
    int step_limit = 100;
    double Ei,Ef;
    double delta_E;
    double boltzmann;
    double tmp_pos[3] = {0,0,0};
    int randatom;
    int step=0;
    writeXYZ(system, system.constants.output_traj, 0, step, 0, 0);
    printf("Step %i :: Energy = %f; diff = %f kcal/mol; \n", 0,stretch_energy(system), 0.0);

    while (!converged) {
        Ei = stretch_energy(system); // + angle_bend_energy(system);
        //Ei = angle_bend_energy(system);

        // select random atom and perturb it.
        randatom = pickRandomAtom(system);
        for (int n=0;n<3;n++) tmp_pos[n] = system.molecules[0].atoms[randatom].pos[n];
        perturbAtom(system, randatom);

        // get new energy
        Ef = 0;
        Ef += stretch_energy(system);
        //Ef += angle_bend_energy(system);
        delta_E = Ef - Ei;

      //  printf("Ef after = %f\n", Ef);

        //boltzmann = exp(-delta_E/0.00001); // just for now..
        //if (getrand() < boltzmann) {
        if (delta_E < 0.001) { 
           // accept
            step++;
            writeXYZ(system, system.constants.output_traj, 0, step, 0, 0);
            printf("Step %i :: Energy = %f; diff = %f kcal/mol; \n", step,Ef, delta_E);
            if (fabs(delta_E) < error_tolerance) {
                printf("Finished with energy = %f kcal/mol \n", Ef);
                converged=1;
            }
        } else {
            // reject
            for (int n=0;n<3;n++) system.molecules[0].atoms[randatom].pos[n] = tmp_pos[n];
            if (step >= step_limit) {
                printf("Finished with energy = %f kcal/mol \n", Ei); 
                converged=1;
            }
        }
    } // end while loop for convergence

} // end optimize
