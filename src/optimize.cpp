#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

using namespace std;

int pickRandomAtom(System &system) {
    return floor(getrand()*(double)system.molecules[0].atoms.size());
}

void perturbAtom(System &system, int i, double mf) {
    double rand[3] = {mf*(getrand()*2-1), mf*(getrand()*2-1), mf*(getrand()*2-1)}; // 3 vals from (-1 -> +1)*movefactor
    for (int n=0;n<3;n++)
        system.molecules[0].atoms[i].pos[n] += rand[n];
    return;
}

double move_factor(double energy, int N) {
    // a simple scalar modifier to the explorable volume for optimization
    // the move_factor scales down to 0 as energy approaches 0
    // thus we save some time optimizing by not trying dumb (big) moves 
    // when energy is very low.
    // normalized by 10 kcal/mol per atom
    return 1.0 -exp(-energy*energy / (10.0*N*N)); // reverse bell-curve
}

void outputEnergies(System &system, int step, double Ef, double delta_E) {
    printf("==============================================================\n");
    printf("Step %i\nEnergy =         %f kcal/mol; \u0394E = %f kcal/mol; \n\n", step,Ef, delta_E);
    printf("Bonds =          %f kcal/mol\nAngle-bends =    %f kcal/mol\nDihedrals =      %f kcal/mol\nIntramolec. LJ = %f kcal/mol\n",
          system.stats.Ustretch.value,
          system.stats.Uangles.value,
          system.stats.Udihedrals.value,
          system.stats.UintraLJ.value); 
}

// Optimize the molecule (ID=0) via MM forcefield(s)
void optimize(System &system) {

    int i;
    // print out all the bonds
    printf("================================================================================\n");
    printf("Dynamically-found Bonds Summary:\n");
    printf("================================================================================\n");
    printf("bond-id :: mol-id :: atom1 :: atom2 :: elements\n");
    for (int n=0; n<system.constants.uniqueBonds.size(); n++) {
        //printf("Atom %i (UFF: %s)\n", i, system.molecules[0].atoms[i].UFFlabel.c_str());
                    int mol=system.constants.uniqueBonds[n].mol;
                    int atom1=system.constants.uniqueBonds[n].atom1;
                    int atom2=system.constants.uniqueBonds[n].atom2;
            printf("%7i :: %6i :: %5i :: %5i :: %4s%1s%4s\n",
                    n,
                    mol,
                    atom1,
                    atom2,
                    system.molecules[mol].atoms[atom1].name.c_str(),
                    "-",
                    system.molecules[mol].atoms[atom2].name.c_str()
                    );
    }
    // and angles
    printf("================================================================================\n");
    printf("Dynamically-found Angles Summary:\n");
    printf("================================================================================\n");
    printf("angle-id :: mol-id :: atom1 :: atom2 :: atom3 :: elements\n");
    for (int n=0;n<system.constants.uniqueAngles.size();n++) {
        int mol = system.constants.uniqueAngles[n].mol;
        int atom1= system.constants.uniqueAngles[n].atom1;
        int atom2 = system.constants.uniqueAngles[n].atom2;
        int atom3 = system.constants.uniqueAngles[n].atom3;
        printf("%8i :: %6i :: %5i :: %5i :: %5i :: %4s%1s%4s%1s%4s\n", n, 
                mol,
                atom1,
                atom2,
                atom3,
                system.molecules[mol].atoms[atom1].name.c_str(),
                "-",
                system.molecules[mol].atoms[atom2].name.c_str(),
                "-",
                system.molecules[mol].atoms[atom3].name.c_str());
    }
    // and Dihedrals
    printf("================================================================================\n");
    printf("Dynamically-found Dihedrals Summary:\n");
    printf("================================================================================\n");
    printf("dihedral-id :: mol-id :: atom1 :: atom2 :: atom3 :: atom4 :: elements\n");
    for (int n=0; n<system.constants.uniqueDihedrals.size();n++) {
        int mol = system.constants.uniqueDihedrals[n].mol;
        int atom1 = system.constants.uniqueDihedrals[n].atom1;
        int atom2 = system.constants.uniqueDihedrals[n].atom2;
        int atom3 = system.constants.uniqueDihedrals[n].atom3;
        int atom4 = system.constants.uniqueDihedrals[n].atom4;
        printf("%11i :: %6i :: %5i :: %5i :: %5i :: %5i :: %4s%1s%4s%1s%4s%1s%4s\n", n,
            mol, atom1,atom2,atom3,atom4,
            system.molecules[mol].atoms[atom1].name.c_str(),
            "-",
            system.molecules[mol].atoms[atom2].name.c_str(),
            "-",
            system.molecules[mol].atoms[atom3].name.c_str(),
            "-",
            system.molecules[mol].atoms[atom4].name.c_str());
    }

    printf("================================================================================\n");
    printf("Dynamically-found UFF atom-types:\n");
    printf("================================================================================\n");
    printf("mol-id :: atom-id :: element :: UFF-label :: num-bonds\n");
    for (int i=0; i<system.molecules.size(); i++) {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
            printf("%6i :: %7i :: %7s :: %9s :: %9i\n",
                    i, j,
                    system.molecules[i].atoms[j].name.c_str(),
                    system.molecules[i].atoms[j].UFFlabel.c_str(),
                    (int)system.molecules[i].atoms[j].bonds.size());
        }
    }
    printf("================================================================================\n");
   

    /* START OPTIMIZATION */
    int converged = 0;
    double error_tolerance = system.constants.opt_error;
    int step_limit = system.constants.opt_step_limit; //100;
    double Ei = totalBondedEnergy(system);
    double Ef;
    double delta_E;
    double boltzmann;
    double tmp_pos[3] = {0,0,0};
    int randatom;
    int step=0;
    writeXYZ(system, system.constants.output_traj, 0, step, 0, 0);
    
    int optmode = system.constants.opt_mode;
    if (optmode == OPTIMIZE_SD)
        printf("STEEPEST DESCENT STRUCTURE OPTIMIZATION\n");
    else if (optmode == OPTIMIZE_MC)
        printf("MONTE CARLO STRUCTURE OPTIMIZATION\n");

    printf("Step %i :: Energy = %f; diff = %f kcal/mol; \n", 0, Ei, 0.0);

    // Monte Carlo sytle opt.
    if (optmode == OPTIMIZE_MC) {
        int N = (int)system.molecules[0].atoms.size();
    while (!converged) {
        Ei = totalBondedEnergy(system);

        // select random atom and perturb it.
        randatom = pickRandomAtom(system);
        for (int n=0;n<3;n++) tmp_pos[n] = system.molecules[0].atoms[randatom].pos[n];
        perturbAtom(system, randatom, move_factor(Ei, N));

        // get new energy
        Ef = totalBondedEnergy(system);
        delta_E = Ef - Ei;

      //  printf("Ef after = %f\n", Ef);

        //boltzmann = exp(-delta_E/0.00001); // just for now..
        //if (getrand() < boltzmann) {
        if (delta_E < 0) { // allow some positives a la Monte Carlo 
           // accept
            step++;
            writeXYZ(system, system.constants.output_traj, 0, step, 0, 0);
            outputEnergies(system, step, Ef, delta_E);
            if (fabs(delta_E) < error_tolerance && delta_E!=0) {
                printf("Finished with energy = %f kcal/mol \n", Ef);
                converged=1;
            }
        } else {
            // reject
            for (int n=0;n<3;n++) system.molecules[0].atoms[randatom].pos[n] = tmp_pos[n];
        }

        // check max-steps convergence
        if (step >= step_limit) {
            printf("Finished with energy = %f kcal/mol \n", Ef); 
            converged=1;
        }

    } // end while loop for convergence
    } // end MC style opt
    
    // steepest desent (follow the negative gradient)
    else if (optmode == OPTIMIZE_SD) {
        const double move_factor = 0.02;
        double grad_mag=0;
        while (!converged) {
            Ei = totalBondedEnergy(system);
            // re-initialize gradient
            for (int i=0; i<system.molecules[0].atoms.size(); i++) 
                for (int n=0;n<3;n++)
                    system.molecules[0].atoms[i].energy_grad[n]=0;
            
            // compute the gradients
            morse_gradient(system);
            angle_bend_gradient(system);
            torsions_gradient(system);
            LJ_intramolec_gradient(system);

            grad_mag = 0;
            // get gradient magnitude
            for (int i=0; i<system.molecules[0].atoms.size(); i++) {
                for (int n=0;n<3;n++)
                    grad_mag += system.molecules[0].atoms[i].energy_grad[n]*system.molecules[0].atoms[i].energy_grad[n];
            }
            grad_mag = sqrt(grad_mag);

            // move the atoms by their (negative!) gradients
            // normalized by the gradient magnitude
            for (int i=0; i<system.molecules[0].atoms.size(); i++) 
                for (int n=0;n<3;n++)
                    system.molecules[0].atoms[i].pos[n] -= move_factor/grad_mag * system.molecules[0].atoms[i].energy_grad[n];

            Ef = totalBondedEnergy(system);
            delta_E = Ef - Ei;

            step++;
            writeXYZ(system, system.constants.output_traj, 0, step, 0, 0);
            outputEnergies(system, step, Ef, delta_E);

            if ((fabs(delta_E) < error_tolerance && delta_E!=0) || step >= step_limit) {
                 printf("Finished with energy = %f kcal/mol \n", Ef);
                 converged=1;
            }
        }
    }

} // end optimize
