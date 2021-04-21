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
    for (int n=0; n<3; n++)
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

void outputEnergies(System &system, int step, double Ef, double delta_E, double sec_per_step) {
    printf("==============================================================\n");
#ifndef WINDOWS
    printf("Optimization Step %i (%.4f sec/step) %s\nEnergy =         %f kcal/mol; \u0394E = %f kcal/mol; \n\n", step, sec_per_step, system.constants.atom_file.c_str(),  Ef * system.constants.kbk, delta_E * system.constants.kbk);
#else
    printf("Optimization Step %i (%.4f sec/step) %s\nEnergy =         %f kcal/mol; dE = %f kcal/mol; \n\n", step, sec_per_step, system.constants.atom_file.c_str(),  Ef * system.constants.kbk, delta_E * system.constants.kbk);
#endif
    printf("Bonds =          %f kcal/mol\nAngle-bends =    %f kcal/mol\nDihedrals =      %f kcal/mol\nIntramolec. LJ = %f kcal/mol\nIntramolec. ES = %f kcal/mol\n",
           system.stats.Ustretch.value * system.constants.kbk,
           system.stats.Uangles.value * system.constants.kbk,
           system.stats.Udihedrals.value * system.constants.kbk,
           system.stats.UintraLJ.value * system.constants.kbk,
           system.stats.UintraES.value * system.constants.kbk);
}

void printBondParameters(System &system) {
    // print out all the bonds
    printf("================================================================================\n");
    printf("Dynamically-found Bonds Summary:\n");
    printf("================================================================================\n");
    printf("bond-id :: mol-id :: atom1 :: atom2 ::  elements :: length   :: r_ij\n");
    for (int n=0; n<system.constants.uniqueBonds.size(); n++) {
        //printf("Atom %i (UFF: %s)\n", i, system.molecules[0].atoms[i].UFFlabel.c_str());
        int mol=system.constants.uniqueBonds[n].mol;
        int atom1=system.constants.uniqueBonds[n].atom1;
        int atom2=system.constants.uniqueBonds[n].atom2;
        double value = system.constants.uniqueBonds[n].value;
        double rij = system.constants.uniqueBonds[n].rij;
        printf("%7i :: %6i :: %5i :: %5i :: %4s%1s%4s :: %5f :: %5f\n",
               n,
               mol,
               atom1,
               atom2,
               system.molecules[mol].atoms[atom1].name.c_str(),
               "-",
               system.molecules[mol].atoms[atom2].name.c_str(),
               value,
               rij
              );
    }
    // and angles
    printf("================================================================================\n");
    printf("Dynamically-found Angles Summary:\n");
    printf("================================================================================\n");
    printf("angle-id :: mol-id :: atom1 :: atom2 :: atom3 ::    elements    :: angle     :: theta_ijk\n");
    for (int n=0; n<system.constants.uniqueAngles.size(); n++) {
        int mol = system.constants.uniqueAngles[n].mol;
        int atom1= system.constants.uniqueAngles[n].atom1;
        int atom2 = system.constants.uniqueAngles[n].atom2;
        int atom3 = system.constants.uniqueAngles[n].atom3;
        double value = system.constants.uniqueAngles[n].value;
        double theta_ijk = system.constants.uniqueAngles[n].theta_ijk*180.0/M_PI;
        printf("%8i :: %6i :: %5i :: %5i :: %5i :: %4s%1s%4s%1s%4s :: %3.5f :: %3.5f\n", n,
               mol,
               atom1,
               atom2,
               atom3,
               system.molecules[mol].atoms[atom1].name.c_str(),
               "-",
               system.molecules[mol].atoms[atom2].name.c_str(),
               "-",
               system.molecules[mol].atoms[atom3].name.c_str(),
               value*180./M_PI,
               theta_ijk);
    }
    // and Dihedrals
    printf("================================================================================\n");
    printf("Dynamically-found Dihedrals Summary:\n");
    printf("================================================================================\n");
    printf("dihedral-id :: mol-id :: atom1 :: atom2 :: atom3 :: atom4 ::       elements      :: angle     :: phi_ijkl\n");
    for (int n=0; n<system.constants.uniqueDihedrals.size(); n++) {
        int mol = system.constants.uniqueDihedrals[n].mol;
        int atom1 = system.constants.uniqueDihedrals[n].atom1;
        int atom2 = system.constants.uniqueDihedrals[n].atom2;
        int atom3 = system.constants.uniqueDihedrals[n].atom3;
        int atom4 = system.constants.uniqueDihedrals[n].atom4;
        double value = system.constants.uniqueDihedrals[n].value;
        double phi_ijkl = system.constants.uniqueDihedrals[n].phi_ijkl*180.0/M_PI;
        printf("%11i :: %6i :: %5i :: %5i :: %5i :: %5i :: %4s%1s%4s%1s%4s%1s%4s :: %3.5f :: %3.5f\n", n,
               mol, atom1,atom2,atom3,atom4,
               system.molecules[mol].atoms[atom1].name.c_str(),
               "-",
               system.molecules[mol].atoms[atom2].name.c_str(),
               "-",
               system.molecules[mol].atoms[atom3].name.c_str(),
               "-",
               system.molecules[mol].atoms[atom4].name.c_str(),
               value*180./M_PI, phi_ijkl);
    }

    // and Impropers
    printf("================================================================================\n");
    printf("Dynamically-found Impropers Summary:\n");
    printf("================================================================================\n");
    printf("improper-id :: mol-id :: atom1 :: atom2 :: atom3 :: atom4 ::       elements      :: angle \n");
    for (int n=0; n<system.constants.uniqueImpropers.size(); n++) {
        int mol = system.constants.uniqueImpropers[n].mol;
        int atom1 = system.constants.uniqueImpropers[n].atom1;
        int atom2 = system.constants.uniqueImpropers[n].atom2;
        int atom3 = system.constants.uniqueImpropers[n].atom3;
        int atom4 = system.constants.uniqueImpropers[n].atom4;
        double value = system.constants.uniqueImpropers[n].value;
        printf("%11i :: %6i :: %5i :: %5i :: %5i :: %5i :: %4s%1s%4s%1s%4s%1s%4s :: %3.5f \n", n,
               mol, atom1,atom2,atom3,atom4,
               system.molecules[mol].atoms[atom1].name.c_str(),
               "-",
               system.molecules[mol].atoms[atom2].name.c_str(),
               "-",
               system.molecules[mol].atoms[atom3].name.c_str(),
               "-",
               system.molecules[mol].atoms[atom4].name.c_str(),
               value*180./M_PI);
    }


    printf("================================================================================\n");
    printf("Dynamically-found UFF atom-types:\n");
    printf("================================================================================\n");
    printf("mol-id :: atom-id :: element :: UFF-label :: num-bonds\n");
    int molecule_limit = 1;
    if (system.constants.md_mode == MD_FLEXIBLE) molecule_limit = system.molecules.size();
    for (int i=0; i<molecule_limit; i++) {
        if (system.molecules[i].frozen && system.constants.mode != "opt" && !system.constants.flexible_frozen && !system.constants.write_lammps) continue; // skip frozens if not optimization mode
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
            printf("%6i :: %7i :: %7s :: %9s :: %9i\n",
                   i, j,
                   system.molecules[i].atoms[j].name.c_str(),
                   system.molecules[i].atoms[j].UFFlabel.c_str(),
                   (int)system.molecules[i].atoms[j].bonds.size());
        }
    }
    printf("================================================================================\n");

    printf("There are:\n %i atoms\n %i bonds\n %i angles\n %i dihedrals\n %i impropers\n %i LJ/ES qualified pairs.\n\n", system.constants.total_atoms, (int)system.constants.uniqueBonds.size(), (int)system.constants.uniqueAngles.size(), (int)system.constants.uniqueDihedrals.size(), (int)system.constants.uniqueImpropers.size(), (int)system.constants.uniqueLJNonBonds.size());
}

// Optimize the molecule (ID=0) via MM forcefield(s)
void optimize(System &system) {

    //int i;
    std::chrono::steady_clock::time_point begin_opt = std::chrono::steady_clock::now();

    printBondParameters(system);

    /* START OPTIMIZATION */
    int converged = 0;
    double error_tolerance = system.constants.opt_error;
    int step_limit = system.constants.opt_step_limit; //100;
    double Ef = totalBondedEnergy(system);
    double Ei;
    double delta_E;
    double tmp_pos[3] = {0,0,0};
    int randatom;
    int step=0;
    writeXYZ(system, system.constants.output_traj, 0, step, 0, 0);

    int optmode = system.constants.opt_mode;
    if (optmode == OPTIMIZE_SD)
        printf("STEEPEST DESCENT STRUCTURE OPTIMIZATION\n");
    else if (optmode == OPTIMIZE_MC)
        printf("MONTE CARLO STRUCTURE OPTIMIZATION\n");

    outputEnergies(system, 0, Ef, 0, 0);

    // Monte Carlo sytle opt.
    if (optmode == OPTIMIZE_MC) {
        int N = (int)system.molecules[0].atoms.size();
        while (!converged) {
            Ei = Ef;

            // select random atom and perturb it.
            randatom = pickRandomAtom(system);
            for (int n=0; n<3; n++) tmp_pos[n] = system.molecules[0].atoms[randatom].pos[n];
            perturbAtom(system, randatom, move_factor(Ei, N));

            // get new energy
            Ef = totalBondedEnergy(system);
            delta_E = Ef - Ei;


            if (delta_E < 0) { // allow some positives a la Monte Carlo
                // accept
                step++;
                writeXYZ(system, system.constants.output_traj, 0, step, 0, 0);

                std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
                double time_elapsed = (std::chrono::duration_cast<std::chrono::microseconds>(end - begin_opt).count()) /1000000.0;
                double sec_per_step = time_elapsed/step;

                outputEnergies(system, step, Ef, delta_E, sec_per_step);
                if (fabs(delta_E) < error_tolerance && delta_E!=0) {
                    printf("Finished with energy = %f kcal/mol \n", Ef * system.constants.kbk);
                    converged=1;
                }
            } else {
                // reject
                for (int n=0; n<3; n++) system.molecules[0].atoms[randatom].pos[n] = tmp_pos[n];
            }

            // check max-steps convergence
            if (step >= step_limit) {
                printf("Finished with energy = %f kcal/mol \n", Ef * system.constants.kbk);
                converged=1;
            }

        } // end while loop for convergence
    } // end MC style opt

    // steepest desent (follow the negative gradient)
    else if (optmode == OPTIMIZE_SD) {
        const double move_factor = 0.02;
        double grad_mag=0;
        while (!converged) {
            Ei = Ef;
            // re-initialize gradient
            for (int i=0; i<system.molecules[0].atoms.size(); i++)
                for (int n=0; n<3; n++)
                    system.molecules[0].atoms[i].force[n]=0;

            //printf("COMPUTING GRADIENTS...\n");
            // compute the gradients
            if (system.constants.opt_bonds)
                morse_gradient(system);
            if (system.constants.opt_angles)
                angle_bend_gradient(system);
            if (system.constants.opt_dihedrals)
                torsions_gradient(system);
            if (system.constants.opt_LJ)
                LJ_intramolec_gradient(system);
            if (system.constants.opt_ES)
                ES_intramolec_gradient(system);

            grad_mag = 0;
            // get gradient magnitude
            for (int i=0; i<system.molecules[0].atoms.size(); i++) {
                for (int n=0; n<3; n++) {
                    //printf("gradient %i[%i] = %f\n", i,n, system.molecules[0].atoms[i].force[n]);
                    grad_mag += system.molecules[0].atoms[i].force[n]*system.molecules[0].atoms[i].force[n];
                }
            }
            grad_mag = sqrt(grad_mag);
            //printf("Gradient magnitude: %f\n", grad_mag);

            // move the atoms by their (negative!) gradients
            // normalized by the gradient magnitude
            for (int i=0; i<system.molecules[0].atoms.size(); i++)
                for (int n=0; n<3; n++)
                    system.molecules[0].atoms[i].pos[n] += move_factor/grad_mag * system.molecules[0].atoms[i].force[n];

            Ef = totalBondedEnergy(system);
            delta_E = Ef - Ei;

            step++;
            writeXYZ(system, system.constants.output_traj, step, step, 0, 0);
            std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
            double time_elapsed = (std::chrono::duration_cast<std::chrono::microseconds>(end - begin_opt).count()) /1000000.0;
            double sec_per_step = time_elapsed/step;
            outputEnergies(system, step, Ef, delta_E, sec_per_step);

            if ((fabs(delta_E) < error_tolerance && delta_E!=0) || step >= step_limit) {
                printf("Finished with energy = %f kcal/mol \n", Ef * system.constants.kbk);
                converged=1;
            }
        }
    }

} // end optimize
