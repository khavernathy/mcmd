#include <string>
#include <algorithm>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

// =================  GET TOTAL ENERGY AND EMERGENT TEMPERATURE FROM SYSTEM STATE ===========================
double * calculateEnergyAndTemp(System &system, double currtime) { // the * is to return an array of doubles as a pointer, not just one double
	double V_total = 0.0;
    double K_total = 0.0, Klin=0, Krot=0, Ek=0.0;
    double v_sum=0.0, avg_v = 0.0;
	double T=0.0, pressure=0;
	

    // KINETIC ENERGIES, VELOCITIES, AND POTENTIALS //
    for (int j=0; j<system.molecules.size(); j++) {
       if (system.constants.md_mode == "molecular") {
            double vsq = 0, wsq = 0;
           for (int n=0; n<3; n++) {
                vsq += system.molecules[j].vel[n] * system.molecules[j].vel[n];
                wsq += system.molecules[j].ang_vel[n] * system.molecules[j].ang_vel[n];
            }
            vsq = sqrt(vsq); wsq = sqrt(wsq);
            v_sum += vsq; // so we're adding up velocities.
            vsq *= vsq;  wsq *= wsq;

            K_total += 0.5 * system.molecules[j].mass * vsq; // linear
            Klin += 0.5 * system.molecules[j].mass * vsq;            

            if (system.constants.md_rotations) {
            K_total += 0.5 * system.molecules[j].inertia * wsq * system.constants.kb / 1e10; // rotational
            Krot += 0.5 * system.molecules[j].inertia * wsq * system.constants.kb / 1e10;
            }
        }
        else if (system.constants.md_mode == "atomic") { 
            for (int i=0; i<system.molecules[j].atoms.size(); i++) {
            double vsq=0;
                for (int n=0; n<3; n++) vsq += system.molecules[j].atoms[i].vel[n] * system.molecules[j].atoms[i].vel[n];
                vsq = sqrt(vsq);
                v_sum += vsq; // sum velocities
                vsq *= vsq;
                K_total += 0.5 * system.molecules[j].atoms[i].m * vsq;
                Klin += 0.5 * system.molecules[j].atoms[i].m * vsq;
            }
        }

        // iteratively sum ROTATIONAL POTENTIAL
        for (int n=0; n<3; n++) {
            V_total += -system.molecules[j].torque[n] * system.molecules[j].d_theta[n];
            //printf("rot pot: %e\n", -system.molecules[j].torque[n] * system.molecules[j].d_theta[n]);
        }

        // iteratively sum LINEAR POTENTIAL
        for (int i=0; i<system.molecules[j].atoms.size(); i++) {
            V_total += system.molecules[j].atoms[i].V;
        }
    }

    avg_v = v_sum / system.molecules.size(); // A/fs
    K_total = K_total / system.constants.kb * 1e10; // convert to K 
    Klin = Klin / system.constants.kb * 1e10;
    Krot = Krot / system.constants.kb * 1e10;
    Ek = (3.0/2.0)*system.constants.temp; // 3/2 NkT, equipartition kinetic.

	// calculate temperature from kinetic energy and number of particles
	// https://en.wikipedia.org/wiki/Thermal_velocity
    T = (avg_v*1e5)*(avg_v*1e5) * system.proto[0].mass * M_PI / 8.0 / system.constants.kb; // NO GOOD FOR MULTISORBATE

   
/*
    if (system.constants.ensemble == "nvt") {
        // grab the sum of F.r's
        double frsum=0;
        for (int i=0; i<system.pairs.size(); i++) {
            frsum += system.pairs[i].fdotr;
        }
        system.stats.fdotrsum.value = frsum;
        system.stats.fdotrsum.calcNewStats(); // makes new average based on new value.

        double volInLiters = system.pbc.volume * system.constants.A32L;
        pressure = (system.stats.count_movables/volInLiters)*
            system.constants.kb * T + (1.0/(3.0*volInLiters))*system.stats.fdotrsum.average*system.constants.kb;
        pressure *= system.constants.JL2ATM; // J/L -> atm
    }
 */

	static double output[8];
	output[0] = K_total;
    output[1] = V_total;
	output[2] = T;
    output[3] = avg_v;
    output[4] = Ek;
    output[5] = Klin;
    output[6] = Krot;
    output[7] = pressure;
	return output;
}


void calculateForces(System &system, string model, double dt) {
	
    // initialize variable for pressure calc in NVT
    // loop through all atoms
	for (int j=0; j <system.molecules.size(); j++) {
	for (int i = 0; i < system.molecules[j].atoms.size(); i++) {
		// initialize stuff
		system.molecules[j].atoms[i].force[0] = 0.0;
		system.molecules[j].atoms[i].force[1] = 0.0;
		system.molecules[j].atoms[i].force[2] = 0.0;
		system.molecules[j].atoms[i].V = 0.0;
	}
	}
	
    if (model == "lj" || model == "ljes" || model == "ljespolar")
        lj_force(system);    
    if (model == "ljes" || model == "ljespolar")
        coulombic_real_force(system);

    
    // atomic forces are done, so now calc molecular values
    for (int i=0; i<system.molecules.size(); i++) {
        system.molecules[i].calc_force();
        if (system.constants.md_rotations && system.molecules[i].atoms.size() > 1) 
            system.molecules[i].calc_torque();
    }

}

// ==================== MOVE ATOMS MD STYLE =========================
/* THIS IS THE MAIN LOOPING FUNCTION. calculateForces() is called within */
void integrate(System &system, double dt) {
    int i,j,n;

    // DEBUG
    int_fast8_t debug=0;
    if (debug == 1) {
        for (j=0; j<system.molecules.size(); j++) {
            if (system.constants.md_mode == "molecular") system.molecules[j].printAll();
            for (i=0; i<system.molecules[j].atoms.size(); i++) {
                if (system.constants.md_mode == "atomic") system.molecules[j].atoms[i].printAll();
            }
        }
    }
    // END IF DEBUG

    // 1a) CHANGE POSITIONS OF PARTICLES
    // save old positions
    for (j=0; j<system.molecules.size(); j++) {
        if (!system.molecules[j].frozen) {
        for (i=0; i<system.molecules[j].atoms.size(); i++) {
            for (n=0; n<3; n++) system.molecules[j].atoms[i].prevpos[n] = system.molecules[j].atoms[i].pos[n];
        } // end for atom i
        } // end if movable
    } // end for molecule j
    // done saving old positions

    // if molecular motion
    if (system.constants.md_mode == "molecular") {
        for (j=0; j<system.molecules.size(); j++) {
            if (!system.molecules[j].frozen) {

            system.molecules[j].calc_pos(dt);
            
              // ROTATION
            if (system.constants.md_rotations && system.molecules[j].atoms.size() > 1) {
            system.molecules[j].calc_ang_pos(dt);

            // rotate molecules
            for (i=0; i<system.molecules[j].atoms.size(); i++) {
                // ROTATE IN X
                double* rotatedx = rotatePoint(system, 
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                0, system.molecules[j].ang_pos[0] * 180.0/M_PI); 
                for (n=0; n<3; n++) 
                    system.molecules[j].atoms[i].pos[n] = rotatedx[n] + system.molecules[j].com[n];

                // ROTATE IN Y
                double* rotatedy = rotatePoint(system, 
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                1, system.molecules[j].ang_pos[1] * 180.0/M_PI); 
                for (n=0; n<3; n++) 
                    system.molecules[j].atoms[i].pos[n] = rotatedy[n] + system.molecules[j].com[n];

                // ROTATE IN Z
                double* rotatedz = rotatePoint(system, 
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                2, system.molecules[j].ang_pos[2] * 180.0/M_PI); 
                for (n=0; n<3; n++) 
                    system.molecules[j].atoms[i].pos[n] = rotatedz[n] + system.molecules[j].com[n];
            } // end loop over atoms i 
            } // end if rotations allowed and >1 atom
            } // end if movable molecule
        } // end for molecules j
    } // end if molecular motion
    // if atomic motion
    else if (system.constants.md_mode == "atomic") {
        for (j=0; j<system.molecules.size(); j++) {
            if (!system.molecules[j].frozen) {
            for (i=0; i<system.molecules[j].atoms.size(); i++) {
                system.molecules[j].atoms[i].calc_pos(dt);
            }
            }
        }
    } // end if atomic motion
    // END POSITION CHANGES

    // 1b) CHECK P.B.C. (move the molecule/atom back in the box if needed)
    if (system.constants.md_pbc) {
        for (j=0; j<system.molecules.size(); j++) {
            if (!system.molecules[j].frozen) {
                checkInTheBox(system,j);
            } // end if movable
	    } // end loop j molecules
    } // end if PBC
   
    // 2) GET NEW C.O.M. for new positions.
    for (j=0; j<system.molecules.size(); j++) system.molecules[j].calc_center_of_mass();
 
    // 3) GET NEW FORCES (AND TORQUES) BASED ON NEW POSITIONS
	calculateForces(system, system.constants.potential_form, dt);

    // 4) GET NEW ACCELERATION AND VELOCITY FOR ALL PARTICLES
	for (j=0; j<system.molecules.size(); j++) {
		if (!system.molecules[j].frozen) { // only movable atoms should move.

            // if atoms allowed to move from molecules
            if (system.constants.md_mode == "atomic") {
                for (i=0; i<system.molecules[j].atoms.size(); i++) {
                    system.molecules[j].atoms[i].calc_acc();
                    system.molecules[j].atoms[i].calc_vel(dt); 
            } // end atomic loop i
            } // end if atomic
            // otherwise handle molecular movement with rigidity.
            else if (system.constants.md_mode == "molecular") {
                // linear
                system.molecules[j].calc_acc();
                system.molecules[j].calc_vel(dt, system.constants.md_vel_goal);
    
                // rotational
                if (system.constants.md_rotations) {
                    system.molecules[j].calc_ang_acc();
                    system.molecules[j].calc_ang_vel(dt);
                }
            }
        } // end if movable
    } // end for j molecules

    // 5) apply heat bath in NVT
    if (system.constants.ensemble == ENSEMBLE_NVT) {
        // loop through all molecules and adjust velocities by Anderson Thermostat method
        // this process makes the NVT MD simulation stochastic/ Markov / MC-like, which is good for equilibration results.
        double probab = system.constants.md_thermostat_probab;
        double ranf;
        if (system.constants.md_mode == "molecular") {
        for (i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen) continue; // skip frozens
            ranf = (double)rand() / (double)RAND_MAX; // 0 -> 1
            if (ranf < probab) {
                // adjust the velocity components of the molecule.
                for (n=0; n<3; n++) {
                    if (system.molecules[i].vel[n] >= 0) system.molecules[i].vel[n] = system.constants.md_vel_goal;
                    else system.molecules[i].vel[n] = -system.constants.md_vel_goal;
                }
            }
        }
        } else if (system.constants.md_mode == "atomic") {
            for (i =0; i<system.molecules.size(); i++) {
                for (j=0; j<system.molecules[i].atoms.size(); j++) {
                    if (system.molecules[i].atoms[j].frozen) continue; // skip frozen atoms
                    ranf = (double)rand() / (double)RAND_MAX; // 0 -> 1
                    if (ranf <probab) {
                        for (n=0; n<3; n++) {
                            if (system.molecules[i].atoms[j].vel[n] >= 0) system.molecules[i].atoms[j].vel[n] = system.constants.md_vel_goal;
                            else system.molecules[i].atoms[j].vel[n] = -system.constants.md_vel_goal;
                        }
                    }            
                }
            }
        }
    } // end if NVT (thermostat)
}// end integrate() function
