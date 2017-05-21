#include <string>
#include <algorithm>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

#define SQRT2  1.414213562373095
// ================ gaussian function ==================
// returns a gaussian-probability-generated velocity basied on mean (temperature goal) velocity and S.D.
double gaussian(double sigma) { // sigma is SD of the gaussian curve

    /* NOTE TO SELF
        the usage of erfInverse here may explode if rand()/RAND_MAX ever yields exactly -1 or +1
        because erf^-1( +-1) = +-inf
        if that ever happens I just need to check for abs(ranf) == 1 and avoid it.
    */

    // assuming mean velocity is zero (+ or - boltzmann velocity for particles yields net 0)
    double ranf = 2*(((double)rand() / (double)RAND_MAX)-0.5); // -1 to +1
    //return abs(sigma*SQRT2*erfInverse(ranf)); // I'm doing absolute value here bc +- is determined by velocity.
    
    // USE THE SIGMA AS MEAN-VALUE, RANDOMIZE DIRECTION. Seems to overinflate the velocities..
    //double pm=1.0;
    //double ranfsign = (double)rand() / (double)RAND_MAX;
    //if (ranfsign < 0.5) pm = -1.0;
    //double displacement = sigma * pm; // this is +/- the goal component velocity (gaussian width is same as mean)

    return sigma*SQRT2*erfInverse(ranf); //  + displacement;
    // if mean was nonzero it would be
    // mean + sigma*SQRT2*erfInverse(ranf);
    // my green notebook ("Space Group Research #2") has notes on this in 2nd divider.
}



// =================  GET TOTAL ENERGY AND EMERGENT TEMPERATURE FROM SYSTEM STATE ===========================
double * calculateEnergyAndTemp(System &system, double currtime) { // the * is to return an array of doubles as a pointer, not just one double
	double V_total = 0.0;
    double K_total = 0.0, Klin=0, Krot=0, Ek=0.0;
    double v_sum=0.0, avg_v = 0.0;
	double T=0.0, pressure=0;
    double vsq=0., wsq=0.;
    double energy_holder=0.;
	int i,j,n;

    // grab fixed potential energy of system
        if (system.constants.md_pbc) {
            V_total += getTotalPotential(system);
            //printf("used getTotalPotential\n");
        } else {
            for (i=0; i<system.molecules.size(); i++) {
                for (j=0; j<system.molecules[i].atoms.size(); j++) {
                    V_total += system.molecules[i].atoms[j].V;
                }
            }
        }

    // KINETIC ENERGIES, VELOCITIES, AND POTENTIALS //
    for (j=0; j<system.molecules.size(); j++) {
       if (system.constants.md_mode == MD_MOLECULAR) {
            vsq = 0; wsq = 0;
           for (n=0; n<3; n++) {
                vsq += system.molecules[j].vel[n] * system.molecules[j].vel[n];
                wsq += system.molecules[j].ang_vel[n] * system.molecules[j].ang_vel[n];
            }
            v_sum += sqrt(vsq); // so we're adding up velocities.

            energy_holder = 0.5 * system.molecules[j].mass * vsq;
            K_total += energy_holder; // linear: kg A^2 / fs^2
            Klin += energy_holder;

            if (system.constants.md_rotations) {
                energy_holder = 0.5 * system.molecules[j].inertia * wsq * system.constants.kb / 1e10;
                K_total += energy_holder; // rotational: kg A^2 / fs^2
                Krot += energy_holder;
            }
        }
        else if (system.constants.md_mode == MD_ATOMIC) {
            for (i=0; i<system.molecules[j].atoms.size(); i++) {
            vsq=0;
                for (n=0; n<3; n++) vsq += system.molecules[j].atoms[i].vel[n] * system.molecules[j].atoms[i].vel[n];
                v_sum += sqrt(vsq); // sum velocities
                energy_holder = 0.5*system.molecules[j].atoms[i].m * vsq;
                K_total += energy_holder;
                Klin += energy_holder;
            }
        }

        /*
         * I'm getting the state potential by the MC function instead of dynamic quantities.
         * It's more foolproof, I do believe.
        // iteratively sum ROTATIONAL POTENTIAL
        for (int n=0; n<3; n++) {
            V_total += -system.molecules[j].torque[n] * system.molecules[j].d_theta[n];
            //printf("rot pot: %e\n", -system.molecules[j].torque[n] * system.molecules[j].d_theta[n]);
        }


        // iteratively sum LINEAR POTENTIAL
        for (int i=0; i<system.molecules[j].atoms.size(); i++) {
            V_total += system.molecules[j].atoms[i].V;
        }
        */
    }

    avg_v = v_sum / system.stats.count_movables; //system.molecules.size(); // A/fs
    K_total = K_total / system.constants.kb * 1e10; // convert to K
    Klin = Klin / system.constants.kb * 1e10; // ""
    Krot = Krot / system.constants.kb * 1e10; // ""
    Ek = (3.0/2.0)*system.constants.temp; // 3/2 NkT, equipartition kinetic.

	// calculate temperature from kinetic energy and number of particles
	// https://en.wikipedia.org/wiki/Thermal_velocity
    T = (avg_v*1e5)*(avg_v*1e5) * system.proto[0].mass * M_PI / 8.0 / system.constants.kb; // NO GOOD FOR MULTISORBATE
    //T = 1e10*(avg_v*avg_v)*system.proto[0].mass / system.constants.kb;

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


void calculateForces(System &system, double dt) {
  int_fast8_t model = system.constants.potential_form;

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

    // GET FORCES
    // CPU style
    if (!system.constants.cuda) {
        // no pbc
        if (!system.constants.md_pbc) {
        if (model == POTENTIAL_LJ || model == POTENTIAL_LJES || model == POTENTIAL_LJESPOLAR || model == POTENTIAL_LJPOLAR)
            lj_force_nopbc(system);
        if (model == POTENTIAL_LJES || model == POTENTIAL_LJESPOLAR)
            coulombic_force_nopbc(system);
        } 
        // pbc
        else {
        if (model == POTENTIAL_LJ || model == POTENTIAL_LJES || model == POTENTIAL_LJESPOLAR || model == POTENTIAL_LJPOLAR)
            lj_force(system);
        if (model == POTENTIAL_LJES || model == POTENTIAL_LJESPOLAR)
            coulombic_real_force(system);
        if (model == POTENTIAL_LJESPOLAR || model == POTENTIAL_LJPOLAR)
            polarization_force(system);
            //int index=0;
            //for (int i=0;i<system.molecules.size();i++)
              //  for (int j=0;j<system.molecules[i].atoms.size();j++){
           //         printf("H[0] force = %f %f %f\n",system.molecules[0].atoms[0].force[0], system.molecules[0].atoms[0].force[1], system.molecules[0].atoms[0].force[2]);
                //    index++;
                //}
            
        }
    // GPU style
    } else {
        #ifdef CUDA
        // CUDA FORCES
        // no pbc    
        if (!system.constants.md_pbc) {
            CUDA_force_nopbc(system);
        // pbc
        } else {
            CUDA_force(system);
        }
        #endif
    }

    // atomic forces are done, so now calc molecular values
    for (int i=0; i<system.molecules.size(); i++) {
        if (system.molecules[i].frozen) continue;
        system.molecules[i].calc_force();
        if (system.constants.md_rotations && system.molecules[i].atoms.size() > 1)
            system.molecules[i].calc_torque();
    }

}

// ==================== MOVE ATOMS MD STYLE =========================
/* THIS IS THE MAIN LOOPING FUNCTION. calculateForces() is called within */
void integrate(System &system, double dt) {
    system.checkpoint("started integrate()");
    int i,j,n;

    // DEBUG
    int_fast8_t debug=0;
    if (debug == 1) {
        for (j=0; j<system.molecules.size(); j++) {
            if (system.constants.md_mode == MD_MOLECULAR) system.molecules[j].printAll();
            for (i=0; i<system.molecules[j].atoms.size(); i++) {
                if (system.constants.md_mode == MD_ATOMIC) system.molecules[j].atoms[i].printAll();
            }
        }
    }
    // END IF DEBUG

    // 1a) CHANGE POSITIONS OF PARTICLES
    /* NO NEED FOR SAVING OLD POS
    // save old positions
    for (j=0; j<system.molecules.size(); j++) {
        if (!system.molecules[j].frozen) {
        for (i=0; i<system.molecules[j].atoms.size(); i++) {
            for (n=0; n<3; n++) system.molecules[j].atoms[i].prevpos[n] = system.molecules[j].atoms[i].pos[n];
        } // end for atom i
        } // end if movable
    } // end for molecule j
    // done saving old positions
    */

    system.checkpoint("moving particles based on forces.");
    // if molecular motion
    if (system.constants.md_mode == MD_MOLECULAR) {
        for (j=0; j<system.molecules.size(); j++) {
            if (!system.molecules[j].frozen) {

            system.molecules[j].calc_pos(dt);

              // ROTATION
            if (system.constants.md_rotations && system.molecules[j].atoms.size() > 1) {
            system.molecules[j].calc_ang_pos(dt);

            // rotate molecules
            for (i=0; i<system.molecules[j].atoms.size(); i++) {
                // ROTATE IN X
                double* rotatedx = rotatePointRadians(system,
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                0, system.molecules[j].ang_pos[0] );
                for (n=0; n<3; n++)
                    system.molecules[j].atoms[i].pos[n] = rotatedx[n] + system.molecules[j].com[n];

                // ROTATE IN Y
                double* rotatedy = rotatePointRadians(system,
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                1, system.molecules[j].ang_pos[1] );
                for (n=0; n<3; n++)
                    system.molecules[j].atoms[i].pos[n] = rotatedy[n] + system.molecules[j].com[n];

                // ROTATE IN Z
                double* rotatedz = rotatePointRadians(system,
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                2, system.molecules[j].ang_pos[2] );
                for (n=0; n<3; n++)
                    system.molecules[j].atoms[i].pos[n] = rotatedz[n] + system.molecules[j].com[n];
            } // end loop over atoms i
            } // end if rotations allowed and >1 atom
            } // end if movable molecule
        } // end for molecules j
    } // end if molecular motion
    // if atomic motion
    else if (system.constants.md_mode == MD_ATOMIC) {
        for (j=0; j<system.molecules.size(); j++) {
            if (!system.molecules[j].frozen) {
            for (i=0; i<system.molecules[j].atoms.size(); i++) {
                system.molecules[j].atoms[i].calc_pos(dt);
            }
            }
        }
    } // end if atomic motion
    // END POSITION CHANGES
    system.checkpoint("done moving particles. Checking PBC for all particles");

    // 1b) CHECK P.B.C. (move the molecule/atom back in the box if needed)
    if (system.constants.md_pbc) {
        for (j=0; j<system.molecules.size(); j++) {
            if (!system.molecules[j].frozen) {
                checkInTheBox(system,j); // also computes COM 
            } // end if movable
	    } // end loop j molecules
    } // end if PBC

    // 3) GET NEW FORCES (AND TORQUES) BASED ON NEW POSITIONS
	system.checkpoint("done checking PBC. Starting calculateForces()");
    calculateForces(system, dt);
    system.checkpoint("Done with calculateForces(). Starting integrator (for a&v)");

    // 4) GET NEW ACCELERATION AND VELOCITY FOR ALL PARTICLES
	// Normal CPU routine
    //if (!system.constants.cuda) {
    for (j=0; j<system.molecules.size(); j++) {
		if (!system.molecules[j].frozen) { // only movable atoms should move.

            // if atoms allowed to move from molecules
            if (system.constants.md_mode == MD_ATOMIC) {
                for (i=0; i<system.molecules[j].atoms.size(); i++) {
                    system.molecules[j].atoms[i].calc_acc();
                    system.molecules[j].atoms[i].calc_vel(dt);
            } // end atomic loop i
            } // end if atomic
            // otherwise handle molecular movement with rigidity.
            else if (system.constants.md_mode == MD_MOLECULAR) {
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
    //}
    // CUDA GPU style
    //else {
     //   #ifdef CUDA
      //  CUDA_verlet(system);
       // #endif
    //}
    system.checkpoint("Done with a,v integration. Starting heat bath (if nvt)");

    // 5) apply heat bath in NVT
    if (system.constants.ensemble == ENSEMBLE_NVT) {
        // loop through all molecules and adjust velocities by Anderson Thermostat method
        // this process makes the NVT MD simulation stochastic/ Markov / MC-like, which is good for equilibration results.
        // the thermostat probability must be recalc'd every new step.
        //double probab = system.constants.md_thermostat_freq *
         //           exp(-system.constants.md_thermostat_freq * system.stats.MDtime);
        //printf("MD thermostat probab = %f\n", probab);
        // or not? 
        double probab = system.constants.md_thermostat_probab;

        double ranf; //, sigma = sqrt(system.constants.kb * system.constants.temp /  system.proto[0].mass) *1e-5; // to A/s
        double sigma = system.constants.md_vel_goal;
        //printf("sigma = %f\n", sigma);
        //double newvel;
        if (system.constants.md_mode == MD_MOLECULAR) {
        for (i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen) continue; // skip frozens
            ranf = (double)rand() / (double)RAND_MAX; // 0 -> 1
            if (ranf < probab) {
                // adjust the velocity components of the molecule.
                for (n=0; n<3; n++) {
                    //printf("gauss(sigma, mean) = gauss(%f, %f) = %f\n", sigma, mean_velocity, gaussian(sigma, mean_velocity));
                    system.molecules[i].vel[n] = gaussian(sigma);
                    //printf("Gaussian molecular vel[%i]: %f\n",n, system.molecules[i].vel[n]);
                    //if (system.molecules[i].vel[n] < 0) system.molecules[i].vel[n] = -newvel;
                    //else system.molecules[i].vel[n] = newvel;
                    /*
                    if (system.molecules[i].vel[n] >= 0) {
                        //system.molecules[i].vel[n] = system.constants.md_vel_goal;
                        system.molecules[i].vel[n] = gaussian(sigma,mean_velocity);
                    }
                    else {
                        //system.molecules[i].vel[n] = -system.constants.md_vel_goal;
                        system.molecules[i].vel[n] = gaussian(sigma, -mean_velocity);
                    }*/

                }
            }
        }
        } else if (system.constants.md_mode == MD_ATOMIC) {
            for (i =0; i<system.molecules.size(); i++) {
                for (j=0; j<system.molecules[i].atoms.size(); j++) {
                    if (system.molecules[i].atoms[j].frozen) continue; // skip frozen atoms
                    ranf = (double)rand() / (double)RAND_MAX; // 0 -> 1
                    if (ranf <probab) {
                        for (n=0; n<3; n++) {
                            system.molecules[i].vel[n] = gaussian(sigma);
                            //if (system.molecules[i].vel[n] < 0) system.molecules[i].vel[n] = -newvel;
                            //else system.molecules[i].vel[n] = newvel;
                    //      if (system.molecules[i].atoms[j].vel[n] >= 0) system.molecules[i].atoms[j].vel[n] = system.constants.md_vel_goal;
                       //     else system.molecules[i].atoms[j].vel[n] = -system.constants.md_vel_goal;
                        }
                    }
                }
            }
        }
    } // end if NVT (thermostat)
    system.checkpoint("Done with heatbath if NVT.");
    system.checkpoint("Done with integrate() function.");
}// end integrate() function


