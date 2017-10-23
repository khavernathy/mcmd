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
        the usage of erfInverse here may explode if getrand() ever yields exactly -1 or +1
        because erf^-1( +-1) = +-inf
        if that ever happens I just need to check for abs(ranf) == 1 and avoid it.
    */

    double ranf = 2*((getrand())-0.5); // -1 to +1

    return sigma*SQRT2*erfInverse(ranf); //  + displacement;
    // my green notebook ("Space Group Research #2") has notes on this in 2nd divider.
}



// =================  GET TOTAL ENERGY AND EMERGENT TEMPERATURE FROM SYSTEM STATE ===========================
double * calculateObservablesMD(System &system, double currtime) { // the * is to return an array of doubles as a pointer, not just one double
	double V_total = 0.0;
    double K_total = 0.0, Klin=0, Krot=0, Ek=0.0;
    double avg_v = 0.0, avg_v2=0.0, v_rms=0.0;
    double avg_v_ALL=0;
	double T=0.0, T_rms=0, pressure=0;
    double vsq=0., wsq=0.;
    double energy_holder=0.;
	int i,j,n,z;

    double v_sum[(int)system.proto.size()];
    double v2_sum[(int)system.proto.size()];
    int N_local[(int)system.proto.size()];


    // grab fixed potential energy of system
        // from PBC (same as Monte Carlo potential)
        if (system.constants.md_pbc) {
            V_total += getTotalPotential(system);
        // from individual atomic contributions (no PBC)
        } else {
            for (i=0; i<system.molecules.size(); i++) {
                for (j=0; j<system.molecules[i].atoms.size(); j++) {
                    V_total += system.molecules[i].atoms[j].V;
                }
            }
        }

    // KINETIC ENERGIES, VELOCITIES, ETC. BY MOLECULE TYPE
    for (z = 0; z<system.proto.size(); z++) {
        v2_sum[z] = 0;
        v_sum[z] = 0;
        N_local[z] = 0;
       for (i=0; i<system.molecules.size(); i++) {
            if (system.proto[z].name == system.molecules[i].name) {
                N_local[z]++;
                if (system.constants.md_mode == MD_MOLECULAR) {
                    vsq=0; wsq=0;
                    for (n=0; n<3; n++) {
                        vsq += system.molecules[i].vel[n]*system.molecules[i].vel[n];
                        wsq += system.molecules[i].ang_vel[n]*system.molecules[i].ang_vel[n];
                    }
                    v2_sum[z] += vsq;
                    v_sum[z] += sqrt(vsq);

                    energy_holder = 0.5 * system.molecules[i].mass * vsq;
                    K_total += energy_holder;
                    Klin += energy_holder;

                    if (system.constants.md_rotations) {
                        // new tensor method.
                        system.molecules[i].calc_inertia_tensor();
                        double wx = system.molecules[i].ang_vel[0];
                        double wy = system.molecules[i].ang_vel[1];
                        double wz = system.molecules[i].ang_vel[2];

                        energy_holder = 0.5 * (system.molecules[i].inertia_tensor[0]*wx*wx +
                            system.molecules[i].inertia_tensor[1]*wy*wy +
                            system.molecules[i].inertia_tensor[2]*wz*wz +
                            2*system.molecules[i].inertia_tensor[3]*wx*wy +
                            2*system.molecules[i].inertia_tensor[4]*wy*wz +
                            2*system.molecules[i].inertia_tensor[5]*wx*wz);

                        energy_holder *= system.constants.kb/1e10;

                        K_total += energy_holder; // rotational: (rad^2)*kg A^2 / fs^2
                        Krot += energy_holder;
                    } // end if rotations
                } // end if molecular motion
                else if (system.constants.md_mode == MD_ATOMIC) {
                    for (j=0; j<system.molecules[i].atoms.size(); j++) {
                        vsq=0;
                        for (n=0; n<3; n++) 
                            vsq += system.molecules[i].atoms[j].vel[n] * system.molecules[i].atoms[j].vel[n];
                        v_sum[z] += sqrt(vsq); // sum velocities
                        v2_sum[z] += vsq;
                        energy_holder = 0.5*system.molecules[i].atoms[j].m * vsq;
                        K_total += energy_holder;
                        Klin += energy_holder;
                    } // end for atoms j in molecule i
                } // end if atomic motion
            } // end if prototype z
        } // end molecule loop 
    } // end prototype loop

    // calculate temperature from kinetic energy and number of particles
	// https://en.wikipedia.org/wiki/Thermal_velocity
    // also McQuarrie Stat. Mech. p358 Elementary Kinetic Theory of Transport in Gases
    for (z=0; z<system.proto.size(); z++) {
        avg_v = v_sum[z] / N_local[z]; //system.stats.count_movables; //system.molecules.size(); // A/fs
        avg_v2 = v2_sum[z] / N_local[z]; //system.stats.count_movables;
        avg_v_ALL += avg_v * (N_local[z] / (double)system.stats.count_movables);
        v_rms = sqrt(avg_v2);
        // contribution to Temperature from this sorbate
        T += 1e10*avg_v*avg_v*system.proto[z].mass*M_PI/8.0/system.constants.kb * (N_local[z] / (double)system.stats.count_movables);
        T_rms += 1e10*v_rms*v_rms*system.proto[z].mass/3.0/system.constants.kb * (N_local[z] / (double)system.stats.count_movables);
        // T_rms is computed here but I'm not using it as the reported T
    }
    //printf("T_rms = %f\n", T_rms);

    // fix units
    K_total = K_total/system.constants.kb * 1e10; // convert to K
    Klin = Klin/system.constants.kb * 1e10; // ""
    Krot = Krot/system.constants.kb * 1e10; // ""
    Ek = (3.0/2.0)*system.constants.temp; // 3/2 NkT, equipartition kinetic.


    // add to partition function
    double tmp=0;
    if (T>0) {
        tmp = -(K_total+V_total)/T; // K/K = unitless
        if (tmp < 10) system.stats.Q.value += exp(tmp);
        //printf("Q += exp(-(%f+%f)/%f) = %e\n", K_total,V_total,T,exp(-(K_total+V_total)/T)); 
    }

	static double output[8];
	output[0] = K_total; 
    output[1] = V_total;
	output[2] = T;
    output[3] = avg_v_ALL;
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
        } // end if PBC
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

    // apply a constant external force if requested
    if (system.constants.md_external_force) {
        // molecular motion
        if (system.constants.md_mode == MD_MOLECULAR) {
            for (int i=0; i<system.molecules.size(); i++) {
                if (!system.molecules[i].frozen) {
                for (int n=0;n<3;n++) {
                    system.molecules[i].force[n] += system.constants.external_force_vector[n];
                }
                }
            }
        } else if (system.constants.md_mode == MD_ATOMIC) {
            for (int i=0; i<system.molecules.size(); i++) {
                for (int j=0; j<system.molecules[i].atoms.size(); j++) {
                    if (!system.molecules[i].atoms[j].frozen) {
                    for (int n=0;n<3;n++) 
                        system.molecules[i].atoms[j].force[n] += system.constants.external_force_vector[n];
                    }
                }
            }
        } // end if molecular else atomic
    } // end if EXTERNAL force
} // end force function.

// ==================== MOVE ATOMS MD STYLE =========================
/* THIS IS THE MAIN INTEGRATOR FUNCTION. calculateForces() is called within */
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

    system.checkpoint("moving particles based on forces.");
    // if molecular motion
    if (system.constants.md_mode == MD_MOLECULAR) {
        double prevangpos[3];
        for (j=0; j<system.molecules.size(); j++) {
            if (!system.molecules[j].frozen) {

            // TRANSLATION
            if (system.constants.md_translations)
                system.molecules[j].calc_pos(dt);

            // ROTATION
            // TESTING ROTATE-BY-DELTA-THETA INSTEAD OF ROTATE-BY-THETA
            // THIS IS THE BEST SETUP THUS FAR. NVE IS ALMOST CONSERVING AND THE SYSTEM IS SPATIALLY STABLE
            if (system.constants.md_rotations && system.molecules[j].atoms.size() > 1) {
            for (int n=0;n<3;n++) prevangpos[n] = system.molecules[j].ang_pos[n];
            system.molecules[j].calc_ang_pos(dt);

            // rotate molecules
            for (i=0; i<system.molecules[j].atoms.size(); i++) {
                // ROTATE IN X
                double* rotatedx = rotatePointRadians(system,
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                0, system.molecules[j].ang_pos[0] - prevangpos[0] );
                for (n=0; n<3; n++)
                    system.molecules[j].atoms[i].pos[n] = rotatedx[n] + system.molecules[j].com[n];

                // ROTATE IN Y
                double* rotatedy = rotatePointRadians(system,
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                1, system.molecules[j].ang_pos[1] - prevangpos[1] );
                for (n=0; n<3; n++)
                    system.molecules[j].atoms[i].pos[n] = rotatedy[n] + system.molecules[j].com[n];

                // ROTATE IN Z
                double* rotatedz = rotatePointRadians(system,
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                2, system.molecules[j].ang_pos[2] - prevangpos[2] );
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
    if (system.constants.md_pbc && system.constants.md_translations) {
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
    for (j=0; j<system.molecules.size(); j++) {
		if (!system.molecules[j].frozen) { // only movable atoms should move.

            // if atoms allowed to move from molecules
            if (system.constants.md_mode == MD_ATOMIC) {
                for (i=0; i<system.molecules[j].atoms.size(); i++) {
                    int nh = (system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER) ? 1 : 0;
                    system.molecules[j].atoms[i].calc_acc(nh, system.constants.lagrange_multiplier);
                    system.molecules[j].atoms[i].calc_vel(dt);
            } // end atomic loop i
            } // end if atomic
            // otherwise handle molecular movement with rigidity.
            else if (system.constants.md_mode == MD_MOLECULAR) {
                // translational
                if (system.constants.md_translations) {
                    int nh = (system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER) ? 1 : 0;
                    system.molecules[j].calc_acc(nh, system.constants.lagrange_multiplier);
                    system.molecules[j].calc_vel(dt);
                }

                // rotational
                if (system.constants.md_rotations) {
                    system.molecules[j].calc_ang_acc();
                    system.molecules[j].calc_ang_vel(dt);
                }
            }
        } // end if movable
    } // end for j molecules
    system.checkpoint("Done with a,v integration. Starting heat bath (if nvt/uvt)");

    // 5) apply heat bath in constant-temp ensembles
    if (system.constants.ensemble == ENSEMBLE_NVT || system.constants.ensemble == ENSEMBLE_UVT) {
        if (system.constants.thermostat_type == THERMOSTAT_ANDERSEN) {
        // loop through all molecules and adjust velocities by Anderson Thermostat method
        // this process makes the NVT MD simulation stochastic/ Markov / MC-like, 
        // which is usually good for obtaining equilibrium quantities.
        double probab = system.constants.md_thermostat_probab;
        double ranf;
        double sigma;
        if (system.constants.md_mode == MD_MOLECULAR && system.constants.md_translations) {
        for (i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen) continue; // skip frozens
            ranf = getrand(); // 0 -> 1
            if (ranf < probab) {
                for (int z=0; z<system.proto.size(); z++) {
                    if (system.proto[z].name == system.molecules[i].name) {
                        sigma = system.proto[z].md_velx_goal;
                        break;
                    }
                }
                // adjust the velocity components of the molecule.
                for (n=0; n<3; n++) {
                    system.molecules[i].vel[n] = gaussian(sigma);
                }
            }
        }
        } else if (system.constants.md_mode == MD_ATOMIC) {
            for (i =0; i<system.molecules.size(); i++) {
                for (j=0; j<system.molecules[i].atoms.size(); j++) {
                    if (system.molecules[i].atoms[j].frozen) continue; // skip frozen atoms
                    ranf = getrand(); // 0 -> 1
                    if (ranf <probab) {
                        for (int z=0; z<system.proto.size(); z++) {
                            if (system.proto[z].name == system.molecules[i].name) {
                                sigma = system.proto[z].md_velx_goal;
                                break;
                            }
                        }   
                        for (n=0; n<3; n++) {
                            system.molecules[i].vel[n] = gaussian(sigma);
                        }
                    }
                }
            }
        }
        }
        // Rapaport p158-159
        else if (system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER) {
            double vdotF_sum = 0;
            double mv2_sum = 0;
            if (system.constants.md_mode == MD_MOLECULAR) {
                for (i=0; i<system.molecules.size(); i++) {
                    if (system.molecules[i].frozen) continue;
                    vdotF_sum += dddotprod(system.molecules[i].vel, system.molecules[i].force);
                    mv2_sum += system.molecules[i].mass * dddotprod(system.molecules[i].vel, system.molecules[i].vel);
                }
                system.constants.lagrange_multiplier = -vdotF_sum / (mv2_sum/system.constants.kb*1e10);
            }
            else if (system.constants.md_mode == MD_ATOMIC) {
                for (i=0; i<system.molecules.size(); i++) {
                    if (system.molecules[i].frozen) continue;
                    for (j=0; j<system.molecules[i].atoms.size(); j++) {
                        vdotF_sum += dddotprod(system.molecules[i].atoms[j].vel, system.molecules[i].atoms[j].force);
                        mv2_sum += system.molecules[i].atoms[j].m * dddotprod(system.molecules[i].atoms[j].vel, system.molecules[i].atoms[j].vel);
                    }
                }
                system.constants.lagrange_multiplier = -vdotF_sum / (mv2_sum/system.constants.kb*1e10);
            }
            // the lagrange multiplier will be applied in the integration calc (for acceleration).
        } 
    } // end if uVT or NVT (thermostat)
    system.checkpoint("Done with heatbath if NVT/uVT.");
    system.checkpoint("Done with integrate() function.");
}// end integrate() function


