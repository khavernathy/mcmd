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
    // highlighted in yellow..
}


void calculateForces(System &system, double dt) {
  int_fast8_t model = system.constants.potential_form;

    //system.checkpoint("Calculating forces.");
    // loop through all atoms
	for (int i=0; i <system.molecules.size(); i++) {
	for (int j = 0; j < system.molecules[i].atoms.size(); j++) {
		// initialize stuff
		system.molecules[i].atoms[j].force[0] = 0.0;
		system.molecules[i].atoms[j].force[1] = 0.0;
		system.molecules[i].atoms[j].force[2] = 0.0;
		system.molecules[i].atoms[j].V = 0.0;
	}
	}
    if (system.constants.calc_pressure_option)
        system.constants.fdotr_sum = 0.0;

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
            /* REPULSION -- DISPERSION */
            if (model == POTENTIAL_LJ || model == POTENTIAL_LJES || model == POTENTIAL_LJESPOLAR || model == POTENTIAL_LJPOLAR) {
                #ifdef OMP
                if (system.constants.openmp_threads > 0)
                    lj_force_omp(system);
                else
                    lj_force(system);
                #else
                lj_force(system);
                #endif  
            } else if (model == POTENTIAL_TT || model == POTENTIAL_TTES || model == POTENTIAL_TTESPOLAR)
                tt_forces(system);

            /* ELECTROSTATICS */
            if (model == POTENTIAL_LJES || model == POTENTIAL_LJESPOLAR || model == POTENTIAL_TTES || model == POTENTIAL_TTESPOLAR) {
                #ifdef OMP
                if (system.constants.openmp_threads > 0)
                    coulombic_real_force_omp(system);
                else
                    coulombic_real_force(system);
                #else
                coulombic_real_force(system);
                #endif
            }

            /* POLARIZATION */    
            if (model == POTENTIAL_LJESPOLAR || model == POTENTIAL_LJPOLAR || model == POTENTIAL_TTESPOLAR)
                #ifdef OMP
                if (system.constants.openmp_threads > 0)
                    polarization_force_omp(system);
                else
                    polarization_force(system);
                #else
                polarization_force(system);
                #endif

            /* BONDING */
            if (system.constants.flexible_frozen || system.constants.md_mode == MD_FLEXIBLE) {
                  if (system.constants.opt_bonds)
                    morse_gradient(system);
                  if (system.constants.opt_angles)
                    angle_bend_gradient(system);
                  if (system.constants.opt_dihedrals)
                    torsions_gradient(system);
                  if (system.constants.opt_LJ)
                    LJ_intramolec_gradient(system);
                  if ((model == POTENTIAL_LJES || model == POTENTIAL_LJESPOLAR) && system.constants.opt_ES) // only if electrostatics is on.
                    ES_intramolec_gradient(system);
            }
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
    if (system.constants.md_external_force && system.stats.MDstep % system.constants.md_external_force_freq == 0) {
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
    //system.checkpoint("Done calculating forces.");
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

    // 1) MOVE ATOMS FROM FORCES/TORQUES
    system.checkpoint("moving particles based on forces.");
    // if molecular motion
    if (system.constants.md_mode == MD_MOLECULAR) {
        double prevangpos[3];
        for (j=0; j<system.molecules.size(); j++) {
            if (!system.molecules[j].frozen) {
            // TRANSLATION
            if (system.constants.md_translations) {
                if (system.constants.integrator == INTEGRATOR_VV)
                    system.molecules[j].calc_pos(dt);
                else if (system.constants.integrator == INTEGRATOR_RK4) {
                    // GOOD FOR SINGLE-ATOM MOLECULES ONLY
                    double k1,k2,k3,k4,dxt,tmp_pos[3];
                    // tmp_pos stores the original position of atom
                    for (int n=0;n<3;n++)
                        tmp_pos[n] = system.molecules[j].atoms[0].pos[n];

                    for (int n=0;n<3;n++) {
                        k1 = dt*system.molecules[j].vel[n];
                        system.molecules[j].atoms[0].pos[n] = tmp_pos[n] + 0.5*k1;
                        singleAtomForceLJ(system,j,0);
                        dxt = system.molecules[j].vel[n] + system.molecules[j].atoms[0].force[n]*1.3806488e-33/system.molecules[j].atoms[0].m*(0.5*dt);
                        
                        k2 = dt*dxt;
                        system.molecules[j].atoms[0].pos[n] = tmp_pos[n] + 0.5*k2;
                        singleAtomForceLJ(system,j,0);
                        dxt = system.molecules[j].vel[n] + system.molecules[j].atoms[0].force[n]*1.3806488e-33/system.molecules[j].atoms[0].m*(0.5*dt);

                        k3 = dt*dxt;
                        system.molecules[j].atoms[0].pos[n] = tmp_pos[n] + k3;
                        singleAtomForceLJ(system,j,0);
                        dxt = system.molecules[j].vel[n] + system.molecules[j].atoms[0].force[n]*1.3806488e-33/system.molecules[j].atoms[0].m*dt;

                        k4 = dt*dxt;
                        system.molecules[j].atoms[0].pos[n] = tmp_pos[n] + (k1 + 2.0*(k2 + k3) + k4)/6.0;
                    }
                }
            }

            // ROTATION
            // TESTING ROTATE-BY-DELTA-THETA INSTEAD OF ROTATE-BY-THETA
            // THIS IS THE BEST SETUP THUS FAR. NVE IS ALMOST CONSERVING AND THE SYSTEM IS SPATIALLY STABLE
            if (system.constants.md_rotations && system.molecules[j].atoms.size() > 1) {
            for (n=0;n<3;n++) prevangpos[n] = system.molecules[j].ang_pos[n];
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
                /*
                 * TODO: apply constraints
                 * within atom class
                 * can call local map of all bonds
                 */
            }
            }
        }
    } // end if atomic motion
    else if (system.constants.md_mode == MD_FLEXIBLE) {
        for (i=0;i<system.molecules.size();i++) {
            if (system.molecules[i].frozen) continue;
            for (j=0;j<system.molecules[i].atoms.size();j++) {
                system.molecules[i].atoms[j].calc_pos(dt);
            }
        }
    }

    // flexible "frozen" (MOF) atoms integration
    if (system.constants.flexible_frozen) {
        // ASSUMES THERE IS ONLY ONE "FROZEN" (MOF) MOLECULE
        for (i=0; i<system.molecules[0].atoms.size(); i++) {
            system.molecules[0].atoms[i].calc_pos(dt);
        }
    }
    // end if frozen (MOF) is flexible
    // END POSITION CHANGES
    system.checkpoint("done moving particles. Checking PBC for all particles");

    // TODO -- MOVE THE BOX PARAMS IF MOF FLEXES NEAR/OUTSIDE THE BOX?
    // 2) CHECK P.B.C. (move the molecule/atom back in the box if needed)
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
    // 4a) MOVABLES
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
            // or flexible molecules
            else if (system.constants.md_mode == MD_FLEXIBLE) {
                int nh = (system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER) ? 1 : 0;
                    for (i=0;i<system.molecules[j].atoms.size();i++) {
                        system.molecules[j].atoms[i].calc_acc(nh, system.constants.lagrange_multiplier);
                        system.molecules[j].atoms[i].calc_vel(dt);
                    }
            }
        } // end if movable
    } // end for j molecules
    // 4b) Flexible "frozen" (MOF) atoms
    if (system.constants.flexible_frozen) {
        int nh = (system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER) ? 1 : 0;
        for (j=0; j<system.molecules[0].atoms.size(); j++) {
          system.molecules[0].atoms[j].calc_acc(nh, system.constants.lagrange_multiplier);
          system.molecules[0].atoms[j].calc_vel(dt);
        }
    } // end flexible MOF dynamics
    system.checkpoint("Done with a,v integration. Starting heat bath (if nvt/uvt)");

    // TODO -- flexible MOF thermostat?
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
            if (system.molecules[i].frozen && !system.constants.flexible_frozen) continue; // skip frozens
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
                    if (system.molecules[i].atoms[j].frozen && !system.constants.flexible_frozen) continue; // skip frozen atoms
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
        // end if atomic mode
        } else if (system.constants.md_mode == MD_FLEXIBLE) {
            for (i=0; i<system.molecules.size(); i++) {
                if (system.molecules[i].frozen && !system.constants.flexible_frozen) continue;
                for (j=0; j<system.molecules[i].atoms.size();j++) {
                    ranf = getrand();
                    if (ranf < probab) {
                        sigma = system.molecules[i].atoms[j].md_velx_goal;
                        break;
                    } 
                    for (n=0;n<3;n++) {
                        system.molecules[i].atoms[j].vel[n] = gaussian(sigma);
                    }
                }
            } 
        }
        } // end Andersen thermostat
        // Rapaport p158-159
        else if (system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER) {
            double vdotF_sum = 0;
            double mv2_sum = 0;
            if (system.constants.md_mode == MD_MOLECULAR) {
                for (i=0; i<system.molecules.size(); i++) {
                    if (system.molecules[i].frozen && !system.constants.flexible_frozen) continue;
                    vdotF_sum += dddotprod(system.molecules[i].vel, system.molecules[i].force);
                    mv2_sum += system.molecules[i].mass * dddotprod(system.molecules[i].vel, system.molecules[i].vel);
                }
                system.constants.lagrange_multiplier = -vdotF_sum / (mv2_sum/system.constants.kb*1e10);
            }
            else if (system.constants.md_mode == MD_ATOMIC || system.constants.md_mode == MD_FLEXIBLE) {
                for (i=0; i<system.molecules.size(); i++) {
                    if (system.molecules[i].frozen && !system.constants.flexible_frozen) continue;
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
