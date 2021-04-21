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


void calculateForces(System &system) {
    double dt = system.constants.md_dt;
    int_fast8_t model = system.constants.potential_form;

    system.checkpoint("Calculating forces.");
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
            system.checkpoint("Computing RD force.");
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
            system.checkpoint("Computing ES force.");
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
            system.checkpoint("Computing Polarization force.");
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
            system.checkpoint("Computing bonded forces.");
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
                    for (int n=0; n<3; n++) {
                        system.molecules[i].force[n] += system.constants.external_force_vector[n];
                    }
                }
            }
        } else if (system.constants.md_mode == MD_ATOMIC || system.constants.md_mode == MD_FLEXIBLE) {
            for (int i=0; i<system.molecules.size(); i++) {
                for (int j=0; j<system.molecules[i].atoms.size(); j++) {
                    if (!system.molecules[i].atoms[j].frozen) {
                        for (int n=0; n<3; n++)
                            system.molecules[i].atoms[j].force[n] += system.constants.external_force_vector[n];
                    }
                }
            }
        } // end if molecular else atomic
    } // end if EXTERNAL force
    system.checkpoint("Done calculating forces.");
} // end force function.



void acceleration_velocity_verlet(System &system) {
    double dt = system.constants.md_dt;
    int i,j;
    for (j=0; j<system.molecules.size(); j++) {
        if (!system.molecules[j].frozen) { // only movable atoms should move.

            // if atoms allowed to move from molecules
            if (system.constants.md_mode == MD_ATOMIC) {
                for (i=0; i<system.molecules[j].atoms.size(); i++) {
                    int nh = (system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER) ? 1 : 0;
                    system.molecules[j].atoms[i].calc_vel_verlet(dt, nh, system.constants.lagrange_multiplier);
                } // end atomic loop i
            } // end if atomic
            // otherwise handle molecular movement with rigidity.
            else if (system.constants.md_mode == MD_MOLECULAR) {
                // translational
                if (system.constants.md_translations) {
                    int nh = (system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER) ? 1 : 0;
                    system.molecules[j].calc_vel_verlet(dt, nh, system.constants.lagrange_multiplier);
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
                for (i=0; i<system.molecules[j].atoms.size(); i++) {
                    system.molecules[j].atoms[i].calc_vel_verlet(dt, nh, system.constants.lagrange_multiplier);
                }
            }
        } // end if movable
    } // end for j molecules
    if (system.constants.flexible_frozen) {
        int nh = (system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER) ? 1 : 0;
        for (j=0; j<system.molecules[0].atoms.size(); j++) {
            system.molecules[0].atoms[j].calc_vel_verlet(dt, nh, system.constants.lagrange_multiplier);
        }
    } // end flexible MOF dynamics
}
// end a and v calculator

void acceleration_velocity_RK4(System &system) {
    double dt = system.constants.md_dt;
    int i,j;
    int nh=0;

    if (system.constants.ensemble==ENSEMBLE_NVT && system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER)
        nh=1;

    for (j=0; j<system.molecules.size(); j++) {
        if (!system.molecules[j].frozen) { // only movable atoms should move.

            // if atoms allowed to move from molecules
            if (system.constants.md_mode == MD_ATOMIC) {
                for (i=0; i<system.molecules[j].atoms.size(); i++) {
                    system.molecules[j].atoms[i].calc_vel(dt,nh,system.constants.lagrange_multiplier);
                } // end atomic loop i
            } // end if atomic
            // otherwise handle molecular movement with rigidity.
            else {
                printf("ERROR: RK4 is not available for molecular or flexible md_mode. Use `md_mode atomic` instead.\n");
                exit(EXIT_FAILURE);
            }
        } // end if movable
    } // end for j molecules
}


void NVT_thermostat_andersen(System &system) {
    int i,j,n;
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
                    for (j=0; j<system.molecules[i].atoms.size(); j++) {
                        ranf = getrand(); // 0 -> 1
                        if (ranf < probab) {
                            sigma = system.molecules[i].atoms[j].md_velx_goal;
                            for (n=0; n<3; n++) {
                                system.molecules[i].atoms[j].vel[n] = gaussian(sigma);
                            }
                        }
                    }
                }
            } // end if flexible
        } // end Andersen thermostat
    } // end if uVT or NVT (thermostat)
}
// end thermostat function

/*
void NVT_NoseHoover_LM(System &system) {
        int i,j;
        // Rapaport p158-159
        if (system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER) {
            double vdotF_sum = 0;
            double mv2_sum = 0;
            if (system.constants.md_mode == MD_MOLECULAR) {
                for (i=0; i<system.molecules.size(); i++) {
                    if (system.molecules[i].frozen && !system.constants.flexible_frozen) continue;
                    vdotF_sum += dddotprod(system.molecules[i].vel, system.molecules[i].force);
                    mv2_sum += system.molecules[i].mass*system.constants.amu2kg * dddotprod(system.molecules[i].vel, system.molecules[i].vel);
                }
                system.constants.lagrange_multiplier = -vdotF_sum / (mv2_sum/system.constants.kb*1e10);
            }
            else if (system.constants.md_mode == MD_ATOMIC || system.constants.md_mode == MD_FLEXIBLE) {
                for (i=0; i<system.molecules.size(); i++) {
                    if (system.molecules[i].frozen && !system.constants.flexible_frozen) continue;
                    for (j=0; j<system.molecules[i].atoms.size(); j++) {
                        vdotF_sum += dddotprod(system.molecules[i].atoms[j].vel, system.molecules[i].atoms[j].force);
                        mv2_sum += system.molecules[i].atoms[j].mass*system.constants.amu2kg * dddotprod(system.molecules[i].atoms[j].vel, system.molecules[i].atoms[j].vel);
                    }
                }
                system.constants.lagrange_multiplier = -vdotF_sum / (mv2_sum/system.constants.kb*1e10);
            }
        }
}
*/

void calculateNHLM_now(System &system) {
    // Rapaport p158-159
    unsigned int i,j;
    double vdotF_sum = 0;
    double mv2_sum = 0;
    if (system.constants.md_mode == MD_MOLECULAR) {
        for (i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen && !system.constants.flexible_frozen) continue;
            vdotF_sum += dddotprod(system.molecules[i].vel, system.molecules[i].force);
            mv2_sum += system.molecules[i].mass*system.constants.amu2kg * dddotprod(system.molecules[i].vel, system.molecules[i].vel);
        }
        system.constants.lagrange_multiplier = -vdotF_sum / (mv2_sum/system.constants.kb*1e10);
    }
    else if (system.constants.md_mode == MD_ATOMIC || system.constants.md_mode == MD_FLEXIBLE) {
        for (i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen && !system.constants.flexible_frozen) continue;
            for (j=0; j<system.molecules[i].atoms.size(); j++) {
                vdotF_sum += dddotprod(system.molecules[i].atoms[j].vel, system.molecules[i].atoms[j].force);
                mv2_sum += system.molecules[i].atoms[j].mass*system.constants.amu2kg * dddotprod(system.molecules[i].atoms[j].vel, system.molecules[i].atoms[j].vel);
            }
        }
        system.constants.lagrange_multiplier = -vdotF_sum / (mv2_sum/system.constants.kb*1e10);
    }
    // units of LM are 1/fs
}

void calculateNH_Q(System &system) {
    double dt=system.constants.md_dt;
    double dof = system.constants.DOF; // computed in main()
    double tau = 20.0*dt; // a minimum to pick is 20dt
    double Q = dof*system.constants.temp*tau*tau;  // units Kfs^2. http://www.courses.physics.helsinki.fi/fys/moldyn/lectures/L5.pdf
    system.constants.NH_Q = Q * system.constants.NH_Q_scale;
    printf("Calculated Nose-Hoover Q parameter = %f K fs^2 with scale factor = %f.\n",system.constants.NH_Q,system.constants.NH_Q_scale);
}

void position_RK4(System &system) {
    // GOOD FOR SINGLE-ATOM MOLECULES ONLY
    double k1,k2,k3,k4,dxt,tmp_pos[3];
    double dt = system.constants.md_dt;
    // tmp_pos stores the original position of atom
    for (int j=0; j<system.molecules.size(); j++) {
        for (int n=0; n<3; n++)
            tmp_pos[n] = system.molecules[j].atoms[0].pos[n];

        for (int n=0; n<3; n++) {
            k1 = dt*system.molecules[j].vel[n];
            system.molecules[j].atoms[0].pos[n] = tmp_pos[n] + 0.5*k1;
            singleAtomForceLJ(system,j,0);
            dxt = system.molecules[j].vel[n] + system.molecules[j].atoms[0].force[n]*system.constants.KA2Afs2/system.molecules[j].atoms[0].mass*(0.5*dt);

            k2 = dt*dxt;
            system.molecules[j].atoms[0].pos[n] = tmp_pos[n] + 0.5*k2;
            singleAtomForceLJ(system,j,0);
            dxt = system.molecules[j].vel[n] + system.molecules[j].atoms[0].force[n]*system.constants.KA2Afs2/system.molecules[j].atoms[0].mass*(0.5*dt);

            k3 = dt*dxt;
            system.molecules[j].atoms[0].pos[n] = tmp_pos[n] + k3;
            singleAtomForceLJ(system,j,0);
            dxt = system.molecules[j].vel[n] + system.molecules[j].atoms[0].force[n]*system.constants.KA2Afs2/system.molecules[j].atoms[0].mass*dt;

            k4 = dt*dxt;
            system.molecules[j].atoms[0].pos[n] = tmp_pos[n] + (k1 + 2.0*(k2 + k3) + k4)/6.0;
        }
    }
}

void position_verlet(System &system) {
    int i,j,n;
    double dt = system.constants.md_dt;
    // if molecular motion
    if (system.constants.md_mode == MD_MOLECULAR) {
        double prevangpos[3];
        for (j=0; j<system.molecules.size(); j++) {
            if (!system.molecules[j].frozen) {
                // TRANSLATION
                if (system.constants.md_translations) {
                    if (system.constants.integrator == INTEGRATOR_VV)
                        system.molecules[j].calc_pos(dt);
                }

                // ROTATION
                // TESTING ROTATE-BY-DELTA-THETA INSTEAD OF ROTATE-BY-THETA
                // THIS IS THE BEST SETUP THUS FAR. NVE IS ALMOST CONSERVING AND THE SYSTEM IS SPATIALLY STABLE
                if (system.constants.md_rotations && system.molecules[j].atoms.size() > 1) {
                    for (n=0; n<3; n++) prevangpos[n] = system.molecules[j].ang_pos[n];
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
    // if atomic/flexible motion
    else if (system.constants.md_mode == MD_ATOMIC || system.constants.md_mode == MD_FLEXIBLE) {
        for (i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen) continue;
            for (j=0; j<system.molecules[i].atoms.size(); j++) {
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
}


void position_VV_NH(System &system) {
    int i,j,n;
    double dt = system.constants.md_dt;
    // if molecular motion
    if (system.constants.md_mode == MD_MOLECULAR) {
        double prevangpos[3];
        for (j=0; j<system.molecules.size(); j++) {
            if (!system.molecules[j].frozen) {
                // TRANSLATION
                if (system.constants.md_translations) {
                    if (system.constants.integrator == INTEGRATOR_VV)
                        system.molecules[j].calc_pos_VV_NH(dt,system.constants.lagrange_multiplier);
                }

                // ROTATION
                // TESTING ROTATE-BY-DELTA-THETA INSTEAD OF ROTATE-BY-THETA
                // THIS IS THE BEST SETUP THUS FAR. NVE IS ALMOST CONSERVING AND THE SYSTEM IS SPATIALLY STABLE
                if (system.constants.md_rotations && system.molecules[j].atoms.size() > 1) {
                    for (n=0; n<3; n++) prevangpos[n] = system.molecules[j].ang_pos[n];
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
    // if atomic/flexible motion
    else if (system.constants.md_mode == MD_ATOMIC || system.constants.md_mode == MD_FLEXIBLE) {
        for (i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen) continue;
            for (j=0; j<system.molecules[i].atoms.size(); j++) {
                system.molecules[i].atoms[j].calc_pos_VV_NH(dt, system.constants.lagrange_multiplier);
            }
        }
    }
    // flexible "frozen" (MOF) atoms integration
    if (system.constants.flexible_frozen) {
        // ASSUMES THERE IS ONLY ONE "FROZEN" (MOF) MOLECULE
        for (i=0; i<system.molecules[0].atoms.size(); i++) {
            system.molecules[0].atoms[i].calc_pos_VV_NH(dt, system.constants.lagrange_multiplier);
        }
    }
    // end if frozen (MOF) is flexible
}

void velocity_VV_NH_final(System &system) {
    // the last step in velocity verlet for Nose Hoover thermostat
    unsigned int i,j;
    double dt = system.constants.md_dt;
    if (system.constants.md_mode==MD_MOLECULAR) { // rigid molecules
        for (i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen) continue;
            system.molecules[i].calc_vel_VV_NH_final(dt,system.constants.lagrange_multiplier);
        }
    } else { // atomic/flexible systems
        for (i=0; i<system.molecules.size(); i++) {
            if (system.molecules[i].frozen && !system.constants.flexible_frozen) continue;
            for (j=0; j<system.molecules[i].atoms.size(); j++) {
                system.molecules[i].atoms[j].calc_vel_VV_NH_final(dt,system.constants.lagrange_multiplier);
            }
        }
    }
}


void doPBCcheck(System &system) {
    unsigned int j;
    if (system.constants.md_pbc && system.constants.md_translations) {
        for (j=0; j<system.molecules.size(); j++) {
            if (!system.molecules[j].frozen) {
                checkInTheBox(system,j); // also computes COM
            } // end if movable
        } // end loop j molecules
    } // end if PBC
}

void updateLM(System &system, int s) {
    unsigned int i,j;
    double dt=system.constants.md_dt;
    double dof = system.constants.DOF; // computed in main()
    double Q = system.constants.NH_Q;  // units Kfs^2. http://www.courses.physics.helsinki.fi/fys/moldyn/lectures/L5.pdf
    double chunk = 0.5*dof*system.constants.temp; // in K. 1/2 per DOF
    double dto2q = 0.5*dt/Q; // units 1/(K*fs)
    double mv2sum=0;
    if (s==0) { // first half step (i.e. use the original v(t))
        if (system.constants.md_mode==MD_MOLECULAR) { // rigid molecules
            for (i=0; i<system.molecules.size(); i++) {
                mv2sum += system.molecules[i].mass * dddotprod(system.molecules[i].ov,system.molecules[i].ov);
            }
        } else { // atomic/flexible systems
            for (i=0; i<system.molecules.size(); i++) {
                for (j=0; j<system.molecules[i].atoms.size(); j++) {
                    mv2sum += system.molecules[i].atoms[j].mass * dddotprod(system.molecules[i].atoms[j].ov, system.molecules[i].atoms[j].ov);
                }
            }
        }
    } else if (s==1) { // next to full dt step (i.e. use the updated v(t+dt/2))
        if (system.constants.md_mode==MD_MOLECULAR) { // rigid molecules
            for (i=0; i<system.molecules.size(); i++) {
                mv2sum += system.molecules[i].mass * dddotprod(system.molecules[i].vel,system.molecules[i].vel);
            }
        } else { // atomic/flexible systems
            for (i=0; i<system.molecules.size(); i++) {
                for (j=0; j<system.molecules[i].atoms.size(); j++) {
                    mv2sum += system.molecules[i].atoms[j].mass * dddotprod(system.molecules[i].atoms[j].vel, system.molecules[i].atoms[j].vel);
                }
            }
        }
    }
    mv2sum *= 0.5 * system.constants.reduced2K; // ,  0.5 mv^2. convert amu*A^2/fs^2 --> K

    system.constants.lagrange_multiplier += dto2q*(mv2sum - chunk); // so units are 1/fs
}


