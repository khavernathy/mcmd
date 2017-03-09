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

            if (system.constants.md_rotations == "on") {
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
	//float dof = 6.0; // degrees of freedom. Not sure what correct value is.
	//T = K_total / ((3.0 * (float)system.constants.total_atoms - dof )); // in kelvin
	// https://en.wikipedia.org/wiki/Thermal_velocity
    T = (avg_v*1e5)*(avg_v*1e5) * system.proto.mass * M_PI / 8.0 / system.constants.kb;

    /*if (system.constants.ensemble == "nvt") {
		// NVT THERMOSTAT: Berendsen :: https://en.wikipedia.org/wiki/Berendsen_thermostat 
		double dTdt = (T - system.constants.prevtemp)/(system.constants.md_dt);
		double tau = system.constants.md_thermostat_constant;
		double T_change = T - system.constants.prevtemp;
		double prod = tau*dTdt;
		// set new velocity according to thermostat
		double scale_v = 8.0*system.constants.kb/system.proto.mass/M_PI;
		scale_v *= (T_change/tau * currtime) - 
		            (T_change/tau * (currtime-system.constants.md_dt)) + 
		            (system.constants.prevtemp);
		scale_v = sqrt(scale_v)*1e-5; // converted to A/fs
		

		//printf("dT/dt: %f;  T_change: %f;  prod: %f; scale_v: %f\n\n", dTdt, T_change,prod,scale_v);
		// END NVT THERMOSTAT
    
		// get emergent PRESSURE for NVT Frenkel p84
        pressure= system.stats.count_movables/system.constants.volume * system.constants.kb * system.constants.temp; // rho * k_b  * T
        double fsfp = system.constants.kb * (system.constants.force_sum_for_pressure / system.constants.MM_interactions); // force sum for pressure, calc'ed during force function in MD  
        pressure += (1.0/(3.0*system.constants.volume)) * fsfp;
        pressure *= 1e27 * system.constants.JL2ATM; // to atm  
    }*/
    
    // reset temperature for NVT thermostat
    //system.constants.prevtemp = T;

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

// ================ GET FORCES ON ATOMS ===============================
void calculateForces(System &system, string model, double dt) {
	
    // initialize variable for pressure calc in NVT
    system.constants.force_sum_for_pressure = 0;
    system.constants.MM_interactions = 0;
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
	
   // int countem=0;
    double eps, sig;
	for (int i = 0; i < (int)system.molecules.size(); i++) {
		for (int j = 0; j < system.molecules[i].atoms.size(); j++) {

		for (int k = 0; k<system.molecules.size(); k++) { // use k=i+1 to neglect intramolec
		for (int l = 0; l<system.molecules[k].atoms.size(); l++) {
        if (k>i || (k==i && l>j)) { // the 1/2 n^2 - n matrix of atomic interaction

                //countem++;
                // check for mixing rules
				eps = system.molecules[i].atoms[j].eps;
                sig = system.molecules[i].atoms[j].sig;

				if (system.molecules[i].atoms[j].eps != system.molecules[k].atoms[l].eps)
					eps = sqrt(eps * system.molecules[k].atoms[l].eps);
				if (system.molecules[i].atoms[j].sig != system.molecules[k].atoms[l].sig)
					sig = 0.5*(sig + system.molecules[k].atoms[l].sig);
            
                //printf("computing interaction of atomID %i and %i\n", system.molecules[i].atoms[j].ID, system.molecules[k].atoms[l].ID);
				// preliminary calculations (distance between atoms, etc.)
				double d[3],di[3],f[3],img[3],rsq,r,rimg,rimg2,u[3],s2,s6,r6, sr, sr2, sr6;
                int n,p,q;
				for (n=0; n<3; n++) 
                    d[n] = system.molecules[i].atoms[j].pos[n] - system.molecules[k].atoms[l].pos[n];
				
                if (system.constants.md_pbc == "on") {
                    // get image r in reciprocal space
                    /*for (p = 0; p <3 ; p++) {
                        img[p]=0;
                        for (q = 0; q < 3; q++) {
                            img[p] += system.pbc.reciprocal_basis[q][p]*d[q];
                        }
                        img[p] = rint(img[p]); // round to integer.
                    }

                    // matrix multiply back into real basis
                    for (p=0; p<3; p++) {
                        di[p]=0;
                        for (q=0; q<3; q++) {
                            di[p] += system.pbc.basis[q][p]*img[q];
                        }
                    }
                    // correct displacement
                    for (n=0; n<3; n++)
                        di[n] = d[n] - di[n];
                    */

                    // pythagorean
	                rsq = dddotprod(d,d);
                    //rimg2 = dddotprod(di, di);
                    r = sqrt(rsq);
                    //rimg = sqrt(rimg2);
                } else {
                    double* distances = getDistanceXYZ(system, i,j,k,l);
                    r = distances[3];
                    for (int n=0; n<3; n++) d[n] = distances[n];
                    rsq = r*r;
                }
    
                // apply 1/2 box cutoff if NVT / NVE:: p.29-30 Computer Simulation of Liquids 1991 Allen Tildesley
                // NOTE: the r_c cutoff is NOT used in MD at all..
                if (d[0] < -system.pbc.x_length/2) d[0] += system.pbc.x_length;
                if (d[0] > system.pbc.x_length/2) d[0] -= system.pbc.x_length;
                if (d[1] < -system.pbc.y_length/2) d[1] += system.pbc.y_length;
                if (d[1] > system.pbc.y_length/2) d[1] -= system.pbc.y_length;
                if (d[2] < -system.pbc.z_length/2) d[2] += system.pbc.z_length;
                if (d[2] > system.pbc.z_length/2) d[2] -= system.pbc.z_length;

                //printf("r = %f; rimg = %f\n", r, rimg);
                //printf("d: %f %f %f;\n", d[0],d[1],d[2]); //,di[0],di[1],di[2]);

				r6 = rsq*rsq*rsq;
				s2 = sig*sig;
				s6 = s2*s2*s2;
    
                if (i != k) { // don't do self-interaction for potential.
                    sr = sig/r;
                    sr2 = sr*sr;    
                    sr6 = sr2*sr2*sr2;
                }
    
                for (int n=0; n<3; n++) u[n] = d[n]/r;
			
	
				// Lennard-Jones force calculations in K/A
				for (int n=0; n<3; n++) f[n] = 24.0*d[n]*eps*(2*(s6*s6)/(r6*r6*rsq) - s6/(r6*rsq));
                for (int n=0; n<3; n++) {
                    system.molecules[i].atoms[j].force[n] += f[n];
				    system.molecules[k].atoms[l].force[n] -= f[n];
                }

				// LJ Potential in K
                if (i != k)
				    system.molecules[i].atoms[j].V += 4.0*eps*(sr6*sr6 - sr6);

				if (model == "ljes" || model == "ljespolar") {
			
				// Coulomb's law electrostatic force. Overwrite fx,fy,fz in K/A
                for (int n=0; n<3; n++) 
                    f[n] = ((system.molecules[i].atoms[j].C*system.constants.E2REDUCED * system.molecules[k].atoms[l].C*system.constants.E2REDUCED)/rsq) * u[n];
                for (int n=0; n<3; n++) {
                    system.molecules[i].atoms[j].force[n] += f[n];
                    system.molecules[k].atoms[l].force[n] -= f[n];
                }

                //Coulombic potential in K
                if (i != k)
				    system.molecules[i].atoms[j].V += ((system.molecules[i].atoms[j].C*system.constants.E2REDUCED * system.molecules[k].atoms[l].C*system.constants.E2REDUCED)/r);
				} // end coulombic addition


            } // end if not self
            } // end loop l
		    } // end loop k 
            } // end loop j
        // atomic forces are done, so now calc molecular values
        system.molecules[i].calc_force();
        if (system.constants.md_rotations == "on" && system.molecules[i].atoms.size() > 1) 
            system.molecules[i].calc_torque();


        // also compute force * r dot products for pressure calc in NVT
        /*if (system.constants.ensemble == "nvt") {
            for (int j=0; j<system.molecules[i].atoms.size(); j++) {
                for (int k=0; k<system.molecules.size(); k++) {
                    for (int l=0; l<system.molecules[k].atoms.size(); l++) {
                        if (k>i || (k==i && l>j)) {
                           if (system.molecules[i].MF == "M" && system.molecules[k].MF == "M") {
                                double* dist = getDistanceXYZ(system, i,j,k,l);
                                for (int n=0; n<3; n++) {
                                    system.constants.force_sum_for_pressure += system.molecules[i].atoms[j].force[n] * dist[n];
                                    system.constants.force_sum_for_pressure += system.molecules[k].atoms[l].force[n] * dist[n];
                                    system.constants.MM_interactions += 2;
                                }
                            }     
                        } 
                    }
                }
            }
        } // end if NVT
        */
    } // end loop i

//printf("count: %i\n",countem);
}

// ==================== MOVE ATOMS MD STYLE =========================
/* THIS IS THE MAIN LOOPING FUNCTION. calculateForces() is called within */
void integrate(System &system, double dt) {

    // DEBUG
    int debug=0;
    if (debug == 1) {
        for (int j=0; j<system.molecules.size(); j++) {
            if (system.constants.md_mode == "molecular") system.molecules[j].printAll();
            for (int i=0; i<system.molecules[j].atoms.size(); i++) {
                if (system.constants.md_mode == "atomic") system.molecules[j].atoms[i].printAll();
            }
        }
    }
    // END IF DEBUG

    // 1a) CHANGE POSITIONS OF PARTICLES
    // save old positions
    for (int j=0; j<system.molecules.size(); j++) {
        if (system.molecules[j].MF == "M") {
        for (int i=0; i<system.molecules[j].atoms.size(); i++) {
            for (int n=0; n<3; n++) system.molecules[j].atoms[i].prevpos[n] = system.molecules[j].atoms[i].pos[n];
        } // end for atom i
        } // end if movable
    } // end for molecule j
    // done saving old positions

    // if molecular motion
    if (system.constants.md_mode == "molecular") {
        for (int j=0; j<system.molecules.size(); j++) {
            if (system.molecules[j].MF == "M") {

            system.molecules[j].calc_pos(dt);
            
              // ROTATION
            if (system.constants.md_rotations == "on" && system.molecules[j].atoms.size() > 1) {
            system.molecules[j].calc_ang_pos(dt);

            // rotate molecules
            for (int i=0; i<system.molecules[j].atoms.size(); i++) {
                // ROTATE IN X
                double* rotatedx = rotatePoint(system, 
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                "x", system.molecules[j].ang_pos[0] * 180.0/M_PI); 
                for (int n=0; n<3; n++) 
                    system.molecules[j].atoms[i].pos[n] = rotatedx[n] + system.molecules[j].com[n];

                // ROTATE IN Y
                double* rotatedy = rotatePoint(system, 
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                "y", system.molecules[j].ang_pos[1] * 180.0/M_PI); 
                for (int n=0; n<3; n++) 
                    system.molecules[j].atoms[i].pos[n] = rotatedy[n] + system.molecules[j].com[n];

                // ROTATE IN Z
                double* rotatedz = rotatePoint(system, 
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                "z", system.molecules[j].ang_pos[2] * 180.0/M_PI); 
                for (int n=0; n<3; n++) 
                    system.molecules[j].atoms[i].pos[n] = rotatedz[n] + system.molecules[j].com[n];
            } // end loop over atoms i 
            } // end if rotations allowed and >1 atom
            } // end if movable molecule
        } // end for molecules j
    } // end if molecular motion
    // if atomic motion
    else if (system.constants.md_mode == "atomic") {
        for (int j=0; j<system.molecules.size(); j++) {
            if (system.molecules[j].MF == "M") {
            for (int i=0; i<system.molecules[j].atoms.size(); i++) {
                system.molecules[j].atoms[i].calc_pos(dt);
            }
            }
        }
    } // end if atomic motion
    // END POSITION CHANGES

    // 1b) CHECK P.B.C. (move the molecule/atom back in the box if needed)
    if (system.constants.md_pbc == "on") {
        for (int j=0; j<system.molecules.size(); j++) {
            if (system.molecules[j].MF == "M") {
                // HEREIN LIES THE PROBLEM.
                checkInTheBox(system,j);
            } // end if movable
	    } // end loop j molecules
    } // end if PBC
   
    // 2) GET NEW C.O.M. for new positions.
    for (int j=0; j<system.molecules.size(); j++) system.molecules[j].calc_center_of_mass();
 
    // 3) GET NEW FORCES (AND TORQUES) BASED ON NEW POSITIONS
	calculateForces(system, system.constants.potential_form, dt);

    // 4) GET NEW ACCELERATION AND VELOCITY FOR ALL PARTICLES
	for (int j=0; j<system.molecules.size(); j++) {
		if (system.molecules[j].MF == "M") { // only movable atoms should move.

            // if atoms allowed to move from molecules
            if (system.constants.md_mode == "atomic") {
                for (int i=0; i<system.molecules[j].atoms.size(); i++) {
                    system.molecules[j].atoms[i].calc_acc();
                    system.molecules[j].atoms[i].calc_vel(dt); 
            } // end atomic loop i
            } // end if atomic
            // otherwise handle molecular movement with rigidity.
            else if (system.constants.md_mode == "molecular") {
                // linear
                system.molecules[j].calc_acc();
                system.molecules[j].calc_vel(dt, system.constants.md_init_vel);
    
                // rotational
                if (system.constants.md_rotations == "on") {
                    system.molecules[j].calc_ang_acc();
                    system.molecules[j].calc_ang_vel(dt);
                }
            }
        } // end if movable
    	} // end for j molecules
}// end integrate() function
