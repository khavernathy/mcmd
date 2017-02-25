#include <string>
#include <algorithm>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

// =================  GET TOTAL ENERGY AND EMERGENT TEMPERATURE FROM SYSTEM STATE ===========================
double * calculateEnergyAndTemp(System &system) { // the * is to return an array of doubles as a pointer, not just one double
	double V_total = 0.0;
    double K_total = 0.0, Klin=0, Krot=0, Ek=0.0;
    double v_sum=0.0, avg_v = 0.0;
	double T=0.0;
	//long unsigned int size = system.atoms.size();
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

        // iteratively sum potential
        for (int i=0; i<system.molecules[j].atoms.size(); i++) V_total += system.molecules[j].atoms[i].V;
		//total energy
		// iteratively sum total energy
	//}
	}

    avg_v = v_sum / system.molecules.size();
    K_total = K_total / system.constants.kb * 1e10; // convert to K 
    Klin = Klin / system.constants.kb * 1e10;
    Krot = Krot / system.constants.kb * 1e10;
    Ek = (3.0/2.0)*system.constants.temp; // 3/2 NkT, equipartition kinetic.

	// calculate temperature from kinetic energy and number of particles
	//float dof = 6.0; // degrees of freedom. Not sure what correct value is.
	//T = K_total / ((3.0 * (float)system.constants.total_atoms - dof )); // in kelvin
	// https://en.wikipedia.org/wiki/Thermal_velocity
    T = (avg_v*1e5)*(avg_v*1e5) * system.proto.mass * M_PI / 8.0 / system.constants.kb;

    //printf("Temperature: %4.2fK ",T); 

	static double output[5];
	output[0] = K_total;
    output[1] = V_total;
	output[2] = T;
    output[3] = avg_v;
    output[4] = Ek;
    output[5] = Klin;
    output[6] = Krot;
	return output;
}

// ================ GET FORCES ON ATOMS ===============================
void calculateForces(System &system, string model, double dt) {
	// loop through all atoms
	for (int j=0; j <system.molecules.size(); j++) {
	for (int i = 0; i < system.molecules[j].atoms.size(); i++) {
		// initialize force and potential to zero.
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
				double dx,dy,dz,rsq,r,ux,uy,uz,s2,s6,r6, sr, sr2, sr6;
				dx = system.molecules[i].atoms[j].pos[0] - system.molecules[k].atoms[l].pos[0];
				dy = system.molecules[i].atoms[j].pos[1] - system.molecules[k].atoms[l].pos[1];
				dz = system.molecules[i].atoms[j].pos[2] - system.molecules[k].atoms[l].pos[2];
				
                // apply 1/2 box cutoff:: p.29-30 Computer Simulation of Liquids 1991 Allen Tildesley
                if (dx > system.constants.x_max) { dx = dx - system.constants.x_length; }
                if (dx < system.constants.x_min) { dx = dx + system.constants.x_length; }
                if (dy > system.constants.y_max) { dy = dy - system.constants.y_length; }
                if (dy < system.constants.y_min) { dy = dy + system.constants.y_length; }
                if (dz > system.constants.z_max) { dz = dz - system.constants.z_length; }
                if (dz < system.constants.z_min) { dz = dz + system.constants.z_length; }
	
				rsq = dx*dx + dy*dy + dz*dz;
				r6 = rsq*rsq*rsq;
				r = sqrt(rsq);
				s2 = sig*sig;
				s6 = s2*s2*s2;
    
                if (i != k) { // don't do self-interaction for potential.
                    sr = sig/r;
                    sr2 = sr*sr;    
                    sr6 = sr2*sr2*sr2;
                }

				ux = dx/r;
				uy = dy/r;
				uz = dz/r;
				
				// Lennard-Jones force calculations
				double fx,fy,fz;
				fx = 24.0*dx*eps*(2*(s6*s6)/(r6*r6*rsq) - s6/(r6*rsq));
				fy = 24.0*dy*eps*(2*(s6*s6)/(r6*r6*rsq) - s6/(r6*rsq));  //   (2*pow(sig,12)*pow(r,-14) - pow(sig,6) * pow(r,-8));
				fz = 24.0*dz*eps*(2*(s6*s6)/(r6*r6*rsq) - s6/(r6*rsq));//    (2*pow(sig,12)*pow(r,-14) - pow(sig,6) * pow(r,-8));
				
				system.molecules[i].atoms[j].force[0] += fx;
				system.molecules[i].atoms[j].force[1] += fy;
				system.molecules[i].atoms[j].force[2] += fz;

				system.molecules[k].atoms[l].force[0] -= fx;
				system.molecules[k].atoms[l].force[1] -= fy;
				system.molecules[k].atoms[l].force[2] -= fz;

				// LJ Potential
                if (i != k)
				    system.molecules[i].atoms[j].V += 4.0*eps*(sr6*sr6 - sr6);

				if (model == "ljes" || model == "ljespolar") {
			
				// Coulomb's law electrostatic force. Overwrite fx,fy,fz
				fx = ((system.molecules[i].atoms[j].C*system.constants.E2REDUCED * system.molecules[k].atoms[l].C*system.constants.E2REDUCED)/rsq) * ux;
				fy = ((system.molecules[i].atoms[j].C*system.constants.E2REDUCED * system.molecules[k].atoms[l].C*system.constants.E2REDUCED)/rsq) * uy;
				fz = ((system.molecules[i].atoms[j].C*system.constants.E2REDUCED * system.molecules[k].atoms[l].C*system.constants.E2REDUCED)/rsq) * uz;	

				system.molecules[i].atoms[j].force[0] += fx;
				system.molecules[i].atoms[j].force[1] += fy;
                system.molecules[i].atoms[j].force[2] += fz;

                system.molecules[k].atoms[l].force[0] -= fx;
                system.molecules[k].atoms[l].force[1] -= fy;
                system.molecules[k].atoms[l].force[2] -= fz;
				
                //Coulombic potential
                if (i != k)
				    system.molecules[i].atoms[j].V += ((system.molecules[i].atoms[j].C*system.constants.E2REDUCED * system.molecules[k].atoms[l].C*system.constants.E2REDUCED)/r);
        
				} // end coulombic addition
            } // end if not self
            } // end loop l
		    } // end loop k 
            } // end loop j
        // atomic forces are done, so now compile molecular values
        system.molecules[i].calc_force();
        if (system.constants.md_rotations == "on" && system.molecules[i].atoms.size() > 1) system.molecules[i].calc_torque();
        //if (system.molecules[i].atoms.size() > 1) system.molecules[i].calc_rotations(dt);
    } // end loop i
//printf("count: %i\n",countem);
}

// ==================== MOVE ATOMS MD STYLE =========================
/* THIS IS THE MAIN LOOPING FUNCTION. calculateForces() is called within */
void integrate(System &system, double dt) {

    // DEBUG
    int debug=1;
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
                system.molecules[j].calc_center_of_mass();

                // ROTATE IN Y
                double* rotatedy = rotatePoint(system, 
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                "y", system.molecules[j].ang_pos[1] * 180.0/M_PI); 
                for (int n=0; n<3; n++) 
                    system.molecules[j].atoms[i].pos[n] = rotatedy[n] + system.molecules[j].com[n];
                system.molecules[j].calc_center_of_mass();

                // ROTATE IN Z
                double* rotatedz = rotatePoint(system, 
                system.molecules[j].atoms[i].pos[0] - system.molecules[j].com[0],
                system.molecules[j].atoms[i].pos[1] - system.molecules[j].com[1],
                system.molecules[j].atoms[i].pos[2] - system.molecules[j].com[2],
                "z", system.molecules[j].ang_pos[2] * 180.0/M_PI); 
                for (int n=0; n<3; n++) 
                    system.molecules[j].atoms[i].pos[n] = rotatedz[n] + system.molecules[j].com[n];
                system.molecules[j].calc_center_of_mass(); 
            } // end loop over atoms i 
            } // end if rotations allowed and >1 atom
        // */
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
    for (int j=0; j<system.molecules.size(); j++) {
    if (system.molecules[j].MF == "M") {
    for (int i=0; i<system.molecules[j].atoms.size(); i++) {
    if (system.constants.md_mode == "molecular") {
        if (system.molecules[j].atoms[i].pos[0]  > system.constants.x_max) {
                for (int k=0; k<system.molecules[j].atoms.size(); k++) {
                    system.molecules[j].atoms[k].pos[0] -= system.constants.x_length;
                }
        }
        else if (system.molecules[j].atoms[i].pos[0] < system.constants.x_min) {
                for (int k=0; k<system.molecules[j].atoms.size(); k++) {
	                system.molecules[j].atoms[k].pos[0] += system.constants.x_length;
                }
        }

        if (system.molecules[j].atoms[i].pos[1] > system.constants.y_max) {
	            for (int k=0; k<system.molecules[j].atoms.size(); k++) {
                    system.molecules[j].atoms[k].pos[1] -= system.constants.y_length;
                }
        }
        else if (system.molecules[j].atoms[i].pos[1] < system.constants.y_min) {
	            for (int k=0; k<system.molecules[j].atoms.size(); k++) {
                    system.molecules[j].atoms[k].pos[1] += system.constants.y_length;
                }
        }

        if (system.molecules[j].atoms[i].pos[2]  > system.constants.z_max) {
	            for (int k=0; k<system.molecules[j].atoms.size(); k++) {
                    system.molecules[j].atoms[k].pos[2] -= system.constants.z_length;
                }
        }
        else if (system.molecules[j].atoms[i].pos[2] < system.constants.z_min) {
	            for (int k=0; k<system.molecules[j].atoms.size(); k++) {
                    system.molecules[j].atoms[k].pos[2] += system.constants.z_length;
                }
        }// end pbc check	
        } // end if molecular
        else if (system.constants.md_mode == "atomic") {
        if (system.molecules[j].atoms[i].pos[0]  > system.constants.x_max) {
                    system.molecules[j].atoms[i].pos[0] -= system.constants.x_length;
        }
        else if (system.molecules[j].atoms[i].pos[0] < system.constants.x_min) {
                    system.molecules[j].atoms[i].pos[0] += system.constants.x_length;
        }

        if (system.molecules[j].atoms[i].pos[1] > system.constants.y_max) {
                    system.molecules[j].atoms[i].pos[1] -= system.constants.y_length;
        }
        else if (system.molecules[j].atoms[i].pos[1] < system.constants.y_min) {
                    system.molecules[j].atoms[i].pos[1] += system.constants.y_length;
        }

        if (system.molecules[j].atoms[i].pos[2]  > system.constants.z_max) {
                    system.molecules[j].atoms[i].pos[2] -= system.constants.z_length;
        }
        else if (system.molecules[j].atoms[i].pos[2] < system.constants.z_min) {
                    system.molecules[j].atoms[i].pos[2] += system.constants.z_length;
        }// end pbc check 
        } // end if atomic
	} // end loop i atoms
    } // end if movable
	} // end loop j molecules
   
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
                system.molecules[j].calc_vel(dt);
    
                // rotational
                if (system.constants.md_rotations == "on") {
                    system.molecules[j].calc_ang_acc();
                    system.molecules[j].calc_ang_vel(dt);
                }
            }
        } // end if movable
    	} // end for j molecules
}// end integrate() function
