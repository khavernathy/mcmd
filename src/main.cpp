/*
 Douglas M. Franz
 University of South Florida
 Space group research
 2017
*/

// needed c/c++ libraries
#include <iostream>
#include <string>
#include <strings.h>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
//#include <float.h>

// c++ code files of this program
#include <constants.cpp>
#include <system.cpp>
#include <mc.cpp> // this will include potential.cpp, which includes lj, coulombic, polar
#include <md.cpp>
#include <system_functions.cpp>
#include <io.cpp>
#include <radial_dist.cpp>

using namespace std;

int main(int argc, char **argv) {

	// start timing
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();	
	double time_elapsed;		
	double sec_per_step;

	srand((unsigned)time(NULL)); // initiate random seed

	// disable output buffering (print everything immediately to output)
	setbuf(stdout, NULL); // makes sure runlog output is fluid on SLURM etc.

    // set up our wonderful system
	System system;
	readInput(system, argv[1]);
	readInAtoms(system, system.constants.atom_file);
	paramOverrideCheck(system);	
	centerCoordinates(system);
	defineBox(system, system.constants.x_length, system.constants.y_length, system.constants.z_length);
    if (system.stats.radial_dist == "on")   
        setupRadialDist(system);

    // CONFIRM ATOMS AND MOLECULES PRINTOUT
/*    
	for (int b=0; b<system.molecules.size(); b++) {
		printf("Molecule PDBID = %i: %s has %i atoms and is %s; The first atom has PDBID = %i\n", system.molecules[b].PDBID, system.molecules[b].name.c_str(), (int)system.molecules[b].atoms.size(), system.molecules[b].MF.c_str(), system.molecules[b].atoms[0].PDBID);
		for (int c=0; c<system.molecules[b].atoms.size(); c++) {
            system.molecules[b].atoms[c].printAll();
            printf("\n");
        }
		
	} 
*/
    if (system.constants.mode == "mc") { // prototype is only used for MC.
	printf("Prototype molecule has PDBID %i ( name %s ) and has %i atoms\n", system.proto.PDBID, system.proto.name.c_str(), (int)system.proto.atoms.size());
    /*for (int i=0; i<system.proto.atoms.size(); i++) {
        system.proto.atoms[i].printAll();
    } */	
    }
    
    // END MOLECULE PRINTOUT

	// clobber files 
	remove( system.constants.output_traj.c_str() ); remove( system.constants.thermo_output.c_str() );
	remove( system.constants.restart_pdb.c_str() );
	remove( system.stats.radial_file.c_str() );

    // Prep thermo output file
    FILE *f = fopen(system.constants.thermo_output.c_str(), "w");
    fprintf(f, "#step #TotalE(K) #LinKE(K)  #RotKE(K)  #PE(K)  #density(g/mL) #temp(K) #pres(atm)\n");
    fclose(f);

	// write initial XYZ for first frame.
   //	writeXYZ(system,system.constants.output_traj,0);

    system.checkpoint("Initial protocols complete. Starting MC or MD.");

	// =========================== MONTE CARLO ===================================
	if (system.constants.mode == "mc") {
	printf("| ================================== |\n");
	printf("|  BEGINNING MONTE CARLO SIMULATION  |\n");
	printf("| ================================== |\n\n");

    //outputCorrtime(system, 0); // do initial output before starting mc

    int frame = 1;
    int stepsize = system.constants.stepsize;
    int finalstep = system.constants.finalstep;
    int corrtime = system.constants.mc_corrtime; // print output every x steps, change volume etc.

    double current_lj_sum=0.0; double current_es_sum=0.0; double lj_average=0.0; double es_average=0.0; double current_polar_sum=0.0; double polar_average=0.0;
    double current_energy_sum=0.0; double current_density_sum=0.0; double current_volume_sum=0.0; double current_z_sum = 0.0;
    double energy_average=0.0; double density_average=0.0; double volume_average=0.0; double z_average = 0.0;
    double current_Nmov_sum=0.0; double Nmov_average=0.0;
    double current_wt_percent_sum=0.0; double wt_percent_average=0.0;
    double current_wt_percent_ME_sum=0.0; double wt_percent_ME_average=0.0;

    if (system.constants.potential_form == "ljespolar") {
        system.constants.old_total_atoms = system.constants.total_atoms;
        int N = 3 * system.constants.total_atoms;
        system.constants.A_matrix= (double **) calloc(N,sizeof(double*));
        for (int i=0; i< N; i++ ) {
            system.constants.A_matrix[i]= (double *) malloc(N*sizeof(double));
        }
    }


    // begin timing for steps "begin_steps"
	std::chrono::steady_clock::time_point begin_steps = std::chrono::steady_clock::now();	
	
	// main MC slice
	int corrtime_iter=1;
	for (int t=0; t <= finalstep; t+=stepsize) {
		system.checkpoint("New step starting."); //printf("Step %i\n",t);	

				// DO MC STEP BRAH
                if (t>0) {
				    runMonteCarloStep(system,system.constants.potential_form);
                    system.constants.old_total_atoms = system.constants.total_atoms;
                    system.checkpoint("...finished runMonteCarloStep");
                }
                         

        // CHECK FOR CORRTIME
        if (t==0 || t % corrtime == 0) { /// write every x steps
         //   outputCorrtime(system, t);				
  
			/* -------------------------------- */	
			// [[[[ CALCULATE OUTPUT VALUES ]]]]		
			/* -------------------------------- */

            // TIMING 
            std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
            time_elapsed = (std::chrono::duration_cast<std::chrono::microseconds>(end - begin_steps).count()) /1000000.0;
            sec_per_step = time_elapsed/t;
            double progress = (((float)t)/finalstep*100);
            double ETA = ((time_elapsed*finalstep/t) - time_elapsed)/60.0;
            double ETA_hrs = ETA/60.0;



			// BOLTZMANN AVERAGES
			double bf_avg = (system.stats.insert_bf_sum + system.stats.remove_bf_sum + system.stats.displace_bf_sum + system.stats.volume_change_bf_sum + system.stats.rotate_bf_sum)/5.0/t;
			double ibf_avg = system.stats.insert_bf_sum/t;
			double rbf_avg = system.stats.remove_bf_sum/t;
			double dbf_avg = system.stats.displace_bf_sum/t;
			double vbf_avg = system.stats.volume_change_bf_sum/t;
			double robf_avg = system.stats.rotate_bf_sum/t; //printf("%f -----------\n",system.stats.rotate_bf_sum);

			// MASS
			double totalmass = 0.0; double movablemass = 0.0; double frozenmass = 0.0;
			for (int c=0; c<system.molecules.size();c++) {
				for (int d=0; d<system.molecules[c].atoms.size(); d++) {
					totalmass += system.molecules[c].atoms[d].m/system.constants.cM/system.constants.NA; // total mass in g
				    if (system.molecules[c].MF == "M")
                        movablemass += system.molecules[c].atoms[d].m/system.constants.cM/system.constants.NA;
                    else if (system.molecules[c].MF == "F")
                        frozenmass += system.molecules[c].atoms[d].m/system.constants.cM/system.constants.NA;
                }
			}

            // N_movables (sorbates, usually)
            current_Nmov_sum += system.stats.count_movables;
            if (t==0) system.constants.initial_sorbates = system.stats.count_movables;
            Nmov_average = current_Nmov_sum/corrtime_iter;
    

			// ENERGY
			double* potentials = getTotalPotential(system,system.constants.potential_form); //new double[2];
                                current_energy_sum+=(potentials[0]+potentials[1]+potentials[2]);
                                if (t==0) system.constants.initial_energy = current_energy_sum;
                                current_lj_sum += potentials[0];
                                current_es_sum += potentials[1];
		                		current_polar_sum += potentials[2];
                                        lj_average = current_lj_sum/corrtime_iter;
                                        es_average = current_es_sum/corrtime_iter;
				                    	polar_average = current_polar_sum/corrtime_iter;
                                        energy_average = current_energy_sum/corrtime_iter;		
            
            // FOR CHEMICAL POTENTIAL dE/dN
            double chemical_potential = (energy_average - system.constants.initial_energy)
                    / (Nmov_average - system.constants.initial_sorbates);		
            
			// VOLUME
			double volume = system.constants.volume; 
			current_volume_sum += volume;
				volume_average = current_volume_sum/corrtime_iter;
	
			// DENSITY
			double density = movablemass/(volume*1e-24); // that's mass in g /mL	
			current_density_sum += density;
				density_average = current_density_sum/corrtime_iter;

            // WT %
            double wt_percent = (movablemass / totalmass)*100;
            double wt_percent_ME = (movablemass / frozenmass)*100;
            current_wt_percent_sum += wt_percent;
            current_wt_percent_ME_sum += wt_percent_ME;
                wt_percent_average = current_wt_percent_sum/corrtime_iter;
                wt_percent_ME_average = current_wt_percent_ME_sum/corrtime_iter;			

			// COMPRESSIBILITY FACTOR Z = PV/nRT  =  atm*L / (mol * J/molK * K)
			// GOOD FOR HOMOGENOUS GASES ONLY!!
            double n_moles_sorb = movablemass/(system.proto.get_mass()*1000*system.constants.NA);
			double Z = (system.constants.pres*(volume*1e-27) * 101.325 ) // PV
                    / ( (n_moles_sorb) * system.constants.R  * system.constants.temp ); // over nRT
			current_z_sum += Z;
				z_average = current_z_sum/corrtime_iter;

			// MC MOVE ACCEPT STATS
			double total_accepts = system.stats.insert_accepts + system.stats.remove_accepts + system.stats.displace_accepts + system.stats.volume_change_accepts + system.stats.rotate_accepts;
			double total_attempts = system.stats.insert_attempts + system.stats.remove_attempts + system.stats.displace_attempts + system.stats.volume_attempts + system.stats.rotate_attempts;

			// PRINT MAIN OUTPUT
			printf("MONTE CARLO\n");
            printf("%s %s\n",system.constants.jobname.c_str(),argv[1]);
			printf("ENSEMBLE: %s; T = %.3f; P = %.3f\n",system.constants.ensemble.c_str(), system.constants.temp, system.constants.pres);
			printf("Input atoms: %s\n",system.constants.atom_file.c_str());
			printf("Step: %i / %i; Progress = %.3f%%\n",t,finalstep,progress);
			printf("Time elapsed = %.2f s = %.3f sec/step; ETA = %.3f min = %.3f hrs\n",time_elapsed,sec_per_step,ETA,ETA_hrs);
			printf("BF avg = %.4f       ( %.3f Ins / %.3f Rem / %.3f Dis / %.3f Rot / %.3f Vol ) \n",
				bf_avg,
				ibf_avg,
				rbf_avg,
				dbf_avg,
				robf_avg,
				vbf_avg);
				
			printf("AR avg = %.4f       ( %.3f Ins / %.3f Rem / %.3f Dis / %.3f Rot / %.3f Vol )  \n",
				total_accepts/total_attempts, 
				(double)system.stats.insert_accepts/system.stats.insert_attempts, 
				(double)system.stats.remove_accepts/system.stats.remove_attempts, 
				(double)system.stats.displace_accepts/system.stats.displace_attempts, 
				(double)system.stats.rotate_accepts/system.stats.rotate_attempts, 
				(double)system.stats.volume_change_accepts/system.stats.volume_attempts); 
			printf("Total accepts: %i ( %.2f%% Ins / %.2f%% Rem / %.2f%% Dis / %.2f%% Rot / %.2f%% Vol )  \n", 
				(int)total_accepts, 
				system.stats.insert_accepts/total_accepts*100, 
				system.stats.remove_accepts/total_accepts*100, 
				system.stats.displace_accepts/total_accepts*100, 
				system.stats.rotate_accepts/total_accepts*100, 
				system.stats.volume_change_accepts/total_accepts*100);
			printf("RD avg =              %.4f K  (LJ = %.4f, LRC = %.4f, LRC_self = %.4f)\n",lj_average,system.constants.lj,system.constants.lj_lrc, system.constants.lj_self_lrc);
			printf("ES avg =              %.4f K  (real = %.4f, recip = %.4f, self = %.4f)\n",es_average,system.constants.coulombic_real, system.constants.coulombic_reciprocal, system.constants.coulombic_self);
			printf("Polar avg =           %.4f K\n",polar_average);
			printf("Total potential avg = %.4f K\n",energy_average);
			printf("Volume avg  = %.2f A^3 = %.2f nm^3\n",volume_average,volume_average/1000.0);
			printf("Density avg = %.6f g/mL = %6f g/L \n",density_average,density_average*1000.0); 
			printf("wt %% = %.4f %%; wt %% ME = %.4f %% \n",wt_percent_average, wt_percent_ME_average);
            printf("Chemical potential avg = %.4f K per sorbate molecule \n", chemical_potential); 
            printf("Compressibility factor Z avg = %.6f (for homogeneous gas %s) \n",z_average,system.proto.name.c_str());
			printf("N_movables avg = %.3f; N_molecules = %i; N_movables = %i; N_sites = %i\n",Nmov_average,(int)system.molecules.size(), system.stats.count_movables, system.constants.total_atoms);
			printf("--------------------\n\n");
			
			// CONSOLIDATE ATOM AND MOLECULE PDBID's
            // quick loop through all atoms to make PDBID's pretty (1->N)
            int molec_counter=1, atom_counter=1;
            for (int i=0; i<system.molecules.size(); i++) {
                system.molecules[i].PDBID = molec_counter;
                for (int j=0; j<system.molecules[i].atoms.size(); j++) {
                    system.molecules[i].atoms[j].PDBID = atom_counter;
                    system.molecules[i].atoms[j].mol_PDBID = system.molecules[i].PDBID;
                    atom_counter++;
                } // end loop j
                molec_counter++;
            } // end loop i

            // WRITE RESTART FILE AND OTHER OUTPUTS
			writeXYZ(system,system.constants.output_traj,frame,t,0);
            frame++;
            writePDB(system, system.constants.restart_pdb);
			writeThermo(system, energy_average, 0.0, 0.0, energy_average, density_average*1000, system.constants.temp, system.constants.pres, t);
            if (system.stats.radial_dist == "on") {
                radialDist(system);
                writeRadialDist(system);		
            }	

            // count the corrtime occurences.
            corrtime_iter++;

		} // END IF CORRTIME	
	} // END MC STEPS LOOP.

	// FINAL EXIT OUTPUT
	printf("Final box length parameters: \n");
	printf("x: %.5f\n",system.constants.x_length);
	printf("y: %.5f\n",system.constants.y_length);
	printf("z: %.5f\n",system.constants.z_length);
	printf("Insert accepts:        %i\n", system.stats.insert_accepts);
	printf("Remove accepts:        %i\n", system.stats.remove_accepts);
	printf("Displace accepts:      %i\n", system.stats.displace_accepts);
	printf("Volume change accepts: %i\n", system.stats.volume_change_accepts);	

	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
        time_elapsed = (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0;
	printf("Total wall time = %f s\n",time_elapsed);

    printf("Freeing data structures... ");
    for (int i=0; i< 3* system.constants.total_atoms; i++) {
        free(system.constants.A_matrix[i]);
    }
    free(system.constants.A_matrix); system.constants.A_matrix = NULL;
    printf("done.\n");
	printf("MC steps completed. Exiting program.\n"); std::exit(0);
	}
	// ===================== END MONTE CARLO ================================================
	

	
	// ===================== MOLECULAR DYNAMICS ==============================================
	else if (system.constants.mode == "md") {
    
        // write initial XYZ
        writeXYZ(system,system.constants.output_traj, 1, 0, 0);	
        int frame = 2; // weird way to initialize but it works for the output file.

	    // assign initial velocities
    double randv; double DEFAULT = 99999.99;
    if (system.constants.md_init_vel != DEFAULT) {
        // if user defined
        for (int i=0; i<system.molecules.size(); i++) {
            for (int n=0; n<3; n++) {
                randv = ((double)rand() / (double)RAND_MAX *2 - 1) * system.constants.md_init_vel;
                system.molecules[i].vel[n] = randv; // that is, +- 0->1 * user param
            }
        }
        printf("Computed initial velocities via user-defined value: %f A/fs\n",system.constants.md_init_vel);
    } else {
        // otherwise use temperature as default via v_RMS
        // default temp is zero so init. vel's will be 0 if no temp is given.
        double v_init = sqrt(8.0 * system.constants.R * system.constants.temp / (M_PI*system.proto.mass*system.constants.NA)) / 1e5; // A/fs
        double v_component = sqrt(v_init*v_init / 3.0);

        double pm = 0;
        for (int i=0; i<system.molecules.size(); i++) {
            for (int n=0; n<3; n++) {
                randv = (double)rand() / (double)RAND_MAX;
                if (randv > 0.5) pm = 1.0;
                else pm = -1.0;

                system.molecules[i].vel[n] = v_component * pm;
            }
        }
        printf("Computed initial velocities from temperature: v_init = %f A/fs\n",v_init);
        system.constants.md_init_vel = v_init;
    } 
     // end initial velocities

    // compute inital COM for all molecules, and moment of inertia
    // (io.cpp handles molecular masses // 
    for (int i=0; i<system.molecules.size(); i++) {
        system.molecules[i].calc_center_of_mass();
        if (system.molecules[i].atoms.size() > 1) system.molecules[i].calc_inertia();
    }

	double dt = system.constants.md_dt; // * 1e-15; //0.1e-15; // 1e-15 is one femptosecond.
	double tf = system.constants.md_ft; // * 1e-15; //100e-15; // 100,000e-15 would be 1e-9 seconds, or 1 nanosecond. 
	int total_steps = floor(tf/dt);
	int count_md_steps = 1;
	
        printf("| ========================================= |\n");
        printf("|  BEGINNING MOLECULAR DYNAMICS SIMULATION  |\n");
        printf("| ========================================= |\n\n");
	

	// begin timing for steps
	std::chrono::steady_clock::time_point begin_steps = std::chrono::steady_clock::now();
	for (double t=dt; t < tf; t=t+dt) {
		//printf("%f\n",t);
		integrate(system,dt);
		
		if (count_md_steps % system.constants.md_corrtime == 0) {  // print every x steps 
				
            // get KE and PE and T at this step.
            double* ETarray = calculateEnergyAndTemp(system, t);
            double KE = ETarray[0];
            double PE = ETarray[1];
            double TE = KE+PE;
            double Temp = ETarray[2];
            double v_avg = ETarray[3];
            double Ek = ETarray[4];
            double Klin = ETarray[5];
            double Krot = ETarray[6];

            // PRESSURE
            double pressure = TE/system.constants.volume * system.constants.kb * 1e30 * 9.86923e-6; // P/V to atm

			// PRINT OUTPUT
			std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
			time_elapsed = (std::chrono::duration_cast<std::chrono::microseconds>(end - begin_steps).count()) /1000000.0;
			sec_per_step = time_elapsed/(double)count_md_steps;				
			double progress = (((float)count_md_steps)/(float)total_steps*100);
			double ETA = ((time_elapsed*(float)total_steps/(float)count_md_steps) - time_elapsed)/60.0;
			double ETA_hrs = ETA/60.0;
       
            printf("MOLECULAR DYNAMICS\n");
            printf("%s %s\n",system.constants.jobname.c_str(),argv[1]);
            printf("Input atoms: %s\n",system.constants.atom_file.c_str());
            printf("ENSEMBLE: %s\n",system.constants.ensemble.c_str());
            printf("Input T: %.3f K; Input P: %.3f atm\n",system.constants.temp, system.constants.pres);
            printf("Step: %i / %i; Progress = %.3f%%; Realtime = %.2f fs\n",count_md_steps,total_steps,progress,t);
			printf("Time elapsed = %.2f s = %.3f sec/step; ETA = %.3f min = %.3f hrs\n",time_elapsed,sec_per_step,ETA,ETA_hrs);
			printf("KE: %.3f K (lin: %.3f , rot: %.3e ); PE: %.3f K; Total E: %.3f K; \nKE(equipartition) = %.3f K\nEmergent T: %.3f K; Average v = %.5f A/fs; v_init = %.5f A/fs\nEmergent Pressure: %.3f atm\n", 
            KE, Klin, Krot, PE, TE, Ek, 
            Temp, v_avg, system.constants.md_init_vel,
            pressure);			

			printf("--------------------\n\n");

            // WRITE OUTPUT FILES 
			writeXYZ(system,system.constants.output_traj,frame,count_md_steps,t);
            frame++;
            writeThermo(system, TE, Klin, Krot, PE, 0.0, Temp, pressure, count_md_steps); 
            writePDB(system,system.constants.restart_pdb);	
            if (system.stats.radial_dist == "on") {
                radialDist(system);
                writeRadialDist(system);
            }
		}
		count_md_steps++;
	} // end MD timestep loop
	}
}        
