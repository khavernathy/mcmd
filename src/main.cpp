/*
 Douglas M. Franz
 University of South Florida
 Space group research
 2017
*/

// needed c/c++ libraries
#include <cmath>
#include <iostream>
#include <ctime>
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
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>
#include <sys/stat.h>
//#include <float.h>

// Windows compilation (c11 instead of c++) doesn't auto-define pi.
#ifdef WINDOWS
#define M_PI 3.14159265359
#endif

// c++ code files of this software
// ORDER MATTERS HERE !!!!!!!!!!!!!
#ifdef MPI
    #include <mpi.h>
#endif
#ifdef OMP
    #include "omp.h"
#endif
#include "usefulmath.cpp"
#include "classes.cpp"
#include "system.cpp"
#include "fugacity.cpp"
#include "distance.cpp"
#include "system_functions.cpp" 
#include "mc.cpp" // this will include potential.cpp, which includes lj, coulombic, polar
#ifdef CUDA
    #include "cudafuncs.cu"  // CUDA STUFF
#endif
#include "md.cpp"
#include "sp.cpp"
#include "optimize.cpp"
#include "io.cpp"
#include "radial_dist.cpp"
#include "averages.cpp"
#include "histogram.cpp"

using namespace std;

int main(int argc, char **argv) {
   
        /*        
        // MPI is on hold for now.
        // set the default MPI params
        int rank=0, size=1;
    
       // MPI setup
        #ifdef MPI
         MPI_Init(&argc, &argv);
         MPI_Comm_rank(MPI_COMM_WORLD, &rank);
         MPI_Comm_size(MPI_COMM_WORLD, &size);
         char processor_name[MPI_MAX_PROCESSOR_NAME];
         int name_len;
         MPI_Get_processor_name(processor_name, &name_len);
         // Print off a hello world message
         printf("Hello world from processor %s, rank %d"
                " out of %d processors\n",
                processor_name, rank, size);
        #endif
         // suffix for filenames for each core running (if MPI)
         string rankstring = to_string(rank);
        // suffix is, e.g. -00004 for 5th core; -01000 for 1001st core
        std::string corenum = std::string(5 - rankstring.length(), '0') + rankstring;
        char suffix[7] = "-";
        std::strcat(suffix,corenum.c_str());
        //printf("suffix: %s\n", suffix); 
         */

    //output current date/time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer,sizeof(buffer),"%m-%d-%Y at %H:%M:%S",timeinfo);
    std::string str(buffer);
    std::cout << "MCMD started on " << str << endl;

    // print system info.
    string hostcom="hostname";
    int zzz=std::system(hostcom.c_str());
    string linuxcheck="/proc/cpuinfo";
    //linux
    if (std::ifstream(linuxcheck.c_str())) {
        string cpucom="cat /proc/cpuinfo  | tail -25 | grep -i 'model name'";
        zzz=std::system(cpucom.c_str());
        zzz=std::system("echo $(mem=$(grep MemTotal /proc/meminfo | awk '{print $2}'); echo $mem | awk {'print $1/1024/1024'})' GB memory available on this node (Linux).'");
    // mac
    } else {
        string cpumac="sysctl -n machdep.cpu.brand_string";
        zzz=std::system(cpumac.c_str());
        zzz=std::system("echo $(mem=$(sysctl hw.memsize | awk {'print $2'}); echo $mem | awk {'print $1/1024/1024/1024.0'})' GB memory available on this node (Mac).'");
    }

   	// start timing for checkpoints
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	double time_elapsed;
	double sec_per_step;

	srand((unsigned)time(NULL)); // initiate random seed, old version
    srand48((unsigned)time(NULL)); // initiate drand48 48-bit integer random seed.

	// disable output buffering (print everything immediately to output)
	setbuf(stdout, NULL); // makes sure runlog output is fluid on SLURM etc.

    // SET UP THE SYSTEM
    System system;
	system.checkpoint("setting up system with main functions...");
    readInput(system, argv[1]); // executable takes the input file as only argument.
	if (system.constants.restart_mode) {
        // restarting an old job. Make a new saved_data folder e.g. "saved_data4" then overwrite the input file.
        int save_folder_number=1;
        struct stat sb;
        string folder_name = "saved_data"+to_string(save_folder_number);
        while (stat(folder_name.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
            save_folder_number++;
            folder_name = "saved_data"+to_string(save_folder_number);
        }
        string command = "mkdir "+folder_name+"; cp * "+folder_name+"; cp restart.pdb "+system.constants.atom_file.c_str();
        int whatever=std::system(command.c_str());
        printf("RESTARTING previous job from input atoms contained in %s/restart.pdb\n",folder_name.c_str());
    }
    if (system.constants.mode == "md") system.constants.auto_reject_option = 0; // FORCE auto-reject off for MD b/c it will give wrong potential energies
    if (system.constants.atom_file != "<none>") {
        readInAtoms(system, system.constants.atom_file);
        consolidatePDBIDs(system);
        // collect movable molecule IDs
        for (int i=0; i<system.molecules.size(); i++)
            if (!system.molecules[i].frozen)
                system.stats.movids.push_back(i);
    }
	paramOverrideCheck(system);
	if (system.constants.autocenter)
        centerCoordinates(system);
    setupBox(system);
    if (system.constants.manual_cutoff) system.pbc.cutoff = system.constants.manual_cutoff_val; // override the cutoff if user-defined.
    if (system.stats.radial_dist) {
        string command = "rm " + system.stats.radial_file + "*";
        int whatever=std::system(command.c_str()); //remove( system.stats.radial_file.c_str() ); 
        setupRadialDist(system);
    }
    if (system.constants.scale_charges)
        scaleCharges(system);

    moleculePrintout(system); // this will confirm the sorbate to the user in the output. Also checks for system.constants.model_name and overrides the prototype sorbate accordingly.
    if (system.constants.write_lammps) 
        writeLAMMPSfiles(system);

    
    if (system.constants.crystalbuild) {
        setupCrystalBuild(system);
        // write an XYZ of the built system by default
        string delit = "rm crystalbuild.xyz";
        int whatever = std::system(delit.c_str());
        writeXYZ(system, "crystalbuild.xyz", 0, 0, 0, 0);
    }

    if (system.constants.histogram_option) {
        if (system.pbc.volume > 100.*100.*100.) {
            std::cout << "ERROR: Histogram cannot be enabled for a box with volume > 10^6 cubic angstroms. Current volume is " << to_string(system.pbc.volume) << " A^3. Use `histogram off`, or use a different box size, e.g. `carbasis 99 99 99 90 90 90`." << endl;
            exit(EXIT_FAILURE);
        }        

        system.grids.histogram = (histogram_t *) calloc(1,sizeof(histogram_t));
        system.grids.avg_histogram = (histogram_t *) calloc(1,sizeof(histogram_t));
        setup_histogram(system);
        allocate_histogram_grid(system);
    }
    setupFugacity(system);
    if (system.constants.bias_uptake != 0 && system.constants.ensemble == ENSEMBLE_UVT)
        setupNBias(system); 
    if (system.constants.fragmaker) {
        string del = "rm fragment-*.xyz";
        int whatev = std::system(del.c_str());
        fragmentMaker(system);
    }

    system.pbc.printBasis();
    initialize(system); // these are just system name sets, nothing more
    printf("SORBATE COUNT: %i\n", (int)system.proto.size());
    printf("VERSION NUMBER: %i\n", 796); // i.e. github commit
    inputValidation(system);
    printf("...input options validated.\n");
    system.checkpoint("...input options validated. Done with system setup functions.");

    // compute inital COM for all molecules, and moment of inertia
    // (io.cpp handles molecular masses //
    for (int i=0; i<system.molecules.size(); i++) {
        system.molecules[i].calc_center_of_mass();
        if (system.molecules[i].atoms.size() > 1) system.molecules[i].calc_inertia();
        for (int n=0;n<3;n++) system.molecules[i].original_com[n] = system.molecules[i].com[n]; // save original molecule COMs for diffusion calculation in MD.
    }

	// *** clobber files
	remove( system.constants.output_traj.c_str() ); remove( system.constants.thermo_output.c_str() );
	remove( system.constants.restart_pdb.c_str() ); remove ( system.constants.output_traj_pdb.c_str() );
	remove( system.constants.output_histogram.c_str() );
	remove( system.constants.dipole_output.c_str() ); remove( system.constants.frozen_pdb.c_str() );
    remove( system.constants.restart_mov_pdb.c_str() ); remove( system.constants.output_traj_movers_pdb.c_str() );
        
    // *** done clobbering files.

    // INITIAL WRITEOUTS
    // Prep thermo output file
    FILE *f = fopen(system.constants.thermo_output.c_str(), "w");
    fprintf(f, "#step #TotalE(K) #LinKE(K)  #RotKE(K)  #PE(K) #RD(K) #ES(K) #POL(K) #density(g/L) #temp(K) #pres(atm) #N\n");
    fclose(f);
    // Prep pdb trajectory if needed
    if (system.constants.pdb_traj_option) {
        if (system.constants.pdb_bigtraj_option) {
            FILE *f = fopen(system.constants.output_traj_pdb.c_str(), "w");
            fclose(f);
        } else {
            // also the movables traj (going to phase-out the old trajectory 
            // which writes frozen atoms every time
            FILE *g = fopen(system.constants.restart_mov_pdb.c_str(), "w");
            fclose(g);
        }
    }
		// Prep histogram if needed
		if (system.constants.histogram_option)
			system.file_pointers.fp_histogram = fopen(system.constants.output_histogram.c_str(), "w");
    // frozen .pdb (just the MOF, usually)
    if (system.stats.count_frozens > 0) {
        writePDBfrozens(system, system.constants.frozen_pdb.c_str());
    }
    // END INTIAL WRITEOUTS

    system.checkpoint("Initial protocols complete. Starting MC or MD.");

    // RESIZE A MATRIX IF POLAR IS ACTIVE (and initialize the dipole file)
    if (system.constants.potential_form == POTENTIAL_LJESPOLAR || system.constants.potential_form == POTENTIAL_LJPOLAR || system.constants.potential_form == POTENTIAL_COMMYESPOLAR) {
				FILE * fp = fopen(system.constants.dipole_output.c_str(), "w");
				fclose(fp);

        double memreqA;
        system.last.total_atoms = system.constants.total_atoms;
        int N = 3 * system.constants.total_atoms;
        system.last.thole_total_atoms = system.constants.total_atoms;
        makeAtomMap(system); // writes the atom indices
        // 1/2 matrix
        if (!system.constants.full_A_matrix_option) {
            system.constants.A_matrix = (double **) calloc(N, sizeof(double*));
            int inc=0, blocksize=3;
            for (int i=0; i<N; i++) {
                system.constants.A_matrix[i] = (double *) malloc(blocksize*sizeof(double));
                inc++;
                if (inc%3==0) blocksize+=3;
            }
            memreqA = (double)sizeof(double)*((3*system.constants.total_atoms*3*system.constants.total_atoms - 3*system.constants.total_atoms)/2.0)/(double)1e6;
        // full matrix
        } else {
            system.constants.A_matrix_full = (double **) calloc(N, sizeof(double*));
            for (int i=0;i<N;i++) {
                system.constants.A_matrix_full[i] = (double *) malloc(N * sizeof(double));
            }
            memreqA = (double)sizeof(double)* ( N*N )/(double)1e6;
        }
        printf("The polarization Thole A-Matrix will require %.2f MB = %.4f GB.\n", memreqA, memreqA/1000.);

}

    // SET UP Ewald k-space if needed
    if (system.constants.mode == "md" && (system.constants.potential_form == POTENTIAL_LJESPOLAR || system.constants.potential_form == POTENTIAL_LJES || system.constants.potential_form == POTENTIAL_COMMYES || system.constants.potential_form == POTENTIAL_COMMYESPOLAR)) {
        // based on k-max, find the number of k-vectors to use in Ewald force.
        int count_ks = 0;
        double kmax = system.constants.ewald_kmax;
        int l[3];
        // define k-space
        for (l[0] = 0; l[0] <= kmax; l[0]++) {
            for (l[1] = (!l[0] ? 0 : -kmax); l[1] <= kmax; l[1]++) {
                for (l[2] = ((!l[0] && !l[1]) ? 1 : -kmax); l[2] <= kmax; l[2]++) {
                    // skip if norm is out of sphere
                    if (l[0]*l[0] + l[1]*l[1] + l[2]*l[2] > kmax*kmax) continue;
                    count_ks++;
                } // end for l[2], n
            } // end for l[1], m
        } // end for l[0], l
        
        system.constants.ewald_num_k = count_ks;
    } // end MD Ewald k-space setup.

    // BEGIN MC OR MD ===========================================================
	// =========================== MONTE CARLO ===================================
	if (system.constants.mode == "mc") {
	printf("\n| ================================== |\n");
	printf("|  BEGINNING MONTE CARLO SIMULATION  |\n");
	printf("| ================================== |\n\n");

    //outputCorrtime(system, 0); // do initial output before starting mc
    int frame = 1;
    int stepsize = system.constants.stepsize;
    int finalstep = system.constants.finalstep;
    int corrtime = system.constants.mc_corrtime; // print output every corrtime steps

    // begin timing for steps "begin_steps"
	std::chrono::steady_clock::time_point begin_steps = std::chrono::steady_clock::now();

	// MAIN MC STEP LOOP
	int corrtime_iter=1;
	for (int t=0; t <= (finalstep-system.constants.step_offset); t+=stepsize) { // 0 is initial step
		system.checkpoint("New MC step starting."); //printf("Step %i\n",t);
        system.stats.MCstep = t;
        system.stats.MCcorrtime_iter = corrtime_iter;

				// DO MC STEP
                if (t!=0) {
                    setCheckpoint(system); // save all the relevant values in case we need to revert something.
                    //make_pairs(system); // establish pair quantities
                    //computeDistances(system);
                    runMonteCarloStep(system);
                    system.checkpoint("...finished runMonteCarloStep");

                    if (system.stats.MCmoveAccepted == false)
                        revertToCheckpoint(system);
                    else if (system.constants.simulated_annealing) { // S.A. only goes when move is accepted.
                        system.constants.temp =
                            system.constants.sa_target +
                            (system.constants.temp - system.constants.sa_target) *
                            system.constants.sa_schedule;
                    }

                    //computeAverages(system);
                } else {
                    computeInitialValues(system);
                }

        // CHECK FOR CORRTIME
        if (t==0 || t % corrtime == 0 || t == finalstep) { /// output every x steps

			// get all observable averages
            if (t>0 || (t==0 && system.stats.count_movables>0)) computeAverages(system);

			// prep histogram for writing.
			if (system.constants.histogram_option) {
			    zero_grid(system.grids.histogram->grid,system);
                population_histogram(system);
                if (t != 0) update_root_histogram(system);
            }
            /* -------------------------------- */
			// [[[[ PRINT OUTPUT VALUES ]]]]
			/* -------------------------------- */

            // TIMING
            std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
            time_elapsed = (std::chrono::duration_cast<std::chrono::microseconds>(end - begin_steps).count()) /1000000.0;
            sec_per_step = time_elapsed/t;
            double progress = (((float)t)/finalstep*100);
            double ETA = ((time_elapsed*finalstep/t) - time_elapsed)/60.0;
            double ETA_hrs = ETA/60.0;
            double efficiency = system.stats.MCeffRsq / time_elapsed;

			// PRINT MAIN OUTPUT
            if (system.constants.ensemble == ENSEMBLE_UVT && system.constants.bias_uptake_switcher)
                printf("MCMD Monte Carlo: %s (%s): Loading bias (%s %i %s)\n", system.constants.jobname.c_str(), argv[1], "on:", (int)ceil(system.constants.bias_uptake), "molecules");
            else if (system.constants.ensemble == ENSEMBLE_UVT && !system.constants.bias_uptake_switcher)
                printf("MCMD Monte Carlo: %s (%s): Loading bias (%s)\n", system.constants.jobname.c_str(), argv[1], "off");
            else
			    printf("MCMD Monte Carlo: %s (%s)\n", system.constants.jobname.c_str(), argv[1]);
            printf("Input atoms: %s\n",system.constants.atom_file.c_str());
            if (!system.constants.simulated_annealing)
			    printf("Ensemble: %s; T = %.3f K\n", system.constants.ensemble_str.c_str(), system.constants.temp);
			else
                printf("Ensemble: %s; T = %.3f K (Simulated annealing on)\n",system.constants.ensemble_str.c_str(), system.constants.temp);
            
            printf("Time elapsed = %.2f s = %.4f sec/step; ETA = %.3f min = %.3f hrs\n",time_elapsed,sec_per_step,ETA,ETA_hrs);
			printf("Step: %i / %i; Progress = %.3f%%; Efficiency = %.3f\n",system.stats.MCstep+system.constants.step_offset,finalstep,progress,efficiency);
			printf("Total accepts: %i ( %.2f%% Ins / %.2f%% Rem / %.2f%% Dis / %.2f%% Vol )  \n",
				(int)system.stats.total_accepts,
			    system.stats.ins_perc,
                system.stats.rem_perc,
                system.stats.dis_perc,
                system.stats.vol_perc);
	
            printf("BF avg = %.4f       ( %.3f Ins / %.3f Rem / %.3f Dis / %.3f Vol ) \n",
				system.stats.bf_avg,
				system.stats.ibf_avg,
				system.stats.rbf_avg,
				system.stats.dbf_avg,
				system.stats.vbf_avg);

			printf("AR avg = %.4f       ( %.3f Ins / %.3f Rem / %.3f Dis / %.3f Vol )  \n",
				(system.stats.ar_tot),
				(system.stats.ar_ins),
				(system.stats.ar_rem),
				(system.stats.ar_dis),
				(system.stats.ar_vol));
		    
            printf("RD avg =              %.5f +- %.5f K (%.2f %%)\n", // (LJ = %.4f, LRC = %.4f, LRC_self = %.4f)\n",
                system.stats.rd.average, system.stats.rd.sd, system.stats.rd.average/system.stats.potential.average *100); //, system.stats.lj.average, system.stats.lj_lrc.average, system.stats.lj_self_lrc.average);
			printf("ES avg =              %.5f +- %.5f K (%.2f %%)\n", //(real = %.4f, recip = %.4f, self = %.4f)\n",
                system.stats.es.average, system.stats.es.sd, system.stats.es.average/system.stats.potential.average *100); //, system.stats.es_real.average, system.stats.es_recip.average, system.stats.es_self.average);
			printf("Polar avg =           %.5f +- %.5f K (%.2f %%); iterations = %.3f +- %.3f\n",
                system.stats.polar.average, system.stats.polar.sd, system.stats.polar.average/system.stats.potential.average*100, system.stats.polar_iterations.average, system.stats.polar_iterations.sd);
			printf("Total potential avg = %.5f +- %.5f K\n",system.stats.potential.average, system.stats.potential.sd);
			printf("Volume avg  = %.2f +- %.2f A^3 = %.2f nm^3\n",system.stats.volume.average, system.stats.volume.sd, system.stats.volume.average/1000.0);
			for (int i=0; i<system.proto.size(); i++) {
                double mmolg = system.stats.wtpME[i].average * 10 / (system.proto[i].mass*1000*system.constants.NA);
                double cm3gSTP = mmolg*22.4;
                double mgg = mmolg * (system.proto[i].mass*1000*system.constants.NA);
                string flspacing = "";
                if (system.proto[i].name.length() == 3) flspacing="           ";// e.g. CO2
                else flspacing="            "; // stuff like H2, O2 (2 chars)
                if (system.stats.count_frozens > 0) {
                    printf("-> %s wt %% =%s   %.5f +- %.5f %%; %.5f cm^3/g (STP)\n", system.proto[i].name.c_str(),flspacing.c_str(), system.stats.wtp[i].average, system.stats.wtp[i].sd, cm3gSTP);
                    printf("      wt %% ME =            %.5f +- %.5f %%; %.5f mmol/g\n",system.stats.wtpME[i].average, system.stats.wtpME[i].sd, mmolg);
                }
                if (system.stats.count_frozens > 0) {
                    printf("      N_movables =         %.5f +- %.5f;   %.5f mg/g\n",
                    system.stats.Nmov[i].average, system.stats.Nmov[i].sd, mgg);
                } else {
                    printf("-> %s N_movables = %.5f +- %.5f;   %.5f mg/g\n",
                    system.proto[i].name.c_str(),system.stats.Nmov[i].average, system.stats.Nmov[i].sd, mgg);
                }
                if (system.stats.excess[i].average > 0 && system.constants.free_volume > 0)
                    printf("      Excess ads. ratio =  %.5f +- %.5f mg/g\n", system.stats.excess[i].average, system.stats.excess[i].sd);
                printf("      Density avg = %.6f +- %.3f g/mL = %6f g/L = kg/m^3\n",system.stats.density[i].average, system.stats.density[i].sd, system.stats.density[i].average*1000.0);
                if (system.proto.size() > 1)
                    printf("      Selectivity = %.3f +- %.3f\n",system.stats.selectivity[i].average, system.stats.selectivity[i].sd);
            } // end prototype molecules loop for uptake data

            if (system.constants.ensemble == ENSEMBLE_UVT) {
                if (system.proto.size() == 1) {
                    if (system.stats.qst.average > 0)
                        printf("Qst = %.5f kJ/mol\n", system.stats.qst.value); //, system.stats.qst.sd);
                    if (system.stats.qst_nvt.average > 0)
                        printf("U/N avg = %.5f kJ/mol\n", system.stats.qst_nvt.value); //, system.stats.qst_nvt.sd);
                }
                printf("N_molecules = %i N_movables = %i N_sites = %i\n", (int)system.molecules.size(), system.stats.count_movables, system.constants.total_atoms);
            }
            if (system.proto.size() == 1 && system.stats.count_frozen_molecules == 0)
                printf("Compressibility factor Z avg = %.6f +- %.6f (for homogeneous gas %s) \n",system.stats.z.average, system.stats.z.sd, system.proto[0].name.c_str());
            if (system.constants.ensemble != ENSEMBLE_NVE && system.proto.size() ==1)
                printf("Heat capacity = %.5f +- %.5f kJ/molK\n", system.stats.heat_capacity.value, system.stats.heat_capacity.sd);
                
            if (system.constants.dist_within_option) {
                printf("N of %s within %.5f A of origin: %.5f +- %.3f (actual: %i)\n", system.constants.dist_within_target.c_str(), system.constants.dist_within_radius, system.stats.dist_within.average, system.stats.dist_within.sd, (int)system.stats.dist_within.value);
            }
            //if (system.stats.Q.value > 0) printf("Q (parition function) = %.5e\n", system.stats.Q.value);
            printf("--------------------\n\n");

            if (system.molecules.size() > 0) {
            consolidatePDBIDs(system);

            // WRITE RESTART FILE AND OTHER OUTPUTS
            if (system.constants.xyz_traj_option)
                writeXYZ(system,system.constants.output_traj,frame,t,0,system.constants.xyz_traj_movers_option);
            frame++;
            writePDB(system, system.constants.restart_pdb); // all atoms
            if (!system.constants.pdb_bigtraj_option) writePDBmovables(system, system.constants.restart_mov_pdb); // only movers
            if (system.constants.pdb_traj_option) {
                if (system.constants.pdb_bigtraj_option)
                    writePDBtraj(system, system.constants.restart_pdb, system.constants.output_traj_pdb, t); // all atoms
                else writePDBtraj(system, system.constants.restart_mov_pdb, system.constants.output_traj_movers_pdb,t); // just movers
            }
            // ONLY WRITES DENSITY FOR FIRST SORBATE
            writeThermo(system, system.stats.potential.value, 0.0, 0.0, system.stats.potential.value, system.stats.rd.value, system.stats.es.value, system.stats.polar.value, system.stats.density[0].value*1000, system.constants.temp, system.constants.pres, t, system.stats.Nmov[0].value);
            if (system.stats.radial_dist) {
                radialDist(system);
                writeRadialDist(system);
            }
						if (t != 0 && system.constants.histogram_option)
							write_histogram(system.file_pointers.fp_histogram, system.grids.avg_histogram->grid, system);

						if ((system.constants.potential_form == POTENTIAL_LJPOLAR || system.constants.potential_form == POTENTIAL_LJESPOLAR) && system.constants.dipole_output_option)
							write_dipole(system);
            } // end if N > 0
            // count the corrtime occurences.
            corrtime_iter++;

		} // END IF CORRTIME
	} // END MC STEPS LOOP.

	// FINAL EXIT OUTPUT
    if (system.constants.ensemble == ENSEMBLE_NPT) {
	    printf("Final basis parameters: \n");
        system.pbc.printBasis();
    }
	printf("Insert accepts:        %i\n", system.stats.insert_accepts);
	printf("Remove accepts:        %i\n", system.stats.remove_accepts);
	printf("Displace accepts:      %i\n", system.stats.displace_accepts);
	printf("Volume change accepts: %i\n", system.stats.volume_change_accepts);
    printf("Auto-rejects (r <= %.5f A): %i\n", system.constants.auto_reject_r, system.constants.rejects);

	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
        time_elapsed = (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0;
	printf("Total wall time = %f s\n",time_elapsed);

    if (system.constants.potential_form == POTENTIAL_LJESPOLAR || system.constants.potential_form == POTENTIAL_LJPOLAR) {
        printf("Freeing data structures... ");
        // 1/2 matrix
        if (!system.constants.full_A_matrix_option) {
            for (int i=0; i< 3* system.constants.total_atoms; i++) {
                free(system.constants.A_matrix[i]);
            }
            free(system.constants.A_matrix); system.constants.A_matrix = NULL;
        // full matrix
        } else {
            for (int i=0; i<3*system.constants.total_atoms; i++) {
                free(system.constants.A_matrix_full[i]);
            }
            free(system.constants.A_matrix_full); system.constants.A_matrix_full = NULL;
        }
    }
    printf("done.\n");
	printf("MC steps completed. Exiting program.\n"); std::exit(0);
	}
	// ===================== END MONTE CARLO ================================================



	// ===================== MOLECULAR DYNAMICS ==============================================
	else if (system.constants.mode == "md") {

        // write initial XYZ
        if (system.constants.xyz_traj_option)
            writeXYZ(system,system.constants.output_traj, 1, 0, 0, system.constants.xyz_traj_movers_option);
        int frame = 2; // weird way to initialize but it works for the output file.
        // and initial PDB
        writePDB(system,system.constants.restart_pdb);
            if (system.constants.pdb_traj_option && system.constants.pdb_bigtraj_option)
                writePDBtraj(system,system.constants.restart_pdb, system.constants.output_traj_pdb, 0);

        initialVelMD(system);
	      // end initial velocities

	double dt = system.constants.md_dt; // * 1e-15; //0.1e-15; // 1e-15 is one femptosecond.
	double tf = system.constants.md_ft; // * 1e-15; //100e-15; // 100,000e-15 would be 1e-9 seconds, or 1 nanosecond.
	int total_steps = floor(tf/dt);
	int count_md_steps = 1;
    double diffusion_d[3] = {0,0,0}, diffusion_sum=0., D[(int)system.proto.size()];
	double KE=0., PE=0., TE=0., Temp=0., v_avg=0., Klin=0., Krot=0., pressure=0.; //, Ek=0.;
    int i,n;
        printf("\n| ========================================= |\n");
        printf("|  BEGINNING MOLECULAR DYNAMICS SIMULATION  |\n");
        printf("| ========================================= |\n\n");

	// begin timing for steps
	std::chrono::steady_clock::time_point begin_steps = std::chrono::steady_clock::now();
    

    computeInitialValues(system);
    // Main MD time loop
	for (double t=dt; t <= tf; t=t+dt) {
        system.stats.MDtime = t;
		
        // Main Molecular Dynamics Loop function (contains forces, movements, etc.)
        if (system.stats.count_movables > 0)
            integrate(system,dt);

        if (system.constants.ensemble == ENSEMBLE_UVT && count_md_steps % system.constants.md_insert_attempt == 0) {
            // try a MC uVT insert/delete
            getTotalPotential(system); // this is needed on-the-spot because of 
                                       // time-evolution of the system. Otherwise, 
                                       // potential is only re-calculated at corrtime.
            double ranf2 = getrand(); // 0->1
            // ADD A MOLECULE
            if (ranf2 < 0.5 || system.constants.bias_uptake_switcher) { // this will force insertions and never removes if the bias loading is activated.
                system.checkpoint("doing molecule add move.");
                addMolecule(system);
                system.checkpoint("done with molecule add move.");
            } // end add
            else { // REMOVE MOLECULE
                system.checkpoint("doing molecule delete move.");
                removeMolecule(system);
                system.checkpoint("done with molecule delete move.");
            } // end add vs. remove
        }		

        if (count_md_steps % system.constants.md_corrtime == 0 || t==dt || t==tf) {  // print every x steps and first and last.
        
            if (system.constants.ensemble == ENSEMBLE_UVT) computeAveragesMDuVT(system); // get averages (uptake etc.) every corrtime. (for uVT MD only)
            if (system.constants.histogram_option) {
				zero_grid(system.grids.histogram->grid,system);
                population_histogram(system);
                if (t != dt) update_root_histogram(system);
            }

            if (system.stats.count_movables > 0) {
            // get KE and PE and T at this step.
            double* ETarray = calculateObservablesMD(system);
            KE = ETarray[0] * system.constants.K2KJMOL;
            PE = ETarray[1] * system.constants.K2KJMOL;
            TE = KE+PE;
            Temp = ETarray[2];
            v_avg = ETarray[3];
            //Ek = ETarray[4]; // Equipartition Kinetic energy (apparently). Not even using.
            Klin = ETarray[5] * system.constants.K2KJMOL;
            Krot = ETarray[6] * system.constants.K2KJMOL;
            pressure = ETarray[7]; // not using this yet. NVT pressure derived from forces/stat mech stuff. Frenkel p84
            system.stats.temperature.value = Temp;
                system.stats.temperature.calcNewStats();

            // calc diffusion
            // R^2 as a function of time should be linear according to physical theory.
            for (int sorbid=0; sorbid < system.proto.size(); sorbid++) {
                int localN = getNlocal(system, sorbid);
                if (localN < 1) continue; // skip N=0 sorbates

                // re-initialize these vars for each sorbate
                for (int h=0;h<3;h++) diffusion_d[h]=0.;
                diffusion_sum=0.; 
                for (i=0; i<system.molecules.size(); i++) {
                    // only consider molecules of this type (for multi-sorb)
                    if (system.molecules[i].name == system.proto[sorbid].name) {
                        for (n=0; n<3; n++)
                            diffusion_d[n] = system.molecules[i].com[n] + system.molecules[i].diffusion_corr[n] - system.molecules[i].original_com[n];
                    
                        diffusion_sum += dddotprod(diffusion_d, diffusion_d); // the net R^2 from start -> now (mean square displacement)
                    }
                } // end all molecules loop
                system.stats.msd[sorbid].value = diffusion_sum;
                system.stats.msd[sorbid].calcNewStats(); // finds and stores average MSD sum
                double avg_msd_sum = system.stats.msd[sorbid].average;

                D[sorbid] = (avg_msd_sum / (localN *6.0*t)); // 6 because 2*dimensionality = 2*3
                D[sorbid] *= 0.1; // A^2 per fs -> cm^2 per sec (CGS units).
            } // end sorbate types loop
            // we've calc'd diffusion coefficients for all sorbates now.

            // PRESSURE (my pathetic nRT/V method)
						// using this until i get a decent NVT frenkel method working.
						// PE since it's I.G. approximation.
						double nmol = (system.proto[0].mass*system.stats.count_movables)/(system.proto[0].mass*system.constants.NA);
						system.stats.pressure.value = nmol*system.constants.R*Temp/(system.pbc.volume*system.constants.A32L) * system.constants.JL2ATM;
						system.stats.pressure.calcNewStats();
            } // end if N>0 (stats calculation) 

			// PRINT OUTPUT
			std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
			time_elapsed = (std::chrono::duration_cast<std::chrono::microseconds>(end - begin_steps).count()) /1000000.0;
			sec_per_step = time_elapsed/(double)count_md_steps;
			double progress = (((float)count_md_steps)/(float)total_steps*100);
			double ETA = ((time_elapsed*(float)total_steps/(float)count_md_steps) - time_elapsed)/60.0;
			double ETA_hrs = ETA/60.0;
            double outputTime; string timeunit;
            if (t > 1e9) {
                outputTime = t/1e9; timeunit="us";
            } else if (t > 1e6) {
                outputTime = t/1e6; timeunit="ns";
            } else if (t > 1e3) {
                outputTime = t/1e3; timeunit="ps";
            } else {
                outputTime = t; timeunit="fs";
            }
            
            if (system.constants.cuda) printf("MCMD (CUDA) Molecular Dynamics: %s (%s)\n", system.constants.jobname.c_str(), argv[1]);
            else printf("MCMD Molecular Dynamics: %s (%s)\n", system.constants.jobname.c_str(), argv[1]);
            printf("Input atoms: %s\n",system.constants.atom_file.c_str());
            printf("Ensemble: %s; N_movables = %i N_atoms = %i\n",system.constants.ensemble_str.c_str(), system.stats.count_movables, system.constants.total_atoms);
            printf("Time elapsed = %.2f s = %.4f sec/step; ETA = %.3f min = %.3f hrs\n",time_elapsed,sec_per_step,ETA,ETA_hrs);
            printf("Step: %i / %i; Progress = %.3f%%; Realtime = %.5f %s\n",count_md_steps,total_steps,progress,outputTime, timeunit.c_str());
            if (system.constants.ensemble == ENSEMBLE_NVT || system.constants.ensemble == ENSEMBLE_UVT) printf("        Input T = %.4f K\n", system.constants.temp); 
            printf("     Emergent T = %.4f +- %.4f K\n", system.stats.temperature.average, system.stats.temperature.sd);
            printf("Instantaneous T = %.4f K\n", Temp);
            printf("     KE = %.3f kJ/mol (lin: %.3f , rot: %.3f )\n",
                  KE, Klin, Krot );
            printf("     PE = %.3f kJ/mol\n",
                  PE
                  );
            if (system.constants.ensemble == ENSEMBLE_NVE)
                printf("Total E = %.3f :: error = %.3f kJ/mol\n", TE, system.constants.md_NVE_err);
            else
                printf("Total E = %.3f kJ/mol\n", TE);
            printf("Average v = %.2f m/s; v_init = %.2f m/s\nEmergent Pressure = %.3f +- %.3f atm (I.G. approx)\n",
                v_avg*1e5, system.constants.md_init_vel*1e5, system.stats.pressure.average, system.stats.pressure.sd );
            if (system.constants.md_pbc || system.constants.ensemble != ENSEMBLE_UVT) { // for now, don't do diffusion unless PBC is on. (checkInTheBox assumes it)
                for (int sorbid=0; sorbid < system.proto.size(); sorbid++) {
                    printf("Diffusion coefficient of %s = %.4e cm^2 / s\n", system.proto[sorbid].name.c_str(), D[sorbid]);
			        if (system.proto.size() == 1)
                        printf("Mean square displacement = %.5f A^2\n", diffusion_sum/system.stats.count_movables);
                }
            }
            //if (system.stats.Q.value > 0) printf("Q (partition function) = %.5e\n", system.stats.Q.value);
            // uptake data if uVT
            if (system.constants.ensemble == ENSEMBLE_UVT) { 
			for (int i=0; i<system.proto.size(); i++) {
                double mmolg = system.stats.wtpME[i].average * 10 / (system.proto[i].mass*1000*system.constants.NA);
                double cm3gSTP = mmolg*22.4;
                double mgg = mmolg * (system.proto[i].mass*1000*system.constants.NA);
                string flspacing = "";
                if (system.proto[i].name.length() == 3) flspacing="           ";// e.g. CO2
                else flspacing="            "; // stuff like H2, O2 (2 chars)
                if (system.stats.count_frozens > 0) {
                    printf("-> %s wt %% =%s   %.5f +- %.5f %%; %.5f cm^3/g (STP)\n", system.proto[i].name.c_str(),flspacing.c_str(), system.stats.wtp[i].average, system.stats.wtp[i].sd, cm3gSTP);
                    printf("      wt %% ME =            %.5f +- %.5f %%; %.5f mmol/g\n",system.stats.wtpME[i].average, system.stats.wtpME[i].sd, mmolg);
                }
                if (system.stats.count_frozens > 0) {
                    printf("      N_movables =         %.5f +- %.5f;   %.5f mg/g\n",
                    system.stats.Nmov[i].average, system.stats.Nmov[i].sd, mgg);
                } else {
                    printf("-> %s N_movables = %.5f +- %.5f;   %.5f mg/g\n",
                    system.proto[i].name.c_str(),system.stats.Nmov[i].average, system.stats.Nmov[i].sd, mgg);
                }
                if (system.stats.excess[i].average > 0 && system.constants.free_volume > 0)
                    printf("      Excess ads. ratio =  %.5f +- %.5f mg/g\n", system.stats.excess[i].average, system.stats.excess[i].sd);
                printf("      Density avg = %.6f +- %.3f g/mL = %6f g/L = kg/m^3\n",system.stats.density[i].average, system.stats.density[i].sd, system.stats.density[i].average*1000.0);
                if (system.proto.size() > 1)
                    printf("      Selectivity = %.3f +- %.3f\n",system.stats.selectivity[i].average, system.stats.selectivity[i].sd);
            } // end prototype molecules loop for uptake data
            
                if (system.proto.size() == 1) {
                    if (system.stats.qst.average > 0)
                        printf("Qst = %.5f kJ/mol\n", system.stats.qst.value); //, system.stats.qst.sd);
                    if (system.stats.qst_nvt.average > 0)
                        printf("U/N avg = %.5f kJ/mol\n", system.stats.qst_nvt.value); //, system.stats.qst_nvt.sd);
                }
            } // end if uVT
            if ((system.constants.ensemble == ENSEMBLE_NVT || system.constants.ensemble == ENSEMBLE_NVE) && system.proto.size() == 1) {
                if (system.stats.heat_capacity.value > 1e-5)
                    printf("Heat capacity = %.5f +- %.5f kJ/molK\n", system.stats.heat_capacity.value, system.stats.heat_capacity.sd);
                else 
                    printf("Heat capacity = %.5e +- %.5e kJ/molK\n", system.stats.heat_capacity.value, system.stats.heat_capacity.sd);
            }
            if (system.constants.potential_form == POTENTIAL_LJESPOLAR || system.constants.potential_form == POTENTIAL_LJPOLAR)
                printf("Polarization dipole iterations = %.3f +- %.3f\n",
                system.stats.polar_iterations.average, system.stats.polar_iterations.sd);

            printf("--------------------\n\n");
            
            if (system.molecules.size() > 0) {
                consolidatePDBIDs(system);
            } // end if  N>0

            // WRITE OUTPUT FILES
            if (system.molecules.size() > 0) { 
            writeThermo(system, TE, Klin, Krot, PE, system.stats.rd.value, system.stats.es.value, system.stats.polar.value, 0.0, system.stats.temperature.average, pressure, count_md_steps, system.stats.Nmov[0].value);
            // restart file.
            writePDB(system, system.constants.restart_pdb); // containing all atoms

            // trajectory file
                if (system.constants.xyz_traj_option)
			        writeXYZ(system,system.constants.output_traj,frame,count_md_steps,t, system.constants.xyz_traj_movers_option);
                    
                if (!system.constants.pdb_bigtraj_option) writePDBmovables(system, system.constants.restart_mov_pdb); // only movers restart frame
                if (system.constants.pdb_traj_option) {
                    if (system.constants.pdb_bigtraj_option)
                        writePDBtraj(system, system.constants.restart_pdb, system.constants.output_traj_pdb, t); // copy all-atoms-restart-PDB to PDB trajectory
                    else writePDBtraj(system, system.constants.restart_mov_pdb, system.constants.output_traj_movers_pdb,t); // just movers to PDB trajectory
                }
            frame++;
            if (system.stats.radial_dist) {
                radialDist(system);
                writeRadialDist(system);
            }
            if (t != dt && system.constants.histogram_option)
				write_histogram(system.file_pointers.fp_histogram, system.grids.avg_histogram->grid, system);
		    
            } // end if N>0, write output files.
        } // end if corrtime (quite sure.)
		count_md_steps++;
	} // end MD timestep loop
	} // end if MD
// ============================= END MOLECULAR DYNAMICS =======================================


// ============================= SINGLE POINT ENERGY ==========================================
    else if (system.constants.mode == "sp") {
        printf("\n| ==================================== |\n");
	    printf("|  BEGINNING SINGLE POINT CALCULATION  |\n");
	    printf("| ==================================== |\n\n");

        singlePointEnergy(system); 
        

    } // end if Single-Point mode (not md or mc)
// ============================ END SINGLE POINT ENERGY =======================================


// ============================ OPTIMIZATION ==================================================
    else if (system.constants.mode == "opt") {
        printf("\n| ==================================== |\n");
        printf("|   BEGINNING STRUCTURE OPTIMIZATION   |\n");
        printf("| ==================================== |\n\n");
    
        optimize(system);
    } // end optimization mode

} // end main()
