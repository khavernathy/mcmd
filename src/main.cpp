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
#ifdef WINDOWS
       #include <string.h>
	#ifndef S_ISDIR
	#define S_ISDIR(mode)  (((mode) & S_IFMT) == S_IFDIR)
	#endif
#else
       #include <strings.h>
#endif
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
//#ifdef WINDOWS
//#define M_PI 3.14159265359
//#endif

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
#include "bonding.cpp"
#include "mc.cpp" // this will include potential.cpp, which includes lj, coulombic, polar
#ifdef CUDA
    #include "cudafuncs.cu"  // CUDA STUFF
#endif
#include "md_functions.cpp"
#include "md.cpp"
#include "sp.cpp"
#include "io.cpp"
#include "optimize.cpp"
#include "radial_dist.cpp"
#include "observables.cpp"
#include "averages.cpp"
#include "histogram.cpp"
#include "main_out.cpp"

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
#ifndef WINDOWS
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
#else
	printf("??? GB memory available on this node (Windows).\n");
#endif
   	// start timing for checkpoints
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  	#ifndef WINDOWS
		srand48((unsigned)time(NULL)); // initiate drand48 48-bit integer random seed.
	#else
		srand((unsigned)time(NULL)); // initiate random seed, old version
	#endif
	// disable output buffering (print everything immediately to output)
	setbuf(stdout, NULL); // makes sure runlog output is fluid on SLURM etc.

    // SET UP THE SYSTEM
    System system;
	system.checkpoint("setting up system with main functions...");
    readInput(system, argv[1]); // executable takes the input file as only argument.
	system.constants.inputfile = argv[1];
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
    if (system.stats.count_frozens > 0) printf("%s\n",getFormulaUnit(system).c_str());
    initialize(system); // these are just system name sets, nothing more
    printf("SORBATE COUNT: %i\n", (int)system.proto.size());
    printf("VERSION NUMBER: %i\n", 1152); // i.e. github commit
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
    remove( system.constants.molec_dipole_output.c_str() );
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

        FILE * fp2 = fopen(system.constants.molec_dipole_output.c_str(), "w");
        fclose(fp2);

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
            memreqA = (double)sizeof(double)*((N*N - N)/2.0)/(double)1e6;
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

    updateMolecularDOFs(system); // input molecules are defaulted to DOF=3 but may need update
    calcDOF(system); // used only for NVT Nose-Hoover thermostat (other functions calculate it on-the-fly)

    #ifdef OMP
        if (system.constants.openmp_threads > 0)
            printf("Running MCMD with OpenMP using %i threads.\n", system.constants.openmp_threads);
    #endif
    

    // BEGIN MC OR MD ===========================================================
	// =========================== MONTE CARLO ===================================
	if (system.constants.mode == "mc") {
	printf("\n| ================================== |\n");
	printf("|  BEGINNING MONTE CARLO SIMULATION  |\n");
	printf("| ================================== |\n\n");

    //outputCorrtime(system, 0); // do initial output before starting mc
    system.constants.frame = 1;
    int stepsize = system.constants.stepsize;
    long int finalstep = system.constants.finalstep;
    int corrtime = system.constants.mc_corrtime; // print output every corrtime steps

    // begin timing for steps "begin_steps"
	system.constants.begin_steps = std::chrono::steady_clock::now();

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
            mc_main_output(system);
                    
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

        system.constants.frame = 1; 
        initialVelMD(system, 1);
        if (system.constants.ensemble==ENSEMBLE_NVT && system.constants.thermostat_type==THERMOSTAT_NOSEHOOVER && system.constants.user_Q==0)
            calculateNH_Q(system); // Q param for Nose Hoover thermostat

        if (system.constants.flexible_frozen || system.constants.md_mode == MD_FLEXIBLE) {
          findBonds(system);
          setBondingParameters(system);
          printBondParameters(system);
        }

	double dt = system.constants.md_dt; // * 1e-15; //0.1e-15; // 1e-15 is one femptosecond.
	double tf = system.constants.md_ft; // * 1e-15; //100e-15; // 100,000e-15 would be 1e-9 seconds, or 1 nanosecond.
	double thing = floor(tf/dt);
    long int total_steps = (long int)thing;
	int count_md_steps = 0;
    int i,j,n;
        printf("\n| ========================================= |\n");
        printf("|  BEGINNING MOLECULAR DYNAMICS SIMULATION  |\n");
        printf("| ========================================= |\n\n");

	// begin timing for steps
	system.constants.begin_steps = std::chrono::steady_clock::now();

    system.checkpoint("Computing initial values for MD.");
    computeInitialValues(system);
    // Main MD time loop
	for (double t=0; t <= tf; t=t+dt) {
        system.checkpoint("Started timestep.");
        system.stats.MDtime = t;
        system.stats.MDstep = count_md_steps;

        // MD integration. the workload is here. First step is unique. Just get F
        if (t==0) {
            system.checkpoint("t=0; calculating forces.");
            calculateForces(system);
            if (system.constants.ensemble==ENSEMBLE_NVT && system.constants.thermostat_type==THERMOSTAT_NOSEHOOVER) {
                calculateNHLM_now(system); // get Lagrange multiplier for initial state
            }
        } else if (system.stats.count_movables > 0 || system.constants.flexible_frozen) {
            integrate(system);
        }

        // first step: update VACF original velocities according to step 1 if the orig's were all 0
        if (t==dt && system.constants.zero_init_vel_flag) {
            for (i=0;i<system.molecules.size();i++) {
                for (n=0;n<3;n++) system.molecules[i].original_vel[n] = system.molecules[i].vel[n];
                for (j=0;j<system.molecules[i].atoms.size();j++) {
                    for (n=0;n<3;n++) system.molecules[i].atoms[j].original_vel[n] = system.molecules[i].atoms[j].vel[n];
                }
            }
        }

        system.checkpoint("check uVT MD.");
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

        if (count_md_steps % system.constants.md_corrtime == 0 || t==0 || t==tf) {  // print every x steps and first and last.

            if (system.constants.ensemble == ENSEMBLE_UVT) computeAveragesMDuVT(system); // get averages (uptake etc.) every corrtime. (for uVT MD only)
            if (system.constants.histogram_option) {
				zero_grid(system.grids.histogram->grid,system);
                population_histogram(system);
                if (t != dt) update_root_histogram(system);
            }

            if (system.stats.count_movables > 0 || system.constants.flexible_frozen) {

                if (system.constants.simulated_annealing) { // S.A. only goes when move is accepted.
                    system.constants.temp =
                        system.constants.sa_target +
                        (system.constants.temp - system.constants.sa_target) *
                        system.constants.sa_schedule;

                    initialVelMD(system, 0); // reset system temperature by velocities
                }

                /* ========================== */
                calculateObservablesMD(system);
                /* ========================== */

            } // end if N>0 (stats calculation)
		
            // MAIN OUTPUT
            md_main_output(system);
        } // end if corrtime
		count_md_steps++;
	} // end MD timestep loop
	} // end if MD
// ============================= END MOLECULAR DYNAMICS =======================================


// ============================= SINGLE POINT ENERGY ==========================================
    else if (system.constants.mode == "sp") {
        printf("\n| ==================================== |\n");
	    printf("|  BEGINNING SINGLE POINT CALCULATION  |\n");
	    printf("| ==================================== |\n\n");

        if (system.pbc.a==0 && system.pbc.b==0 && system.pbc.c==0 && system.pbc.alpha==0 && system.pbc.beta==0 && system.pbc.gamma==0)
            system.constants.all_pbc=0; // force no PBC if no box given

        singlePointEnergy(system);


    } // end if Single-Point mode (not md or mc)
// ============================ END SINGLE POINT ENERGY =======================================


// ============================ OPTIMIZATION ==================================================
    else if (system.constants.mode == "opt") {
        printf("\n| ==================================== |\n");
        printf("|   BEGINNING STRUCTURE OPTIMIZATION   |\n");
        printf("| ==================================== |\n\n");

        if (system.pbc.a==0 && system.pbc.b==0 && system.pbc.c==0 && system.pbc.alpha==0 && system.pbc.beta==0 && system.pbc.gamma==0) {
            system.constants.all_pbc=0; // force no PBC if no box given
            system.pbc.cutoff = 25.0; // default cutoff
        }
        findBonds(system);
        setBondingParameters(system);
        optimize(system);
    } // end optimization mode



    // Final timing stats.
	std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
    system.constants.time_elapsed = (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()) /1000000.0;
	printf("Total wall time = %f s\n",system.constants.time_elapsed);


} // end main()
