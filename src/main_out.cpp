#include <time.h>
#include <chrono>
#include <sys/stat.h>

using namespace std;

void md_main_output(System &system) {
// PRINT OUTPUT
            double t = system.stats.MDtime;
			int count_md_steps = system.stats.MDstep;
            double dt = system.constants.md_dt;
            double thing = floor(system.constants.md_ft / system.constants.md_dt);
            long int total_steps = (long int)thing;
            system.constants.end= std::chrono::steady_clock::now();
			system.constants.time_elapsed = (std::chrono::duration_cast<std::chrono::microseconds>(system.constants.end - system.constants.begin_steps).count()) /1000000.0;
			system.constants.sec_per_step = system.constants.time_elapsed/(double)count_md_steps;
			double progress = (((float)count_md_steps)/(float)total_steps*100);
			double ETA = ((system.constants.time_elapsed*(float)total_steps/(float)count_md_steps) - system.constants.time_elapsed)/60.0;
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

            if (system.constants.cuda) printf("MCMD (CUDA) Molecular Dynamics: %s (%s)\n", system.constants.jobname.c_str(), system.constants.inputfile.c_str());
            else printf("MCMD Molecular Dynamics: %s (%s)\n", system.constants.jobname.c_str(), system.constants.inputfile.c_str());
            printf("Input atoms: %s\n",system.constants.atom_file.c_str());
            printf("Ensemble: %s; N_movables = %i N_atoms = %i\n",system.constants.ensemble_str.c_str(), system.stats.count_movables, system.constants.total_atoms);
            printf("Time elapsed = %.2f s = %.4f sec/step; ETA = %.3f min = %.3f hrs\n",system.constants.time_elapsed,system.constants.sec_per_step,ETA,ETA_hrs);
            printf("Step: %i / %li; Progress = %.3f%%; Realtime = %.5f %s\n",count_md_steps,total_steps,progress,outputTime, timeunit.c_str());
            if (system.constants.ensemble == ENSEMBLE_NVT || system.constants.ensemble == ENSEMBLE_UVT) {
                if (system.constants.simulated_annealing)
                    printf("        Input T = %.4f K | simulated annealing (on)\n", system.constants.temp);
                else 
                    printf("        Input T = %.4f K | simulated annealing (off)\n", system.constants.temp);
            }
            printf("      Average T = %.4f +- %.4f K\n", system.stats.temperature.average, system.stats.temperature.sd);
            printf("Instantaneous T = %.4f K\n", system.stats.temperature.value);
            if (system.constants.ensemble == ENSEMBLE_NVT && system.constants.calc_pressure_option) {
                printf("      Average P = %.4f +- %.4f atm\n", system.stats.pressure.average, system.stats.pressure.sd);
                printf("Instantaneous P = %.4f atm\n", system.stats.pressure.value);
            }
            printf("     KE = %.3f kJ/mol (lin: %.3f , rot: %.3f )\n",
                  system.stats.kinetic.value*system.constants.K2KJMOL, system.stats.Klin.value*system.constants.K2KJMOL, system.stats.Krot.value*system.constants.K2KJMOL );
            printf("     PE = %.3f kJ/mol\n",
                  system.stats.potential.value*system.constants.K2KJMOL
                  );
            printf("          RD = %.3f kJ/mol\n", system.stats.rd.value*system.constants.K2KJMOL);
            printf("          ES = %.3f kJ/mol\n", system.stats.es.value*system.constants.K2KJMOL);
            printf("         Pol = %.3f kJ/mol\n", system.stats.polar.value*system.constants.K2KJMOL);
            printf("      Bonded = %.3f kJ/mol\n", system.stats.bonded.value*system.constants.K2KJMOL);
        
            if (system.constants.ensemble == ENSEMBLE_NVE)
                printf("Total E = %.3f :: error = %.3f kJ/mol ( %.3f %% )\n", system.stats.totalE.value*system.constants.K2KJMOL, system.constants.md_NVE_err, system.constants.md_NVE_err/(fabs(system.constants.md_initial_energy_NVE)*system.constants.K2KJMOL)*100.);
            else
                printf("Total E = %.3f kJ/mol\n", system.stats.totalE.value*system.constants.K2KJMOL);
            printf("Average v = %.2f m/s; v_init = %.2f m/s\n",
                system.stats.avg_v.value*1e5, system.constants.md_init_vel*1e5);
            if (system.constants.md_pbc || system.constants.ensemble != ENSEMBLE_UVT) { // for now, don't do diffusion unless PBC is on. (checkInTheBox assumes it)
                for (int sorbid=0; sorbid < system.proto.size(); sorbid++) {
                    printf("Diffusion coefficient of %s = %.4e cm^2 / s\n", system.proto[sorbid].name.c_str(), system.stats.diffusion[sorbid].value);
                    printf("    %s MSD = %.5f A^2\n", system.proto[sorbid].name.c_str(), system.stats.msd[sorbid].value);
		            printf("VACF of %s = %f\n", system.proto[sorbid].name.c_str(), system.stats.vacf[sorbid].value);
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
                    printf("      avg N/f.u. =         %.5f +- %.5f\n", system.stats.Nmov[i].average / (double)system.constants.num_fu, system.stats.Nmov[i].sd / (double)system.constants.num_fu);
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
                printf("Heat capacity = %.5f +- %.5f kJ/molK\n", system.stats.heat_capacity.value, system.stats.heat_capacity.sd);
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
            writeThermo(system, system.stats.totalE.value, system.stats.Klin.value, system.stats.Krot.value, system.stats.potential.value, system.stats.rd.value, system.stats.es.value, system.stats.polar.value, 0.0, system.stats.temperature.value, system.stats.pressure.value, count_md_steps, system.stats.Nmov[0].value);
            // restart file.
            writePDB(system, system.constants.restart_pdb); // containing all atoms
            writePDBrestartBak(system, system.constants.restart_pdb, system.constants.restart_pdb_bak);
            // trajectory file
                if (system.constants.xyz_traj_option)
			        writeXYZ(system,system.constants.output_traj,system.constants.frame,count_md_steps,t, system.constants.xyz_traj_movers_option);

                if (!system.constants.pdb_bigtraj_option) writePDBmovables(system, system.constants.restart_mov_pdb); // only movers restart frame
                if (system.constants.pdb_traj_option) {
                    if (system.constants.pdb_bigtraj_option)
                        writePDBtraj(system, system.constants.restart_pdb, system.constants.output_traj_pdb, t); // copy all-atoms-restart-PDB to PDB trajectory
                    else writePDBtraj(system, system.constants.restart_mov_pdb, system.constants.output_traj_movers_pdb,t); // just movers to PDB trajectory
                }
            system.constants.frame++;
            if (system.stats.radial_dist) {
                radialDist(system);
                writeRadialDist(system);
            }
            if (t != dt && system.constants.histogram_option)
				write_histogram(system.file_pointers.fp_histogram, system.grids.avg_histogram->grid, system);
            if ((system.constants.potential_form == POTENTIAL_LJPOLAR || system.constants.potential_form == POTENTIAL_LJESPOLAR) && system.constants.dipole_output_option) {
                write_dipole(system, count_md_steps);
                write_molec_dipole(system, count_md_steps);
            }
            } // end if N>0, write output files.
}
// end MD main output
