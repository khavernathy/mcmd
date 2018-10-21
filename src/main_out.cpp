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
            if (system.constants.ensemble == ENSEMBLE_NVT && system.constants.calc_pressure_option && system.constants.openmp_threads > 0) {
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
                    printf("Diffusion coefficient of %s = %.4e cm^2/s\n", system.proto[sorbid].name.c_str(), system.stats.diffusion[sorbid].value);
                    printf("    %s MSD = %.5e A^2\n", system.proto[sorbid].name.c_str(), system.stats.msd[sorbid].value);
		            printf("VACF of %s = %.5e A^2/fs^2\n", system.proto[sorbid].name.c_str(), system.stats.vacf[sorbid].value);
                }
            }
            //if (system.stats.Q.value > 0) printf("Q (partition function) = %.5e\n", system.stats.Q.value);
            // uptake data if uVT
            if (system.constants.ensemble == ENSEMBLE_UVT) {
			for (int i=0; i<system.proto.size(); i++) {
                double mmolg = system.stats.wtpME[i].average * 10 / (system.proto[i].mass*system.constants.amu2kg*1000*system.constants.NA);
                double cm3gSTP = mmolg*22.4;
                double mgg = mmolg * (system.proto[i].mass*system.constants.amu2kg*1000*system.constants.NA);
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


void mc_main_output(System &system) {
            // TIMING
            system.constants.end= std::chrono::steady_clock::now();
            system.constants.time_elapsed = (std::chrono::duration_cast<std::chrono::microseconds>(system.constants.end - system.constants.begin_steps).count()) /1000000.0;
            int t = system.stats.MCstep;
            system.constants.sec_per_step = system.constants.time_elapsed/t;
            double progress = (((float)t)/system.constants.finalstep*100);
            double ETA = ((system.constants.time_elapsed*system.constants.finalstep/t) - system.constants.time_elapsed)/60.0;
            double ETA_hrs = ETA/60.0;
            double efficiency = system.stats.MCeffRsq / system.constants.time_elapsed;

			// PRINT MAIN OUTPUT
            if (system.constants.ensemble == ENSEMBLE_UVT && system.constants.bias_uptake_switcher)
                printf("MCMD Monte Carlo: %s (%s): Loading bias (%s %i %s)\n", system.constants.jobname.c_str(), system.constants.inputfile.c_str(), "on:", (int)ceil(system.constants.bias_uptake), "molecules");
            else if (system.constants.ensemble == ENSEMBLE_UVT && !system.constants.bias_uptake_switcher)
                printf("MCMD Monte Carlo: %s (%s): Loading bias (%s)\n", system.constants.jobname.c_str(), system.constants.inputfile.c_str(), "off");
            else
			    printf("MCMD Monte Carlo: %s (%s)\n", system.constants.jobname.c_str(), system.constants.inputfile.c_str());
            printf("Input atoms: %s\n",system.constants.atom_file.c_str());
            if (!system.constants.simulated_annealing)
			    printf("Ensemble: %s; T = %.3f K\n", system.constants.ensemble_str.c_str(), system.constants.temp);
			else
                printf("Ensemble: %s; T = %.3f K (Simulated annealing on)\n",system.constants.ensemble_str.c_str(), system.constants.temp);

            printf("Time elapsed = %.2f s = %.4f sec/step; ETA = %.3f min = %.3f hrs\n",system.constants.time_elapsed,system.constants.sec_per_step,ETA,ETA_hrs);
			printf("Step: %i / %li; Progress = %.3f%%; Efficiency = %.3f\n",system.stats.MCstep+system.constants.step_offset,system.constants.finalstep,progress,efficiency);
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
			printf("Total potential avg = %.5f +- %.5f K ( %.3f +- %.3f kJ/mol )\n",system.stats.potential.average, system.stats.potential.sd, system.stats.potential.average*system.constants.kb/1000.0*system.constants.NA, system.stats.potential.sd*system.constants.kb/1000.0*system.constants.NA);
			printf("Volume avg  = %.2f +- %.2f A^3 = %.2f nm^3\n",system.stats.volume.average, system.stats.volume.sd, system.stats.volume.average/1000.0);
			for (int i=0; i<system.proto.size(); i++) {
                double mmolg = system.stats.wtpME[i].average * 10 / (system.proto[i].mass*system.constants.amu2kg*1000*system.constants.NA);
                double cm3gSTP = mmolg*22.4;
                double mgg = mmolg * (system.proto[i].mass*system.constants.amu2kg*1000*system.constants.NA);
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
                writeXYZ(system,system.constants.output_traj,system.constants.frame,t,0,system.constants.xyz_traj_movers_option);
            system.constants.frame++;
            writePDB(system, system.constants.restart_pdb); // all atoms
            writePDBrestartBak(system, system.constants.restart_pdb, system.constants.restart_pdb_bak);
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

			if ((system.constants.potential_form == POTENTIAL_LJPOLAR || system.constants.potential_form == POTENTIAL_LJESPOLAR) && system.constants.dipole_output_option) {
				write_dipole(system, t);
                write_molec_dipole(system, t);
            }
            } // end if N > 0

}
