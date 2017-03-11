#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <stdlib.h>

using namespace std;


void computeInitialValues(System &system) {
	
    double t=1.0;

    // MASS OF SYSTEM
    system.stats.totalmass = 0.0; system.stats.movablemass = 0.0; system.stats.frozenmass = 0.0;
	for (int c=0; c<system.molecules.size();c++) {
		for (int d=0; d<system.molecules[c].atoms.size(); d++) {
            double thismass = system.molecules[c].atoms[d].m/system.constants.cM/system.constants.NA;
			system.stats.totalmass += thismass; // total mass in g
		    if (system.molecules[c].MF == "M")
                system.stats.movablemass += thismass;
            else if (system.molecules[c].MF == "F")
                system.stats.frozenmass += thismass;
        }
	}

    // N_movables (sorbates, usually)
    system.stats.current_Nmov_sum += system.stats.count_movables;
    system.constants.initial_sorbates = system.stats.count_movables;
    system.stats.Nmov_average = system.stats.current_Nmov_sum/t;

	// ENERGY
    double* placeholder = getTotalPotential(system, system.constants.potential_form); // getTotPot sends values to stat vars
            system.stats.rd_average = system.stats.current_rd_sum/t;
                system.stats.lj_lrc_avg = system.stats.current_lj_lrc_sum/t;
                system.stats.lj_self_lrc_avg = system.stats.current_lj_lrc_sum/t;
                system.stats.lj_avg = system.stats.current_lj_sum/t;

            system.stats.es_average = system.stats.current_es_sum/t;
                system.stats.coulombic_self_avg = system.stats.current_coulombic_self_sum/t;
                system.stats.coulombic_real_avg = system.stats.current_coulombic_real_sum/t;
                system.stats.coulombic_reciprocal_avg = system.stats.current_coulombic_reciprocal_sum/t;

        	system.stats.polar_average = system.stats.current_polar_sum/t;

            system.stats.energy_average = system.stats.current_energy_sum/t;		
    
    // CHEMICAL POTENTIAL dE/dN
    system.stats.chemical_potential = (system.stats.energy_average - system.constants.initial_energy)
            / (system.stats.Nmov_average - system.constants.initial_sorbates);		

	// VOLUME
	system.stats.volume = system.pbc.volume; 
	system.stats.current_volume_sum += system.stats.volume;
		system.stats.volume_average = system.stats.current_volume_sum/t;

	// DENSITY
	system.stats.density = system.stats.movablemass/(system.stats.volume*1e-24); // that's mass in g /mL	
	system.stats.current_density_sum += system.stats.density;
		system.stats.density_average = system.stats.current_density_sum/t;

    // WT %
    system.stats.wt_percent = (system.stats.movablemass / system.stats.totalmass)*100;
    system.stats.wt_percent_ME = (system.stats.movablemass / system.stats.frozenmass)*100;
    system.stats.current_wt_percent_sum += system.stats.wt_percent;
    system.stats.current_wt_percent_ME_sum += system.stats.wt_percent_ME;
        system.stats.wt_percent_average = system.stats.current_wt_percent_sum/t;
        system.stats.wt_percent_ME_average = system.stats.current_wt_percent_ME_sum/t;			

	// COMPRESSIBILITY FACTOR Z = PV/nRT  =  atm*L / (mol * J/molK * K)
	// GOOD FOR HOMOGENOUS GASES ONLY!!
    double n_moles_sorb = system.stats.movablemass/(system.proto.get_mass()*1000*system.constants.NA);
	system.stats.Z = (system.constants.pres*(system.stats.volume*1e-27) * 101.325 ) // PV
            / ( (n_moles_sorb) * system.constants.R  * system.constants.temp ); // over nRT
	system.stats.current_z_sum += system.stats.Z;
		system.stats.z_average = system.stats.current_z_sum/t;


}



void computeAverages(System &system) {
    system.checkpoint("started computeAverages()");
    double t = (double)system.stats.MCstep;
    double t1 = t+1.0; // for values which already got initial-value treatment

    // MC MOVE ACCEPT STATS
	system.stats.total_accepts = system.stats.insert_accepts + system.stats.remove_accepts + system.stats.displace_accepts + system.stats.volume_change_accepts;
	system.stats.total_attempts = system.stats.insert_attempts + system.stats.remove_attempts + system.stats.displace_attempts + system.stats.volume_attempts;
    //printf("total attempts = %i\n",system.stats.total_attempts);
    
    // BOLTZMANN AVERAGES
	system.stats.bf_avg = (system.stats.insert_bf_sum + system.stats.remove_bf_sum + system.stats.displace_bf_sum + system.stats.volume_change_bf_sum)/4.0/t;
	system.stats.ibf_avg = system.stats.insert_bf_sum/t;
	system.stats.rbf_avg = system.stats.remove_bf_sum/t;
	system.stats.dbf_avg = system.stats.displace_bf_sum/t;
	system.stats.vbf_avg = system.stats.volume_change_bf_sum/t;
    system.checkpoint("got BF averages");

    // ACCEPTANCE RATIO AVERAGES
    system.stats.ar_tot = (system.stats.total_attempts == 0 ) ? 0 : system.stats.total_accepts/(double)system.stats.total_attempts;
    system.stats.ar_ins = (system.stats.insert_attempts == 0 ) ? 0 : system.stats.insert_accepts/(double)system.stats.insert_attempts;
    system.stats.ar_rem = (system.stats.remove_attempts == 0 ) ? 0 :system.stats.remove_accepts/(double)system.stats.remove_attempts;
    system.stats.ar_dis = (system.stats.displace_attempts == 0) ? 0 :system.stats.displace_accepts/(double)system.stats.displace_attempts;
    system.stats.ar_vol = (system.stats.volume_attempts == 0) ? 0 : system.stats.volume_change_accepts/(double)system.stats.volume_attempts;
    system.checkpoint("got AR averages");

    // PERCENT MOVES
    system.stats.ins_perc = (system.stats.total_accepts == 0) ? 0 : system.stats.insert_accepts/(double)system.stats.total_accepts*100;
    system.stats.rem_perc = (system.stats.total_accepts == 0) ? 0 : system.stats.remove_accepts/(double)system.stats.total_accepts*100;
    system.stats.dis_perc = (system.stats.total_accepts == 0) ? 0 : system.stats.displace_accepts/(double)system.stats.total_accepts*100;
    system.stats.vol_perc = (system.stats.total_accepts == 0) ? 0 : system.stats.volume_change_accepts/(double)system.stats.total_accepts*100;

    system.checkpoint("done with boltzmann stuff.");
	// MASS OF SYSTEM
    if (system.constants.ensemble == "uvt") { // only uvt changes mass.
	system.stats.totalmass = 0.0; system.stats.movablemass = 0.0; system.stats.frozenmass = 0.0;
	for (int c=0; c<system.molecules.size();c++) {
		for (int d=0; d<system.molecules[c].atoms.size(); d++) {
            double thismass = system.molecules[c].atoms[d].m/system.constants.cM/system.constants.NA;
			system.stats.totalmass += thismass; // total mass in g
		    if (system.molecules[c].MF == "M")
                system.stats.movablemass += thismass;
            else if (system.molecules[c].MF == "F")
                system.stats.frozenmass += thismass;
        }
	}
    }

    // N_movables (sorbates, usually)
    system.stats.current_Nmov_sum += system.stats.count_movables;
    system.stats.Nmov_average = system.stats.current_Nmov_sum/t1;

	// ENERGY
    // 1) RD
    system.stats.rd_average = system.stats.current_rd_sum/t1;
        system.stats.lj_lrc_avg = system.stats.current_lj_lrc_sum/t1;
        system.stats.lj_self_lrc_avg = system.stats.current_lj_lrc_sum/t1;
        system.stats.lj_avg = system.stats.current_lj_sum/t1;
    // 2) ES
    system.stats.es_average = system.stats.current_es_sum/t1;
        system.stats.coulombic_self_avg = system.stats.current_coulombic_self_sum/t1;
        system.stats.coulombic_real_avg = system.stats.current_coulombic_real_sum/t1;
        system.stats.coulombic_reciprocal_avg = system.stats.current_coulombic_reciprocal_sum/t1;
    // 3) POLAR
	system.stats.polar_average = system.stats.current_polar_sum/t1;
    // 4) TOTAL
    system.stats.energy_average = system.stats.current_energy_sum/t1;		
    
    // CHEMICAL POTENTIAL dE/dN
    system.stats.chemical_potential = (system.stats.energy_average - system.constants.initial_energy)
            / (system.stats.Nmov_average - system.constants.initial_sorbates);		
    
    // QST
    if (system.constants.ensemble != "nve") { // T must be fixed for Qst
        // NU (for qst)
        system.stats.NU = system.stats.totalU*system.stats.count_movables;
        system.stats.current_NU_sum += system.stats.NU;
        system.stats.NU_average = system.stats.current_NU_sum / t;

        // Nsq (for qst)
        system.stats.Nsq = system.stats.count_movables * system.stats.count_movables;
        system.stats.current_Nsq_sum += system.stats.Nsq;
        system.stats.Nsq_average = system.stats.current_Nsq_sum / t;

        // Qst
            system.stats.qst = -(system.stats.NU_average - system.stats.Nmov_average * system.stats.energy_average);
            system.stats.qst /= (system.stats.Nsq_average - system.stats.Nmov_average * system.stats.Nmov_average);
            system.stats.qst += system.constants.temp; 
            system.stats.qst *= system.constants.kb * system.constants.NA * 1e-3; // to kJ/mol
            if (0 != system.stats.Nsq_average - system.stats.Nmov_average * system.stats.Nmov_average) {
                system.stats.current_qst_sum += system.stats.qst;
                system.stats.qst_average = system.stats.current_qst_sum / system.stats.qst_counter;
                system.stats.qst_counter++;
                //printf("qst = %f; sum = %f, avg = %f, denom = %f\n", system.stats.qst, system.stats.current_qst_sum,system.stats.qst_average, (system.stats.Nsq_average - system.stats.Nmov_average * system.stats.Nmov_average));

            }
    }

	// VOLUME
	system.stats.volume = system.pbc.volume; 
	system.stats.current_volume_sum += system.stats.volume;
		system.stats.volume_average = system.stats.current_volume_sum/t1;

	// DENSITY
	system.stats.density = system.stats.movablemass/(system.stats.volume*1e-24); // that's mass in g /mL	
	system.stats.current_density_sum += system.stats.density;
		system.stats.density_average = system.stats.current_density_sum/t1;

    // WT %
    system.stats.wt_percent = (system.stats.movablemass / system.stats.totalmass)*100;
    system.stats.wt_percent_ME = (system.stats.movablemass / system.stats.frozenmass)*100;
    system.stats.current_wt_percent_sum += system.stats.wt_percent;
    system.stats.current_wt_percent_ME_sum += system.stats.wt_percent_ME;
        system.stats.wt_percent_average = system.stats.current_wt_percent_sum/t1;
        system.stats.wt_percent_ME_average = system.stats.current_wt_percent_ME_sum/t1;			

	// COMPRESSIBILITY FACTOR Z = PV/nRT  =  atm*L / (mol * J/molK * K)
	// GOOD FOR HOMOGENOUS GASES ONLY!!
    double n_moles_sorb = system.stats.movablemass/(system.proto.get_mass()*1000*system.constants.NA);
	system.stats.Z = (system.constants.pres*(system.stats.volume*1e-27) * 101.325 ) // PV
            / ( (n_moles_sorb) * system.constants.R  * system.constants.temp ); // over nRT
	system.stats.current_z_sum += system.stats.Z;
		system.stats.z_average = system.stats.current_z_sum/t1;

	
}
