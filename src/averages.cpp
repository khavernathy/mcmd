#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <stdlib.h>

using namespace std;


void computeInitialValues(System &system) {
	
    // MASS OF SYSTEM
    system.stats.totalmass.value = 0.0; system.stats.movablemass.value = 0.0; system.stats.frozenmass.value = 0.0;
	for (int c=0; c<system.molecules.size();c++) {
		for (int d=0; d<system.molecules[c].atoms.size(); d++) {
            double thismass = system.molecules[c].atoms[d].m/system.constants.cM/system.constants.NA;
			system.stats.totalmass.value += thismass; // total mass in g
		    if (system.molecules[c].MF == "M")
                system.stats.movablemass.value += thismass;
            else if (system.molecules[c].MF == "F")
                system.stats.frozenmass.value += thismass;
        }
	}

    // N_movables (sorbates, usually)
    system.constants.initial_sorbates = system.stats.count_movables;
    system.stats.Nmov.value = system.stats.count_movables;
    system.stats.Nmov.average = system.stats.count_movables;

	// ENERGY
    double* placeholder = getTotalPotential(system, system.constants.potential_form); // getTotPot sends values to stat vars
            system.stats.rd.average = system.stats.rd.value;
                system.stats.lj_lrc.average = system.stats.lj_lrc.value;
                system.stats.lj_self_lrc.average = system.stats.lj_lrc.value;
                system.stats.lj.average = system.stats.lj.value;

            system.stats.es.average = system.stats.es.value;
                system.stats.es_self.average = system.stats.es_self.value;
                system.stats.es_real.average = system.stats.es_real.value;
                system.stats.es_recip.average = system.stats.es_recip.value;

        	system.stats.polar.average = system.stats.polar.value;

            system.stats.potential.average = system.stats.potential.value;		
            system.constants.initial_energy = system.stats.potential.value;    

    // CHEMICAL POTENTIAL dE/dN
    if (system.constants.ensemble != "npt") {
    system.stats.chempot.value = (system.stats.potential.average - system.constants.initial_energy)
            / (system.stats.Nmov.average - system.constants.initial_sorbates);		
    }
    
	// VOLUME
	system.stats.volume.average = system.pbc.volume; 
    system.stats.volume.value = system.pbc.volume;

	// DENSITY
	system.stats.density.value = system.stats.movablemass.value/(system.stats.volume.value*1e-24); // that's mass in g /mL	
		system.stats.density.average = system.stats.density.value;

    // WT %
    system.stats.wtp.value = (system.stats.movablemass.value / system.stats.totalmass.value)*100;
    system.stats.wtpME.value = (system.stats.movablemass.value / system.stats.frozenmass.value)*100;
        system.stats.wtp.average = system.stats.wtp.value;
        system.stats.wtpME.average = system.stats.wtpME.value;

	// COMPRESSIBILITY FACTOR Z = PV/nRT  =  atm*L / (mol * J/molK * K)
	// GOOD FOR HOMOGENOUS GASES ONLY!!
    double n_moles_sorb = system.stats.movablemass.value/(system.proto.get_mass()*1000*system.constants.NA);
	system.stats.z.value = (system.constants.pres*(system.stats.volume.value*1e-27) * 101.325 ) // PV
            / ( (n_moles_sorb) * system.constants.R  * system.constants.temp ); // over nRT
		system.stats.z.average = system.stats.z.value;


}



void computeAverages(System &system) {
    system.checkpoint("started computeAverages()");
    double t = system.stats.MCstep;

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
    if (system.constants.ensemble == "uvt" && system.stats.Nmov.value != system.last.Nmov) { // only uvt changes mass
	system.stats.totalmass.value = 0.0; system.stats.movablemass.value = 0.0; system.stats.frozenmass.value = 0.0;
	for (int c=0; c<system.molecules.size();c++) {
		for (int d=0; d<system.molecules[c].atoms.size(); d++) {
            double thismass = system.molecules[c].atoms[d].m/system.constants.cM/system.constants.NA;
			system.stats.totalmass.value += thismass; // total mass in g
		    if (system.molecules[c].MF == "M")
                system.stats.movablemass.value += thismass;
            else if (system.molecules[c].MF == "F")
                system.stats.frozenmass.value += thismass;
        }
	}
    }

    // N_movables (sorbates, usually)
    system.stats.Nmov.value = system.stats.count_movables;
    system.stats.Nmov.calcNewStats();

	// ENERGY
    system.stats.rd.calcNewStats();
        system.stats.lj.calcNewStats();
        system.stats.lj_lrc.calcNewStats();
        system.stats.lj_self_lrc.calcNewStats();
    system.stats.es.calcNewStats();
        system.stats.es_real.calcNewStats();
        system.stats.es_self.calcNewStats();
        system.stats.es_recip.calcNewStats();
    system.stats.polar.calcNewStats();
    system.stats.potential.calcNewStats();

    // CHEMICAL POTENTIAL dE/dN
    if (system.constants.ensemble != "npt") {
        system.stats.chempot.value = (system.stats.potential.average - system.constants.initial_energy)
                / (system.stats.Nmov.average - system.constants.initial_sorbates);		
        system.stats.chempot.calcNewStats();
    }
    // QST
    if (system.constants.ensemble == "uvt") { // T must be fixed for Qst
        // NU (for qst)
        system.stats.NU.value = system.stats.potential.value*system.stats.count_movables;
        system.stats.NU.calcNewStats();

        // Nsq (for qst)
        system.stats.Nsq.value = system.stats.count_movables * system.stats.count_movables;
        system.stats.Nsq.calcNewStats();

        // Qst
            if (0 != system.stats.Nsq.average - system.stats.Nmov.average * system.stats.Nmov.average) {
            double qst = -(system.stats.NU.average - system.stats.Nmov.average * system.stats.potential.average);
            qst /= (system.stats.Nsq.average - system.stats.Nmov.average * system.stats.Nmov.average);
            qst += system.constants.temp; 
            qst *= system.constants.kb * system.constants.NA * 1e-3; // to kJ/mol
                system.stats.qst.value = qst;
                system.stats.qst.calcNewStats();

            }
    }

	// VOLUME
	system.stats.volume.value = system.pbc.volume; 
        system.stats.volume.calcNewStats();

	// DENSITY
	system.stats.density.value = system.stats.movablemass.value/(system.stats.volume.value*1e-24); // that's mass in g /mL	
        system.stats.density.calcNewStats();

    // WT %
    system.stats.wtp.value = (system.stats.movablemass.value / system.stats.totalmass.value)*100;
    system.stats.wtpME.value = (system.stats.movablemass.value / system.stats.frozenmass.value)*100;
        system.stats.wtp.calcNewStats();
        system.stats.wtpME.calcNewStats();

	// COMPRESSIBILITY FACTOR Z = PV/nRT  =  atm*L / (mol * J/molK * K)
	// GOOD FOR HOMOGENOUS GASES ONLY!!
    double n_moles_sorb = system.stats.movablemass.value/(system.proto.get_mass()*1000*system.constants.NA);
	system.stats.z.value = (system.constants.pres*(system.stats.volume.value*1e-27) * 101.325 ) // PV
            / ( (n_moles_sorb) * system.constants.R  * system.constants.temp ); // over nRT
        system.stats.z.calcNewStats();

    system.checkpoint("finished computeAverages()");	
}
