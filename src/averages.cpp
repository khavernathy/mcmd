#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <string.h>

using namespace std;


void computeInitialValues(System &system) {

    // MASS OF SYSTEM
    system.stats.totalmass.value = 0.0; system.stats.frozenmass.value = 0.0;
    for (int i=0; i<system.proto.size(); i++)
        system.stats.movablemass[i].value = 0.0;
	for (int c=0; c<system.molecules.size();c++) {
        string thismolname = system.molecules[c].name;
		for (int d=0; d<system.molecules[c].atoms.size(); d++) {
            double thismass = system.molecules[c].atoms[d].m/system.constants.cM/system.constants.NA;
			system.stats.totalmass.value += thismass; // total mass in g
		    if (!system.molecules[c].frozen) {
                for (int i=0; i<system.proto.size(); i++) {
                    if (system.proto[i].name == thismolname)
                        system.stats.movablemass[i].value += thismass;
                }

            }
            else if (system.molecules[c].frozen)
                system.stats.frozenmass.value += thismass;
        }
	}

    // N_movables (sorbates, usually)
    system.constants.initial_sorbates = system.stats.count_movables;
    for (int i=0; i<system.molecules.size(); i++) {
        if (system.molecules[i].frozen) continue;
        string thismolname = system.molecules[i].name;
        for (int j=0; j<system.proto.size(); j++)
            if (thismolname == system.proto[j].name)
                system.stats.Nmov[j].value++;
    }

    for (int i=0; i<system.proto.size(); i++)
        system.stats.Nmov[i].average = system.stats.Nmov[i].value;

    if (system.constants.dist_within_option) {
        countAtomInRadius(system, system.constants.dist_within_target, system.constants.dist_within_radius);
        system.stats.dist_within.average = system.stats.dist_within.value;
    }

	// ENERGY
    double placeholder = getTotalPotential(system); // getTotPot sends values to stat vars
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

	// VOLUME
	system.stats.volume.average = system.pbc.volume;
    system.stats.volume.value = system.pbc.volume;

	// DENSITY
    for (int i=0; i<system.proto.size(); i++) {
	    system.stats.density[i].value = system.stats.movablemass[i].value/(system.stats.volume.value*1e-24); // that's mass in g /mL
		system.stats.density[i].average = system.stats.density[i].value;
    }

    // WT % 
    for (int i=0; i<system.proto.size(); i++) {
        system.stats.wtp[i].value = (system.stats.movablemass[i].value / system.stats.totalmass.value)*100;
        system.stats.wtpME[i].value = (system.stats.movablemass[i].value / system.stats.frozenmass.value)*100;
        system.stats.wtp[i].average = system.stats.wtp[i].value;
        system.stats.wtpME[i].average = system.stats.wtpME[i].value;

    }



	// COMPRESSIBILITY FACTOR Z = PV/nRT  =  atm*L / (mol * J/molK * K)
	// GOOD FOR HOMOGENOUS GASES ONLY!!
    double n_moles_sorb = system.stats.movablemass[0].value/(system.proto[0].get_mass()*1000*system.constants.NA);
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
    system.stats.totalmass.value = 0.0; system.stats.frozenmass.value = 0.0;
    for (int i=0; i<system.proto.size(); i++)
        system.stats.movablemass[i].value = 0.0;
	for (int c=0; c<system.molecules.size();c++) {
        string thismolname = system.molecules[c].name;
		for (int d=0; d<system.molecules[c].atoms.size(); d++) {
            double thismass = system.molecules[c].atoms[d].m/system.constants.cM/system.constants.NA;
			system.stats.totalmass.value += thismass; // total mass in g
		    if (!system.molecules[c].frozen) {
                for (int i=0; i<system.proto.size(); i++) {
                    if (system.proto[i].name == thismolname)
                        system.stats.movablemass[i].value += thismass;
                }

            }
            else if (system.molecules[c].frozen)
                system.stats.frozenmass.value += thismass;
        }
	}

    // N_movables (sorbates, usually)
    for (int i=0; i<system.proto.size(); i++) system.stats.Nmov[i].value = 0; // initialize b4 counting.
    for (int i=0; i<system.molecules.size(); i++) {
        if (system.molecules[i].frozen) continue;
        string thismolname = system.molecules[i].name;
        for (int j=0; j<system.proto.size(); j++)
            if (thismolname == system.proto[j].name)
                system.stats.Nmov[j].value++;
    }

    for (int i=0; i<system.proto.size(); i++)
        system.stats.Nmov[i].calcNewStats();

    if (system.constants.dist_within_option) {
        countAtomInRadius(system, system.constants.dist_within_target, system.constants.dist_within_radius);
        system.stats.dist_within.calcNewStats();
    }

    // SELECTIVITY :: N / (other Ns)
    double num, denom;
    for (int i=0; i<system.proto.size(); i++) {
        num = system.stats.Nmov[i].average;
        denom=0;
        for (int j=0; j<system.proto.size(); j++) {
            if (i == j) continue;
            denom += system.stats.Nmov[j].average;
        }
        if (denom != 0.0) system.stats.selectivity[i].value = num/denom; // so selec. will remain unchanged if zero sorbates in system.
        system.stats.selectivity[i].calcNewStats();
    }

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

    // Q (partition function)
    system.stats.Q.value += exp(-system.stats.potential.value / system.constants.temp); // K/K = unitless

    // QST
    if (system.constants.ensemble == ENSEMBLE_UVT && system.proto.size() == 1) { // T must be fixed for Qst

        // NU (for qst)
        system.stats.NU.value = system.stats.potential.value*system.stats.count_movables;
        system.stats.NU.calcNewStats();

        // Nsq (for qst)
        system.stats.Nsq.value = system.stats.count_movables * system.stats.count_movables;
        system.stats.Nsq.calcNewStats();

        // Qst
            if (0 != system.stats.Nsq.average - system.stats.Nmov[0].average * system.stats.Nmov[0].average) {
            double qst = -(system.stats.NU.average - system.stats.Nmov[0].average * system.stats.potential.average);
            qst /= (system.stats.Nsq.average - system.stats.Nmov[0].average * system.stats.Nmov[0].average);
            qst += system.constants.temp;
            qst *= system.constants.kb * system.constants.NA * 1e-3; // to kJ/mol
                system.stats.qst.value = qst;
                system.stats.qst.calcNewStats();

            if (0 != system.stats.Nmov[0].average) {
            double qst_nvt = -system.stats.potential.average * system.constants.NA * system.constants.kb * 1e-3 / system.stats.Nmov[0].average;
                system.stats.qst_nvt.value = qst_nvt;
                system.stats.qst_nvt.calcNewStats();
            }

        }
    }

	// VOLUME
	system.stats.volume.value = system.pbc.volume;
        system.stats.volume.calcNewStats();

	// DENSITY
    for (int i=0; i<system.proto.size(); i++) {
	    system.stats.density[i].value = system.stats.movablemass[i].value/(system.stats.volume.value*1e-24); // that's mass in g /mL
        system.stats.density[i].calcNewStats();
    }

    // WT % / excess adsorption
    for (int i=0; i<system.proto.size(); i++) {
        system.stats.wtp[i].value = (system.stats.movablemass[i].value / system.stats.totalmass.value)*100;
        system.stats.wtpME[i].value = (system.stats.movablemass[i].value / system.stats.frozenmass.value)*100;
            system.stats.wtp[i].calcNewStats();
            system.stats.wtpME[i].calcNewStats();
            
            double mm = system.proto[i].mass * 1000 * system.constants.NA; // molar mass
            double frozmm = system.stats.frozenmass.value * system.constants.NA;// ""
            system.stats.excess[i].value = 1e3*(mm*system.stats.Nmov[i].average - (mm * system.constants.free_volume * system.proto[i].fugacity * system.constants.ATM2REDUCED) / system.constants.temp) / 
            frozmm;  // to mg/g

        system.stats.excess[i].calcNewStats();

    }


	// COMPRESSIBILITY FACTOR Z = PV/nRT  =  atm*L / (mol * J/molK * K)
	// GOOD FOR HOMOGENOUS GASES ONLY!!
    double n_moles_sorb = system.stats.movablemass[0].value/(system.proto[0].get_mass()*1000*system.constants.NA);
	system.stats.z.value = (system.constants.pres*(system.stats.volume.value*1e-27) * 101.325 ) // PV
            / ( (n_moles_sorb) * system.constants.R  * system.constants.temp ); // over nRT
        system.stats.z.calcNewStats();

    system.checkpoint("finished computeAverages()");
}


// special variation for uVT Molecular Dynamics only
void computeAveragesMDuVT(System &system) {
    system.checkpoint("started computeAverages()");

    // MASS OF SYSTEM
    system.stats.totalmass.value = 0.0; system.stats.frozenmass.value = 0.0;
    for (int i=0; i<system.proto.size(); i++)
        system.stats.movablemass[i].value = 0.0;
	for (int c=0; c<system.molecules.size();c++) {
        string thismolname = system.molecules[c].name;
		for (int d=0; d<system.molecules[c].atoms.size(); d++) {
            double thismass = system.molecules[c].atoms[d].m/system.constants.cM/system.constants.NA;
			system.stats.totalmass.value += thismass; // total mass in g
		    if (!system.molecules[c].frozen) {
                for (int i=0; i<system.proto.size(); i++) {
                    if (system.proto[i].name == thismolname)
                        system.stats.movablemass[i].value += thismass;
                }

            }
            else if (system.molecules[c].frozen)
                system.stats.frozenmass.value += thismass;
        }
	}

    // N_movables (sorbates, usually)
    for (int i=0; i<system.proto.size(); i++) system.stats.Nmov[i].value = 0; // initialize b4 counting.
    for (int i=0; i<system.molecules.size(); i++) {
        if (system.molecules[i].frozen) continue;
        string thismolname = system.molecules[i].name;
        for (int j=0; j<system.proto.size(); j++)
            if (thismolname == system.proto[j].name)
                system.stats.Nmov[j].value++;
    }

    for (int i=0; i<system.proto.size(); i++)
        system.stats.Nmov[i].calcNewStats();

    if (system.constants.dist_within_option) {
        countAtomInRadius(system, system.constants.dist_within_target, system.constants.dist_within_radius);
        system.stats.dist_within.calcNewStats();
    }

    // SELECTIVITY :: N / (other Ns)
    double num, denom;
    for (int i=0; i<system.proto.size(); i++) {
        num = system.stats.Nmov[i].average;
        denom=0;
        for (int j=0; j<system.proto.size(); j++) {
            if (i == j) continue;
            denom += system.stats.Nmov[j].average;
        }
        if (denom != 0.0) system.stats.selectivity[i].value = num/denom; // so selec. will remain unchanged if zero sorbates in system.
        system.stats.selectivity[i].calcNewStats();
    }

    // QST
    if (system.constants.ensemble == ENSEMBLE_UVT && system.proto.size() == 1) { // T must be fixed for Qst

        // NU (for qst)
        system.stats.NU.value = system.stats.potential.value*system.stats.count_movables;
        system.stats.NU.calcNewStats();

        // Nsq (for qst)
        system.stats.Nsq.value = system.stats.count_movables * system.stats.count_movables;
        system.stats.Nsq.calcNewStats();

        // Qst
            if (0 != system.stats.Nsq.average - system.stats.Nmov[0].average * system.stats.Nmov[0].average) {
            double qst = -(system.stats.NU.average - system.stats.Nmov[0].average * system.stats.potential.average);
            qst /= (system.stats.Nsq.average - system.stats.Nmov[0].average * system.stats.Nmov[0].average);
            qst += system.constants.temp;
            qst *= system.constants.kb * system.constants.NA * 1e-3; // to kJ/mol
                system.stats.qst.value = qst;
                system.stats.qst.calcNewStats();

            if (0 != system.stats.Nmov[0].average) {
            double qst_nvt = -system.stats.potential.average * system.constants.NA * system.constants.kb * 1e-3 / system.stats.Nmov[0].average;
                system.stats.qst_nvt.value = qst_nvt;
                system.stats.qst_nvt.calcNewStats();
            }

        }
    }

	// VOLUME
	system.stats.volume.value = system.pbc.volume;
        system.stats.volume.calcNewStats();

	// DENSITY
    for (int i=0; i<system.proto.size(); i++) {
	    system.stats.density[i].value = system.stats.movablemass[i].value/(system.stats.volume.value*1e-24); // that's mass in g /mL
        system.stats.density[i].calcNewStats();
    }

    // WT % / excess adsorption
    for (int i=0; i<system.proto.size(); i++) {
        system.stats.wtp[i].value = (system.stats.movablemass[i].value / system.stats.totalmass.value)*100;
        system.stats.wtpME[i].value = (system.stats.movablemass[i].value / system.stats.frozenmass.value)*100;
            system.stats.wtp[i].calcNewStats();
            system.stats.wtpME[i].calcNewStats();
            
            double mm = system.proto[i].mass * 1000 * system.constants.NA; // molar mass
            double frozmm = system.stats.frozenmass.value * system.constants.NA;// ""
            system.stats.excess[i].value = 1e3*(mm*system.stats.Nmov[i].average - (mm * system.constants.free_volume * system.proto[i].fugacity * system.constants.ATM2REDUCED) / system.constants.temp) / 
            frozmm;  // to mg/g

        system.stats.excess[i].calcNewStats();

    }


	// COMPRESSIBILITY FACTOR Z = PV/nRT  =  atm*L / (mol * J/molK * K)
	// GOOD FOR HOMOGENOUS GASES ONLY!!
    double n_moles_sorb = system.stats.movablemass[0].value/(system.proto[0].get_mass()*1000*system.constants.NA);
	system.stats.z.value = (system.constants.pres*(system.stats.volume.value*1e-27) * 101.325 ) // PV
            / ( (n_moles_sorb) * system.constants.R  * system.constants.temp ); // over nRT
        system.stats.z.calcNewStats();

    system.checkpoint("finished computeAverages()");
}
