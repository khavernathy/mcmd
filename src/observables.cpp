#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <string.h>

using namespace std;

/* these calculations are only a select few
   most calculations are performed within averages.cpp
   e.g. kinetic energy, diffusion, Cv, etc.
*/

double calcTemperature(System &system, int * N_local, double * v2_sum) {
    // returns the whole system temperature
    // get Temperature. Frenkel p64 and p84
    // each sorbate has its own contribution to total T,
    // scaled by representation N_sorb/N_total
    // TODO --- integrate the movable MOF temperature here by its D.O.F.
    double T=0;

    if (system.constants.md_mode == MD_MOLECULAR) {
        double dof_total = 0;
        for (unsigned int z=0; z<system.proto.size(); z++) {
            if (N_local[z] < 1) continue; // no contribution from N=0 sorbates
            double Tcontrib = 1e10*v2_sum[z]*system.proto[z].mass*system.constants.amu2kg/N_local[z]/system.proto[z].dof/system.constants.kb;
            Tcontrib *= ((double)N_local[z]/system.stats.count_movables); // ratio of this type N to total N
            T += Tcontrib;
            dof_total += system.proto[z].dof * N_local[z];
        }
        if (system.constants.ensemble==ENSEMBLE_NVT && system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER)
            T *= (dof_total)/(dof_total+1.); // there is 1 extra DOF for NH thermostat.
    }
    else if (system.constants.md_mode == MD_FLEXIBLE || system.constants.md_mode == MD_ATOMIC) {
        double mv2_sum=0;
        double v2=0;
        unsigned int dof=0;
        double N = (system.constants.flexible_frozen ? system.constants.total_atoms : system.constants.total_atoms - system.stats.count_frozens);
        for (unsigned int i=0;i<system.molecules.size();i++) {
            if (system.molecules[i].frozen && !system.constants.flexible_frozen) continue;
            for (unsigned int j=0;j<system.molecules[i].atoms.size();j++) {
                v2 = dddotprod(system.molecules[i].atoms[j].vel, system.molecules[i].atoms[j].vel);
                mv2_sum += system.molecules[i].atoms[j].mass*system.constants.amu2kg*v2;
            }
        }
        dof = N*3.0 - 3.0; // - (int)system.constants.uniqueBonds.size();
        T = 1e10*mv2_sum/(double)dof/system.constants.kb;
        if (system.constants.ensemble == ENSEMBLE_NVT && system.constants.thermostat_type==THERMOSTAT_NOSEHOOVER)
            T *= dof/(dof+1.0); // there is 1 extra DOF for NH thermostat.
    }
    // temperature from the MOF itself (???)
    /*
    if (system.constants.flexible_frozen) {
        double v2=0.; // vx^2 + vy^2 + vz^2
        for (int j=0;j<system.molecules[0].atoms.size();j++) {
          for (int n=0;n<3;n++)
            v2 += system.molecules[0].atoms[j].vel[n]*system.molecules[0].atoms[j].vel[n];
          // assume DOF of the atom is 3.0
          // TODO :::: address DOF exactly for the whole system
          double Tcontrib = 1e10*v2*system.molecules[0].atoms[j].m/1.0/3.0/system.constants.kb;
          Tcontrib /= system.stats.count_frozens;
          T += Tcontrib;
        }
    }
    */
    return T;
}

double calcPressureNVT(System &system) {
    double P=0;
    system.stats.fdotr_sum.value = system.constants.fdotr_sum;
    system.stats.fdotr_sum.calcNewStats(); 
    double V = system.stats.volume.value;
    int N = (system.constants.flexible_frozen ? system.constants.total_atoms : system.constants.total_atoms - system.stats.count_frozens);
    double rho = N/V;
    double T = system.constants.temp;

    P = rho*T + 1./(3.0*V)*system.stats.fdotr_sum.average;  // in K/A^3
    P *= system.constants.KA32ATM; // to atm
    return P;    
}

double calcDOF(System &system) {
    int i;
    double dof=0;
    // molecular/rigid
    if (system.constants.md_mode == MD_MOLECULAR) {
        for (i=0;i<system.molecules.size();i++) {
            dof += system.molecules[i].dof;
        }
        if (system.constants.ensemble == ENSEMBLE_NVT && system.constants.thermostat_type==THERMOSTAT_NOSEHOOVER)
            dof += 1.0;
    // atomic/flexible
    } else {
        double N = (system.constants.flexible_frozen ? system.constants.total_atoms : system.constants.total_atoms - system.stats.count_frozens);
        dof = 3.0*N - 3.0;
        if (system.constants.ensemble == ENSEMBLE_NVT && system.constants.thermostat_type==THERMOSTAT_NOSEHOOVER)
            dof += 1.0;
    }
    system.constants.DOF = dof;
    return dof;
}
