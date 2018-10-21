#include <string>
#include <algorithm>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>

// ==================== MOVE ATOMS MD STYLE =========================
/* THIS IS THE MAIN INTEGRATOR FUNCTION. calculateForces() is called within */
void integrate(System &system) {
    system.checkpoint("started integrate()");
    int i,j,n;
    // DEBUG
    int_fast8_t debug=0;
    if (debug == 1) {
        for (j=0; j<system.molecules.size(); j++) {
            if (system.constants.md_mode == MD_MOLECULAR) system.molecules[j].printAll();
            for (i=0; i<system.molecules[j].atoms.size(); i++) {
                if (system.constants.md_mode == MD_ATOMIC) system.molecules[j].atoms[i].printAll();
            }
        }
    }
    // END IF DEBUG

    // Remember, on step 0 the force is pre-calculated, before the integration happens.

    /* ------------- NVE --------------- */
    // NVE velocity verlet
    if (system.constants.ensemble == ENSEMBLE_NVE && system.constants.integrator == INTEGRATOR_VV) {
        acceleration_velocity_verlet(system); // 1/2 step init
        position_verlet(system);
        doPBCcheck(system);
        calculateForces(system);
        acceleration_velocity_verlet(system); // 1/2 step final
    }
    // NVE RK4
    else if (system.constants.ensemble == ENSEMBLE_NVE && system.constants.integrator== INTEGRATOR_RK4) {
       acceleration_velocity_RK4(system); 
       position_RK4(system); 
       doPBCcheck(system);
       calculateForces(system);
    }

    /* ------------- NVT --------------- */
    // NVT velocity verlet, andersen
    else if (system.constants.ensemble == ENSEMBLE_NVT && system.constants.thermostat_type == THERMOSTAT_ANDERSEN && system.constants.integrator == INTEGRATOR_VV) {
        acceleration_velocity_verlet(system); // 1/2 step init
        NVT_thermostat_andersen(system);
        position_verlet(system);
        doPBCcheck(system);
        calculateForces(system);
        acceleration_velocity_verlet(system); // 1/2 step final
    }
    // NVT velocity verlet, nose hoover
    else if (system.constants.ensemble == ENSEMBLE_NVT && system.constants.thermostat_type==THERMOSTAT_NOSEHOOVER && system.constants.integrator == INTEGRATOR_VV) {
        // https://www2.ph.ed.ac.uk/~dmarendu/MVP/MVP03.pdf
        position_VV_NH(system); 
        acceleration_velocity_verlet(system);
        calculateForces(system);
        updateLM(system,0); // first 1/2 step, using v(t)
        updateLM(system,1); // last 1/2 step, using v(t+dt/2)
        velocity_VV_NH_final(system);
    }     
    // NVT RK4, andersen
    else if (system.constants.ensemble == ENSEMBLE_NVT && system.constants.thermostat_type == THERMOSTAT_ANDERSEN && system.constants.integrator == INTEGRATOR_RK4) {
        acceleration_velocity_RK4(system);
        NVT_thermostat_andersen(system); // apply stochastic boltzmann velocities
        position_RK4(system);
        doPBCcheck(system);
        calculateForces(system);
    }
    // NVT RK4, nose hoover 
    else if (system.constants.ensemble == ENSEMBLE_NVT && system.constants.thermostat_type == THERMOSTAT_NOSEHOOVER && system.constants.integrator == INTEGRATOR_RK4) {
        acceleration_velocity_RK4(system); // will pick up NH flag in integrator
        calculateNHLM_now(system); // calculate the friction term
        position_RK4(system);
        doPBCcheck(system);
        calculateForces(system);
    }

    /* ------------- uVT --------------- */
    // uVT velocity verlet, andersen
    else if (system.constants.ensemble == ENSEMBLE_UVT && system.constants.thermostat_type == THERMOSTAT_ANDERSEN && system.constants.integrator == INTEGRATOR_VV) {
        acceleration_velocity_verlet(system); // 1/2 step init
        NVT_thermostat_andersen(system);
        position_verlet(system);
        doPBCcheck(system);
        calculateForces(system);
        acceleration_velocity_verlet(system); // 1/2 step final
    }
    else if (system.constants.ensemble == ENSEMBLE_UVT) {
        printf("ERROR: uVT molecular dynamics only available with Andersen thermostat and velocity verlet integrator. Use `thermostat andersen` and `integrator vv`.\n");
        exit(EXIT_FAILURE);
    }

    system.checkpoint("Done with integrate() function.");
}// end integrate() function
