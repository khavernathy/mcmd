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

    // 1) MOVE ATOMS FROM FORCES/TORQUES
    system.checkpoint("moving particles based on forces.");
    position(system);
    // END POSITION CHANGES
    system.checkpoint("done moving particles. Checking PBC for all particles");

    // 2) CHECK P.B.C. (move the molecule/atom back in the box if needed)
    if (system.constants.md_pbc && system.constants.md_translations) {
        for (j=0; j<system.molecules.size(); j++) {
            if (!system.molecules[j].frozen) {
                checkInTheBox(system,j); // also computes COM
            } // end if movable
	    } // end loop j molecules
    } // end if PBC

    // 3) GET NEW FORCES (AND TORQUES) BASED ON NEW POSITIONS
	system.checkpoint("done checking PBC. Starting calculateForces()");
    calculateForces(system);
    system.checkpoint("Done with calculateForces(). Starting integrator (for a&v)");

    // 4) GET NEW ACCELERATION AND VELOCITY FOR ALL PARTICLES
    acceleration_velocity(system); 
    system.checkpoint("Done with a,v integration. Starting heat bath (if nvt/uvt)");

    // TODO -- flexible MOF thermostat?
    // 5) apply heat bath in constant-temp ensembles
    NVT_thermostat(system);
    system.checkpoint("Done with heatbath if NVT/uVT.");
    





system.checkpoint("Done with integrate() function.");
}// end integrate() function
