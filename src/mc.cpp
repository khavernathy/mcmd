#include <string>
#include <algorithm>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>
#include <potential.cpp>
#include <rotatepoint.cpp>
#include <boltzmann.cpp>
#include <moves.cpp>

// PHAST2 NOT INCLUDED YET

// ================== MAIN MC FUNCTION. HANDLES MOVE TYPES AND BOLTZMANN ACCEPTANCE =================
// ACCORDING TO DIFFERENT ENSEMBLES
void runMonteCarloStep(System &system, string model) {
    system.checkpoint("Entered runMonteCarloStep().");
    system.stats.MCmoveAccepted = false; // reset acceptance checker

	// VOLUME MOVE (only NPT)
	if (system.constants.ensemble == ENSEMBLE_NPT) {
		double VCP = system.constants.vcp_factor/(double)system.stats.count_movables; // Volume Change Probability
		double ranf = (double)rand() / (double)RAND_MAX; // between 0 and 1
        	if (ranf < VCP) {
                    system.checkpoint("doing a volume change move.");
                	changeVolumeMove(system);
                    system.checkpoint("done with volume change move.");
			return; // we tried a volume change, so exit MC step.
		}
	}

	// ADD / REMOVE (only uVT)
	// We'll choose 0-0.5 for add; 0.5-1 for remove (equal prob.)
	if (system.constants.ensemble == ENSEMBLE_UVT) {
		double IRP = system.constants.insert_factor;  
		double ranf = (double)rand() / (double)RAND_MAX; // between 0 and 1
		if (ranf < IRP) {
            //system.checkpoint("doing an add or remove.");
			// we're going to do an add/remove now.
			double ranf2 = (double)rand() / (double)RAND_MAX; // 0->1
			// ADD A MOLECULE
			if (ranf2 < 0.5) {
                system.checkpoint("doing molecule add move.");
				addMolecule(system, model);			
                system.checkpoint("done with molecule add move.");
			} // end add
 
			else { // REMOVE MOLECULE
                system.checkpoint("doing molecule delete move.");
				removeMolecule(system, model);
                system.checkpoint("done with molecule delete move.");
			} // end add vs. remove		
		return; // we did the add or remove so exit MC step.
		} // end doing an add/remove. 
	} // end if uVT		

	
	// DISPLACE / ROTATE :: final default (for all: NPT, uVT, NVT, NVE); NVE has special BoltzFact tho.
	// make sure it's a movable molecule
	system.checkpoint("NOT volume/add/remove :: Starting displace or rotate..");
    displaceMolecule(system, model);
    system.checkpoint("done with displace/rotate");
    return; // done with move, so exit MC step	
}
