#include <iostream>
#include <string>
#include <strings.h>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <chrono>

/* READ IN THE STARTING COORDINATES, CHARGES, MOLECULE ID'S FROM PDB */
void readInAtoms(System &system, string filename) {
	string line;
	ifstream myfile (filename); //("Neons.xyz"); // ("rhtMOF9Zn.xyz"); //("water_1527.dat"); // test2.dat
	if (myfile.is_open())
	{
        //int master_index=-1;
		std::string::size_type sz;     // alias of size_t
		// loop through each line
                int current_mol_id=-1; // Initializer. Will be changed.
                int mol_counter=-1;
		bool prototype_made = false;
		bool first_mover_passed = false;
		int first_mover_id = -1;
        Molecule whatev;
        system.proto.push_back(whatev); // the first prototype.
		Molecule current_molecule; // initializer. Will be overwritten
		while ( getline (myfile,line) )
		{
			vector<string> myvector;
      			istringstream iss(line);
			//ostream_iterator<string> out_it (cout,",");	
			copy(
				istream_iterator<string>(iss),
				istream_iterator<string>(),
				back_inserter(myvector) // "normally" out_it goes here.
			);
	
            // skip blank lines
            if (myvector.size() != 0) {
                if (myvector[0] != "ATOM") continue; // skip anything in file that isn't an atom
                if (myvector[2] == "X" && myvector[3] == "BOX") continue; // skip box vertices

			//temporary class instance current_atom
			Atom current_atom;
			current_atom.name = myvector[2];
            // I have a database of defaults in classes.cpp
            // Those defaults will load unless the input column is there
            if (9 < myvector.size() && myvector[9] != "default") current_atom.m = stod(myvector[9])*system.constants.cM;
            else current_atom.m = system.constants.masses[current_atom.name]; 
			
            if (12 < myvector.size() && myvector[12] != "default") current_atom.eps = stod(myvector[12]);
            else current_atom.eps = system.constants.eps[current_atom.name];
			
            if (13 < myvector.size() && myvector[13] != "default") current_atom.sig = stod(myvector[13]);
            else current_atom.sig = system.constants.sigs[current_atom.name];

            if (11 < myvector.size() && myvector[11] != "default") current_atom.polar = stod(myvector[11]);			
            else current_atom.polar = system.constants.polars[current_atom.name];
            //==============================================================
            current_atom.V = 0.0;
			current_atom.K = 0.0;
			current_atom.E = 0.0;
			current_atom.PDBID = stoi(myvector[1]); // pulled from input pdb column 2
			current_atom.mol_name = myvector[3];
			current_atom.MF = myvector[4];
			// flag the first moving molecule as prototype sorbate
			// this will only return true once.
			if (myvector[4] == "M" && first_mover_passed == false) {
				first_mover_passed = true;
				first_mover_id = stoi(myvector[5]);
				system.proto[0].name = myvector[3];
				system.proto[0].MF = myvector[4];
				system.proto[0].PDBID = stoi(myvector[5]);
			}
			current_atom.mol_PDBID = stoi(myvector[5]);
			//current_mol_id = stoi(myvector[5]);
			current_atom.pos[0] = stod(myvector[6]); 
			current_atom.pos[1] = stod(myvector[7]);
			current_atom.pos[2] = stod(myvector[8]);
			current_atom.prevpos[0] = current_atom.pos[0];
			current_atom.prevpos[1] = current_atom.pos[1];
			current_atom.prevpos[2] = current_atom.pos[2];
            current_atom.C = stod(myvector[10]) * system.constants.E2REDUCED;
			//system.atoms.push_back(current_atom); // to master atom list
			
			// create new molecule if needed.
			if (current_mol_id != stoi(myvector[5])) {
				// first save the previous molecule only if not the first instance
				//if (current_mol_id != -1)
				//	system.molecules.push_back(current_molecule);	
				
				// make a new molecule.
				mol_counter++;
				current_mol_id = stoi(myvector[5]);
				Molecule current_molecule;
				current_molecule.PDBID = current_mol_id;
				current_molecule.name = myvector[3];
				current_molecule.MF = myvector[4];
				//current_molecule.init_ang_vel(); // for rotations in MD
                system.molecules.push_back(current_molecule); // make the molecule

				if (myvector[4] == "M")
					system.stats.count_movables++;
                else if (myvector[4] == "F") {
                    system.stats.count_frozen_molecules++; // add +1 frozen molecule
                }

			}
					
			// add atom to current molecule by default
			system.molecules[mol_counter].atoms.push_back(current_atom);
            system.molecules[mol_counter].mass += current_atom.m;

			// and add current atom to prototype only if its in the first mover
			if (current_mol_id == first_mover_id) {
				system.proto[0].atoms.push_back(current_atom);
                system.proto[0].mass += current_atom.m;	
			}
			system.constants.total_atoms++;	// add +1 to master atoms count

            if (myvector[4] == "F")
                system.stats.count_frozens++; // add +1 to frozen atoms count	
			// this would print the whole line.
			//cout << line << '\n';
            } // end if vector size nonzero
		}
		myfile.close();
		
	}
	else {
        if (system.constants.sorbate_name.size() > 0) return;
  
        printf("ERROR: Unable to open %s. Exiting.\n",filename.c_str()); std::exit(0); 
    } 
}


/* WRITE FRAME COORDINATE FOR TRAJECTORY FILE */
void writeXYZ(System &system, string filename, int frame, int step, double realtime) {
	
	ofstream myfile;
	myfile.open (filename, ios_base::app);
	//long unsigned int size = system.atoms.size();

	if (system.constants.com_option == "on")
 		myfile << to_string(system.constants.total_atoms + 1) + "\nFrame " + to_string(frame) + "; Step count: " + to_string(step) + "; Realtime (MD) = " + to_string(realtime) + "fs\n";
	else
		myfile << to_string(system.constants.total_atoms) + "\nFrame " + to_string(frame) + "; Step count: " + to_string(step) + "; Realtime (MD) = " + to_string(realtime) + " fs\n";

	for (int j = 0; j < system.molecules.size(); j++) {
		for (int i = 0; i < system.molecules[j].atoms.size(); i++) {
		myfile << system.molecules[j].atoms[i].name;
		myfile <<  "   ";
		myfile << system.molecules[j].atoms[i].pos[0]; 
		myfile <<  "   ";
		myfile << system.molecules[j].atoms[i].pos[1];
		myfile <<  "   "; 
		myfile << system.molecules[j].atoms[i].pos[2];
		myfile << "\n";
		}
	}

	if (system.constants.com_option == "on") {	
		// get center of mass
        if (system.stats.count_movables > 0) {
		    double* comfinal = centerOfMass(system);
		    myfile << "COM   " << to_string(comfinal[0]) << "   " << to_string(comfinal[1]) << "   " << to_string(comfinal[2]) << "\n";
        	delete[] comfinal;
        } else {
            myfile << "COM   " << "0.0" << "   " << "0.0" << "   " << "0.0" << "\n";
        }
	}

  	myfile.close();
}

/* WRITE PDB TRAJECTORY (TAKES restart.pdb and appends to trajectory file */
void writePDBtraj(System &system, string restartfile, string trajfile, int step) {
    std::ifstream ifile(restartfile.c_str(), std::ios::in);
    std::ofstream ofile(trajfile.c_str(), std::ios::out | std::ios::app);

    ofile << "REMARK step=" << step << "\n";
    ofile << "REMARK total_molecules=" << system.molecules.size() << ", total_atoms=" << system.constants.total_atoms << "\n";
    ofile << "REMARK frozen_molecules=" << system.stats.count_frozen_molecules << ", movable_molecules=" << system.stats.count_movables << "\n";
    ofile << "REMARK frozen_atoms=" << system.stats.count_frozens << ", movable_atoms=" << (system.constants.total_atoms - system.stats.count_frozens) << "\n";

    if (!ifile.is_open()) {
        printf("Error opening restart file!\n");
        exit(1);
    } else {
        ofile << ifile.rdbuf(); // append contents of restartfile into trajfile
    }

    ofile << "ENDMDL\n";
}


/* WRITE PDB RESTART FILE EVERY CORRTIME */
void writePDB(System &system, string filename) {
	remove ( filename.c_str() );
	
FILE *f = fopen(filename.c_str(), "w");
if (f == NULL)
{
    printf("Error opening PDB restart file!\n");
    exit(1);
}
	for (int j=0; j<system.molecules.size(); j++) {
		for (int i=0; i<system.molecules[j].atoms.size(); i++) {
        if (system.constants.pdb_long == "off") {
        // this is default. VMD requires the "true" %8.3f
		fprintf(f, "ATOM  %5i %4s %3s %1s %3i    %8.3f%8.3f%8.3f %3.5f %3.5f %f %f %f\n",
            system.molecules[j].atoms[i].PDBID, // col 2
            system.molecules[j].atoms[i].name.c_str(), // 3
            system.molecules[j].atoms[i].mol_name.c_str(), // 4
            system.molecules[j].atoms[i].MF.c_str(), // 5
            system.molecules[j].atoms[i].mol_PDBID, // 6
            system.molecules[j].atoms[i].pos[0], // 7
            system.molecules[j].atoms[i].pos[1],  // 8
            system.molecules[j].atoms[i].pos[2], //9
            system.molecules[j].atoms[i].m/system.constants.cM, // 10 
            system.molecules[j].atoms[i].C/system.constants.E2REDUCED,  // 11
            system.molecules[j].atoms[i].polar, // 12
            system.molecules[j].atoms[i].eps,  //13
            system.molecules[j].atoms[i].sig); //14
		}
        else if (system.constants.pdb_long == "on") {
        fprintf(f, "ATOM  %5i %4s %3s %1s %3i %8.6f %8.6f %8.6f %3.6f %3.6f %f %f %f\n",
            system.molecules[j].atoms[i].PDBID, // col 2
            system.molecules[j].atoms[i].name.c_str(), // 3
            system.molecules[j].atoms[i].mol_name.c_str(), // 4
            system.molecules[j].atoms[i].MF.c_str(), // 5
            system.molecules[j].atoms[i].mol_PDBID, // 6
            system.molecules[j].atoms[i].pos[0], // 7
            system.molecules[j].atoms[i].pos[1],  // 8
            system.molecules[j].atoms[i].pos[2], //9
            system.molecules[j].atoms[i].m/system.constants.cM, // 10 
            system.molecules[j].atoms[i].C/system.constants.E2REDUCED,  // 11
            system.molecules[j].atoms[i].polar, // 12
            system.molecules[j].atoms[i].eps,  //13
            system.molecules[j].atoms[i].sig); //14
        }
        } // end for atoms
	} // end for molecules



    // and draw the box if user desires
    if (system.constants.draw_box_option == "on") {

        int i,j,k,p,q,diff,l,m,n;
        int box_labels[2][2][2];
        double box_occupancy[3];
        double box_pos[3];
        int last_mol_index = system.molecules.size() - 1;
        int last_mol_pdbid = system.molecules[last_mol_index].PDBID;
        int last_atom_pdbid = system.molecules[last_mol_index].atoms[system.molecules[last_mol_index].atoms.size() - 1].PDBID;
        int atom_box = last_atom_pdbid + 1;
        int molecule_box = last_mol_pdbid + 1;

        // draw the box points
        for(i = 0; i < 2; i++) {
            for(j = 0; j < 2; j++) {
                for(k = 0; k < 2; k++) {

                // make this frozen 
                fprintf(f, "ATOM  ");
                fprintf(f, "%5d", atom_box);
                fprintf(f, " %-4.45s", "X");
                fprintf(f, " %-3.3s ", "BOX");
                fprintf(f, "%-1.1s", "F");
                fprintf(f, " %4d   ", molecule_box);

                // box coords 
                box_occupancy[0] = ((double)i) - 0.5;
                box_occupancy[1] = ((double)j) - 0.5;
                box_occupancy[2] = ((double)k) - 0.5;

                
                for(p = 0; p < 3; p++)
                    for(q = 0, box_pos[p] = 0; q < 3; q++)
                        box_pos[p] += system.pbc.basis[q][p]*box_occupancy[q];

                for(p = 0; p < 3; p++)
                    if(system.constants.pdb_long != "on")
                        fprintf(f, "%8.3f", box_pos[p]);
                    else
                        fprintf(f, "%11.6f ", box_pos[p]);

                // null interactions 
                fprintf(f, " %8.4f", 0.0);
                fprintf(f, " %8.4f", 0.0);
                fprintf(f, " %8.5f", 0.0);
                fprintf(f, " %8.5f", 0.0);
                fprintf(f, " %8.5f", 0.0);
                fprintf(f, "\n");

                box_labels[i][j][k] = atom_box;
                ++atom_box;

                } // for k
            } // for j
        } // for i

        // and draw the connecting lines
        for(i = 0; i < 2; i++) {
            for(j = 0; j < 2; j++) {
                for(k = 0; k < 2; k++) {

                    for(l = 0; l < 2; l++) {
                        for(m = 0; m < 2; m++) {
                            for(n = 0; n < 2; n++) {

                                    diff = fabs(i - l) + fabs(j - m) + fabs(k - n);
                                    if(diff == 1)
                                        fprintf(f, "CONECT %4d %4d\n", box_labels[i][j][k], box_labels[l][m][n]);

                            } // n 
                        } // m 
                    } // l 


                } // k
            } // j 
        } // i 

    } // if draw box is on
    // (end drawing the box in .pbd restart)

fclose(f);
}

/* WRITE RUNNING ENERGY AVERAGE EVERY CORRTIME */
void writeThermo(System &system, double TE, double LKE, double RKE, double PE, double density, double temp, double pressure, int step) {
    FILE *f = fopen(system.constants.thermo_output.c_str(), "a");
    if (f == NULL) { printf("Error opening thermo data file!\n"); exit(1); }

    fprintf(f, "%i  %f  %f  %f  %f  %f  %f  %f\n", 
        step, TE, LKE, RKE, PE, density, temp, pressure);

    fclose(f);
}



/* READ INPUT FILE PARAMETERS AND OPTIONS */
void readInput(System &system, char* filename) {
	printf("Reading input parameters from %s.\n",filename);

	string line;
	ifstream myfile (filename);
	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{
			vector<string> lc;
			istringstream iss(line);
			copy(
				istream_iterator<string>(iss),
				istream_iterator<string>(),
				back_inserter(lc)
			);

			//std::cout << lc[0] << ' ';
			//printf("%\n",lc[0].c_str());			
			if (!lc.empty()) { // ignore blank lines	
	
			if (!strncasecmp(lc[0].c_str(), "!", 1) || (!strncasecmp(lc[0].c_str(), "#", 1))) {
				continue; // treat ! and # as comments
			
			} else if (!strcasecmp(lc[0].c_str(),"name")) {
				system.constants.jobname = lc[1].c_str();
				std::cout << "Got job name = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "mode")) {
				system.constants.mode = lc[1].c_str();
				std::cout << "Got mode = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "ensemble")) {
				system.constants.ensemble = lc[1].c_str();
				std::cout << "Got ensemble = " << lc[1].c_str(); printf("\n");
			
            } else if (!strcasecmp(lc[0].c_str(), "sorbate_name")) {
                system.constants.sorbate_name.push_back(lc[1].c_str());
                std::cout << "Got sorbate model name 1 = " << lc[1].c_str(); printf("\n");

                for (int i=2; i<=10; i++) { // so max sorbates is 10.
                    if (lc.size() >= (i+1)) {
                        system.constants.sorbate_name.push_back(lc[i].c_str());
                        std::cout << "Got sorbate model name " << i << " = " << lc[i].c_str(); printf("\n");
                    }
                }
            
            } else if (!strcasecmp(lc[0].c_str(), "sorbate_fugacities")) {
                system.constants.sorbate_fugacity.push_back(atof(lc[1].c_str()));
                std::cout << "Got fugacity for sorbate 1 = " << lc[1].c_str(); printf("\n");

                for (int i=2; i<=10; i++) {
                    if (lc.size() >= (i+1)) {
                        system.constants.sorbate_fugacity.push_back(atof(lc[i].c_str()));
                        std::cout << "Got fugacity for sorbate " << i << " = " << lc[i].c_str(); printf("\n");
                    }
                }

            // BASIS STUFF. 
            // If user inputs x_length, y_length, z_length, assume 90deg. angles
            } else if (!strcasecmp(lc[0].c_str(), "x_length")) {
                system.pbc.x_length = atof(lc[1].c_str());
                system.pbc.basis[0][0] = atof(lc[1].c_str());
                system.pbc.basis[0][1] = 0;
                system.pbc.basis[0][2] = 0;
                std::cout << "Got x_length = " << lc[1].c_str() << " A"; printf("\n");
            
            } else if (!strcasecmp(lc[0].c_str(), "y_length")) {
                system.pbc.y_length = atof(lc[1].c_str());
                system.pbc.basis[1][1] = atof(lc[1].c_str());
                system.pbc.basis[1][0] = 0;
                system.pbc.basis[1][2] = 0;
                std::cout << "Got y_length = " << lc[1].c_str() << " A"; printf("\n");
            
            } else if (!strcasecmp(lc[0].c_str(), "z_length")) {
                system.pbc.z_length = atof(lc[1].c_str());
                system.pbc.basis[2][2] = atof(lc[1].c_str());
                system.pbc.basis[2][0] = 0;
                system.pbc.basis[2][1] = 0;
                std::cout << "Got z_length = " << lc[1].c_str() << " A"; printf("\n");
        
                system.pbc.calcCarBasis();

            // OR EXACT BASIS INPUT (by vectors)
            } else if (!strcasecmp(lc[0].c_str(), "basis1")) {
                for (int n=0; n<3; n++)
                    system.pbc.basis[0][n] = atof(lc[n+1].c_str());
                system.pbc.x_length = system.pbc.basis[0][0];

                std::cout << "Got basis1 = " << lc[1].c_str() << " " << lc[2].c_str() << " " << lc[3].c_str(); printf("\n");
            } else if (!strcasecmp(lc[0].c_str(), "basis2")) {
                for (int n=0; n<3; n++)
                    system.pbc.basis[1][n] = atof(lc[n+1].c_str());
                system.pbc.y_length = system.pbc.basis[1][1];

                std:: cout << "Got basis2 = " << lc[1].c_str() << " " << lc[2].c_str() << " " << lc[3].c_str(); printf("\n");
            } else if (!strcasecmp(lc[0].c_str(), "basis3")) {
                for (int n=0; n<3; n++)
                    system.pbc.basis[2][n] = atof(lc[n+1].c_str());
                system.pbc.z_length = system.pbc.basis[2][2];

                std:: cout << "Got basis3 = " << lc[1].c_str() << " " << lc[2].c_str() << " " << lc[3].c_str(); printf("\n");

                system.pbc.calcCarBasis();

            } else if (!strcasecmp(lc[0].c_str(), "carbasis")) {
                double a = atof(lc[1].c_str());
                double b = atof(lc[2].c_str());
                double c = atof(lc[3].c_str());
                double alpha = atof(lc[4].c_str());
                double beta = atof(lc[5].c_str());
                double gamma = atof(lc[6].c_str());

                system.pbc.a = a;
                system.pbc.b = b;
                system.pbc.c = c;
                system.pbc.alpha = alpha;
                system.pbc.beta = beta;
                system.pbc.gamma = gamma;
                system.pbc.x_length = a;
                system.pbc.y_length = b;
                system.pbc.z_length = c;

                system.pbc.calcNormalBasis();

                std::cout << "Got .car basis: a,b,c = " << lc[1].c_str() << ", " << lc[2].c_str() << ", " << lc[3].c_str(); printf("\n");
                std::cout << "Got .car basis alpha,beta,gamma = " << lc[4].c_str() << ", " << lc[5].c_str() << ", " << lc[6].c_str(); printf("\n");

	
			} else if (!strcasecmp(lc[0].c_str(), "input_atoms")) {
				system.constants.atom_file = lc[1].c_str();
				std::cout << "Got input atoms file name = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "restart_pdb")) {
				system.constants.restart_pdb = lc[1].c_str();
				std::cout << "Got restart output file = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "thermo_output")) {
				system.constants.thermo_output = lc[1].c_str();
				std::cout << "Got thermo output file = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "stepsize")) {
				system.constants.stepsize = atoi(lc[1].c_str());
				std::cout << "Got step size = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "finalstep")) {
				system.constants.finalstep = atoi(lc[1].c_str());
				std::cout << "Got total steps = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "dist_within")) {
                system.constants.dist_within_option = lc[1].c_str(); 
                std::cout << "Got dist_within option = " << lc[1].c_str(); printf("\n");
        		
            } else if (!strcasecmp(lc[0].c_str(), "dist_within_target")) {
                system.constants.dist_within_target = lc[1].c_str();
                std::cout << "Got dist_within_target atom = " << lc[1].c_str(); printf("\n");

	        } else if (!strcasecmp(lc[0].c_str(), "dist_within_radius")) {
                system.constants.dist_within_radius = atof(lc[1].c_str());
                std::cout << "Got dist_within_radius = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "auto_center")) {
                system.constants.autocenter = lc[1].c_str();
                std::cout << "Got auto-center-atoms-to-origin option = " << lc[1].c_str(); printf("\n");

			} else if (!strcasecmp(lc[0].c_str(), "md_corrtime")) {
				system.constants.md_corrtime = atoi(lc[1].c_str());
				std::cout << "Got MD corrtime = " << lc[1].c_str(); printf("\n");
			
            } else if (!strcasecmp(lc[0].c_str(), "md_init_vel")) {
                system.constants.md_init_vel = atof(lc[1].c_str());
                std::cout << "Got MD initial velocity for all molecules = " << lc[1].c_str(); printf("\n");        
			} else if (!strcasecmp(lc[0].c_str(), "md_mode")) {
                system.constants.md_mode = lc[1].c_str();
                std::cout << "Got MD mode = " << lc[1].c_str(); printf("\n");
            
            } else if (!strcasecmp(lc[0].c_str(), "md_pbc")) {
                system.constants.md_pbc = lc[1].c_str();
                std::cout << "Got MD PBC option = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "mc_pbc")) {
                system.constants.mc_pbc = lc[1].c_str();
                std::cout << "Got MC PBC option = " << lc[1].c_str(); printf("\n");
                if (system.constants.mc_pbc == "off") {
                    system.constants.ewald_es = "off";
                    system.constants.rd_lrc = "off";
                    system.constants.polar_pbc = "off";
                }

            } else if (!strcasecmp(lc[0].c_str(), "simulated_annealing")) {
                system.constants.simulated_annealing = lc[1].c_str();
                std::cout << "Got simulated annealing option = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "simulated_annealing_target")) {
                system.constants.sa_target = atof(lc[1].c_str());
                std::cout << "Got simulated annealing target temperature = " << lc[1].c_str() << " K"; printf("\n");
            
            } else if (!strcasecmp(lc[0].c_str(), "simulated_annealing_schedule")) {
                system.constants.sa_schedule = atof(lc[1].c_str());
                std::cout << "Got simulated annealing schedule = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "draw_box_option")) {
                system.constants.draw_box_option = lc[1].c_str();
                std::cout << "Got draw-box-option for PDB output = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "mc_corrtime")) {
				system.constants.mc_corrtime = atoi(lc[1].c_str());
				std::cout << "Got MC corrtime = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "temperature")) {
				system.constants.temp = atof(lc[1].c_str());
				std::cout << "Got temperature = " << lc[1].c_str() << " K"; printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "pressure")) {
				system.constants.pres = atof(lc[1].c_str());
				std::cout << "Got pressure = " << lc[1].c_str() << " atm"; printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "insert_factor")) {
				system.constants.insert_factor = atof(lc[1].c_str());
				std::cout << "Got insert/delete factor = " << lc[1].c_str(); printf("\n");

			} else if (!strcasecmp(lc[0].c_str(), "volume_change")) {
				system.constants.volume_change = atof(lc[1].c_str());
				std::cout << "Got volume change factor = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "rotate_option")) {
				system.constants.rotate_option = lc[1].c_str();
				std::cout << "Got rotate option = " << lc[1].c_str(); printf("\n");

			} else if (!strcasecmp(lc[0].c_str(), "rotate_angle_factor")) {
				system.constants.rotate_angle_factor = atof(lc[1].c_str());
				std::cout << "Got rotate angle factor = " << lc[1].c_str(); printf("\n");

			} else if (!strcasecmp(lc[0].c_str(), "output_traj")) {
				system.constants.output_traj = lc[1].c_str();
				std::cout << "Got output trajectory XYZ filename = " << lc[1].c_str(); printf("\n");
			
            } else if (!strcasecmp(lc[0].c_str(), "output_traj_pdb")) {
                system.constants.output_traj_pdb = lc[1].c_str();
                std::cout << "Got output trajectory PDB filename = " << lc[1].c_str(); printf("\n");

			} else if (!strcasecmp(lc[0].c_str(), "vcp_factor")) {
				system.constants.vcp_factor = atof(lc[1].c_str());
				std::cout << "Got volume change probability factor = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "displace_factor")) {
				system.constants.displace_factor = atof(lc[1].c_str());
				std::cout << "Got displace factor = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "potential_form")) {
				system.constants.potential_form = lc[1].c_str();
				std::cout << "Got potential form = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "polar_iter")) {
				system.constants.polar_max_iter = atoi(lc[1].c_str());
				std::cout << "Got polarization iterations = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "com_option")) {
				system.constants.com_option = lc[1].c_str();
				std::cout << "Got center-of-mass option = " << lc[1].c_str(); printf("\n");
			
            /* DEPRECATED
			} else if (!strcasecmp(lc[0].c_str(),  "rotate_prob")) {
				system.constants.rotate_prob = atof(lc[1].c_str());
				std::cout << "Got rotate probability = " << lc[1].c_str(); printf("\n");
			*/
			} else if (!strcasecmp(lc[0].c_str(), "md_dt")) {
				system.constants.md_dt = atof(lc[1].c_str());
				std::cout << "Got MD timestep = " << lc[1].c_str() << " fs"; printf("\n");
			
                system.constants.md_thermostat_probab = system.constants.md_thermostat_freq * 
                    exp(-system.constants.md_thermostat_freq * system.constants.md_dt);

			} else if (!strcasecmp(lc[0].c_str(), "md_ft")) {
				system.constants.md_ft = atof(lc[1].c_str());
				std::cout << "Got MD final step = " << lc[1].c_str() << " fs"; printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "md_rotations")) {
                system.constants.md_rotations = lc[1].c_str();
                std::cout << "Got MD rotations option = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "sig_override")) {		
				system.constants.sig_override[lc[1]] = atof(lc[2].c_str());
				std::cout << "Got LJ sigma override for " << lc[1].c_str() << " = " << lc[2].c_str(); printf("\n");	
			
			} else if (!strcasecmp(lc[0].c_str(), "eps_override")) {
				system.constants.eps_override[lc[1]] = atof(lc[2].c_str());
				std::cout << "Got LJ epsilon override for " << lc[1].c_str() << " = " << lc[2].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "radial_dist")) { 
                system.stats.radial_dist = lc[1].c_str();
                std::cout << "Got radial distribution option = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "radial_bin_size")) {
                system.stats.radial_bin_size = atof(lc[1].c_str());
                std::cout << "Got radial bin size = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "radial_max_dist")) {
                system.stats.radial_max_dist = atof(lc[1].c_str());
                std::cout << "Got radial maximum distance = " << lc[1].c_str(); printf("\n");        
    
            } else if (!strcasecmp(lc[0].c_str(), "radial_centroid")) {
                system.stats.radial_centroid = lc[1].c_str();
                std::cout << "Got radial centroid = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "radial_counterpart")) { 
                system.stats.radial_counterpart = lc[1].c_str();
                std::cout << "Got radial counterpart = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "radial_file")) {
                system.stats.radial_file = lc[1].c_str();
                std::cout << "Got radial dist. file = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "checkpoints_option")) {
                system.constants.checkpoints_option = lc[1].c_str();
                std::cout << "Got checkpoints option = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "total_energy")) {
                system.constants.total_energy = atof(lc[1].c_str());
                std::cout << "Got NVE total energy constant = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "rd_lrc")) {
                system.constants.rd_lrc = lc[1].c_str();
                std::cout << "Got RD long range correction option = " << lc[1].c_str(); printf("\n");

            } else if (!strcasecmp(lc[0].c_str(), "ewald_es")) {
                system.constants.ewald_es = lc[1].c_str();
                std::cout << "Got Ewald electrostatics option = " << lc[1].c_str(); printf("\n");        
        
            } else if (!strcasecmp(lc[0].c_str(), "pdb_long")) {
                system.constants.pdb_long = lc[1].c_str();
                std::cout << "Got option for PDB long float output = " << lc[1].c_str(); printf("\n");

            } else { std::cout << "WARNING: INPUT '" << lc[0].c_str() << "' UNRECOGNIZED."; printf("\n");}
			} // end if line not blank
		} // end while reading lines
	} // end if file open
} // end read input function

void paramOverrideCheck(System &system) {
	// LJ sigma/eps override if needed. 
	if ((int)system.constants.sig_override.size() > 0) {	
       map<string, double>::iterator it;
       for ( it = system.constants.sig_override.begin(); it != system.constants.sig_override.end(); it++ ) {
            for (int i=0; i<system.molecules.size(); i++) {
            for (int j=0; j<system.molecules[i].atoms.size(); j++) {
                if (system.molecules[i].atoms[j].name == it->first)
                    system.molecules[i].atoms[j].sig = it->second;
            }
            }            
        } // end map loop     
    } // end if sigma overrides

	if ((int)system.constants.eps_override.size() > 0) {
        map<string, double>::iterator it;
        for (it = system.constants.eps_override.begin(); it != system.constants.eps_override.end(); it++) {
            for (int i=0; i<system.molecules.size(); i++) {
            for (int j=0; j<system.molecules[i].atoms.size(); j++) {
                if (system.molecules[i].atoms[j].name == it->first)
                    system.molecules[i].atoms[j].eps = it->second;
            }
            }
        } // end map loop
    } // end if epsilon override
}

