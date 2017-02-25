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
	
            if (myvector.size() != 0) {

			//temporary class instance current_atom
			Atom current_atom;
            //master_index++;
            //current_atom.master_index = master_index;
			current_atom.name = myvector[2];
            // I have a database of defaults in constants.cpp
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
				system.proto.name = myvector[3];
				system.proto.MF = myvector[4];
				system.proto.PDBID = stoi(myvector[5]);
			}
			current_atom.mol_PDBID = stoi(myvector[5]);
			//current_mol_id = stoi(myvector[5]);
			current_atom.pos[0] = stod(myvector[6]); 
			current_atom.pos[1] = stod(myvector[7]);
			current_atom.pos[2] = stod(myvector[8]);
			current_atom.prevpos[0] = current_atom.pos[0];
			current_atom.prevpos[1] = current_atom.pos[1];
			current_atom.prevpos[2] = current_atom.pos[2];
            current_atom.C = stod(myvector[10]);
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
                else if (myvector[4] == "F")
                    system.stats.count_frozens++;

			}
					
			// add atom to current molecule by default
			system.molecules[mol_counter].atoms.push_back(current_atom);
            system.molecules[mol_counter].mass += current_atom.m;

			// and add current atom to prototype only if its in the first mover
			if (current_mol_id == first_mover_id) {
				system.proto.atoms.push_back(current_atom);
                system.proto.mass += current_atom.m;	
			}
			system.constants.total_atoms++;		
	
			// this would print the whole line.
			//cout << line << '\n';
            } // end if vector size nonzero
		}
		myfile.close();
		
	}
	else { printf("ERROR: Unable to open %s. Exiting.\n",filename.c_str()); std::exit(0); } 
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
		double* comfinal = centerOfMass(system);
		myfile << "COM   " << to_string(comfinal[0]) << "   " << to_string(comfinal[1]) << "   " << to_string(comfinal[2]) << "\n";
        	delete[] comfinal;
	}

  	myfile.close();
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
            system.molecules[j].atoms[i].C,  // 11
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
            system.molecules[j].atoms[i].C,  // 11
            system.molecules[j].atoms[i].polar, // 12
            system.molecules[j].atoms[i].eps,  //13
            system.molecules[j].atoms[i].sig); //14
        }
        } // end for atoms
	} // end for molecules

fclose(f);
}

/* WRITE RUNNING ENERGY AVERAGE EVERY CORRTIME */
void writeEnergy(System &system, double energy, int step) {
	ofstream myfile;
	myfile.open (system.constants.energy_output, ios_base::app);
	myfile << to_string(step) + " " + to_string(energy) + "\n";
	myfile.close();
}


/* WRITE RUNNING DENSITY AVERAGE EVERY CORRTIME */
void writeDensity(System &system, double density, int step) {
	ofstream myfile;
	myfile.open (system.constants.density_output, ios_base::app);
	myfile << to_string(step) + " " + to_string(density) + "\n";
	myfile.close();
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
			
			} else if (!strcasecmp(lc[0].c_str(), "x_length")) {
				system.constants.x_length = atof(lc[1].c_str());
				std::cout << "Got x_length = " << lc[1].c_str() << " A"; printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "y_length")) {
				system.constants.y_length = atof(lc[1].c_str());
				std::cout << "Got y_length = " << lc[1].c_str() << " A"; printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "z_length")) {
				system.constants.z_length = atof(lc[1].c_str());
				std::cout << "Got z_length = " << lc[1].c_str() << " A"; printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "input_atoms")) {
				system.constants.atom_file = lc[1].c_str();
				std::cout << "Got input atoms file name = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "restart_pdb")) {
				system.constants.restart_pdb = lc[1].c_str();
				std::cout << "Got restart output file = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "energy_output")) {
				system.constants.energy_output = lc[1].c_str();
				std::cout << "Got energy output file = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "density_output")) {
				system.constants.density_output = lc[1].c_str();
				std::cout << "Got density output file = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "stepsize")) {
				system.constants.stepsize = atoi(lc[1].c_str());
				std::cout << "Got step size = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "finalstep")) {
				system.constants.finalstep = atoi(lc[1].c_str());
				std::cout << "Got total steps = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "md_corrtime")) {
				system.constants.md_corrtime = atoi(lc[1].c_str());
				std::cout << "Got MD corrtime = " << lc[1].c_str(); printf("\n");
			
            } else if (!strcasecmp(lc[0].c_str(), "md_init_vel")) {
                system.constants.md_init_vel = atof(lc[1].c_str());
                std::cout << "Got MD initial velocity for all molecules = " << lc[1].c_str(); printf("\n");        

			} else if (!strcasecmp(lc[0].c_str(), "md_mode")) {
                system.constants.md_mode = lc[1].c_str();
                std::cout << "Got MD mode = " << lc[1].c_str(); printf("\n");

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
				std::cout << "Got output filename = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "volume_change_option")) {
				system.constants.volume_change_option = lc[1].c_str();
				std::cout << "Got volume change option = " << lc[1].c_str(); printf("\n");
			
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
			
			} else if (!strcasecmp(lc[0].c_str(),  "rotate_prob")) {
				system.constants.rotate_prob = atof(lc[1].c_str());
				std::cout << "Got rotate probability = " << lc[1].c_str(); printf("\n");
			
			} else if (!strcasecmp(lc[0].c_str(), "md_dt")) {
				system.constants.md_dt = atof(lc[1].c_str());
				std::cout << "Got MD timestep = " << lc[1].c_str() << " fs"; printf("\n");
			
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

