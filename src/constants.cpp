#include <stdio.h>
#include <math.h>
#include <string>
#include <map>
#include <vector>

using namespace std;

class Constants {
	public:
		Constants();
		double e,kb,kbk,fs,cC,keSI,ke,eV,cM,cA,cJ,NA,cV,R,mpmc2uff,uff2mpmc,ATM2REDUCED,kg2em, E2REDUCED, TORQUE2REDUCED, FORCE2REDUCED, DEBYE2SKA; // all defined below.
		string jobname="default_jobname";
        string mode; // "mc" or "md" 
        string checkpoints_option="off"; // enables checkpoints for debuggin
        string ensemble; // "nve", "nvt", "nvp", "uvt"
        string atom_file; // input atoms .pdb file
        string output_traj="traj.xyz"; // system trajectory in xyz 
        string restart_pdb="restart.pdb"; // a file to re-kick an unfinished run
        string energy_output="energy.dat"; // logs energy average as f'n of step
        string density_output="density.dat"; // logs avg. density as f'n of step
        string volume_change_option; // kind of useless now but used for NPT,
        string potential_form="lj"; // "lj", "ljes", "ljespolar", "phast2" models for potential
        string com_option="off"; // enables computation of center-of-mass and logs in output_traj
        string rotate_option="on"; // MC ONLY: can deactivate rotates if wanted. 
        string md_rotations="on"; // MD only.
		string rd_lrc="on"; // long range corrections for LJ RD
        string ewald_es="off"; // ewald method for electrostatic potential calculation.
        string pdb_long="off"; // on would force long coordinate/charge output
        map <string,double> masses; // mass database for defaults.
		map <string,double> sigs; // LJ sigmas database for defaults. NOT r_m (as in UFF). Defined below
		map <string,double> eps; // LJ epsions database for defaults. Defined lated
		map <string,double> phast2_c6; map <string,double> phast2_c8; map <string,double> phast2_c10; map <string,double> phast2_sigs; map <string,double> phast2_eps; map <string,double> phast2_polar; // phast2 defaults
		map <string,double> polars; // polarizability database for defaults. Defined below
		double volume; // in A^3
        double temp=0.0; //273.15; // in K
        double pres=1.0; // in atm
        double volume_change=2.5; // a factor used in volume change BF, mpmc default 2.5
        double vcp_factor=1.0; // a factor used in probability for volume change. this / num_sorbates is good per Frenkel
        double displace_factor=2.5; // in A, amount to displace * +-0 to 1 
        double insert_factor=0.5; // probability to do insert or delete (instead of displace/rotate) in uvt
        double rotate_prob=0.5; // prob to rotate instead of displace when displace/rotate is selected
        double rotate_angle_factor=360.0; // angle +- to rotate if rotate selected 
		int stepsize=1; // obvi
        int finalstep; // user defined for MC
        int  mc_corrtime=1000; // default 1k cuz I used that a lot for mpmc research
        int  md_corrtime; // user defined for MD
		double x_length,y_length,z_length,x_max,y_max,z_max,x_min,y_min,z_min; // box parameters, in A
		double cutoff;

        double md_init_vel=99999.99; // placeholder value. Will be overwritten. A / fs^2. User can set. Will be random +- up to this num.
        double v_init = 0.0; // initial velocity for homogenous gases, based on T and molar mass.
        double md_dt=0.1, md_ft; // MD timestep and final time, in fs
        string md_mode = "molecular"; // default is to keep molecules rigid (bonded)
		map <string,double> sig_override;
        map <string,double> eps_override; // feature for overriding preset LJ params (for developing LJ models). 0.0 are defaults which will be overwritten if option is used. sig=A; eps=K
	
    	int total_atoms=0; // actual sites, not "atoms" persay
        int old_total_atoms =0; // for use in resizing thole matrix in uvt

        double total_energy=0.0; // for NVE, in K, user defined

        int initial_sorbates=0.0; // for safekeeping to calculate chemical potential in uVT
        double initial_energy=0.0; // "" ""

        // Ewald (for ES)
        double ewald_alpha=3.5/cutoff;
        double ewald_kmax = 7;
        
        // Wolf (for polarization)
        double polar_wolf_alpha = 0.13;
        double polar_damp = 2.1304;
        double polar_gamma = 1.03;
        int polar_max_iter = 4;
        double **A_matrix, **B_matrix, C_matrix[3][3];
        double polar_precision=0.0;
        int iter_success;
        double polar_rrms;
        double dipole_rrms = 0.0;

        // SYSTEM ENERGIES
        double coulombic_real = 0;
        double coulombic_reciprocal = 0;
        double coulombic_self = 0;
        double total_coulombic = 0;
        double lj_lrc =0;
        double lj_self_lrc = 0;
        double lj = 0;
        double total_rd = 0;
};

class Stats {
	public:
		Stats();
        
        string radial_dist = "off"; // default is no radial distribution
        string radial_file = "radial_distribution.dat"; // default filename for output.
        double radial_bin_size = 0.1; // bin counts will be considered for this range in A
        double radial_max_dist = 10.0; // maximum r to consider in rad. dist.
        vector<long unsigned int> radial_bins; // holds the counters       
        string radial_centroid, radial_counterpart; // the two atoms to get distance between, user def.
		
		double insert_bf_sum = 0; // insert boltzmanns added up
		double remove_bf_sum = 0; // ...
		double displace_bf_sum = 0;
		double volume_change_bf_sum = 0;
		double rotate_bf_sum = 0;

		int insert_accepts = 0; int insert_attempts=0;// Counters for successful moves. uVT only
		int remove_accepts = 0; int remove_attempts=0;// uVT only
		int displace_accepts = 0; int displace_attempts=0;
		int volume_change_accepts = 0; int volume_attempts=0;// NPT only
		int rotate_accepts = 0; int rotate_attempts=0;// mpmc doesn't have this but i do
		int count_movables = 0; // this is the SORBATE MOLECULES, e.g. 2H2 means 2, not 4
        int count_frozens = 0; // frozen ATOMS, not molecules (which is normally just 1)

        double polar_iterations=0;
};

Stats::Stats() {}


class Atom {
	public:
		Atom();
		string name; // element or label, e.g. H or H2G
        string mol_name; // molecule name that the atom belongs to
        string MF; // movable/frozen (M or F)
		int mol_PDBID; // the molecule's PDBID that this atom belongs to
		double m=0.0; // mass, kg. This is the only one I'm keeping SI as of now.
        double eps=0.0; // LJ param in K
        double sig=0.0; // LJ param in A -- the real sigma, not r_m (as in UFF) 
        double polar=0.0; // polarizability in A^3
        double C=0.0; // charge in e
        double V=0.0; // potential energy in K
        double K=0.0; // kinetic energy in K
        double E=0.0; // total energy in K
		int PDBID; // the atom's PDBID (from input)
        double rank_metric;       
 
		vector<double> pos = vector<double>(3); // 0=x; 1=y; 2=z
		vector<double> prevpos = vector<double>(3);
		vector<double> force = vector<double>(3); // force, K / A
		vector<double> vel = vector<double>(3); // velocity, A/fs. Not using b/c I switched to rigid molecular motion
		vector<double> acc = vector<double>(3); // acceleration, A/fs^2
        vector<double> old_acc = vector<double>(3);
		vector<double> torque = vector<double>(3);
        vector<double> dip = vector<double>(3); // dipole in e*A ... I think
		vector<double> newdip = vector<double>(3);
        vector<double> olddip = vector<double>(3);
        double dipole_rrms=0;    
        vector<double> efield = vector<double>(3); // external electric field
        vector<double> efield_self = vector<double>(3); // self electric field
        vector<double> efield_induced = vector<double>(3); // induced e field
        vector<double> efield_induced_change = vector<double>(3); 

        double * get_acc() {
            static double output[3];
            for (int n=0; n<3; n++) output[n] = acc[n];
            return output;
        }

        void calc_acc() {
            for (int n=0; n<3; n++) {
                old_acc[n] = acc[n];
                acc[n] = force[n]*1.3806488e-33 / m; // to A / fs^2
            }
        } 

        void calc_vel(double dt) {
            for (int n=0; n<3; n++) vel[n] = vel[n] + 0.5*(acc[n] + old_acc[n])*dt; // that's where VV comes into play. 1/2 * (a - prev_a)
        }

        void calc_pos(double dt) {
            for (int n=0; n<3; n++) pos[n] = pos[n] + vel[n] * dt + 0.5*acc[n] * dt * dt;
        }

        /* for debugging */
        void printAll() {
            printf("====================\natom (PDBID %i) %s on molecule %s (PBDID %i) is %s \n -----> m = %e; eps = %f; sig = %f; C = %f\nforce: %f %f %f; \nacc: %f %f %f; \nold_acc: %f %f %f; \nvel: %f %f %f\n", PDBID, name.c_str(), mol_name.c_str(), mol_PDBID, MF.c_str(), m, eps, sig,C, force[0], force[1], force[2], acc[0], acc[1], acc[2], old_acc[0], old_acc[1], old_acc[2], vel[0], vel[1], vel[2]);
        }
};

Atom::Atom() {}

class Molecule {
	public:
		Molecule();
		vector<Atom> atoms; // vector that holds this molecule's atoms
		int PDBID; // the molecule's PDBID (from input)
		string name; // the molecule name/label (from input), e.g. H2 or MOF
        string MF; // movable/frozen: M or F
        vector<double> force = vector<double>(3); // K / A
        vector<double> torque = vector<double>(3); // t = r x F
        vector<double> com = vector<double>(3); // center of mass for molecule. Using for MD rotations
        //vector<double> prev_com = vector<double>(3);
        vector<double> acc = vector<double>(3); // A / fs^2
        vector<double> old_acc = vector<double>(3);
        vector<double> vel = vector<double>(3); // A / fs
        vector<double> ang_vel = vector<double>(3); // w, angular velocity
        vector<double> ang_acc = vector<double>(3); // alpha, angular acceleration
        vector<double> old_ang_acc = vector<double>(3); // old ang. acceleration
        vector<double> ang_pos = vector<double>(3); // angular rotation distance in rad
        double mass=0.0;
        double inertia=0.0; //moment of inertia

        double get_mass() {
            // ( mass is generated at input in io.cpp )
            return mass;
        }
            
        void calc_inertia() {
            for (int i=0; i<atoms.size(); i++) {
                double rsq = (atoms[i].pos[0] - com[0])*(atoms[i].pos[0] - com[0]) + (atoms[i].pos[1] - com[1])*(atoms[i].pos[1] - com[1]) + (atoms[i].pos[2] - com[2])*(atoms[i].pos[2] - com[2]);
                inertia += atoms[i].m * rsq; // kg * A^2 
            }
            inertia = inertia/1.3806488e-23/1e20*1e30; // to K fs^2
        }

        // angular acceleration
        void calc_ang_acc() {
            for (int n=0; n<3; n++) {
                old_ang_acc[n] = ang_acc[n];
                ang_acc[n] = torque[n] / inertia; // in rad / fs^2
            }
        }

        // linear acceleration
        void calc_acc() {
            for (int n=0; n<3; n++) {
                old_acc[n] = acc[n];
                acc[n] = force[n] * 1.3806488e-33 /mass; // convert from K/A/kg to A/fs^2
            }
        }

        // angular velocity
        void calc_ang_vel(double dt) {
            for (int n=0; n<3; n++) {
                ang_vel[n] = ang_vel[n] + 0.5*(ang_acc[n] * old_ang_acc[n])*dt;
            }
        }

        // linear velocity
        void calc_vel(double dt) {
            for (int n=0; n<3; n++) vel[n] = vel[n] + 0.5*(acc[n] + old_acc[n])*dt; // in A/fs. vel. verlet
        }
   
        // angular position // in rad
        void calc_ang_pos(double dt) {
            for (int n=0; n<3; n++) {
                ang_pos[n] = ang_pos[n] + ang_vel[n] * dt + 0.5*ang_acc[n] * dt * dt;
                if (ang_pos[n] > 0.008) ang_pos[n] = 0.008; // SET THE ROTATION CAP -- rad/fs
            }
        }

        // linear position 
        void calc_pos(double dt) {
            for (int i=0; i<atoms.size(); i++) {
              for (int n=0; n<3; n++) atoms[i].pos[n] = atoms[i].pos[n] + vel[n] * dt + 0.5*acc[n] * dt * dt;
            }
        }        

        void calc_force() {
            // A molecule's force IS the sum of its atoms' forces
            // this is external force on the molecule.
            force[0]=0; force[1]=0; force[2]=0;
            for (int i=0; i<atoms.size(); i++) {
                for (int n=0; n<3; n++) force[n] += atoms[i].force[n]; // in K/A 
            }
        }

        void calc_center_of_mass() {
            // assigns the current center of mass of the molecule based on positions of atoms
            double x_mass_sum=0.0; double y_mass_sum=0.0; double z_mass_sum=0.0; double mass_sum=0.0;

            for (int i=0; i<atoms.size(); i++) {
                double atom_mass = atoms[i].m;
                mass_sum += atom_mass;

                x_mass_sum += atoms[i].pos[0]*atom_mass;
                y_mass_sum += atoms[i].pos[1]*atom_mass;
                z_mass_sum += atoms[i].pos[2]*atom_mass;
            }

            double comx = x_mass_sum/mass_sum;
            double comy = y_mass_sum/mass_sum;
            double comz = z_mass_sum/mass_sum;
            
            com[0] = comx; com[1] = comy; com[2] = comz;
        }

        void calc_torque() {
            for (int n=0; n<3; n++) torque[n] = 0.0;
            // torque is the cross product rxF NOT Fxr, the sum of all atoms relative to molecule's com.
            for (int i=0; i<atoms.size(); i++) {
                atoms[i].torque[0] = (atoms[i].pos[1]-com[1]) * atoms[i].force[2] - (atoms[i].pos[2]-com[2]) * atoms[i].force[1];
                atoms[i].torque[1] = (atoms[i].pos[2]-com[2]) * atoms[i].force[0] - (atoms[i].pos[0]-com[0]) * atoms[i].force[2];
                atoms[i].torque[2] = (atoms[i].pos[0]-com[0]) * atoms[i].force[1] - (atoms[i].pos[1]-com[1]) * atoms[i].force[0];
                // molecular torque = sum of atomic torques
                for (int n=0; n<3; n++) torque[n] += atoms[i].torque[n]; // in K            
            } // end atomic loop
        } // end calc_torque()

        // for debugging
        void printAll() {
            printf("====================\nmolecule PDBID=%i :: mass: %e; inertia: %e; \nforce: %f %f %f; \nacc: %f %f %f; \nold_acc: %f %f %f; \nvel: %f %f %f; \ncom: %f %f %f; \ntorque: %f %f %f \nang_acc: %f %f %f \nold_ang_acc: %f %f %f \nang_vel: %f %f %f; \nang_pos: %f %f %f (in degrees) \n",
            PDBID,mass,inertia,
            force[0], force[1], force[2], acc[0], acc[1], acc[2], 
            old_acc[0], old_acc[1], old_acc[2], vel[0], vel[1], vel[2], com[0], com[1], com[2],
            torque[0], torque[1], torque[2], ang_acc[0], ang_acc[1], ang_acc[2], 
            old_ang_acc[0], old_ang_acc[1], old_ang_acc[2], ang_vel[0], ang_vel[1], ang_vel[2],
            ang_pos[0]*180.0/M_PI, ang_pos[1]*180.0/M_PI, ang_pos[2]*180.0/M_PI); //com[0], com[1], com[2]);
        
        }
};

Molecule::Molecule() {}

Constants::Constants() {
	e = 2.71828183; // ya boi Euler
	kb = 1.3806488e-23; // Boltzmann's in J/K
	kbk = 0.0019872041; // Boltzmann's in kcal/(mol K)
    fs = 1.0e-15; // fs -> second
	cC = 1.60217662e-19; //  e -> coulombs 
	keSI = 8.9875517873681764e9; // ke, Coulomb's constant, Nm^2/C^2 or Jm/C^2.     
	ke = keSI/kb*1e10*cC*cC; // ke in KA / e^2
    eV = 6.242e18; // 1J = eV electron volts
	cM = 1.660578e-27; // kg / particle from g/mol
	cA = 1.0e-10; // 1 angstroem = cA meters
	cJ = 6.94786e-21; // 1 kcal/mol = cJ Joules
	NA = 6.022140857e23; //  particles per mol
	cV = 10.0e-30; // alpha * cV = alpha in m^3
	R = 8.3144598; // J / mol K
    mpmc2uff = pow(2.0,(1.0/6.0)); // mpmc sig * mpmc2uff = UFF sig
    uff2mpmc = 1.0/mpmc2uff; // UFF sig * uff2mpmc = mpmc sig (the RIGHT sigma for LJ)
    ATM2REDUCED = 0.0073389366; // atm -> K/A^3
    kg2em = 9.10938291e-31; // kg -> electron mass in kg
    E2REDUCED = 408.7816; // e -> sqrt(K*A)
    TORQUE2REDUCED = kb * 1e-30 * 1e20; // K -> kg A^2 / fs^2
    FORCE2REDUCED = kb * 1e-30 * 1e20; // K/A -> kg A/fs^2
    DEBYE2SKA = 85.10597636; // debye to ? MPMC reduced

	// MASS VALUES g/mol -> kg/particle
	masses["HB"] = 2.016*cM; // buch model	h2
	masses["H2G"] = 0.0*cM;
    masses["H2E"] = 1.008*cM;
    masses["H2N"] = 0.0*cM;
    masses["HW"] = 1.008*cM; // H in water ( my model)
    masses["HT"] = 1.008*cM; // H in TIP3P
    masses["H"] = 1.008*cM;
	masses["He"] = 4.002602*cM;
	masses["Li"] = 6.941*cM;
	masses["Be"] = 9.012182*cM;
	masses["B"] = 10.811*cM;
	masses["C"] = 12.011*cM;
	masses["N"] = 14.007*cM;
	masses["O"] = 15.9998*cM;
	masses["O2"] = 32.0*cM; // my O2 model
	masses["OW"] = 15.9998*cM; // O in water (my model)
    masses["OT"] = 15.9998*cM; // O in TIP3P
    masses["F"] = 18.998*cM;
	masses["Ne"] = 20.1797*cM;
	masses["Na"] = 22.98976928*cM;
//mg
//al
	masses["Si"] = 28.085*cM;
	masses["P"] = 30.973*cM;
	masses["S"] = 32.06*cM;
	masses["Cl"] = 35.45*cM;
	masses["Ar"] = 39.948*cM;
// per4
	masses["Fe"] = 55.845*cM;
	masses["Cu"] = 63.546*cM;
    masses["Zn"] = 65.39*cM;
// per4
	masses["Kr"] = 83.798*cM;
	masses["Ru"] = 101.07*cM;
//per5
	masses["Xe"] = 131.293*cM;
//per6
	masses["Rn"] = 222.0176*cM;


	// LJ SIGMA VALUES (A)
	sigs["HB"] = 2.96; // buch model h2 (from mpmc sig)
	sigs["H2G"] = 3.2293; // mpmc -> meters
    sigs["H2E"] = 0.0;
    sigs["H2N"] = 2.3406;
    sigs["H"] = 2.886*uff2mpmc; // Rappe=2.886; I changed to 0.5 for MD water model; 0.3 for MOF5.
	sigs["HW"] = 0.5; // H in water, my model; old (sticky MD) 0.4
    sigs["HT"] = 0.0; // H in TIP3P
    sigs["He"] = 2.362*uff2mpmc;
	sigs["Li"] = 2.451*uff2mpmc;
	sigs["Be"] = 2.745*uff2mpmc;
	sigs["B"] = 4.083*uff2mpmc;
	sigs["C"] = 3.851*uff2mpmc; // Rappe = 3.851; I changed to 1.3 for MD MOF5
	sigs["N"] = 3.66*uff2mpmc;
	sigs["O"] = 3.5*uff2mpmc; // Rappe=3.5; I changed to 1.5 for MD water model (SPC = 3.166); 1.3 for MOF5
	sigs["O2"] = 4.0*uff2mpmc; // my O2 model // roughly accurate for densities -> 10atm
	sigs["OW"] = 1.5; // O in water (my model) -- old (sticky MD) 1.4
    sigs["OT"] = 3.15061; // O in TIP3P
    sigs["F"] = 3.364*uff2mpmc;
	sigs["Ne"] = 3.243*uff2mpmc;
    sigs["Cu"] = 3.495*uff2mpmc;
	sigs["Na"] = 2.983*uff2mpmc;
// mg al ...
	sigs["Cl"] = 3.947*uff2mpmc;
	sigs["Ar"] = 3.868*uff2mpmc;
// k ca ..
	sigs["Fe"] = 2.912*uff2mpmc;
// co ni..
	sigs["Zn"] = 2.763*uff2mpmc;
// ga ge ..
	sigs["Kr"] = 4.141*uff2mpmc;
	sigs["Ru"] = 2.963*uff2mpmc;
//per5
	sigs["Xe"] = 4.404*uff2mpmc;
//per6	
	sigs["Rn"] = 4.765*uff2mpmc;


	// LJ EPSILON VALUES (kcal/mol) -> K
	eps["HB"] = 0.06796; // buch model h2
    eps["H2G"] = 8.8516; // bss model h2
    eps["H2E"] = 0.0; // bss
    eps["H2N"] = 4.0659; // bss
	eps["H"] = 0.044/kbk; // rappe verbatim
	eps["HW"] = 0.044/kbk; // H in water (my model)
    eps["HT"] = 0.0; // H in TIP3P
    eps["He"] = 0.056/kbk;
	eps["Li"] = 0.025/kbk;
	eps["Be"] = 0.085/kbk;
	eps["B"] = 0.18/kbk;
	eps["C"] = 0.105/kbk;
	eps["N"] = 0.069/kbk;
	eps["O"] = 0.06/kbk;
	eps["O2"] = 0.06/kbk; // my O2 model
    eps["OW"] = 0.06/kbk; // O in water (my model)
    eps["OT"] = 0.6364/0.0083144621; // O in TIP3P 
	eps["F"] = 0.05/kbk;
	eps["Ne"] = 0.042/kbk;
    eps["Na"] = 0.03/kbk;
//mg al si ...
	eps["Cl"] = 0.227/kbk;
	eps["Ar"] = 0.185/kbk;
// k ca sc ..
	eps["Fe"] = 0.013/kbk;
// co ni..
    eps["Cu"] = 0.005/kbk;
	eps["Zn"] = 0.124/kbk;
// ga ge ..
	eps["Kr"] = 0.220/kbk;
//per5
	eps["Ru"] = 0.056/kbk;
	eps["Xe"] = 0.332/kbk;
//per6
	eps["Rn"] = 0.248/kbk;

	// POLARIZABILITIES  // in A^3 
	polars["H"] = 0.41380;//*cV/ke;
	polars["HW"] = 0.41380;//*cV/ke; // H in water (my model)
    polars["B"] = 0.6634;//*cV/ke;
	polars["C"] = 1.2866;//*cV/ke;
	polars["N"] = 0.97157;//*cV/ke;
	polars["O"] = 0.852;//*cV/ke;
    polars["OW"] = 0.852;//*cV/ke; // O in water (my model)
	polars["Na"] = 24.11;//*cV/ke; // from paper https://www.researchgate.net/publication/45896756_Absolute_and_ratio_measurements_of_the_polarizability_of_Na_K_and_Rb_with_an_atom_interferometer
	polars["P"] = 3.35;//*cV/ke;
	polars["Cl"] = 2.40028;//*cV/ke;
	polars["Cu"] = 2.19630;//*cV/ke;
	polars["Zn"] = 1.98870;//*cV/ke;
	//polars["Br"]
	polars["Ru"] = 1.98870;//*cV/ke; // THIS IS THE PARAM FOR Zn. Borrowed for now.
	polars["Pd"] = 5.25926;//*cV/ke;
	polars["Pt"] = 8.56281;//*cV/ke;


	// He-PHAST2
    // NEEDS FIXIN
	phast2_c6["He"] = 1.407164;
	phast2_c8["He"] = 11.136350;
	phast2_c10["He"] = 107.964;
	phast2_sigs["He"] = 2.18205*cA; // A -> m
	phast2_eps["He"] = 4.49880*kb; // K -> J
	phast2_polar["He"] = 0.20494*cV/ke; // A^3 -> m^3 -> C^2 m^2 / J
}
