#include <stdio.h>
#include <math.h>
#include <string>
#include <map>
#include <vector>
#include <stdint.h>

using namespace std;

enum {
    POTENTIAL_LJ, // the lennard-jones potential
    POTENTIAL_LJES,
    POTENTIAL_LJPOLAR,
    POTENTIAL_LJESPOLAR,
    POTENTIAL_COMMY, // the communist potential 
    POTENTIAL_COMMYES,
    POTENTIAL_COMMYESPOLAR,
    POTENTIAL_TT, // the tang-toennies potential
    POTENTIAL_TTES,
    POTENTIAL_TTESPOLAR
};
enum {
    ENSEMBLE_UVT,
    ENSEMBLE_NVT,
    ENSEMBLE_NPT,
    ENSEMBLE_NVE
};
enum {
    MOVETYPE_DISPLACE,
    MOVETYPE_INSERT,
    MOVETYPE_REMOVE,
    MOVETYPE_VOLUME
};
enum {
    MD_ATOMIC,
    MD_MOLECULAR
};
enum {
    THERMOSTAT_ANDERSEN,
    THERMOSTAT_NOSEHOOVER
};
enum {
    OPTIMIZE_MC,  // Monte Carlo sytle opt.
    OPTIMIZE_SD   // steepest descent opt.
};

/* the below stuff was more-or-less adopted from mpmc code */
typedef struct _histogram {
    int ***grid;
    int x_dim=0, y_dim=0, z_dim=0;
    double origin[3] = {0,0,0};
    double delta[3][3] = { {0,0,0},{0,0,0},{0,0,0} };
    int count[3];
    int n_data_points=0;
    int norm_total=0;
} histogram_t;

class FilePointer {
  public:
    FilePointer();
    FILE *fp_histogram;
};
FilePointer::FilePointer() {}

// grids (for histogram)
class Grid {
    public:
        Grid();

        histogram_t *histogram;
        histogram_t *avg_histogram;

};
Grid::Grid() {}
/* end stuff for histogram */

// Constants is sort-of a misnomer for some things in this class but you get the idea.
class Constants {
	public:
		Constants();
        const double e = 2.71828183; // ya boi Euler
        const double kb = 1.3806488e-23; // Boltzmann's in J/K
        const double kbk = 0.0019872041; // Boltzmann's in kcal/(mol K)
        const double fs = 1.0e-15; // fs -> second
        const double cC = 1.60217662e-19; //  e -> coulombs
        const double keSI = 8.9875517873681764e9; // ke, Coulomb's constant, Nm^2/C^2 or Jm/C^2.
        const double ke = keSI/kb*1e10*cC*cC; // ke in KA / e^2
        const double eV = 6.242e18; // 1J = eV electron volts
        const double cM = 1.660578e-27; // kg / particle from g/mol
        const double cA = 1.0e-10; // 1 angstroem = cA meters
        const double cJ = 6.94786e-21; // 1 kcal/mol = cJ Joules
        const double NA = 6.022140857e23; //  particles per mol
        const double kek = keSI * cC*cC * 0.00239006 * 1e10 * NA; // ke in kcalA / mol e^2
        const double cV = 10.0e-30; // alpha * cV = alpha in m^3
        const double R = 8.3144598; // J / mol K
        const double mpmc2uff = pow(2.0,(1.0/6.0)); // mpmc sig * mpmc2uff = UFF sig
        const double uff2mpmc = 1.0/mpmc2uff; // UFF sig * uff2mpmc = mpmc sig (the RIGHT sigma for LJ)
        const double ATM2REDUCED = 0.0073389366; // atm -> K/A^3
        const double kg2em = 9.10938291e-31; // kg -> electron mass in kg
        const double E2REDUCED = 408.7816; // e -> sqrt(K*A)
        const double TORQUE2REDUCED = kb * 1e-30 * 1e20; // K -> kg A^2 / fs^2
        const double FORCE2REDUCED = kb * 1e-30 * 1e20; // K/A -> kg A/fs^2
        const double DEBYE2SKA = 85.10597636; // debye to ? MPMC reduced
        const double JL2ATM = 0.00986923297; // J/L to atm
        const double A32L = 1e-27; // A^3 to liters.
        const double K2KJMOL = kb*NA/1000; // K -> kJ/mol
        const double HBARC = 22898848.135746032; // in K*A
        const double vand2mpmc = 0.14818471127642288; // au^3 * this = A^3
        const double eA2D = 1./0.20819434; // eA * this = Debye; or eA^2 * this = Debye*A, etc.
        const double au2D = 2.5411968777103207; // au * this = Debye
        const double K2Eh = 1.0/3.1577513e5;//   0.0000032; // K * this = Hartrees
        const double bohr = 0.529177; // Bohr radius a_B in Angstroms


		string jobname="default_jobname";
        string mode; // "mc" or "md"
        int_fast8_t checkpoints_option=0; // enables checkpoints for debuggin
        int_fast8_t ensemble;
        string ensemble_str;
        //int_fast8_t movetype;
        string atom_file = "<none>"; // input atoms .pdb file
        string output_traj="traj.xyz"; // system trajectory in xyz
        string output_traj_pdb="traj.pdb"; // system trajectory in pdb
        string output_traj_movers_pdb="traj_movers.pdb"; // system traj (only movers) in pdb
        string restart_pdb="restart.pdb"; // a file to re-kick an unfinished run
        string thermo_output="thermo.dat"; // a file for all thermodynamic info
        string output_histogram="histogram.dx"; // histogram information, viewable in VMD
        string dipole_output="dipoles.dat"; // only used when polarization is on.
        string restart_mov_pdb="restart_movables.pdb"; // a restart file with only movable molecules to save i/o
        string frozen_pdb="frozen.pdb"; // a pdb of frozen atoms that is made at startup
        int_fast8_t potential_form = POTENTIAL_LJ; // "lj", "ljes", "ljespolar", "phast2" models for potential
        vector<string> sorbate_name; // e.g. h2_bssp, h2_bss, co2*, co2, co2_trappe, c2h2, etc.
        vector<double> sorbate_fugacity; // holds the fugacities for multi-sorbate gases.
        vector<int> sorbate_dof; // holds degrees of freedom for sorbates
        int_fast8_t xyz_traj_option=1; // option for xyz trajectory, default off
        int_fast8_t xyz_traj_movers_option=1; // option for smaller xyz traj (with only movers)
        int_fast8_t pdb_traj_option=0; // option to write PDB trajectory . default on
        int_fast8_t pdb_bigtraj_option=0;// option to write trajectory WITH frozen atoms (though frozen.pdb gets written at startup)
        int_fast8_t dipole_output_option=1; // for dipole output (polar only)
        int_fast8_t com_option=0; // enables computation of center-of-mass and logs in output_traj
        int_fast8_t rotate_option=1; // MC ONLY: can deactivate rotates if wanted.
        int_fast8_t draw_box_option=1; // option to draw the box for visualization in restart.pdb
        int_fast8_t rd_lrc=1; // long range corrections for LJ RD
        int_fast8_t ewald_es=1; // ewald method for electrostatic potential calculation.
        int_fast8_t pdb_long=0; // on would force long coordinate/charge output
        int_fast8_t dist_within_option=0; // a function to calculate atom distances within a certain radius of origin
        string dist_within_target; // the atom to find in above option
        double dist_within_radius; // the radius within which to search from origin
        int_fast8_t histogram_option = 0; // output histogram data which can be plotted in VMD
        int_fast8_t autocenter = 1; // center all atoms about origin automatically. can opt out of it.
        int_fast8_t restart_mode = 0; // option to restart job automatically (by searching for restart.pdb)
        int_fast8_t no_zero_option = 0; // option to disallow zero sorbates in the simulation. default off.
        int_fast8_t simulated_annealing = 0; // sim. ann.
        double sa_target = 0.0; // target temperature for annealing.
        double sa_schedule = 0.9999; // T-change factor for annealing.
        double free_volume=0; // for excess adsorption calculation, A^3. must be user-input
        unsigned int fugacity_single=0; // set for single-sorbate fugacity calculation at startup.
        string fugacity_single_sorbate; // the sorbate molecule to get fugacity of. h2/n2/co2/ch4

        int_fast8_t feynman_hibbs = 0;
        unsigned int fh_order = 4;

        // DEFAULT ELEMENT/SITE PARAMETERS
        map <string,double> masses; // mass database for defaults.
		map <string,double> sigs; // LJ sigmas database for defaults. (mostly UFF). Defined below
		map <string,double> eps; // LJ epsions database for defaults. (mostly UFF). Defined below
        map <string,double>UFF4MOFsigs; // UFF4MOF sigma params. http://pubs.acs.org.ezproxy.lib.usf.edu/doi/pdf/10.1021/acs.jctc.6b00664
		map <string,double> phast2_c6; map <string,double> phast2_c8; map <string,double> phast2_c10; map <string,double> phast2_sigs; map <string,double> phast2_eps; map <string,double> phast2_polar; // phast2 defaults
		map <string,double> polars; // polarizability database for defaults. Mostly vanD. Defined below
        map <string,double> elements; // H:1 ; He:2 ; etc.
        map <string,double> radii; // VDW radii of elements
        map <string,double> UFF_bonds; // UFF equilibrium bond-lengths in A
        map <string,double> UFF_angles; // UFF theta equilibrium angles, in degrees
        map <string,double> UFF_Z; // UFF Effective charge Z
        map <string,double> UFF_torsions; // UFF torsion angle barriers, kcal/mol
        map <string,double> UFF_electroneg; // UFF electronegativity, in eV

        unsigned int lj_uff=0; // 1 would default all atoms to UFF LJ parameters (override the input)
        unsigned int polars_vand=0; // 1 would defaul all atoms to van Duijnen polarizablitiy parameters

        double total_energy=0.0; // for MC NVE, in K, user defined
        double volume; // in A^3
        double temp=0.0; //273.15; // in K
        double prevtemp = 0.0; // previous temp for NVT MD thermostat
        double pres=1.0; // in atm
        double volume_change=0.25; // a factor used in volume change BF, mpmc default 0.25
        double vcp_factor=1.0; // a factor used in probability for volume change. this / num_sorbates is good per Frenkel
        double displace_factor=2.5; // up to +- this number in A
        double insert_factor=0.5; // probability to do insert or delete (instead of displace/rotate) in uvt
// DEPRECATED double rotate_prob=0.5; // prob to rotate instead of displace when displace/rotate is selected
        
        double bias_uptake =0; // to accelerate uVT MC. User-defined. This will be converted to N before running MC
        string bias_uptake_unit="N"; // units for bias. Default N_movables
        int_fast8_t bias_uptake_switcher=0;
        unsigned int scale_charges=0; // option to scale charges
        double scale_charges_factor; // multiply by this to get new charges for system. 


        double rotate_angle_factor=360; // 0 -> this number to rotate if rotate selected
		unsigned int stepsize=1; // obvi
        unsigned int finalstep=-1; // user defined for MC. Will error-out if not given in MC
        unsigned int  mc_corrtime=1000; // default 1k cuz I used that a lot for mpmc research
        int_fast8_t mc_pbc=1; // PBC in monte carlo, default on
        unsigned int currentprotoid=0; // for getting fugacity for the boltzmann factor.
        unsigned int step_offset=0; // a parameter used to change the step output in the output files (e.g. after a restart) 
        unsigned int readinxyz=0; // option to read an XYZ file for input instead of PDB

        // MD STUFF
        unsigned int  md_corrtime=50; // user defined for MD
        int_fast8_t md_pbc=1; // PBC in molecular dynamics. default on.
        int_fast8_t md_rotations=1; // MD only, rotational motion (of molecules) on
        int_fast8_t md_translations=1; // MD only, translation motion on
        int_fast8_t md_manual_init_vel=0; // flag for user-defined initial velocities; default off
        double md_init_vel=0.; //  A / fs. User can set.
        double md_dt=0.1, md_ft=10000; // MD timestep and final time, in fs
        int_fast8_t md_mode = MD_MOLECULAR; // default is to keep molecules rigid (bonded)
		double md_initial_energy_NVE = 0; // NVE first Newtonian energy, for measuring integrator error over time; in K.
        double md_NVE_err = 0; // delta(total E) in kJ/mol
        double md_thermostat_freq = 0.05; // a value used to calculate probability of a heat-bath collision with molecule i; Frenkel uses 0.01 and 0.001 as examples; but no matter what, a boltzmann distribution is generated
        double md_thermostat_probab = md_thermostat_freq * exp(-md_thermostat_freq * md_dt);
        unsigned int md_insert_attempt=20; // uVT MD. Number of timesteps to try insert/delete. Default every 20 steps.
        unsigned int md_external_force = 0; // option for constant external force in MD
        double external_force_vector[3] = {0,0,0}; // Fx,Fy,Fz stored in K/A.
        double lagrange_multiplier = 0; // used for Nose-Hoover NVT thermostat.
        unsigned int thermostat_type = THERMOSTAT_NOSEHOOVER; // thermostat type for NVT temperature fixture.


        map <string,double> sig_override;
        map <string,double> eps_override; // feature for overriding preset LJ params (for developing LJ models). 0.0 are defaults which will be overwritten if option is used. sig=A; eps=K

    	unsigned int total_atoms=0; // actual sites, not "atoms" persay
                unsigned int initial_sorbates=0.0; // for safekeeping to calculate chemical potential in uVT
        double initial_energy=0.0; // "" ""

        // Ewald (for ES)
        double ewald_alpha; // =3.5/cutoff; Really, sqrt(alpha), by Ewald formula.
                            // i have also seen 2.5 / r_c for this quantity (Rapaport, Art of M.D.)
        double ewald_kmax = 7; // suitable for most cases.
        //double** ewald_k; // holds 3D k-space vectors for Ewald summation in Force calc for MD.
        // actually its faster to make on-the-fly
        unsigned int ewald_num_k; // number of Ewald k vectors stored in ewald_k
        unsigned int kspace_option=1; // default on; option to include k-space contribution to FORCES in md

        // Wolf (for polarization)
        //unsigned int polar_iterative=1; // turn iterative on. If off, will just do one iteration of dipole calc and get polar energy
        double polar_wolf_alpha = 0.13;
        double polar_damp = 2.1304;
        double polar_gamma = 1.03;
        unsigned int polar_max_iter = 4;
        double polar_rmin = 0; // minimum polarizable distance between atoms, for ranking.
        double **B_matrix, C_matrix[3][3];
        double **A_matrix;
        double **A_matrix_full; // takes more memory but computations are faster
        unsigned int full_A_matrix_option = 1; // option to enable the full A matrix for polarization
        double polar_precision=0; // dipole precision required during iterative routine
        unsigned int iter_success=0; // flag for polarization iteration failure (importance for acceptance of moves!)
        int_fast8_t polar_rrms =0;
        double dipole_rrms = 0.0;
        int_fast8_t polar_gs_ranked = 1;
        int_fast8_t polar_gs = 0;
        int_fast8_t polar_palmo = 1;
        int_fast8_t polar_pbc = 1; // default periodic polar
        //unsigned int thole_total_atoms = 0;

        int_fast8_t all_pbc=1;
        int_fast8_t auto_reject_option=1; // enables/disables
        double auto_reject_r=0.76; // Angstroms. If r < this value for any relevant pair, MC will auto-reject a move immediately
        int_fast8_t auto_reject=0; // on or off (for an individual step!). will go on if auto_reject_r is triggered
        unsigned int rejects=0; // counter

        unsigned int manual_cutoff=0; // on/off for user-defined pair-interaction cutoff in A.
        double manual_cutoff_val=0; // in A
    
        unsigned int cuda=0; // CUDA OPTION FOR GPU CALCULATIONS (MD only)
        unsigned int cuda_block_size = 256; // this was the best of a test of 32,64,128,256 on a project from class I took. Can play with this to see how it changes perf.

        unsigned int crystalbuild=0; // option to dynamically build a crystal box to a supercell
        unsigned int crystalbuild_x=1, crystalbuild_y=1, crystalbuild_z = 1; // duplication # in each dim. 
        unsigned int crystalbuild_includemovers=0; // option to include movable molecules in the crystal builder. default off.
        unsigned int charge_sum_check = 1; // option to check the total system charge before simulation. Default on.

        unsigned int fragmaker=0; // option to create fragments at startup.
        unsigned int numfrags=0; // number of fragments to create in fragmentMaker function
        vector<int> fragsize = {250}; // num. of atoms in a frag, default
        double frag_bondlength = 2.1; // Angstroms, default.

        unsigned int write_lammps = 0; // option to write out LAMMPS input files 

        /* SINGLE-POINT OPTIONS */
        unsigned int user_charge=0; // molecular charge for single pounsigned int calc's

        /* OPTIMIZATION OPTIONS */
        double bondlength = 1.9; // for bonding MD, default
        double opt_error = 0.0001; // error, in kcal/mol for convergence
        unsigned int opt_step_limit = 10000; // limit for convergence steps
        unsigned int opt_mode = OPTIMIZE_SD; // steepest descent style optimization as default
        unsigned int opt_bonds = 1; // option to include bonding energy/gradients
        unsigned int opt_angles = 1; // "" angles ""
        unsigned int opt_dihedrals = 1; // "" dihedrals ""
        unsigned int opt_LJ = 1; // "" non-bond LJ ""
        unsigned int opt_ES = 1; // "" non-bond electrostatics

        struct UniqueBond {
            unsigned int mol, atom1, atom2; double  value; // indices of molecule + atoms of bond
            double BO, rij, kij, alpha, Dij;
        };
        vector<UniqueBond> uniqueBonds; // holds all unique bonds
        struct UniqueAngle {
            unsigned int mol, atom1, atom2, atom3; double value;
            double rij, rjk, theta_ijk, C2, C1, C0; // rik, K_ijk;
            double t1,t2,t3;
        };
        vector<UniqueAngle> uniqueAngles; // holds all unique 3-atom angles
        struct UniqueDihedral {
            unsigned int mol, atom1, atom2, atom3, atom4; double value;
            double phi_ijkl, vjk, n;
        };
        vector<UniqueDihedral> uniqueDihedrals; // holds all unique 4-atom dihedrals.
        struct UniqueLJNonBond {
            unsigned int mol, atom1, atom2;         
            double sig, eps;
        };
        vector<UniqueLJNonBond> uniqueLJNonBonds; // holds all unique non-bonded pair interactions (no cutoff)
        struct UniqueChargeNonBond {
            unsigned int mol, atom1, atom2;
            double chargeprod;
        };  
        vector <UniqueChargeNonBond> uniqueChargeNonBonds;

};

class Pbc {
    public:
        Pbc();

        double x_length, y_length, z_length; // for moving molecules back into box on other side, computing in calcBoxVertices
        double x_min,x_max,y_min,y_max,z_min,z_max;
        double basis[3][3];
        double reciprocal_basis[3][3];
		double cutoff=0.; // to truncate long-range interactions
        double mincutoff=2.55; // to neglect close-range interactions (in Single-Pounsigned int calc)
        double volume, inverse_volume, old_volume;
        double a=0, b=0, c=0, alpha=0, beta=0, gamma=0;
        double box_vertices[8][3];
        double A[6], B[6], C[6], D[6]; // these are coefficients for plane equations for PBC
            /* structure of box_points
                0 : -x, -y, -z
                1 : -x, -y, +z
                2 : -x, +y, -z
                3 : -x, +y, +z
                4 : +x, -y, -z
                5 : +x, -y, +z
                6 : +x, +y, -z
                7 : +x, +y, +z
            */

        double maxx=0,maxy=0,maxz=0,minx=0,miny=0,minz=0;
        double lengthx=0, lengthy=0, lengthz=0;

        void printBasis() {
            printf("\n:: --- Box (basis) information --- ::\n"); 
            printf(":: basis1 %.5f %.5f %.5f\n", basis[0][0], basis[0][1], basis[0][2]);
            printf(":: basis2 %.5f %.5f %.5f\n", basis[1][0], basis[1][1], basis[1][2]);
            printf(":: basis3 %.5f %.5f %.5f\n", basis[2][0], basis[2][1], basis[2][2]);
            printf(":: Basis vectors: { a = %9.5f; b = %9.5f; c = %9.5f }\n", a, b,c);
            printf(":: Basis angles:  { α = %9.5f; β = %9.5f; γ = %9.5f }\n", alpha,beta,gamma);
            printf(":: Box vertices ::\n");
            for (unsigned int n=0; n<8; n++)
                printf("   -> %i : %9.5f %9.5f %9.5f\n", n, box_vertices[n][0], box_vertices[n][1], box_vertices[n][2]);
            printf(":: PBC Cutoff = %.5f\n", cutoff);
            for (unsigned int n=0; n<6; n++) {
                printf(":: Plane %i equation :: %.5fx + %.5fy + %.5fz + %.5f = 0\n",
                    n, A[n], B[n], C[n], D[n]);
            }
            printf(":: x_length = %.5f; y_length = %.5f; z_length = %.5f\n", x_length, y_length, z_length);
            if (alpha == 90 && beta == 90 && gamma == 90) { // these don't get calculated or used for weird (non 90-90-90) unit cells
              printf(":: x_max = %.5f; y_max = %.5f; z_max = %.5f\n", x_max, y_max, z_max);
              printf(":: x_min = %.5f; y_min = %.5f; z_min = %.5f\n", x_min, y_min, z_min);
            }
            printf(":: --- End box information --- ::\n\n");
        }

        void calcPlane(unsigned int p1index, unsigned int p2index, unsigned int p3index, unsigned int planeIndex) { // 3 points define a plane.
            double vector1[3], vector2[3];

            // 1) get 3 points (indexes for box vertices provided in arguments)
            // 2) make 2 planar vectors AB, AC
            for (unsigned int n=0; n<3; n++) {
                vector1[n] = box_vertices[p2index][n] - box_vertices[p1index][n];
                vector2[n] = box_vertices[p3index][n] - box_vertices[p1index][n];
            }

            // 3) calculate normal vector to the plane
            double* normal = crossprod(vector1, vector2);

            // 4) plane equation is thus defined
            A[planeIndex] = normal[0];
            B[planeIndex] = normal[1];
            C[planeIndex] = normal[2];
            D[planeIndex] = -dddotprod(normal,box_vertices[p1index]);
            // Thus the plane equation is Ax + By + Cz + D = 0
        }

        void calcPlanes() {
            /* i drew a cube :-)
                                    The A[6],B[6],C[6],D[6] arrays will be used to make plane equations
             2 /------------/ 6     p0 :   0123 plane (-x)
              /|   p3      /|       p1 :   4567 plane (+x)
             / |          / |       p2 :   0145 plane (-y)
          3 |------------|7 |   p1  p3 :   2367 plane (+y)
            |  |         |  |       p4 :   0246 plane (-z) (not shown)
      p0    |  |---------|--| 4     p5 :   1357 plane (+z) (not shown)
            | / 0  ___   |  /
            |/    /p2/   | /        The vertices are defined in box_vertices[8][3].
          1 |____________|/ 5       3 points define a plane, so I'll use the first 3 for the above planes


            */

            // 3 points and plane index
            // quite sure that the plane NEVER differs if a different set of 3 points is used (even for triclinic cells)
            // so it's safe to just pick the first 3 points of the plane.
            calcPlane(0,1,2,0);
            calcPlane(4,5,6,1);
            calcPlane(0,1,4,2);
            calcPlane(2,3,6,3);
            calcPlane(0,2,4,4);
            calcPlane(1,3,5,5);
        }

        void calcVolume() {
            double newvolume;
            newvolume =  basis[0][0]*(basis[1][1]*basis[2][2] - basis[1][2]*basis[2][1]);
            newvolume += basis[0][1]*(basis[1][2]*basis[2][0] - basis[1][0]*basis[2][2]);
            newvolume += basis[0][2]*(basis[1][0]*basis[2][1] - basis[2][1]*basis[2][0]);
            volume = newvolume;
            inverse_volume = 1.0/volume;
        }

        void calcCutoff() {
            if (cutoff != 0.) return; // mpmc only changes the cutoff if it's nonzero
            double MAXVALUE = 1e40; unsigned int MAX_VECT_COEF = 5;
			unsigned int i, j, k, p;
			double curr_mag;
			double short_mag = MAXVALUE;
			double curr_vec[3];
			if ( volume <= 0 ) cutoff = MAXVALUE;

            // smallest vector problem
			for ( i=-MAX_VECT_COEF; i<=MAX_VECT_COEF; i++ ) {
				for ( j=-MAX_VECT_COEF; j<=MAX_VECT_COEF; j++ ) {
				    for ( k=-MAX_VECT_COEF; k<=MAX_VECT_COEF; k++ ) {
				        if ( i == 0 && j == 0 && k == 0 ) continue;
				        for ( p = 0; p < 3; p++ )
				            curr_vec[p] = i*basis[0][p] + j*basis[1][p] + k*basis[2][p];
				            curr_mag = sqrt(
                                (curr_vec[0] * curr_vec[0]) +
                                (curr_vec[1] * curr_vec[1]) +
                                (curr_vec[2] * curr_vec[2])
                            );
				        if ( curr_mag < short_mag ) short_mag = curr_mag;
				    }
				}
			}
			cutoff = 0.5*short_mag;
        }

        void calcRecip() {
			// assumes volume and inverse_volume are already calc'd
            reciprocal_basis[0][0] = inverse_volume*(basis[1][1]*basis[2][2] - basis[1][2]*basis[2][1]);
			reciprocal_basis[0][1] = inverse_volume*(basis[0][2]*basis[2][1] - basis[0][1]*basis[2][2]);
			reciprocal_basis[0][2] = inverse_volume*(basis[0][1]*basis[1][2] - basis[0][2]*basis[1][1]);

			reciprocal_basis[1][0] = inverse_volume*(basis[1][2]*basis[2][0] - basis[1][0]*basis[2][2]);
			reciprocal_basis[1][1] = inverse_volume*(basis[0][0]*basis[2][2] - basis[0][2]*basis[2][0]);
			reciprocal_basis[1][2] = inverse_volume*(basis[0][2]*basis[1][0] - basis[0][0]*basis[1][2]);

			reciprocal_basis[2][0] = inverse_volume*(basis[1][0]*basis[2][1] - basis[1][1]*basis[2][0]);
			reciprocal_basis[2][1] = inverse_volume*(basis[0][1]*basis[2][0] - basis[0][0]*basis[2][1]);
			reciprocal_basis[2][2] = inverse_volume*(basis[0][0]*basis[1][1] - basis[0][1]*basis[1][0]);
        }


        void calcCarBasis() {
            // this function is called if normal basis is supplied by user
            a = sqrt(dddotprod(basis[0], basis[0]));
            b = sqrt(dddotprod(basis[1], basis[1]));
            c = sqrt(dddotprod(basis[2], basis[2]));
            alpha = 180.0/M_PI*acos( dddotprod(basis[1],basis[2]) / sqrt( dddotprod(basis[1], basis[1]) * dddotprod(basis[2], basis[2]) ));
            beta = 180.0/M_PI*acos( dddotprod(basis[2],basis[0]) / sqrt( dddotprod(basis[0], basis[0]) * dddotprod(basis[2], basis[2]) ) );
            gamma = 180.0/M_PI*acos( dddotprod(basis[0],basis[1]) / sqrt( dddotprod(basis[1], basis[1]) * dddotprod(basis[0], basis[0]) ) );
        }

        void calcNormalBasis() {
                double b0[3] = {0,0,0};
                double b1[3] = {0,0,0};
                double b2[3] = {0,0,0};

                b0[0] = a;
                b0[1] = b*cos(M_PI/180.0 * gamma);
                b0[2] = c*cos(M_PI/180.0 * beta);

                b1[0] = 0;
                b1[1] = b*sin(M_PI/180.0 * gamma);
                b1[2] = ( (0*0 + b*0 + 0*c) - (b0[1]*b0[2]) )/b1[1];

                b2[0] = 0;
                b2[1] = 0;
                b2[2] = sqrt( c*c - b0[2]*b0[2] - b1[2]*b1[2] );

                // I'm transposing it manually
                basis[0][0] = b0[0];
                basis[0][1] = b1[0];
                basis[0][2] = b2[0];

                basis[1][0] = b0[1];
                basis[1][1] = b1[1];
                basis[1][2] = b2[1];

                basis[2][0] = b0[2];
                basis[2][1] = b1[2];
                basis[2][2] = b2[2];

        }

        void calcBoxVertices() {
			// calculates the 3D points that encompass the crystalline simulation box.
		    unsigned int i,j,k,p,q,count=0;
		    //unsigned int box_labels[2][2][2];
		    double box_occupancy[3];
		    double box_pos[3];

		    // draw the box points
		    for(i = 0; i < 2; i++) {
		        for(j = 0; j < 2; j++) {
		            for(k = 0; k < 2; k++) {

		                /* box coords */
		                box_occupancy[0] = ((double)i) - 0.5;
		                box_occupancy[1] = ((double)j) - 0.5;
		                box_occupancy[2] = ((double)k) - 0.5;

		                for(p = 0; p < 3; p++) {
		                    for(q = 0, box_pos[p] = 0; q < 3; q++) {
		                        box_pos[p] += basis[q][p]*box_occupancy[q];
                            }
                        }

                        for (unsigned int n=0; n<3; n++)
                            box_vertices[count][n] = box_pos[n]; // box_points[0 -> 7] will be defined.

                        count++;
		            } // for k
		        } // for j
		    } // for i

            x_length = box_vertices[5][0] - box_vertices[1][0]; // box lengths based on the front-left corner of box
            y_length = box_vertices[3][1] - box_vertices[1][1]; // kind of crude way but I think it's foolproof..
            z_length = box_vertices[1][2] - box_vertices[0][2];
        }

        void calcMaxMin() {
            // re-initialize maximums and minimums
            maxx=maxy=maxz=minx=miny=minz=0;
            for (unsigned int n=4; n<=7; n++) if (box_vertices[n][0] > maxx) maxx = box_vertices[n][0];
            for (unsigned int n=0; n<=3; n++) if (box_vertices[n][0] < minx) minx = box_vertices[n][0];

            for (unsigned int n=2; n<=3; n++) if (box_vertices[n][1] > maxy) maxy=box_vertices[n][1];
            for (unsigned int n=6; n<=7; n++) if (box_vertices[n][1] > maxy) maxy=box_vertices[n][1];
            for (unsigned int n=0;n<=1; n++) if (box_vertices[n][1] < miny) miny = box_vertices[n][1];
            for (unsigned int n=4;n<=5; n++) if (box_vertices[n][1] < miny) miny = box_vertices[n][1];

            for (unsigned int n=1; n<=7; n+=2) if (box_vertices[n][2] > maxz) maxz = box_vertices[n][2];
            for (unsigned int n=0; n<=6; n+=2) if (box_vertices[n][2] < minz) minz = box_vertices[n][2];

            lengthx = maxx-minx;
            lengthy = maxy-miny;
            lengthz = maxz-minz;

        }

};

Pbc::Pbc() {}

class Stats {
	public:
		Stats();

        unsigned int MCstep=0, MCcorrtime_iter; // keeps track of steps and coortimes for averages.
        double MDtime=0; // the time in fs of MD simulation
        bool MCmoveAccepted;
        double MCeffRsq; // for calculating Monte Carlo efficiency, roughly, based on successful displaces

        int_fast8_t radial_dist = 0; // default is no radial distribution
        string radial_file = "radial_distribution.dat"; // default filename for output.
        double radial_bin_size = 0.1; // bin counts will be considered for this range in A
        double radial_max_dist = 10.0; // maximum r to consider in rad. dist.
        vector<vector<long unsigned int>> radial_bins; // holds the counters for each g(r)
        vector<string> radial_centroid, radial_counterpart; // the two atoms to get distance between, user def.

		double insert_bf_sum = 0; // insert boltzmanns added up in boltzmann.cpp
		double remove_bf_sum = 0; // ...
		double displace_bf_sum = 0;
		double volume_change_bf_sum = 0;

		unsigned int insert_accepts = 0; unsigned int insert_attempts=0;// Counters for successful moves. uVT only
		unsigned int remove_accepts = 0; unsigned int remove_attempts=0;// uVT only
		unsigned int displace_accepts = 0; unsigned int displace_attempts=0;
		unsigned int volume_change_accepts = 0; unsigned int volume_attempts=0;// NPT only
		unsigned int total_accepts=0; unsigned int total_attempts=0;

        double ar_tot=0, ar_ins=0, ar_rem=0, ar_dis=0, ar_vol=0; // BF acceptance ratios
        double ins_perc=0, rem_perc=0, dis_perc=0, vol_perc=0; // BF percentage of moves

        double bf_avg, ibf_avg, rbf_avg, dbf_avg, vbf_avg;

        unsigned int count_movables = 0; // this is the SORBATE MOLECULES, e.g. 2H2 means 2, not 4
        vector<int> movids; // holds the internal IDs of movable molecules (in system.molecules)
        unsigned int count_frozens = 0; // frozen ATOMS, not molecules (which is normally just 1)
        unsigned int count_frozen_molecules=0; // frozen MOLECULES; normally 1

        struct obs_t { // observable types
            string name;
            double counter=0.0;
            double value=0;
            double average=0;
            double sd=0;

            void calcNewStats() { // gets new stats based on new val
                // the assumption here is that the new val was already calculated and provided (in "value")
                double x = value;
                double prevavg = average;
                double prevsd = sd; //printf("counter %f\n",obs.counter);
                counter = counter+1.0; //printf("counter %f\n",obs.counter);
                average = ((counter-1.0)*average + x)/counter;
                double operand =  prevsd*prevsd + prevavg*prevavg - average*average +((x*x - prevsd*prevsd - prevavg*prevavg)/counter);
                (operand > 0) ? sd = sqrt( operand ) : sd = 0;
            }

        // these are the instances of the "observables type" obs_t.
        // i.e. all the quantities of interest which can be averaged and dealt with.
        } Nsq,NU,qst,qst_nvt, 
          rd,es,polar, // all energies in K unless otherwise stated
          potential,potential_sq,
          kinetic, kinetic_sq,
          totalE, totalE_sq,
          Klin, Krot, Ulin, Urot, // linear (translational) and rotational energies
          lj_lrc, lj_self_lrc, lj, // lj_self_lrc is really just "self" energy for LJ potential.
          es_self, es_real, es_recip,
          tt, tt_lrc, tt_self,
          Ustretch, Uangles, Udihedrals, UintraLJ, UintraES, Ubonded_tot, // bonded energies in kcal/mol
          chempot,totalmass, frozenmass, 
          pressure,temperature, fdotrsum, 
          dist_within, diffusion, volume, z, // z = PV/nRT
          Q, polar_iterations, heat_capacity; // Q = partition func.

        unsigned int max_sorbs=10;
        vector<obs_t> msd = vector<obs_t>(max_sorbs); 
        vector<obs_t> wtp = vector<obs_t>(max_sorbs);
        vector<obs_t> wtpME = vector<obs_t>(max_sorbs);
        vector<obs_t> Nmov = vector<obs_t>(max_sorbs);
        vector<obs_t> movablemass = vector<obs_t>(max_sorbs);
        vector<obs_t> density = vector<obs_t>(max_sorbs); 
        vector<obs_t> selectivity = vector<obs_t>(max_sorbs);
        vector<obs_t> excess = vector<obs_t>(max_sorbs);

};

Stats::Stats() {}

// contains information about an individual pair marked by i,j,k,l
class Pair {
    public:
        Pair();
        double r; double prev_r;
        double d[3]; double prev_d[3];
        unsigned int recalculate = 1; // a flag to queue recalculation of energies
        double eps = 0; // LJ param, mixing rule'd
        double sig = 0; // LJ param, mixing rule'd
        double rd_energy=0;
            double lj=0;
            double lj_lrc=0;
        double es_energy=0;
        double pol_energy=0;
        double fdotr=0; // F.r dot prod. Needed to get emergent pressure in MD NVT

};
Pair::Pair() {}


// stores variables to return to, if move rejected, or for checkpointing.
class Last {
    public:
        Last();
        double Nsq,NU,qst,qst_nvt,rd,es,polar,potential,potential_sq,volume,z,
            lj_lrc,lj_self_lrc,lj,es_self,es_real,es_recip,chempot,totalmass,
            frozenmass,pressure,temperature, fdotrsum, dist_within, csp, diffusion;

        unsigned int total_atoms, thole_total_atoms;

        unsigned int max_sorbs = 10;
        vector<double> wtp = vector<double>(max_sorbs);
        vector<double> wtpME = vector<double>(max_sorbs);
        vector<double> Nmov = vector<double>(max_sorbs);
        vector<double> movablemass = vector<double>(max_sorbs);
        vector<double> density = vector<double>(max_sorbs); 
        vector<double> selectivity = vector<double>(max_sorbs);
        vector<double> excess = vector<double>(max_sorbs);
};

Last::Last() {}

class Atom {
	public:
		Atom();
        //Atom(const Atom& rhs) { /* for cloning */ }
        //Atom& operator=(const Atom& rhs) {};
        string name; // element or label, e.g. H or H2G
        string mol_name; // molecule name that the atom belongs to
        int_fast8_t frozen; // movable/frozen (0 or 1)
		unsigned int mol_PDBID; // the molecule's PDBID that this atom belongs to
		double m=0.0; // mass, kg. This is the only one I'm keeping SI as of now.
        double eps=0.0; // LJ param in K
        double sig=0.0; // LJ param in A -- the real sigma, not r_m (as in UFF)
        double polar=0.0; // polarizability in A^3
        double C=0.0; // charge in e
        double V=0.0; // potential energy in K
        //double K=0.0; // kinetic energy in K
        //double E=0.0; // total energy in K
		// Tang-Toennies params.
        double c6=0,c8=0,c10=0;
        
        unsigned int PDBID; // the atom's PDBID (from input)
        double rank_metric;  // for polarization sorting

        double pos[3] = {0,0,0};
		//double prevpos[3] = {0,0,0};
        double force[3] = {0,0,0};
        double vel[3] = {0,0,0};
        double acc[3] = {0,0,0};
        double old_acc[3] = {0,0,0};
        double torque[3] = {0,0,0};
        double dip[3] = {0,0,0};
        double newdip[3] = {0,0,0};
        double olddip[3] = {0,0,0};
        double efield[3] = {0,0,0};
        double efield_self[3] = {0,0,0};
        double efield_induced[3] = {0,0,0};
        double efield_induced_change[3] = {0,0,0};
        double energy_grad[3] = {0,0,0};
        double dipole_rrms=0;
		
        /*vector<double> force = vector<double>(3); // force, K / A */ // old, slower way to store 3-value vectors.
       
        string UFFlabel; // the UFF-style label for this atom. 
        vector<int> bonds; // <IDs of bonded atoms>
        
        double * get_acc() {
            static double output[3];
            for (unsigned int n=0; n<3; n++) output[n] = acc[n];
            return output;
        }

        void calc_acc(unsigned int nh, double lm) {
            for (unsigned int n=0; n<3; n++) {
                old_acc[n] = acc[n];
                acc[n] = force[n]*1.3806488e-33 / m; // to A / fs^2
                if (nh)
                    acc[n] += lm*vel[n]; // apply Nose-Hoover via Lagrange Multiplier; Rapaport p158
            }
        }

        void calc_vel(double dt) {
            for (unsigned int n=0; n<3; n++) vel[n] = vel[n] + 0.5*(acc[n] + old_acc[n])*dt; // that's where VV comes into play. 1/2 * (a - prev_a)
        }

        void calc_pos(double dt) {
            for (unsigned int n=0; n<3; n++) pos[n] = pos[n] + vel[n] * dt + 0.5*acc[n] * dt * dt;
        }

        /* for debugging */
        void printAll() {
            printf("atom (PDBID %i) %s on molecule %s (PBDID %i) frozen= %i \n -----> m = %f amu; eps = %f K; sig = %f A; alpha = %f A^3; q = %f e\n", PDBID, name.c_str(), mol_name.c_str(), mol_PDBID, frozen, m/(1.660578e-27), eps, sig, polar, C/408.7816);
        }
};

Atom::Atom() {}


class Molecule {
	public:
		Molecule();
		vector<Atom> atoms; // vector that holds this molecule's atoms
		unsigned int PDBID; // the molecule's PDBID (from input)
		string name; // the molecule name/label (from input), e.g. H2 or MOF
        int_fast8_t frozen; //0 or 1
        // arrays are way faster than vectors.
        double force[3] = {0,0,0};
        double torque[3] = {0,0,0};
        double com[3] = {0,0,0};
        double original_com[3] = {0,0,0}; // for diffision calc
        double diffusion_corr[3] = {0,0,0}; // for diffusion calc (accounts for PBC)
        double acc[3] = {0,0,0};
        double old_acc[3] = {0,0,0};
        double vel[3] = {0,0,0};
        double ang_vel[3] = {0,0,0};
        double ang_acc[3] = {0,0,0};
        double old_ang_acc[3] = {0,0,0};
        double ang_pos[3] = {0,0,0};
        double mass=0.0;
        double inertia=0.0; //moment of inertia. stored in K fs^2
        double inertia_tensor[6] = {0,0,0,0,0,0}; // xx,yy,zz,xy,yz,xz
        double fugacity=0.0;
        double md_velx_goal = 0.0; // for Andersen thermostat, unique for unique molecules.
        unsigned int dof=3; // degrees of freedom.

        void reInitialize() {
            // if there are no atoms, don't bother
            if (atoms.size() > 0) {
            while (!atoms.empty()) atoms.pop_back();
            mass=0;
            inertia=0;
            for (unsigned int n=0; n<3; n++) {
                com[n] = 0;
                force[n]=0;
                torque[n]=0;
                acc[n]=0;
                old_acc[n]=0;
                vel[n]=0;
                ang_vel[n]=0;
                ang_acc[n]=0;
                old_ang_acc[n]=0;
                ang_pos[n]=0;
            }
            name = "";
            PDBID=0;
            frozen = 0; // movable
            dof=0;
            }
        }

        double get_mass() {
            // ( mass is generated at input in io.cpp )
            return mass;
        }

        void calc_inertia() {
            for (unsigned int i=0; i<atoms.size(); i++) {
                double rsq = (atoms[i].pos[0] - com[0])*(atoms[i].pos[0] - com[0]) + (atoms[i].pos[1] - com[1])*(atoms[i].pos[1] - com[1]) + (atoms[i].pos[2] - com[2])*(atoms[i].pos[2] - com[2]);
                inertia += atoms[i].m * rsq; // kg * A^2
            }
            inertia = inertia/1.3806488e-23/1e20*1e30; // to K fs^2
        }

        void calc_inertia_tensor() {
            // xx,yy,zz,xy,yz,xz
            for (unsigned int n=0;n<6;n++) inertia_tensor[n]=0; // reset to 0
            for (unsigned int i=0; i<atoms.size(); i++) {
                double x = atoms[i].pos[0]-com[0];
                double y = atoms[i].pos[1]-com[1];
                double z = atoms[i].pos[2]-com[2];
                double x2 = x*x, y2=y*y, z2=z*z;
                double m = atoms[i].m;

                inertia_tensor[0] += m*(y2+z2);
                inertia_tensor[1] += m*(x2+z2);
                inertia_tensor[2] += m*(x2+y2);
                inertia_tensor[3] -= m*x*y;
                inertia_tensor[4] -= m*y*z;
                inertia_tensor[5] -= m*x*z; // all in kg*A^2
            }
            for (unsigned int n=0;n<6;n++) inertia_tensor[n] = inertia_tensor[n]/1.3806488e-23/1e20*1e30; // to K fs^2
        }

        // angular acceleration
        void calc_ang_acc() {
            for (unsigned int n=0; n<3; n++) {
                old_ang_acc[n] = ang_acc[n];
                ang_acc[n] = torque[n] / inertia; // in rad / fs^2
            }
        }

        // linear acceleration
        void calc_acc(unsigned int nh, double lm) {
            for (unsigned int n=0; n<3; n++) {
                old_acc[n] = acc[n];
                acc[n] = force[n]*1.3806488e-33 / mass; // to A / fs^2
                if (nh)
                    acc[n] += lm*vel[n]; // apply Nose-Hoover via Lagrange Multiplier; Rapaport p158
            }
        }

        // angular velocity
        void calc_ang_vel(double dt) {
            for (unsigned int n=0; n<3; n++) {
                ang_vel[n] = ang_vel[n] + 0.5*(ang_acc[n] + old_ang_acc[n])*dt; // in rad/fs
            }
        }

        // linear velocity
        void calc_vel(double dt) {
            for (unsigned int n=0; n<3; n++) {
                vel[n] = vel[n] + 0.5*(acc[n] + old_acc[n])*dt; // in A/fs. vel. verlet
            }
        }

        // angular position // in rad
        void calc_ang_pos(double dt) {
            for (unsigned int n=0; n<3; n++) {
                ang_pos[n] = ang_pos[n] + ang_vel[n] * dt + 0.5*ang_acc[n] * dt * dt;
            }
        }

        // linear position
        void calc_pos(double dt) {
            for (unsigned int i=0; i<atoms.size(); i++) {
              for (unsigned int n=0; n<3; n++) atoms[i].pos[n] = atoms[i].pos[n] + vel[n] * dt + 0.5*acc[n] * dt * dt;
            }
        }

        void calc_force() {
            // A molecule's force IS the sum of its atoms' forces
            // this is external force on the molecule.
            force[0]=0; force[1]=0; force[2]=0;
            for (unsigned int i=0; i<atoms.size(); i++) {
                for (unsigned int n=0; n<3; n++) force[n] += atoms[i].force[n]; // in K/A
            }
        }

        void calc_center_of_mass() {
            // assigns the current center of mass of the molecule based on positions of atoms
            double x_mass_sum=0.0; double y_mass_sum=0.0; double z_mass_sum=0.0;// double mass_sum=0.0;

            for (unsigned int i=0; i<atoms.size(); i++) {
                double atom_mass = atoms[i].m;

                x_mass_sum += atoms[i].pos[0]*atom_mass;
                y_mass_sum += atoms[i].pos[1]*atom_mass;
                z_mass_sum += atoms[i].pos[2]*atom_mass;
            }

            com[0] = x_mass_sum/mass;//_sum;
            com[1] = y_mass_sum/mass;//_sum;
            com[2] = z_mass_sum/mass;//_sum;

        }

        void calc_torque() {
            for (unsigned int n=0; n<3; n++) torque[n] = 0.0;
            // torque is the cross product rxF NOT Fxr, the sum of all atoms relative to molecule's com.
            for (unsigned int i=0; i<atoms.size(); i++) {
                atoms[i].torque[0] = (atoms[i].pos[1]-com[1]) * atoms[i].force[2] - (atoms[i].pos[2]-com[2]) * atoms[i].force[1];
                atoms[i].torque[1] = (atoms[i].pos[2]-com[2]) * atoms[i].force[0] - (atoms[i].pos[0]-com[0]) * atoms[i].force[2];
                atoms[i].torque[2] = (atoms[i].pos[0]-com[0]) * atoms[i].force[1] - (atoms[i].pos[1]-com[1]) * atoms[i].force[0];
                // molecular torque = sum of atomic torques
                for (unsigned int n=0; n<3; n++) torque[n] += atoms[i].torque[n]; // in K
            } // end atomic loop
        } // end calc_torque()

        // for debugging
        void printAll() {
            printf("====================\nmolecule PDBID=%i :: mass: %f amu; inertia: %e; \nname = %s; frozen = %i; \nforce: %f %f %f; \nacc: %f %f %f; \nold_acc: %f %f %f; \nvel: %f %f %f; \ncom: %f %f %f; \ntorque: %f %f %f \nang_acc: %f %f %f \nold_ang_acc: %f %f %f \nang_vel: %f %f %f; \nang_pos: %f %f %f (in degrees) \n",
            PDBID,mass/1.660578e-27,inertia,name.c_str(),frozen,
            force[0], force[1], force[2], acc[0], acc[1], acc[2],
            old_acc[0], old_acc[1], old_acc[2], vel[0], vel[1], vel[2], com[0], com[1], com[2],
            torque[0], torque[1], torque[2], ang_acc[0], ang_acc[1], ang_acc[2],
            old_ang_acc[0], old_ang_acc[1], old_ang_acc[2], ang_vel[0], ang_vel[1], ang_vel[2],
            ang_pos[0]*180.0/M_PI, ang_pos[1]*180.0/M_PI, ang_pos[2]*180.0/M_PI); //com[0], com[1], com[2]);

        }
};

Molecule::Molecule() {}

Constants::Constants() {
	// ATOM DEFAULTS LIBRARY
	// MASS VALUES g/mol -> kg/particle
	masses["HB"] = 2.016*cM; // buch model	h2
	masses["H2G"] = 0.0*cM;
    masses["H2E"] = 1.008*cM;
    masses["H2N"] = 0.0*cM;
    masses["HW"] = 1.008*cM; // H in water ( my model)
    masses["HT"] = 1.008*cM; // H in TIP3P
    masses["H"] = 1.0079*cM;
    masses["HS"] = 1.0079*cM;
	masses["He"] = 4.002602*cM;
	masses["Li"] = 6.941*cM;
	masses["Be"] = 9.012182*cM;
	masses["B"] = 10.811*cM;
	masses["C"] = 12.011*cM;
	masses["C_p"] = 12.0107*cM; // C SAPT
    masses["C_s"] = 12.0107*cM;
    masses["C_t"] = 12.0107*cM;
    masses["C_a"] = 12.0107*cM;
    masses["C_en"] = 12.0107*cM;
    masses["C_yn"] = 12.0107*cM;
    masses["C_ony"] = 12.0107*cM;
    masses["N_sp3"] = 14.0067*cM; // N SAPT
    masses["N_sp2"] = 14.0067*cM;
    masses["N"] = 14.007*cM;
    masses["O_sp3"] = 15.9994*cM; // O SAPT
    masses["O_sp2"] = 15.9994*cM;
	masses["O"] = 15.9998*cM;
	masses["O2"] = 32.0*cM; // my O2 model
	masses["OW"] = 15.9998*cM; // O in water (my model)
    masses["OT"] = 15.9998*cM; // O in TIP3P
    masses["F"] = 18.998*cM;
	masses["Ne"] = 20.1797*cM;
	masses["Na"] = 22.98976928*cM;
    masses["Mg"] = 24.305*cM;
    masses["Al"] = 26.981538*cM;
	masses["Si"] = 28.085*cM;
	masses["P"] = 30.973*cM;
	masses["S"] = 32.06*cM;
    masses["SS"] = 32.065*cM; /// S SAPT
	masses["Cl"] = 35.45*cM;
	masses["Ar"] = 39.948*cM;
    masses["ArS"] = 39.948*cM; // Ar SAPT
masses["K"] = 39.0983*cM;
masses["Ca"] = 40.078*cM;
masses["Sc"] = 44.955910*cM;
masses["Ti"] = 47.867*cM;
masses["V"] = 50.9415*cM;
masses["Cr"] = 51.9961*cM;
masses["Mn"] = 54.938049*cM;
masses["Fe"] = 55.845*cM;
masses["Co"] = 58.933200*cM;
masses["Ni"] = 58.6934*cM;
masses["Cu"] = 63.546*cM;
masses["Zn"] = 65.39*cM;
masses["Ga"] = 69.723*cM;
masses["Ge"] = 72.61*cM;
masses["As"] = 74.92160*cM;
masses["Se"] = 78.96*cM;
masses["Br"] = 79.904*cM;
	masses["Kr"] = 83.798*cM;

masses["Rb"] = 85.4678*cM;
masses["Sr"] = 87.62*cM;
masses["Y"] = 88.90585*cM;
masses["Zr"] = 91.224*cM;
masses["Nb"] = 92.90638*cM;
masses["Mo"] = 95.94*cM;
masses["Tc"] = 98.0*cM;
masses["Ru"] = 101.07*cM;
masses["Rh"] = 102.90550*cM;
masses["Pd"] = 106.42*cM;
masses["Ag"] = 107.8682*cM;
masses["Cd"] = 112.411*cM;
masses["In"] = 114.818*cM;
masses["Sn"] = 118.710*cM;
masses["Sb"] = 121.760*cM;
masses["Te"] = 127.60*cM;
masses["I"] = 126.90447*cM;
masses["Xe"] = 131.29*cM;

masses["Cs"] = 132.90545*cM;
masses["Ba"] = 137.327*cM;
masses["La"] = 138.9055*cM;
masses["Hf"] = 178.49*cM;
masses["Ta"] = 180.9479*cM;
masses["W"] = 183.84*cM;
masses["Re"] = 186.207*cM;
masses["Os"] = 190.23*cM;
masses["Ir"] = 192.217*cM;
masses["Pt"] = 195.078*cM;
masses["Au"] = 196.96655*cM;
masses["Hg"] = 200.59*cM;
masses["Tl"] = 204.3833*cM;
masses["Pb"] = 207.2*cM;
masses["Bi"] = 208.98038*cM;
masses["Po"] = 209.0*cM;
masses["At"] = 210.*cM;
masses["Rn"] = 222.0176*cM;

masses["Fr"] = 223.*cM;
masses["Ra"] = 226.*cM;
masses["Ac"] = 227.*cM;
masses["Rf"] = 261.*cM;
masses["Db"] = 262.*cM;
masses["Sg"] = 263.*cM;
masses["Bh"] = 262.*cM;
masses["Hs"] = 265.*cM;
masses["Mt"] = 266.*cM;

masses["Ce"] = 140.116*cM;
masses["Pr"] = 140.90765*cM;
masses["Nd"] = 144.24*cM;
masses["Pm"] = 145.*cM;
masses["Sm"] = 150.36*cM;
masses["Eu"] = 151.964*cM;
masses["Gd"] = 157.25*cM;
masses["Tb"] = 158.92534*cM;
masses["Dy"] = 162.5*cM;
masses["Ho"] = 164.93032*cM;
masses["Er"] = 167.26*cM;
masses["Tm"] = 168.93421*cM;
masses["Yb"] = 173.04*cM;
masses["Lu"] = 174.967*cM;

masses["Th"] = 232.0381*cM;
masses["Pa"] = 231.03588*cM;
masses["U"] = 238.0289*cM;
masses["Np"] = 237.*cM;
masses["Pu"] = 244.*cM;
masses["Am"] = 243.*cM;
masses["Cm"] = 247.*cM;
masses["Bk"] = 247.*cM;
masses["Cf"] = 251.*cM;
masses["Es"] = 252.*cM;
masses["Fm"] = 257.*cM;
masses["Md"] = 258.*cM;
masses["No"] = 259.*cM;
masses["Lr"] = 262.*cM;

	// LJ SIGMA VALUES (A)
	sigs["HB"] = 2.96; // buch model h2 (from mpmc sig)
	sigs["H2G"] = 3.2293; // mpmc -> meters
    sigs["H2E"] = 0.0;
    sigs["H2N"] = 2.3406;
	sigs["HW"] = 0.5; // H in water, my model; old (sticky MD) 0.4
    sigs["HT"] = 0.0; // H in TIP3P
    sigs["HS"] = 3.09061; // H SAPT (Adam)
    sigs["C_s"] = 3.5786; // C SAPT
    sigs["C_p"] = 3.63956;
    sigs["C_t"] = 3.37707;
    sigs["C_a"] = 3.65947;
    sigs["C_en"] = 3.74892;
    sigs["C_yn"] = 3.79838;
    sigs["C_ony"] = 3.56023;
    sigs["N_sp3"] = 3.32588; // N SAPT
    sigs["N_sp2"] = 3.45133;
	sigs["OW"] = 1.5; // O in water (my model) -- old (sticky MD) 1.4
    sigs["OT"] = 3.15061; // O in TIP3P
    sigs["O_sp3"] = 3.15611; // O SAPT
    sigs["O_sp2"] = 3.34161;
    sigs["SS"] = 3.8899; // S SAPT
    sigs["ArS"] = 3.37191; // Ar SAPT
    // UFF
   sigs["H"] = 2.886*uff2mpmc;
sigs["He"] = 2.362*uff2mpmc;
sigs["Li"] = 2.451*uff2mpmc;
sigs["Be"] = 2.745*uff2mpmc;
sigs["B"] = 4.083*uff2mpmc;
sigs["C"] = 3.851*uff2mpmc;
sigs["N"] = 3.660*uff2mpmc;
sigs["O"] = 3.50*uff2mpmc;
sigs["F"] = 3.364*uff2mpmc;
sigs["Ne"] = 3.243*uff2mpmc;
sigs["Na"] = 2.983*uff2mpmc;
sigs["Mg"] = 3.021*uff2mpmc;
sigs["Al"] = 4.499*uff2mpmc;
sigs["Si"] = 4.295*uff2mpmc;
sigs["P"] = 4.147*uff2mpmc;
sigs["S"] = 4.035*uff2mpmc;
sigs["Cl"] = 3.947*uff2mpmc;
sigs["Ar"] = 3.868*uff2mpmc;
sigs["K"] = 3.812*uff2mpmc;
sigs["Ca"] = 3.399*uff2mpmc;
sigs["Sc"] = 3.295*uff2mpmc;
sigs["Ti"] = 3.175*uff2mpmc;
sigs["V"] = 3.144*uff2mpmc;
sigs["Cr"] = 3.023*uff2mpmc;
sigs["Mn"] = 2.961*uff2mpmc;
sigs["Fe"] = 2.912*uff2mpmc;
sigs["Co"] = 2.872*uff2mpmc;
sigs["Ni"] = 2.834*uff2mpmc;
sigs["Cu"] = 3.495*uff2mpmc;
sigs["Zn"] = 2.763*uff2mpmc;
sigs["Ga"] = 4.383*uff2mpmc;
sigs["Ge"] = 4.280*uff2mpmc;
sigs["As"] = 4.230*uff2mpmc;
sigs["Se"] = 4.205*uff2mpmc;
sigs["Br"] = 4.189*uff2mpmc;
sigs["Kr"] = 4.141*uff2mpmc;
sigs["Rb"] = 4.114*uff2mpmc;
sigs["Sr"] = 3.641*uff2mpmc;
sigs["Y"] = 3.345*uff2mpmc;
sigs["Zr"] = 3.124*uff2mpmc;
sigs["Nb"] = 3.165*uff2mpmc;
sigs["Mo"] = 3.052*uff2mpmc;
sigs["Tc"] = 2.998*uff2mpmc;
sigs["Ru"] = 2.963*uff2mpmc;
sigs["Rh"] = 2.929*uff2mpmc;
sigs["Pd"] = 2.899*uff2mpmc;
sigs["Ag"] = 3.148*uff2mpmc;
sigs["Cd"] = 2.848*uff2mpmc;
sigs["In"] = 4.463*uff2mpmc;
sigs["Sn"] = 4.392*uff2mpmc;
sigs["Sb"] = 4.420*uff2mpmc;
sigs["Te"] = 4.470*uff2mpmc;
sigs["I"] = 4.5*uff2mpmc;
sigs["Xe"] = 4.404*uff2mpmc;
sigs["Cs"] = 4.517*uff2mpmc;
sigs["Ba"] = 3.703*uff2mpmc;
sigs["La"] = 3.522*uff2mpmc;
sigs["Ce"] = 3.556*uff2mpmc;
sigs["Pr"] = 3.606*uff2mpmc;
sigs["Nd"] = 3.575*uff2mpmc;
sigs["Pm"] = 3.547*uff2mpmc;
sigs["Sm"] = 3.52*uff2mpmc;
sigs["Eu"] = 3.493*uff2mpmc;
sigs["Gd"] = 3.368*uff2mpmc;
sigs["Tb"] = 3.451*uff2mpmc;
sigs["Dy"] = 3.428*uff2mpmc;
sigs["Ho"] = 3.409*uff2mpmc;
sigs["Er"] = 3.391*uff2mpmc;
sigs["Tm"] = 3.374*uff2mpmc;
sigs["Yb"] = 3.355*uff2mpmc;
sigs["Lu"] = 3.64*uff2mpmc;
sigs["Hf"] = 3.141*uff2mpmc;
sigs["Ta"] = 3.17*uff2mpmc;
sigs["W"] = 3.069*uff2mpmc;
sigs["Re"] = 2.954*uff2mpmc;
sigs["Os"] = 3.12*uff2mpmc;
sigs["Ir"] = 2.84*uff2mpmc;
sigs["Pt"] = 2.754*uff2mpmc;
sigs["Au"] = 3.293*uff2mpmc;
sigs["Hg"] = 2.705*uff2mpmc;
sigs["Tl"] = 4.347*uff2mpmc;
sigs["Pb"] = 4.297*uff2mpmc;
sigs["Bi"] = 4.37*uff2mpmc;
sigs["Po"] = 4.709*uff2mpmc;
sigs["At"] = 4.750*uff2mpmc;
sigs["Rn"] = 4.765*uff2mpmc;
sigs["Fr"] = 4.90*uff2mpmc;
sigs["Ra"] = 3.677*uff2mpmc;
sigs["Ac"] = 3.478*uff2mpmc;
sigs["Th"] = 3.396*uff2mpmc;
sigs["Pa"] = 3.424*uff2mpmc;
sigs["U"] = 3.395*uff2mpmc;
sigs["Np"] = 3.424*uff2mpmc;
sigs["Pu"] = 3.424*uff2mpmc;
sigs["Am"] = 3.381*uff2mpmc;
sigs["Cm"] = 3.326*uff2mpmc;
sigs["Bk"] = 3.339*uff2mpmc;
sigs["Cf"] = 3.313*uff2mpmc;
sigs["Es"] = 3.299*uff2mpmc;
sigs["Fm"] = 3.286*uff2mpmc;
sigs["Md"] = 3.274*uff2mpmc;
sigs["No"] = 3.248*uff2mpmc;
sigs["Lr"] = 3.236*uff2mpmc; 

    // UFF4MOF sigs
    // nevermind, they only give bonding parameters.

	// LJ EPSILON VALUES ( /kbk means kcal/mol -> K)
	eps["HB"] = 34.20; // buch model h2
    eps["H2G"] = 8.8516; // bss model h2
    eps["H2E"] = 0.0; // bss
    eps["H2N"] = 4.0659; // bss
    eps["HT"] = 0.0; // H in TIP3P
    eps["HS"] = 0.66563; // H SAPT
    eps["C_p"] = 36.692; // C SAPT
    eps["C_s"] = 31.35824;
    eps["C_t"] = 41.45435;
    eps["C_a"] = 22.30908;
    eps["C_en"] = 26.88878;
    eps["C_yn"] = 22.40343;
    eps["C_ony"] = 18.09254;
    eps["N_sp3"] = 36.97995; // N SAPT
    eps["N_sp2"] = 24.25732;
    eps["OT"] = 0.6364/0.0083144621; // O in TIP3P
	eps["O_sp3"] = 30.01345; // O SAPT
    eps["O_sp2"] = 21.81177;
    eps["SS"] = 53.02994; // S SAPT
    eps["ArS"] = 128.32680; // Ar SAPT

// UFF eps
eps["H"] = 0.044/kbk;
eps["He"] = 0.056/kbk;
eps["Li"] = 0.025/kbk;
eps["Be"] = 0.085/kbk;
eps["B"] = 0.180/kbk;
eps["C"] = 0.105/kbk;
eps["N"] = 0.069/kbk;
eps["O"] = 0.06/kbk;
eps["F"] = 0.05/kbk;
eps["Ne"] = 0.042/kbk;
eps["Na"] = 0.03/kbk;
eps["Mg"] = 0.111/kbk;
eps["Al"] = 0.505/kbk;
eps["Si"] = 0.402/kbk;
eps["P"] = 0.305/kbk;
eps["S"] = 0.274/kbk;
eps["Cl"] = 0.227/kbk;
eps["Ar"] = 0.185/kbk;
eps["K"] = 0.035/kbk;
eps["Ca"] = 0.238/kbk;
eps["Sc"] = 0.019/kbk;
eps["Ti"] = 0.017/kbk;
eps["V"] = 0.016/kbk;
eps["Cr"] = 0.015/kbk;
eps["Mn"] = 0.013/kbk;
eps["Fe"] = 0.013/kbk;
eps["Co"] = 0.014/kbk;
eps["Ni"] = 0.015/kbk;
eps["Cu"] = 0.005/kbk;
eps["Zn"] = 0.124/kbk;
eps["Ga"] = 0.415/kbk;
eps["Ge"] = 0.379/kbk;
eps["As"] = 0.309/kbk;
eps["Se"] = 0.291/kbk;
eps["Br"] = 0.251/kbk;
eps["Kr"] = 0.220/kbk;
eps["Rb"] = 0.04/kbk;
eps["Sr"] = 0.235/kbk;
eps["Y"] = 0.072/kbk;
eps["Zr"] = 0.069/kbk;
eps["Nb"] = 0.059/kbk;
eps["Mo"] = 0.056/kbk;
eps["Tc"] = 0.048/kbk;
eps["Ru"] = 0.056/kbk;
eps["Rh"] = 0.053/kbk;
eps["Pd"] = 0.048/kbk;
eps["Ag"] = 0.036/kbk;
eps["Cd"] = 0.228/kbk;
eps["In"] = 0.599/kbk;
eps["Sn"] = 0.567/kbk;
eps["Sb"] = 0.449/kbk;
eps["Te"] = 0.398/kbk;
eps["I"] = 0.339/kbk;
eps["Xe"] = 0.332/kbk;
eps["Cs"] = 0.045/kbk;
eps["Ba"] = 0.364/kbk;
eps["La"] = 0.017/kbk;
eps["Ce"] = 0.013/kbk;
eps["Pr"] = 0.010/kbk;
eps["Nd"] = 0.009/kbk;
eps["Pm"] = 0.008/kbk;
eps["Sm"] = 0.008/kbk;
eps["Eu"] = 0.008/kbk;
eps["Gd"] = 0.009/kbk;
eps["Tb"] = 0.007/kbk;
eps["Dy"] = 0.007/kbk;
eps["Ho"] = 0.007/kbk;
eps["Er"] = 0.007/kbk;
eps["Tm"] = 0.006/kbk;
eps["Yb"] = 0.228/kbk;
eps["Lu"] = 0.041/kbk;
eps["Hf"] = 0.072/kbk;
eps["Ta"] = 0.081/kbk;
eps["W"] = 0.067/kbk;
eps["Re"] = 0.066/kbk;
eps["Os"] = 0.037/kbk;
eps["Ir"] = 0.073/kbk;
eps["Pt"] = 0.080/kbk;
eps["Au"] = 0.039/kbk;
eps["Hg"] = 0.385/kbk;
eps["Tl"] = 0.680/kbk;
eps["Pb"] = 0.663/kbk;
eps["Bi"] = 0.518/kbk;
eps["Po"] = 0.325/kbk;
eps["At"] = 0.284/kbk;
eps["Rn"] = 0.248/kbk;
eps["Fr"] = 0.050/kbk;
eps["Ra"] = 0.404/kbk;
eps["Ac"] = 0.033/kbk;
eps["Th"] = 0.026/kbk;
eps["Pa"] = 0.022/kbk;
eps["U"] = 0.022/kbk;
eps["Np"] = 0.019/kbk;
eps["Pu"] = 0.016/kbk;
eps["Am"] = 0.014/kbk;
eps["Cm"] = 0.013/kbk;
eps["Bk"] = 0.013/kbk;
eps["Cf"] = 0.013/kbk;
eps["Es"] = 0.012/kbk;
eps["Fm"] = 0.012/kbk;
eps["Md"] = 0.011/kbk;
eps["No"] = 0.011/kbk;
eps["Lr"] = 0.011/kbk;

	// POLARIZABILITIES  // in A^3 
    // these are VAN DUIJNEN EXPONENTIAL DAMPING POLARIZABILITIES
    // IT WOULD BE DIFFERENT FOR LINEAR DAMPING
    polars["H"] = 2.7927*vand2mpmc;    
    polars["C"] = 8.6959*vand2mpmc;
    polars["N"] = 6.5565*vand2mpmc;
    polars["O"] = 5.7494*vand2mpmc;
    polars["F"] = 3.0013*vand2mpmc;
    polars["S"] = 16.6984*vand2mpmc;
    polars["Cl"] = 16.1979*vand2mpmc;
    polars["Br"] = 23.5714*vand2mpmc;
    polars["I"] = 36.9880*vand2mpmc;
	//polars["H"] = 0.41380;//*cV/ke;
	polars["HW"] = 0.41380;//*cV/ke; // H in water (my model)
    polars["HS"] = 0.41380; // H SAPT
    polars["B"] = 0.6634;//*cV/ke;
	//polars["C"] = 1.2866;//*cV/ke;
    polars["C_p"] = polars["C_s"] = polars["C_t"] = polars["C_a"] = polars["C_en"] = polars["C_yn"] = polars["C_ony"] = 1.2866; // C SAPT
	//polars["N"] = 0.97157;//*cV/ke;
    polars["N_sp3"] = polars["N_sp2"] = 0.97157; // N SAPT
	//polars["O"] = 0.852;//*cV/ke;
    polars["OW"] = 0.852;//*cV/ke; // O in water (my model)
    polars["O_sp3"] = polars["O_sp2"] = 0.852; // O SAPT
	polars["Na"] = 24.11;//*cV/ke; // from paper https://www.researchgate.net/publication/45896756_Absolute_and_ratio_measurements_of_the_polarizability_of_Na_K_and_Rb_with_an_atom_interferometer
	polars["P"] = 3.35;//*cV/ke;
    polars["SS"] = 2.474; // S SAPT
	polars["Cl"] = 2.40028;//*cV/ke;
    polars["ArS"] = 1.63922; // Ar SAPT
	polars["Co"] = 3.26440; // from Kay/Shannelle.
    polars["Cu"] = 2.19630;//*cV/ke;
	polars["Zn"] = 1.98870;//*cV/ke;
	polars["Br"] = 3.49300; // used in MPM-1-Br work
	polars["Ru"] = 5.191; // I calculated this by Adam's fitting method.
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

// ATOMIC NUMBERS Z (proton count)
elements["H"] = 1.0;
elements["He"] = 2.0;
elements["Li"] = 3.0;
elements["Be"] = 4.0;
elements["B"] = 5.0;
elements["C"] = 6.0;
elements["N"] = 7.0;
elements["O"] = 8.0;
elements["F"] = 9.0;
elements["Ne"] = 10.0;
elements["Na"] = 11.0;
elements["Mg"] = 12.0;
elements["Al"] = 13.0;
elements["Si"] = 14.0;
elements["P"] = 15.0;
elements["S"] = 16.0;
elements["Cl"] = 17.0;
elements["Ar"] = 18.0;
elements["K"] = 19.0;
elements["Ca"] = 20.0;
elements["Sc"] = 21.0;
elements["Ti"] = 22.0;
elements["V"] = 23.0;
elements["Cr"] = 24.0;
elements["Mn"] = 25.0;
elements["Fe"] = 26.0;
elements["Co"] = 27.0;
elements["Ni"] = 28.0;
elements["Cu"] = 29.0;
elements["Zn"] = 30.0;
elements["Ga"] = 31.0;
elements["Ge"] = 32.0;
elements["As"] = 33.0;
elements["Se"] = 34.0;
elements["Br"] = 35.0;
elements["Kr"] = 36.0;
elements["Rb"] = 37.0;
elements["Sr"] = 38.0;
elements["Y"] = 39.0;
elements["Zr"] = 40.0;
elements["Nb"] = 41.0;
elements["Mo"] = 42.0;
elements["Tc"] = 43.0;
elements["Ru"] = 44.0;
elements["Rh"] = 45.0;
elements["Pd"] = 46.0;
elements["Ag"] = 47.0;
elements["Cd"] = 48.0;
elements["In"] = 49.0;
elements["Sn"] = 50.0;
elements["Sb"] = 51.0;
elements["Te"] = 52.0;
elements["I"] = 53.0;
elements["Xe"] = 54.0;
elements["Cs"] = 55.0;
elements["Ba"] = 56.0;
elements["La"] = 57.0;
elements["Ce"] = 58.0;
elements["Pr"] = 59.0;
elements["Nd"] = 60.0;
elements["Pm"] = 61.0;
elements["Sm"] = 62.0;
elements["Eu"] = 63.0;
elements["Gd"] = 64.0;
elements["Tb"] = 65.0;
elements["Dy"] = 66.0;
elements["Ho"] = 67.0;
elements["Er"] = 68.0;
elements["Tm"] = 69.0;
elements["Yb"] = 70.0;
elements["Lu"] = 71.0;
elements["Hf"] = 72.0;
elements["Ta"] = 73.0;
elements["W"] = 74.0;
elements["Re"] = 75.0;
elements["Os"] = 76.0;
elements["Ir"] = 77.0;
elements["Pt"] = 78.0;
elements["Au"] = 79.0;
elements["Hg"] = 80.0;
elements["Tl"] = 81.0;
elements["Pb"] = 82.0;
elements["Bi"] = 83.0;
elements["Po"] = 84.0;
elements["At"] = 85.0;
elements["Rn"] = 86.0;
elements["Fr"] = 87.0;
elements["Ra"] = 88.0;
elements["Ac"] = 89.0;
elements["Th"] = 90.0;
elements["Pa"] = 91.0;
elements["U"] = 92.0;
elements["Np"] = 93.0;
elements["Pu"] = 94.0;
elements["Am"] = 95.0;
elements["Cm"] = 96.0;
elements["Bk"] = 97.0;
elements["Cf"] = 98.0;
elements["Es"] = 99.0;
elements["Fm"] = 100.0;
elements["Md"] = 101.0;
elements["No"] = 102.0;
elements["Lr"] = 103.0;
elements["Rf"] = 104.0;
elements["Db"] = 105.0;
elements["Sg"] = 106.0;
elements["Bh"] = 107.0;
elements["Hs"] = 108.0;
elements["Mt"] = 109.0;

// ATOMIC RADII (VDW RADII) in Angstroms
// http://www.crystalmaker.com/support/tutorials/crystalmaker/atomic-radii/index.html
// E Clementi, D L Raimondi, W P Reinhardt (1963) J Chem Phys. 38:2686.
radii["H"] = 0.53;
radii["He"] = 0.31;
radii["Li"] = 1.67;
radii["Be"] = 1.12;
radii["B"] = 0.87;
radii["C"] = 0.67;
radii["N"] = 0.56;
radii["O"] = 0.48;
radii["F"] = 0.42;
radii["Ne"] = 0.38;
radii["Na"] = 1.9;
radii["Mg"] = 1.45;
radii["Al"] = 1.18;
radii["Si"] = 1.11;
radii["P"] = 0.98;
radii["S"] = 0.88;
radii["Cl"] = 0.79;
radii["Ar"] = 0.71;
radii["K"] = 2.43;
radii["Ca"] = 1.94;
radii["Sc"] = 1.84;
radii["Ti"] = 1.76;
radii["V"] = 1.71;
radii["Cr"] = 1.66;
radii["Mn"] = 1.61;
radii["Fe"] = 1.56;
radii["Co"] = 1.52;
radii["Ni"] = 1.49;
radii["Cu"] = 1.45;
radii["Zn"] = 1.42;
radii["Ga"] = 1.36;
radii["Ge"] = 1.25;
radii["As"] = 1.14;
radii["Se"] = 1.03;
radii["Br"] = 0.94;
radii["Kr"] = 0.88;
radii["Rb"] = 2.65;
radii["Sr"] = 2.19;
radii["Y"] = 2.12;
radii["Zr"] = 2.06;
radii["Nb"] = 1.98;
radii["Mo"] = 1.9;
radii["Tc"] = 1.83;
radii["Ru"] = 1.78;
radii["Rh"] = 1.73;
radii["Pd"] = 1.69;
radii["Ag"] = 1.65;
radii["Cd"] = 1.61;
radii["In"] = 1.56;
radii["Sn"] = 1.45;
radii["Sb"] = 1.33;
radii["Te"] = 1.23;
radii["I"] = 1.15;
radii["Xe"] = 1.08;
radii["Cs"] = 2.98;
radii["Ba"] = 2.53;
radii["La"] = 1.95;
radii["Ce"] = 1.85;
radii["Pr"] = 2.47;
radii["Nd"] = 2.06;
radii["Pm"] = 2.05;
radii["Sm"] = 2.38;
radii["Eu"] = 2.31;
radii["Gd"] = 2.33;
radii["Tb"] = 2.25;
radii["Dy"] = 2.28;
radii["Ho"] = 2.26;
radii["Er"] = 2.26;
radii["Tm"] = 2.22;
radii["Yb"] = 2.22;
radii["Lu"] = 2.17;
radii["Hf"] = 2.08;
radii["Ta"] = 2;
radii["W"] = 1.93;
radii["Re"] = 1.88;
radii["Os"] = 1.85;
radii["Ir"] = 1.8;
radii["Pt"] = 1.77;
radii["Au"] = 1.74;
radii["Hg"] = 1.71;
radii["Tl"] = 1.56;
radii["Pb"] = 1.54;
radii["Bi"] = 1.43;
radii["Po"] = 1.35;
radii["At"] = 1.27;
radii["Rn"] = 1.2;
radii["Ac"] = 1.95;
radii["Th"] = 1.8;
radii["Pa"] = 1.8;
radii["U"] = 1.75;
radii["Np"] = 1.75;
radii["Pu"] = 1.75;
radii["Am"] = 1.75;


// =============
// UFF BONDING PARAMETERS
// bonds in Angstrom; Angles in degrees; Eff. Charge (Z) in e(lementary charge) 
UFF_bonds["H_"] = 0.354; UFF_angles["H_"] = 180; UFF_Z["H_"] = 0.712;
UFF_bonds["H_b"] = 0.46; UFF_angles["H_b"] = 83.5; UFF_Z["H_b"] = 0.712;
UFF_bonds["He4+4"] = 0.849; UFF_angles["He4+4"] = 90; UFF_Z["He4+4"] = 0.098;
UFF_bonds["Li"] = 1.336; UFF_angles["Li"] = 180; UFF_Z["Li"] = 1.026;
UFF_bonds["Be3+2"] = 1.074; UFF_angles["Be3+2"] = 109.47; UFF_Z["Be3+2"] = 1.565;
UFF_bonds["B_3"] = 0.838; UFF_angles["B_3"] = 109.47; UFF_Z["B_3"] = 1.755;
UFF_bonds["B_2"] = 0.828; UFF_angles["B_2"] = 120; UFF_Z["B_2"] = 1.755;
UFF_bonds["C_3"] = 0.757; UFF_angles["C_3"] = 109.47; UFF_Z["C_3"] = 1.912;
UFF_bonds["C_R"] = 0.729; UFF_angles["C_R"] = 120; UFF_Z["C_R"] = 1.912;
UFF_bonds["C_2"] = 0.732; UFF_angles["C_2"] = 120; UFF_Z["C_2"] = 1.912;
UFF_bonds["C_1"] = 0.706; UFF_angles["C_1"] = 180; UFF_Z["C_1"] = 1.912;
UFF_bonds["N_3"] = 0.7; UFF_angles["N_3"] = 106.7; UFF_Z["N_3"] = 2.544;
UFF_bonds["N_R"] = 0.699; UFF_angles["N_R"] = 120; UFF_Z["N_R"] = 2.544;
UFF_bonds["N_2"] = 0.685; UFF_angles["N_2"] = 111.2; UFF_Z["N_2"] = 2.544;
UFF_bonds["N_1"] = 0.656; UFF_angles["N_1"] = 180; UFF_Z["N_1"] = 2.544;
UFF_bonds["O_3"] = 0.658; UFF_angles["O_3"] = 104.51; UFF_Z["O_3"] = 2.3;
UFF_bonds["O_3_z"] = 0.528; UFF_angles["O_3_z"] = 146; UFF_Z["O_3_z"] = 2.3;
UFF_bonds["O_R"] = 0.68; UFF_angles["O_R"] = 110; UFF_Z["O_R"] = 2.3;
UFF_bonds["O_2"] = 0.634; UFF_angles["O_2"] = 120; UFF_Z["O_2"] = 2.3;
UFF_bonds["O_1"] = 0.639; UFF_angles["O_1"] = 180; UFF_Z["O_1"] = 2.3;
UFF_bonds["F_"] = 0.668; UFF_angles["F_"] = 180; UFF_Z["F_"] = 1.735;
UFF_bonds["Ne4+4"] = 0.92; UFF_angles["Ne4+4"] = 90; UFF_Z["Ne4+4"] = 0.194;
UFF_bonds["Na"] = 1.539; UFF_angles["Na"] = 180; UFF_Z["Na"] = 1.081;
UFF_bonds["Mg3+2"] = 1.421; UFF_angles["Mg3+2"] = 109.47; UFF_Z["Mg3+2"] = 1.787;
UFF_bonds["Al3"] = 1.244; UFF_angles["Al3"] = 109.47; UFF_Z["Al3"] = 1.792;
UFF_bonds["Si3"] = 1.117; UFF_angles["Si3"] = 109.47; UFF_Z["Si3"] = 2.323;
UFF_bonds["P_3+3"] = 1.101; UFF_angles["P_3+3"] = 93.8; UFF_Z["P_3+3"] = 2.863;
UFF_bonds["P_3+5"] = 1.056; UFF_angles["P_3+5"] = 109.47; UFF_Z["P_3+5"] = 2.863;
UFF_bonds["P_3+q"] = 1.056; UFF_angles["P_3+q"] = 109.47; UFF_Z["P_3+q"] = 2.863;
UFF_bonds["S_3+2"] = 1.064; UFF_angles["S_3+2"] = 92.1; UFF_Z["S_3+2"] = 2.703;
UFF_bonds["S_3+4"] = 1.049; UFF_angles["S_3+4"] = 103.2; UFF_Z["S_3+4"] = 2.703;
UFF_bonds["S_3+6"] = 1.027; UFF_angles["S_3+6"] = 109.47; UFF_Z["S_3+6"] = 2.703;
UFF_bonds["S_R"] = 1.077; UFF_angles["S_R"] = 92.2; UFF_Z["S_R"] = 2.703;
UFF_bonds["S_2"] = 0.854; UFF_angles["S_2"] = 120; UFF_Z["S_2"] = 2.703;
UFF_bonds["Cl"] = 1.044; UFF_angles["Cl"] = 180; UFF_Z["Cl"] = 2.348;
UFF_bonds["Ar4+4"] = 1.032; UFF_angles["Ar4+4"] = 90; UFF_Z["Ar4+4"] = 0.3;
UFF_bonds["K_"] = 1.953; UFF_angles["K_"] = 180; UFF_Z["K_"] = 1.165;
UFF_bonds["Ca6+2"] = 1.761; UFF_angles["Ca6+2"] = 90; UFF_Z["Ca6+2"] = 2.141;
UFF_bonds["Sc3+3"] = 1.513; UFF_angles["Sc3+3"] = 109.47; UFF_Z["Sc3+3"] = 2.592;
UFF_bonds["Ti3+4"] = 1.412; UFF_angles["Ti3+4"] = 109.47; UFF_Z["Ti3+4"] = 2.659;
UFF_bonds["Ti6+4"] = 1.412; UFF_angles["Ti6+4"] = 90; UFF_Z["Ti6+4"] = 2.659;
UFF_bonds["V_3+5"] = 1.402; UFF_angles["V_3+5"] = 109.47; UFF_Z["V_3+5"] = 2.679;
UFF_bonds["Cr6+3"] = 1.345; UFF_angles["Cr6+3"] = 90; UFF_Z["Cr6+3"] = 2.463;
UFF_bonds["Mn6+2"] = 1.382; UFF_angles["Mn6+2"] = 90; UFF_Z["Mn6+2"] = 2.43;
UFF_bonds["Fe3+2"] = 1.27; UFF_angles["Fe3+2"] = 109.47; UFF_Z["Fe3+2"] = 2.43;
UFF_bonds["Fe6+2"] = 1.335; UFF_angles["Fe6+2"] = 90; UFF_Z["Fe6+2"] = 2.43;
UFF_bonds["Co6+3"] = 1.241; UFF_angles["Co6+3"] = 90; UFF_Z["Co6+3"] = 2.43;
UFF_bonds["Ni4+2"] = 1.164; UFF_angles["Ni4+2"] = 90; UFF_Z["Ni4+2"] = 2.43;
UFF_bonds["Cu3+1"] = 1.302; UFF_angles["Cu3+1"] = 109.47; UFF_Z["Cu3+1"] = 1.756;
UFF_bonds["Zn3+2"] = 1.193; UFF_angles["Zn3+2"] = 109.47; UFF_Z["Zn3+2"] = 1.308;
UFF_bonds["Ga3+3"] = 1.26; UFF_angles["Ga3+3"] = 109.47; UFF_Z["Ga3+3"] = 1.821;
UFF_bonds["Ge3"] = 1.197; UFF_angles["Ge3"] = 109.47; UFF_Z["Ge3"] = 2.789;
UFF_bonds["As3+3"] = 1.211; UFF_angles["As3+3"] = 92.1; UFF_Z["As3+3"] = 2.864;
UFF_bonds["Se3+2"] = 1.19; UFF_angles["Se3+2"] = 90.6; UFF_Z["Se3+2"] = 2.764;
UFF_bonds["Br"] = 1.192; UFF_angles["Br"] = 180; UFF_Z["Br"] = 2.519;
UFF_bonds["Kr4+4"] = 1.147; UFF_angles["Kr4+4"] = 90; UFF_Z["Kr4+4"] = 0.452;
UFF_bonds["Rb"] = 2.26; UFF_angles["Rb"] = 180; UFF_Z["Rb"] = 1.592;
UFF_bonds["Sr6+2"] = 2.052; UFF_angles["Sr6+2"] = 90; UFF_Z["Sr6+2"] = 2.449;
UFF_bonds["Y_3+3"] = 1.698; UFF_angles["Y_3+3"] = 109.47; UFF_Z["Y_3+3"] = 3.257;
UFF_bonds["Zr3+4"] = 1.564; UFF_angles["Zr3+4"] = 109.47; UFF_Z["Zr3+4"] = 3.667;
UFF_bonds["Nb3+5"] = 1.473; UFF_angles["Nb3+5"] = 109.47; UFF_Z["Nb3+5"] = 3.618;
UFF_bonds["Mo6+6"] = 1.467; UFF_angles["Mo6+6"] = 90; UFF_Z["Mo6+6"] = 3.4;
UFF_bonds["Mo3+6"] = 1.484; UFF_angles["Mo3+6"] = 109.47; UFF_Z["Mo3+6"] = 3.4;
UFF_bonds["Tc6+5"] = 1.322; UFF_angles["Tc6+5"] = 90; UFF_Z["Tc6+5"] = 3.4;
UFF_bonds["Ru6+2"] = 1.478; UFF_angles["Ru6+2"] = 90; UFF_Z["Ru6+2"] = 3.4;
UFF_bonds["Rh6+3"] = 1.332; UFF_angles["Rh6+3"] = 90; UFF_Z["Rh6+3"] = 3.508;
UFF_bonds["Pd4+2"] = 1.338; UFF_angles["Pd4+2"] = 90; UFF_Z["Pd4+2"] = 3.21;
UFF_bonds["Ag1+1"] = 1.386; UFF_angles["Ag1+1"] = 180; UFF_Z["Ag1+1"] = 1.956;
UFF_bonds["Cd3+2"] = 1.403; UFF_angles["Cd3+2"] = 109.47; UFF_Z["Cd3+2"] = 1.65;
UFF_bonds["In3+3"] = 1.459; UFF_angles["In3+3"] = 109.47; UFF_Z["In3+3"] = 2.07;
UFF_bonds["Sn3"] = 1.398; UFF_angles["Sn3"] = 109.47; UFF_Z["Sn3"] = 2.961;
UFF_bonds["Sb3+3"] = 1.407; UFF_angles["Sb3+3"] = 91.6; UFF_Z["Sb3+3"] = 2.704;
UFF_bonds["Te3+2"] = 1.386; UFF_angles["Te3+2"] = 90.25; UFF_Z["Te3+2"] = 2.882;
UFF_bonds["I_"] = 1.382; UFF_angles["I_"] = 180; UFF_Z["I_"] = 2.65;
UFF_bonds["Xe4+4"] = 1.267; UFF_angles["Xe4+4"] = 90; UFF_Z["Xe4+4"] = 0.556;
UFF_bonds["Cs"] = 2.57; UFF_angles["Cs"] = 180; UFF_Z["Cs"] = 1.573;
UFF_bonds["Ba6+2"] = 2.277; UFF_angles["Ba6+2"] = 90; UFF_Z["Ba6+2"] = 2.727;
UFF_bonds["La3+3"] = 1.943; UFF_angles["La3+3"] = 109.47; UFF_Z["La3+3"] = 3.3;
UFF_bonds["Ce6+3"] = 1.841; UFF_angles["Ce6+3"] = 90; UFF_Z["Ce6+3"] = 3.3;
UFF_bonds["Pr6+3"] = 1.823; UFF_angles["Pr6+3"] = 90; UFF_Z["Pr6+3"] = 3.3;
UFF_bonds["Nd6+3"] = 1.816; UFF_angles["Nd6+3"] = 90; UFF_Z["Nd6+3"] = 3.3;
UFF_bonds["Pm6+3"] = 1.801; UFF_angles["Pm6+3"] = 90; UFF_Z["Pm6+3"] = 3.3;
UFF_bonds["Sm6+3"] = 1.78; UFF_angles["Sm6+3"] = 90; UFF_Z["Sm6+3"] = 3.3;
UFF_bonds["Eu6+3"] = 1.771; UFF_angles["Eu6+3"] = 90; UFF_Z["Eu6+3"] = 3.3;
UFF_bonds["Gd6+3"] = 1.735; UFF_angles["Gd6+3"] = 90; UFF_Z["Gd6+3"] = 3.3;
UFF_bonds["Tb6+3"] = 1.732; UFF_angles["Tb6+3"] = 90; UFF_Z["Tb6+3"] = 3.3;
UFF_bonds["Dy6+3"] = 1.71; UFF_angles["Dy6+3"] = 90; UFF_Z["Dy6+3"] = 3.3;
UFF_bonds["Ho6+3"] = 1.696; UFF_angles["Ho6+3"] = 90; UFF_Z["Ho6+3"] = 3.416;
UFF_bonds["Er6+3"] = 1.673; UFF_angles["Er6+3"] = 90; UFF_Z["Er6+3"] = 3.3;
UFF_bonds["Tm6+3"] = 1.66; UFF_angles["Tm6+3"] = 90; UFF_Z["Tm6+3"] = 3.3;
UFF_bonds["Yb6+3"] = 1.637; UFF_angles["Yb6+3"] = 90; UFF_Z["Yb6+3"] = 2.618;
UFF_bonds["Lu6+3"] = 1.671; UFF_angles["Lu6+3"] = 90; UFF_Z["Lu6+3"] = 3.271;
UFF_bonds["Hf3+4"] = 1.611; UFF_angles["Hf3+4"] = 109.47; UFF_Z["Hf3+4"] = 3.921;
UFF_bonds["Ta3+5"] = 1.511; UFF_angles["Ta3+5"] = 109.47; UFF_Z["Ta3+5"] = 4.075;
UFF_bonds["W_6+6"] = 1.392; UFF_angles["W_6+6"] = 90; UFF_Z["W_6+6"] = 3.7;
UFF_bonds["W_3+4"] = 1.526; UFF_angles["W_3+4"] = 109.47; UFF_Z["W_3+4"] = 3.7;
UFF_bonds["W_3+6"] = 1.38; UFF_angles["W_3+6"] = 109.47; UFF_Z["W_3+6"] = 3.7;
UFF_bonds["Re6+5"] = 1.372; UFF_angles["Re6+5"] = 90; UFF_Z["Re6+5"] = 3.7;
UFF_bonds["Re3+7"] = 1.314; UFF_angles["Re3+7"] = 109.47; UFF_Z["Re3+7"] = 3.7;
UFF_bonds["Os6+6"] = 1.372; UFF_angles["Os6+6"] = 90; UFF_Z["Os6+6"] = 3.7;
UFF_bonds["Ir6+3"] = 1.371; UFF_angles["Ir6+3"] = 90; UFF_Z["Ir6+3"] = 3.731;
UFF_bonds["Pt4+2"] = 1.364; UFF_angles["Pt4+2"] = 90; UFF_Z["Pt4+2"] = 3.382;
UFF_bonds["Au4+3"] = 1.262; UFF_angles["Au4+3"] = 90; UFF_Z["Au4+3"] = 2.625;
UFF_bonds["Hg1+2"] = 1.34; UFF_angles["Hg1+2"] = 180; UFF_Z["Hg1+2"] = 1.75;
UFF_bonds["Tl3+3"] = 1.518; UFF_angles["Tl3+3"] = 120; UFF_Z["Tl3+3"] = 2.068;
UFF_bonds["Pb3"] = 1.459; UFF_angles["Pb3"] = 109.47; UFF_Z["Pb3"] = 2.846;
UFF_bonds["Bi3+3"] = 1.512; UFF_angles["Bi3+3"] = 90; UFF_Z["Bi3+3"] = 2.47;
UFF_bonds["Po3+2"] = 1.5; UFF_angles["Po3+2"] = 90; UFF_Z["Po3+2"] = 2.33;
UFF_bonds["At"] = 1.545; UFF_angles["At"] = 180; UFF_Z["At"] = 2.24;
UFF_bonds["Rn4+4"] = 1.42; UFF_angles["Rn4+4"] = 90; UFF_Z["Rn4+4"] = 0.583;
UFF_bonds["Fr"] = 2.88; UFF_angles["Fr"] = 180; UFF_Z["Fr"] = 1.847;
UFF_bonds["Ra6+2"] = 2.512; UFF_angles["Ra6+2"] = 90; UFF_Z["Ra6+2"] = 2.92;
UFF_bonds["Ac6+3"] = 1.983; UFF_angles["Ac6+3"] = 90; UFF_Z["Ac6+3"] = 3.9;
UFF_bonds["Th6+4"] = 1.721; UFF_angles["Th6+4"] = 90; UFF_Z["Th6+4"] = 4.202;
UFF_bonds["Pa6+4"] = 1.711; UFF_angles["Pa6+4"] = 90; UFF_Z["Pa6+4"] = 3.9;
UFF_bonds["U_6+4"] = 1.684; UFF_angles["U_6+4"] = 90; UFF_Z["U_6+4"] = 3.9;
UFF_bonds["Np6+4"] = 1.666; UFF_angles["Np6+4"] = 90; UFF_Z["Np6+4"] = 3.9;
UFF_bonds["Pu6+4"] = 1.657; UFF_angles["Pu6+4"] = 90; UFF_Z["Pu6+4"] = 3.9;
UFF_bonds["Am6+4"] = 1.66; UFF_angles["Am6+4"] = 90; UFF_Z["Am6+4"] = 3.9;
UFF_bonds["Cm6+3"] = 1.801; UFF_angles["Cm6+3"] = 90; UFF_Z["Cm6+3"] = 3.9;
UFF_bonds["Bk6+3"] = 1.761; UFF_angles["Bk6+3"] = 90; UFF_Z["Bk6+3"] = 3.9;
UFF_bonds["Cf6+3"] = 1.75; UFF_angles["Cf6+3"] = 90; UFF_Z["Cf6+3"] = 3.9;
UFF_bonds["Es6+3"] = 1.724; UFF_angles["Es6+3"] = 90; UFF_Z["Es6+3"] = 3.9;
UFF_bonds["Fm6+3"] = 1.712; UFF_angles["Fm6+3"] = 90; UFF_Z["Fm6+3"] = 3.9;
UFF_bonds["Md6+3"] = 1.689; UFF_angles["Md6+3"] = 90; UFF_Z["Md6+3"] = 3.9;
UFF_bonds["No6+3"] = 1.679; UFF_angles["No6+3"] = 90; UFF_Z["No6+3"] = 3.9;
UFF_bonds["Lw6+3"] = 1.698; UFF_angles["Lw6+3"] = 90; UFF_Z["Lw6+3"] = 3.9;

UFF_torsions["C_3"] = 2.119;
UFF_torsions["N_3"] = 0.45;
UFF_torsions["O_3"] = 0.018;
UFF_torsions["Si3"] = 1.225;
UFF_torsions["P_3"] = 2.4;
UFF_torsions["S_3"] = 0.484;
UFF_torsions["Ge3"] = 0.701;
UFF_torsions["As3"] = 1.5;
UFF_torsions["Se3"] = 0.335;
UFF_torsions["Sn3"] = 0.199;
UFF_torsions["Sb3"] = 1.1;
UFF_torsions["Te3"] = 0.3;
UFF_torsions["Pb3"] = 0.1;
UFF_torsions["Bi3"] = 1;
UFF_torsions["Po3"] = 0.3;

// separate paper
// http://pubs.acs.org/doi/pdf/10.1021/j100161a070
// Rappe et al. The Journal of Physical Chemistry, Vol. 95, No. 8, 1991
UFF_electroneg["Li"] = 3.006;
UFF_electroneg["C"] = 5.343;
UFF_electroneg["N"] = 6.899;
UFF_electroneg["O"] = 8.741;
UFF_electroneg["F"] = 10.874;
UFF_electroneg["Na"] = 2.843;
UFF_electroneg["Si"] = 4.168;
UFF_electroneg["P"] = 5.463;
UFF_electroneg["S"] = 6.928;
UFF_electroneg["Cl"] = 8.564;
UFF_electroneg["K"] = 2.421;
UFF_electroneg["Br"] = 7.79;
UFF_electroneg["Rb"] = 2.331;
UFF_electroneg["I"] = 6.822;
UFF_electroneg["Cs"] = 2.183;
UFF_electroneg["H"] = 4.528;
// extrapolated from ptable.com ratio with (F electroneg=3.98)
// for MOFs and such.
UFF_electroneg["Zn"] = 4.508065;
UFF_electroneg["Co"] = 5.136462;
UFF_electroneg["Fe"] = 4.999854;
UFF_electroneg["Cu"] = 5.191106;
UFF_electroneg["Ni"] = 5.218427;
UFF_electroneg["Mg"] = 3.579131;
UFF_electroneg["Mn"] = 4.234849;
UFF_electroneg["Cr"] = 4.535387;
UFF_electroneg["Ti"] = 4.207528;
UFF_electroneg["Ca"] = 2.732161;
UFF_electroneg["V"] = 4.453422;
UFF_electroneg["Al"] = 4.398779;
UFF_electroneg["In"] = 4.863246;
UFF_electroneg["Ru"] = 6.010754;
UFF_electroneg["Cd"] = 4.617352;


// UFF4MOF additions
// Heine et al. J. Chem. Theory Comput. 2014, 10, 880−891
UFF_bonds["O_3_f"] = 0.634;
UFF_bonds["O_2_z"] = 0.528;
UFF_bonds["Al6+3"] = 1.22;
UFF_bonds["Sc6+3"] = 1.44;
UFF_bonds["Ti4+2"] = 1.38;
UFF_bonds["V_4+2"] = 1.18;
UFF_bonds["V_6+3"] = 1.3;
UFF_bonds["Cr4+2"] = 1.1;
UFF_bonds["Cr6f3"] = 1.28;
UFF_bonds["Mn6+3"] = 1.34;
UFF_bonds["Mn4+2"] = 1.26;
UFF_bonds["Fe6+3"] = 1.32;
UFF_bonds["Fe4+2"] = 1.1;
UFF_bonds["Co3+2"] = 1.24;
UFF_bonds["Co4+2"] = 1.16;
UFF_bonds["Cu4+2"] = 1.28;
UFF_bonds["Zn4+2"] = 1.34;
UFF_bonds["Zn3f2"] = 1.24;

UFF_angles["O_3_f"] = 109.47;
UFF_angles["O_2_z"] = 120;
UFF_angles["Al6+3"] = 90;
UFF_angles["Sc6+3"] = 90;
UFF_angles["Ti4+2"] = 90;
UFF_angles["V_4+2"] = 90;
UFF_angles["V_6+3"] = 90;
UFF_angles["Cr4+2"] = 90;
UFF_angles["Cr6f3"] = 90;
UFF_angles["Mn6+3"] = 90;
UFF_angles["Mn4+2"] = 90;
UFF_angles["Fe6+3"] = 90;
UFF_angles["Fe4+2"] = 90;
UFF_angles["Co3+2"] = 109.47;
UFF_angles["Co4+2"] = 90;
UFF_angles["Cu4+2"] = 90;
UFF_angles["Zn4+2"] = 90;
UFF_angles["Zn3f2"] = 109.47;

UFF_Z["O_3_f"] = 2.3;
UFF_Z["O_2_z"] = 2.3;
UFF_Z["Al6+3"] = 1.792;
UFF_Z["Sc6+3"] = 2.595;
UFF_Z["Ti4+2"] = 2.659;
UFF_Z["V_4+2"] = 2.679;
UFF_Z["V_6+3"] = 2.679;
UFF_Z["Cr4+2"] = 2.463;
UFF_Z["Cr6f3"] = 2.463;
UFF_Z["Mn6+3"] = 2.43;
UFF_Z["Mn4+2"] = 2.43;
UFF_Z["Fe6+3"] = 2.43;
UFF_Z["Fe4+2"] = 2.43;
UFF_Z["Co3+2"] = 1.308;
UFF_Z["Co4+2"] = 1.308;
UFF_Z["Cu4+2"] = 2.43;
UFF_Z["Zn4+2"] = 1.308;
UFF_Z["Zn3f2"] = 1.308;

}
