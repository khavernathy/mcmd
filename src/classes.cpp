#include <stdio.h>
#include <math.h>
#include <string>
#include <map>
#include <vector>
#include <stdint.h>

using namespace std;

enum {
    POTENTIAL_LJ,
    POTENTIAL_LJES,
    POTENTIAL_LJESPOLAR,
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

// Constants is sort-of a misnomer for some things in this class but you get the idea.
class Constants {
	public:
		Constants();
		double e,kb,kbk,fs,cC,keSI,ke,eV,cM,cA,cJ,NA,cV,R,mpmc2uff,uff2mpmc,
            ATM2REDUCED,kg2em, E2REDUCED, TORQUE2REDUCED, FORCE2REDUCED,
            DEBYE2SKA, JL2ATM, A32L, K2KJMOL; // all defined below.
		string jobname="default_jobname";
        string mode; // "mc" or "md"
        int_fast8_t checkpoints_option=0; // enables checkpoints for debuggin
        int_fast8_t ensemble;
        string ensemble_str;
        //int_fast8_t movetype;
        string atom_file; // input atoms .pdb file
        string output_traj="traj.xyz"; // system trajectory in xyz
        string output_traj_pdb="traj.pdb"; // system trajectory in pdb
        string restart_pdb="restart.pdb"; // a file to re-kick an unfinished run
        string thermo_output="thermo.dat"; // a file for all thermodynamic info
        int_fast8_t potential_form = POTENTIAL_LJ; // "lj", "ljes", "ljespolar", "phast2" models for potential
        //string xyz_traj_option="off"; // default no xyz traj because we'll use the PDB traj instead.
        vector<string> sorbate_name; // e.g. h2_bssp, h2_bss, co2*, co2, co2_trappe, c2h2, etc.
        vector<double> sorbate_fugacity; // holds the fugacities for multi-sorbate gases.
        int_fast8_t pdb_traj_option=1; // option to write PDB trajectory (in addition to xyz). default on
        int_fast8_t xyz_traj_option=0; // option for xyz trajectory, default off
        int_fast8_t com_option=0; // enables computation of center-of-mass and logs in output_traj
        int_fast8_t rotate_option=1; // MC ONLY: can deactivate rotates if wanted.
        int_fast8_t draw_box_option=1; // option to draw the box for visualization in restart.pdb
        int_fast8_t rd_lrc=1; // long range corrections for LJ RD
        int_fast8_t ewald_es=1; // ewald method for electrostatic potential calculation.
        int_fast8_t pdb_long=0; // on would force long coordinate/charge output
        int_fast8_t dist_within_option=0; // a function to calculate atom distances within a certain radius of origin
        string dist_within_target; // the atom to find in above option
        double dist_within_radius; // the radius within which to search from origin
        int_fast8_t autocenter = 1; // center all atoms about origin automatically. can opt out of it.
        int_fast8_t simulated_annealing = 0; // sim. ann.
        double sa_target = 0.0; // target temperature for annealing.
        double sa_schedule = 0.9999; // T-change factor for annealing.

        map <string,double> masses; // mass database for defaults.
		map <string,double> sigs; // LJ sigmas database for defaults. NOT r_m (as in UFF). Defined below
		map <string,double> eps; // LJ epsions database for defaults. Defined lated
		map <string,double> phast2_c6; map <string,double> phast2_c8; map <string,double> phast2_c10; map <string,double> phast2_sigs; map <string,double> phast2_eps; map <string,double> phast2_polar; // phast2 defaults
		map <string,double> polars; // polarizability database for defaults. Defined below
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
        double rotate_angle_factor=360; // 0 -> this number to rotate if rotate selected
		int stepsize=1; // obvi
        int finalstep; // user defined for MC
        int  mc_corrtime=1000; // default 1k cuz I used that a lot for mpmc research
        int_fast8_t mc_pbc=1; // PBC in monte carlo, default on
        int currentprotoid=0; // for getting fugacity for the boltzmann factor.

        // MD STUFF
        int  md_corrtime=50; // user defined for MD
        int_fast8_t md_pbc=1; // PBC in molecular dynamics. default on.
        int_fast8_t md_rotations=1; // MD only.
        double md_init_vel=99999.99; // placeholder value. Will be overwritten. A / fs. User can set. Will be random +- up to this num.
        double md_vel_goal=0; // 1D vector component of the velocity (which has magnitude md_init_vel)
        double md_dt=0.1, md_ft=10000; // MD timestep and final time, in fs
        int_fast8_t md_mode = MD_MOLECULAR; // default is to keep molecules rigid (bonded)
		double md_thermostat_freq = 0.001; // a value used to calculate probability of a heat-bath collision with molecule i; Frenkel uses 0.01 and 0.001 as examples; but no matter what, a boltzmann distribution is generated
        double md_thermostat_probab = md_thermostat_freq * exp(-md_thermostat_freq * md_dt);

        map <string,double> sig_override;
        map <string,double> eps_override; // feature for overriding preset LJ params (for developing LJ models). 0.0 are defaults which will be overwritten if option is used. sig=A; eps=K

    	int total_atoms=0; // actual sites, not "atoms" persay
                int initial_sorbates=0.0; // for safekeeping to calculate chemical potential in uVT
        double initial_energy=0.0; // "" ""

        // Ewald (for ES)
        double ewald_alpha; // =3.5/cutoff;
        double ewald_kmax = 7;

        // Wolf (for polarization)
        double polar_wolf_alpha = 0.13;
        double polar_damp = 2.1304;
        double polar_gamma = 1.03;
        int polar_max_iter = 4;
        double **A_matrix, **B_matrix, C_matrix[3][3];
        int polar_precision=0;
        int iter_success;
        int_fast8_t polar_rrms =0;
        double dipole_rrms = 0.0;
        int_fast8_t polar_gs_ranked = 1;
        int_fast8_t polar_gs = 0;
        int_fast8_t polar_palmo = 1;
        int_fast8_t polar_pbc = 1; // default periodic polar
        //int thole_total_atoms = 0;

        int all_pbc=1;

};

class Pbc {
    public:
        Pbc();

        double x_length, y_length, z_length; // for moving molecules back into box on other side, computing in calcBoxVertices
        double x_min,x_max,y_min,y_max,z_min,z_max;
        double basis[3][3];
        double reciprocal_basis[3][3];
		double cutoff=0.;
        double volume, inverse_volume, old_volume;
        double a, b, c, alpha, beta, gamma;
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
            printf("basis1 %.5f %.5f %.5f\n", basis[0][0], basis[0][1], basis[0][2]);
            printf("basis2 %.5f %.5f %.5f\n", basis[1][0], basis[1][1], basis[1][2]);
            printf("basis3 %.5f %.5f %.5f\n", basis[2][0], basis[2][1], basis[2][2]);
            printf("Basis vectors: { a = %9.5f; b = %9.5f; c = %9.5f }\n", a, b,c);
            printf("Basis angles:  { α = %9.5f; β = %9.5f; γ = %9.5f }\n", alpha,beta,gamma);
            printf("Box vertices ::\n");
            for (int n=0; n<8; n++)
                printf("   -> %i : %9.5f %9.5f %9.5f\n", n, box_vertices[n][0], box_vertices[n][1], box_vertices[n][2]);
            printf("PBC Cutoff = %.5f\n", cutoff);
            for (int n=0; n<6; n++) {
                printf("Plane %i equation :: %.5fx + %.5fy + %.5fz + %.5f = 0\n",
                    n, A[n], B[n], C[n], D[n]);
            }
            printf("x_length = %.5f; y_length = %.5f; z_length = %.5f\n", x_length, y_length, z_length);
            if (alpha == 90 && beta == 90 && gamma == 90) { // these don't get calculated or used for weird (non 90-90-90) unit cells
              printf("x_max = %.5f; y_max = %.5f; z_max = %.5f\n", x_max, y_max, z_max);
              printf("x_min = %.5f; y_min = %.5f; z_min = %.5f\n", x_min, y_min, z_min);
            }
        }

        void calcPlane(int p1index, int p2index, int p3index, int planeIndex) { // 3 points define a plane.
            double vector1[3], vector2[3];

            // 1) get 3 points (indexes for box vertices provided in arguments)
            // 2) make 2 planar vectors AB, AC
            for (int n=0; n<3; n++) {
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
             2 /------------/ 6     0 :   0123 plane (-x)
              /|           /|       1 :   4567 plane (+x)
             / |          / |       2 :   0145 plane (-y)
          3 |------------|7 |       3 :   2367 plane (+y)
            |  |         |  |       4 :   0246 plane (-z)
            |  |---------|--| 4     5 :   1357 plane (+z)
            | / 0        |  /
            |/           | /        The vertices are defined in box_vertices[8][3].
          1 |____________|/ 5       3 points define a plane, so I'll use the first 3 for the above planes

            */

            // 3 points and plane index
            calcPlane(0,1,2,0); // I'm not sure if I'll ever have a system where the 4th point
            calcPlane(4,5,6,1); // results in a different plane than the first 3...
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
            double MAXVALUE = 1e40; int MAX_VECT_COEF = 5;
			int i, j, k, p;
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
		    int i,j,k,p,q,count=0;
		    int box_labels[2][2][2];
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

                        for (int n=0; n<3; n++)
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
            for (int n=4; n<=7; n++) if (box_vertices[n][0] > maxx) maxx = box_vertices[n][0];
            for (int n=0; n<=3; n++) if (box_vertices[n][0] < minx) minx = box_vertices[n][0];

            for (int n=2; n<=3; n++) if (box_vertices[n][1] > maxy) maxy=box_vertices[n][1];
            for (int n=6; n<=7; n++) if (box_vertices[n][1] > maxy) maxy=box_vertices[n][1];
            for (int n=0;n<=1; n++) if (box_vertices[n][1] < miny) miny = box_vertices[n][1];
            for (int n=4;n<=5; n++) if (box_vertices[n][1] < miny) miny = box_vertices[n][1];

            for (int n=1; n<=7; n+=2) if (box_vertices[n][2] > maxz) maxz = box_vertices[n][2];
            for (int n=0; n<=6; n+=2) if (box_vertices[n][2] < minz) minz = box_vertices[n][2];

            lengthx = maxx-minx;
            lengthy = maxy-miny;
            lengthz = maxz-minz;

        }

};

Pbc::Pbc() {}

class Stats {
	public:
		Stats();

        int MCstep=0, MCcorrtime_iter; // keeps track of steps and coortimes for averages.
        bool MCmoveAccepted;
        double MCeffRsq; // for calculating Monte Carlo efficiency, roughly, based on successful displaces

        int_fast8_t radial_dist = 0; // default is no radial distribution
        string radial_file = "radial_distribution.dat"; // default filename for output.
        double radial_bin_size = 0.1; // bin counts will be considered for this range in A
        double radial_max_dist = 10.0; // maximum r to consider in rad. dist.
        vector<long unsigned int> radial_bins; // holds the counters
        string radial_centroid, radial_counterpart; // the two atoms to get distance between, user def.

		double insert_bf_sum = 0; // insert boltzmanns added up in boltzmann.cpp
		double remove_bf_sum = 0; // ...
		double displace_bf_sum = 0;
		double volume_change_bf_sum = 0;

		int insert_accepts = 0; int insert_attempts=0;// Counters for successful moves. uVT only
		int remove_accepts = 0; int remove_attempts=0;// uVT only
		int displace_accepts = 0; int displace_attempts=0;
		int volume_change_accepts = 0; int volume_attempts=0;// NPT only
		int total_accepts=0; int total_attempts=0;

        double ar_tot=0, ar_ins=0, ar_rem=0, ar_dis=0, ar_vol=0; // BF acceptance ratios
        double ins_perc=0, rem_perc=0, dis_perc=0, vol_perc=0; // BF percentage of moves

        double bf_avg, ibf_avg, rbf_avg, dbf_avg, vbf_avg;

        int count_movables = 0; // this is the SORBATE MOLECULES, e.g. 2H2 means 2, not 4
        int count_frozens = 0; // frozen ATOMS, not molecules (which is normally just 1)
        int count_frozen_molecules=0; // frozen MOLECULES; normally 1

        double polar_iterations=0;

        struct obs_t {
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
                //sd = sqrt( ((counter-2.0)*prevsd + (x - average)*(x - prevavg) ) / (counter-1.0));
                double operand =  prevsd*prevsd + prevavg*prevavg - average*average +((x*x - prevsd*prevsd - prevavg*prevavg)/counter);
                (operand > 0) ? sd = sqrt( operand ) : sd = 0;

                //if (name == "es") {
                //if (name == "es" || name == "es_self" || name == "es_real" || name == "es_recip") {
                //printf("observable %14s :: counter = %5f; value = %-10.5f; prevavg = %-10.5f; average = %-10.5f; prevsd = %-10.5f; sd = %-10.5f\n",
                //    name.c_str(), counter, value, prevavg, average, prevsd, sd);
                //}
            }


        } Nsq,NU,qst,qst_nvt,rd,es,polar,potential,volume,z,
            lj_lrc,lj_self_lrc,lj,es_self,es_real,es_recip,chempot,totalmass,
            frozenmass, pressure,temperature, fdotrsum, dist_within, csp, diffusion;

        vector<obs_t> wtp = vector<obs_t>(10);
        vector<obs_t> wtpME = vector<obs_t>(10);
        vector<obs_t> Nmov = vector<obs_t>(10);
        vector<obs_t> movablemass = vector<obs_t>(10);
        vector<obs_t> density = vector<obs_t>(10); // gotta have multiples for multi-sorbate.
        vector<obs_t> selectivity = vector<obs_t>(10);

};

Stats::Stats() {}

// contains information about an individual pair marked by i,j,k,l
class Pair {
    public:
        Pair();
        double r; double prev_r;
        double d[3]; double prev_d[3];
        int recalculate = 1; // a flag to queue recalculation of energies
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
        double Nsq,NU,qst,qst_nvt,rd,es,polar,potential,volume,z,
            lj_lrc,lj_self_lrc,lj,es_self,es_real,es_recip,chempot,totalmass,
            frozenmass,pressure,temperature, fdotrsum, dist_within, csp, diffusion;

        int total_atoms, thole_total_atoms;

        vector<double> wtp = vector<double>(10);
        vector<double> wtpME = vector<double>(10);
        vector<double> Nmov = vector<double>(10);
        vector<double> movablemass = vector<double>(10);
        vector<double> density = vector<double>(10); // for multi-sorbate
        vector<double> selectivity = vector<double>(10);
};

Last::Last() {}

class Atom {
	public:
		Atom();
        string name; // element or label, e.g. H or H2G
        string mol_name; // molecule name that the atom belongs to
        int_fast8_t frozen; // movable/frozen (0 or 1)
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
        double rank_metric;  // for polarization sorting

        double pos[3] = {0,0,0};
		double prevpos[3] = {0,0,0};
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
        double dipole_rrms=0;
		/*vector<double> force = vector<double>(3); // force, K / A
*/
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
            printf("====================\natom (PDBID %i) %s on molecule %s (PBDID %i) frozen= %i \n -----> m = %e; eps = %f; sig = %f; C = %f\nforce: %f %f %f; \nacc: %f %f %f; \nold_acc: %f %f %f; \nvel: %f %f %f\n", PDBID, name.c_str(), mol_name.c_str(), mol_PDBID, frozen, m, eps, sig,C, force[0], force[1], force[2], acc[0], acc[1], acc[2], old_acc[0], old_acc[1], old_acc[2], vel[0], vel[1], vel[2]);
        }
};

Atom::Atom() {}


class Molecule {
	public:
		Molecule();
		vector<Atom> atoms; // vector that holds this molecule's atoms
		int PDBID; // the molecule's PDBID (from input)
		string name; // the molecule name/label (from input), e.g. H2 or MOF
        int_fast8_t frozen; //0 or 1
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
        double d_theta[3] = {0,0,0};
        //vector<double> com = vector<double>(3); // center of mass for molecule. Using for MD rotations
        double mass=0.0;
        double inertia=0.0; //moment of inertia. stored in K fs^2
        double fugacity;

        void reInitialize() {
            // if there are no atoms, don't bother
            if (atoms.size() > 0) {
            while (!atoms.empty()) atoms.pop_back();
            mass=0;
            inertia=0;
            for (int n=0; n<3; n++) {
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
                d_theta[n]=0;
            }
            name = "";
            PDBID=0;
            frozen = 0; // movable
            }
        }

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
            double cap = 0.0000015;
            for (int n=0; n<3; n++) {
                ang_vel[n] = ang_vel[n] + 0.5*(ang_acc[n] * old_ang_acc[n])*dt;
                if (ang_vel[n] > cap) ang_vel[n] = cap;
                else if (ang_vel[n] < -cap) ang_vel[n] = -cap;
            }
        }

        // linear velocity
        void calc_vel(double dt, double goal) {
            //double booster=0.0005; // for NVT thermostat.
            for (int n=0; n<3; n++) {
                vel[n] = vel[n] + 0.5*(acc[n] + old_acc[n])*dt; // in A/fs. vel. verlet
            }
            /*
           double vmag = sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
            for (int n=0; n<3; n++) {
                if (vel[n] < 0) {  // THE NVT THERMOSTAT :: it changes velocities to try to approach the initial v
                    if (vel[n] < -goal) vel[n] += booster;
                    else if (vel[n] > -goal) vel[n] -= booster;
                } else if (vel[n] > 0) {
                    if (vel[n] < goal) vel[n] += booster;
                    else if (vel[n] > goal) vel[n] -= booster;
                }
            }
            */
        }

        // angular position // in rad
        void calc_ang_pos(double dt) {
            double theta[3];
            //double cap = 0.0005;
            for (int n=0; n<3; n++) {
                theta[n] = ang_pos[n];
                ang_pos[n] = ang_pos[n] + ang_vel[n] * dt + 0.5*ang_acc[n] * dt * dt;
                //if (ang_pos[n] > cap) ang_pos[n] = cap; // SET THE ROTATION CAP -- rad/fs
                //else if (ang_pos[n] < -cap) ang_pos[n] = -cap;
                d_theta[n] = ang_pos[n] - theta[n];
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
            printf("====================\nmolecule PDBID=%i :: mass: %e; inertia: %e; frozen = %i; \nforce: %f %f %f; \nacc: %f %f %f; \nold_acc: %f %f %f; \nvel: %f %f %f; \ncom: %f %f %f; \ntorque: %f %f %f \nang_acc: %f %f %f \nold_ang_acc: %f %f %f \nang_vel: %f %f %f; \nang_pos: %f %f %f (in degrees) \n",
            PDBID,mass,inertia,frozen,
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
    JL2ATM = 0.00986923297; // J/L to atm
    A32L = 1e-27; // A^3 to liters.
    K2KJMOL = kb*NA/1000; // K -> kJ/mol

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
//mg
//al
	masses["Si"] = 28.085*cM;
	masses["P"] = 30.973*cM;
	masses["S"] = 32.06*cM;
    masses["SS"] = 32.065*cM; /// S SAPT
	masses["Cl"] = 35.45*cM;
	masses["Ar"] = 39.948*cM;
    masses["ArS"] = 39.948*cM; // Ar SAPT
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
    sigs["HS"] = 3.09061; // H SAPT (Adam)
    sigs["He"] = 2.362*uff2mpmc;
	sigs["Li"] = 2.451*uff2mpmc;
	sigs["Be"] = 2.745*uff2mpmc;
	sigs["B"] = 4.083*uff2mpmc;
	sigs["C"] = 3.851*uff2mpmc; // Rappe = 3.851; I changed to 1.3 for MD MOF5
    sigs["C_s"] = 3.5786; // C SAPT
    sigs["C_p"] = 3.63956;
    sigs["C_t"] = 3.37707;
    sigs["C_a"] = 3.65947;
    sigs["C_en"] = 3.74892;
    sigs["C_yn"] = 3.79838;
    sigs["C_ony"] = 3.56023;
	sigs["N"] = 3.66*uff2mpmc;
    sigs["N_sp3"] = 3.32588; // N SAPT
    sigs["N_sp2"] = 3.45133;
	sigs["O"] = 3.5*uff2mpmc; // Rappe=3.5; I changed to 1.5 for MD water model (SPC = 3.166); 1.3 for MOF5
	sigs["O2"] = 4.0*uff2mpmc; // my O2 model // roughly accurate for densities -> 10atm
	sigs["OW"] = 1.5; // O in water (my model) -- old (sticky MD) 1.4
    sigs["OT"] = 3.15061; // O in TIP3P
    sigs["O_sp3"] = 3.15611; // O SAPT
    sigs["O_sp2"] = 3.34161;
    sigs["F"] = 3.364*uff2mpmc;
	sigs["Ne"] = 3.243*uff2mpmc;
    sigs["Cu"] = 3.495*uff2mpmc;
	sigs["Na"] = 2.983*uff2mpmc;
// mg al ...
    sigs["SS"] = 3.8899; // S SAPT
	sigs["Cl"] = 3.947*uff2mpmc;
	sigs["Ar"] = 3.868*uff2mpmc;
    sigs["ArS"] = 3.37191; // Ar SAPT
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
	eps["HB"] = 34.20; // buch model h2
    eps["H2G"] = 8.8516; // bss model h2
    eps["H2E"] = 0.0; // bss
    eps["H2N"] = 4.0659; // bss
	eps["H"] = 0.044/kbk; // rappe verbatim
	eps["HW"] = 0.044/kbk; // H in water (my model)
    eps["HT"] = 0.0; // H in TIP3P
    eps["HS"] = 0.66563; // H SAPT
    eps["He"] = 0.056/kbk;
	eps["Li"] = 0.025/kbk;
	eps["Be"] = 0.085/kbk;
	eps["B"] = 0.18/kbk;
	eps["C"] = 0.105/kbk;
    eps["C_p"] = 36.692; // C SAPT
    eps["C_s"] = 31.35824;
    eps["C_t"] = 41.45435;
    eps["C_a"] = 22.30908;
    eps["C_en"] = 26.88878;
    eps["C_yn"] = 22.40343;
    eps["C_ony"] = 18.09254;
	eps["N"] = 0.069/kbk;
    eps["N_sp3"] = 36.97995; // N SAPT
    eps["N_sp2"] = 24.25732;
	eps["O"] = 0.06/kbk;
	eps["O2"] = 0.06/kbk; // my O2 model
    eps["OW"] = 0.06/kbk; // O in water (my model)
    eps["OT"] = 0.6364/0.0083144621; // O in TIP3P
	eps["O_sp3"] = 30.01345; // O SAPT
    eps["O_sp2"] = 21.81177;
    eps["F"] = 0.05/kbk;
	eps["Ne"] = 0.042/kbk;
    eps["Na"] = 0.03/kbk;
//mg al si ...
    eps["SS"] = 53.02994; // S SAPT
	eps["Cl"] = 0.227/kbk;
	eps["Ar"] = 0.185/kbk;
    eps["ArS"] = 128.32680; // Ar SAPT
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
    polars["HS"] = 0.41380; // H SAPT
    polars["B"] = 0.6634;//*cV/ke;
	polars["C"] = 1.2866;//*cV/ke;
    polars["C_p"] = polars["C_s"] = polars["C_t"] = polars["C_a"] = polars["C_en"] = polars["C_yn"] = polars["C_ony"] = 1.2866; // C SAPT
	polars["N"] = 0.97157;//*cV/ke;
    polars["N_sp3"] = polars["N_sp2"] = 0.97157; // N SAPT
	polars["O"] = 0.852;//*cV/ke;
    polars["OW"] = 0.852;//*cV/ke; // O in water (my model)
    polars["O_sp3"] = polars["O_sp2"] = 0.852; // O SAPT
	polars["Na"] = 24.11;//*cV/ke; // from paper https://www.researchgate.net/publication/45896756_Absolute_and_ratio_measurements_of_the_polarizability_of_Na_K_and_Rb_with_an_atom_interferometer
	polars["P"] = 3.35;//*cV/ke;
    polars["SS"] = 2.474; // S SAPT
	polars["Cl"] = 2.40028;//*cV/ke;
    polars["ArS"] = 1.63922; // Ar SAPT
	polars["Cu"] = 2.19630;//*cV/ke;
	polars["Zn"] = 1.98870;//*cV/ke;
	//polars["Br"]
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
}
