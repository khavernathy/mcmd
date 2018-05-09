/* Douglas M Franz
 * Space group, USF, 2017
 *
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>

using namespace std;

string convertElement(System &system, string label) {
    string mystring = label;
    string my_sub_string;

    //printf("input %s\n",mystring.c_str());
    for(int n=0; n<10; n++) {
            //printf("searching for %i\n",n);
            string my_sub_string = to_string(n);

            std::size_t found = mystring.find(my_sub_string);

            if(found != std::string::npos) {
                    std::cout << found << endl;
    	//	if((int)found == 1) mystring = mystring.substr(0,1);
                    mystring = mystring.substr(0,(int)found);
    		//printf("mystring after cut %s\n",mystring.c_str());
            }
    }
  //  printf("before check 2  %s\n",mystring.c_str());
    if((int)(mystring.length()) == 2) {
            string first = mystring.substr(0,1);
            string second = mystring.substr(1,1);
            std::transform(second.begin(), second.end(), second.begin(), ::tolower);
            mystring = first + second;
    }
  //  printf("final %s\n",mystring.c_str());
    return mystring;
}

bool find_cycle(System &system, unsigned int mol, unsigned int i) {

    unsigned int b1,b2,b3,b4,b5,b6;
    unsigned int a2,a3,a4,a5,a6,a7;

    // see if atom i on molecule mol is on a 6-membered ring
    for (b1=0; b1<system.molecules[mol].atoms[i].bonds.size(); b1++) {
        a2 = system.molecules[mol].atoms[i].bonds[b1];
        if (system.molecules[mol].atoms[a2].name == "Zn") continue; // avoid false ring in MOF-5 type
        for (b2=0; b2<system.molecules[mol].atoms[a2].bonds.size(); b2++) {
            a3 = system.molecules[mol].atoms[a2].bonds[b2];
            if (a3 == i) continue; // don't go backwards..
            for (b3=0; b3<system.molecules[mol].atoms[a3].bonds.size();b3++) {
                a4 = system.molecules[mol].atoms[a3].bonds[b3];
                if (a4 == a2) continue;
                for (b4=0; b4<system.molecules[mol].atoms[a4].bonds.size();b4++) {
                    a5 = system.molecules[mol].atoms[a4].bonds[b4];
                    if (a5 == a3) continue;
                    for (b5=0; b5<system.molecules[mol].atoms[a5].bonds.size();b5++) {
                        a6 = system.molecules[mol].atoms[a5].bonds[b5];
                        if (a6 == a4) continue;
                        for (b6=0; b6<system.molecules[mol].atoms[a6].bonds.size();b6++) {
                            a7 = system.molecules[mol].atoms[a6].bonds[b6];
                            if (a7 == a5) continue;
                            if (a7 == i)
                                return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}


// function to determine UFF atom-type based on
// name of element and number of bonds
string getUFFlabel(System &system, string name, int num_bonds, int mol, int i) {
    name = convertElement(system,name); // convert weird names like C12 to C
    // starting with just the organic-y atom-types.
    if (name == "H") {
        return "H_"; // assume it's never H_b (borane hydrogen)
    } else if (name == "B") {
        if (num_bonds == 3) return "B_2";
        else if (num_bonds == 4) return "B_3";
    } else if (name == "C") {
        if (find_cycle(system, mol, i) && num_bonds != 4) return "C_R";
        else if (num_bonds == 2) return "C_1";
        else if (num_bonds == 3) return "C_2";
        else if (num_bonds == 4) return "C_3";
        // need to dynamically account for resonant C_R too...
    } else if (name == "N") {
        if (find_cycle(system,mol,i) && num_bonds != 4) return "N_R";
        else if (num_bonds == 1) return "N_1";
        else if (num_bonds == 2) return "N_3";
        else if (num_bonds == 3) return "N_2";
        else if (num_bonds == 4) return "N_3";
    } else if (name == "O") {
        int ZnCount=0;
        for (unsigned int z=0; z<system.molecules[mol].atoms[i].bonds.size(); z++) {
            if (system.molecules[mol].atoms[system.molecules[mol].atoms[i].bonds[z]].name == "Zn") {
                ZnCount++;
            }
        }
        if (ZnCount == 1) return "O_2";
        else if (ZnCount == 4) return "O_3_f"; // MOF-5 Zn cluster center O
        else if (find_cycle(system,mol,i) && num_bonds != 4) return "O_R";
        else if (num_bonds == 1) return "O_1";
        else if (num_bonds == 2) return "O_3";
        else if (num_bonds == 3) return "O_2";
        else if (num_bonds == 4) return "O_3";
        else return "O_3";
    } else if (name == "F") {
        return "F_";
    } else if (name == "Al") {
        return "Al6+3"; // UFF4MOF
    } else if (name == "P") {
        return "P_3+3";
        // + other weird geometries
    } else if (name == "S") {
        if (find_cycle(system,mol,i)) return "S_R";
        else return "S_2";
        // + other weird geometries...
    } else if (name == "Cl") {
        return "Cl";
    } else if (name == "Sc") {
        return "Sc6+3"; // UFF4MOF
    } else if (name == "Ti") {
        return "Ti4+2"; // UFF4MOF
    } else if (name == "V") {
        return "V_4+2"; // UFF4MOF
    } else if (name == "Cr") {
        return "Cr4+2"; // UFF4MOF
    } else if (name == "Co") {
        return "Co6+3";
    } else if (name == "Cu") {
        return "Cu4+2"; // UFF4MOF default
    } else if (name == "Zn") {
        //return "Zn4+2"; // UFF4MOF paddlewheel type
        return "Zn3f2"; // this is the MOF-5 type of Zn
    } else if (name == "Br") {
        return "Br";
    } else if (name == "I") {
        return "I_";
    } else if (name == "Ru") {
        return "Ru6+2";
    } else if (name == "Cd") {
        return "Cd1f1";
    } else if (name == "Si") {
        return "SiF6"; // default sifsix Si
    }


    return "NOTFOUND";
}


// return true if this should count as a bond.
bool qualify_bond(System &system, double r, unsigned int mol, unsigned int i, unsigned int j) {
    string a1 = system.molecules[mol].atoms[i].name;
    string a2 = system.molecules[mol].atoms[j].name;

    double bondlength = system.constants.bondlength;

    if (a1=="H" && a2=="H" && system.molecules[mol].atoms.size() != 2)
        return false;
    else if ((a1=="H" || a2=="H") && r > 1.3)
        return false;
    else if ((a1=="Zn" || a2=="Zn") && r <= 2.1) // Zn4O group
        return true;
    else if ((a1=="Cu" || a2=="Cu") && r <= 2.1) // Cu paddlewheel - O
        return true;
    else if ((a1=="Cu" && a2=="Cu") && r <= 2.9) // Cu paddlewheel
        return true;
    else if ((a1=="Ru" || a2=="Ru") && r <= 2.1) // Ru +2  complexes
        return true;
    else if ((a1=="Co" || a2=="Co") && r <= 2.25) // Co complexes
        return true;
    else if (((a1=="Cd" && a2=="F") || (a1=="F" && a2=="Cd")) && r <= 2.4) // Cd SIFSIX type
        return true;
    else if (r > bondlength)
        return false;

    else return true;
}

/*
 * Essentially everything below comes from the parameters and
 * functional forms prescribed in UFF, via
 * Rappe et al.
 * J. Am. Chem. Soc. Vol. 114, No. 25, 1992
 * */

// provide rij parameter given needed information (atom IDs)
double get_rij(System &system, unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
    const double ri = system.constants.UFF_bonds[system.molecules[i].atoms[j].UFFlabel.c_str()];
    const double rj = system.constants.UFF_bonds[system.molecules[k].atoms[l].UFFlabel.c_str()];
    const double Xi = system.constants.UFF_electroneg[system.molecules[i].atoms[j].name.c_str()];
    const double Xj = system.constants.UFF_electroneg[system.molecules[k].atoms[l].name.c_str()];
    const double BO = 1.0; // assume all single bonds for now.
    const double lambda = 0.1332; // Bond-order correction parameter

    const double rBO = -lambda*(ri + rj)*log(BO);
    const double Xij = sqrt(Xi) - sqrt(Xj);
    const double rEN = ri*rj*(Xij*Xij)/(Xi*ri + Xj*rj);
    return ri + rj + rBO + rEN; // in Angstroms
}

// get Force constant kij for a bond.
double get_kij(System &system, unsigned int i, unsigned int j, unsigned int k, unsigned int l, double rij) {
    const double Zi = system.constants.UFF_Z[system.molecules[i].atoms[j].UFFlabel.c_str()];
    const double Zj = system.constants.UFF_Z[system.molecules[k].atoms[l].UFFlabel.c_str()];
    return 664.12*Zi*Zj/(rij*rij*rij); // in kcal/molA^2, as-is
}

double get_BO(string a1, string a2) {
    if ((a1.find("_R") != std::string::npos && a2.find("_R") != std::string::npos)
     || (a1.find("_R") != std::string::npos && a2.find("_2") != std::string::npos)
     || (a1.find("_2") != std::string::npos && a2.find("_R") != std::string::npos))
       return 1.5; // resonant, carboxyl, etc.
    else if (a1.find("_2") != std::string::npos && a2.find("_2") != std::string::npos)
        return 2.0; // sp2
    else if (a1.find("_1") != std::string::npos && a2.find("_1") != std::string::npos)
        return 3.0; // sp
    else if ((a1.find("O_") != std::string::npos && a2.find("Zn") != std::string::npos)
         || (a1.find("Zn") != std::string::npos && a2.find("O_") != std::string::npos))
        return 0.5; // Zn--O cluster
    else if ((a1.find("Cu") != std::string::npos && a2.find("O_") != std::string::npos)
         ||  (a1.find("O_") != std::string::npos && a2.find("Cu") != std::string::npos))
        return 0.5; // Cu--O paddlewheel
    else if (a1.find("Cu") != std::string::npos && a2.find("Cu") != std::string::npos)
        return 0.25; // Cu paddlewheel

    else return 1.0; // sp3, other single bonds.
}

// get the total potential from bond stretches
// via the Morse potential
double stretch_energy(System &system) {
    double potential = 0;
    double alpha,Dij,rij; // kij; bond params
    double mainterm; // main chunk of potential, to be squared..
    unsigned int i,j,l; // atom indices
    double r; // actual, current distance for pair.
    /* ============================ */
    //double BO;
    /* ============================ */

    // loop through bonds of this atom.
    for (unsigned int it=0; it<system.constants.uniqueBonds.size(); it++) {
        i = system.constants.uniqueBonds[it].mol;
        j = system.constants.uniqueBonds[it].atom1;
        l = system.constants.uniqueBonds[it].atom2;

        rij = system.constants.uniqueBonds[it].rij; // in A
        //kij = system.constants.uniqueBonds[it].kij; // in kcal/molA^2
//          printf("rij = %f; kij = %f\n", rij, kij);

    //    BO = system.constants.uniqueBonds[it].BO;
        Dij = system.constants.uniqueBonds[it].Dij;// in kcal/mol
        alpha = system.constants.uniqueBonds[it].alpha; // in 1/A

        double* distances = getDistanceXYZ(system, i,j,i,l);
        r = distances[3];
        mainterm = exp(-alpha*(r-rij)) - 1.0; // unitless
        if (mainterm==0) continue; // skip 0-energy
        potential += Dij*(mainterm*mainterm); // in kcal/mol
        //printf("bond %i %i energy = %f\n", j,l, Dij*(mainterm*mainterm));
    }

    potential /= system.constants.kbk; // kcal/mol -> K
    system.stats.Ustretch.value = potential;
    return potential; // in kcal/mol
}

// get angle-force parameter for IJK triplet
double get_Kijk(System &system, double rik, double t1, double t2, double t3) {
    return  t1/(rik*rik*rik*rik*rik)*(3.0*t3*(1-t2*t2) - rik*rik*t2);
}

// get the angle ABC where B is center atom, on molecule i
double get_angle(System &system, unsigned int i, unsigned int A, unsigned int B, unsigned int C) {
    // https://stackoverflow.com/questions/19729831/angle-between-3-points-in-3d-space
    double AB[3] = {0,0,0};
    double BC[3] = {0,0,0};

    double* ABdistances = getDistanceXYZ(system, i, A, i, B);
    for (unsigned int n=0;n<3;n++) {
        AB[n] = ABdistances[n];
    }
    double* BCdistances = getDistanceXYZ(system, i, C, i, B);
    for (unsigned int n=0;n<3;n++) {
        BC[n] = BCdistances[n];
    }

    const double dotprod = dddotprod(AB,BC);
    const double ABm = sqrt(dddotprod(AB,AB));
    const double BCm = sqrt(dddotprod(BC,BC));

    double arg = dotprod/(ABm*BCm);
    if (arg > 1.0) arg=1.0;
    else if (arg < -1.0) arg=-1.0;

    return acos(arg); // returns in radians
}

// get r_ik, a parameter for angle bends, ** different from r_ij (and r_jk) **
double get_rik(System &system, double rij, double rjk, double t3, double angle) {
    // angle is in radians
    return sqrt(rij*rij + rjk*rjk - 2*t3*cos(angle));
}

// get the total potential from angle bends
// via simple Fourier small cosine expansion
double angle_bend_energy(System &system) {
    double potential=0;
    //const double deg2rad = M_PI/180.0;
    unsigned int i,j,l,m;
    double rij, rjk, rik, K_ijk, C0, C1, C2; //theta_ijk; // angle-bend params
    double t1,t2,t3;
    double angle; // the actual angle IJK

    for (unsigned int it=0; it<system.constants.uniqueAngles.size(); it++) {
        i = system.constants.uniqueAngles[it].mol;
        j = system.constants.uniqueAngles[it].atom1;
        l = system.constants.uniqueAngles[it].atom2;
        m = system.constants.uniqueAngles[it].atom3;
        rij = system.constants.uniqueAngles[it].rij; // in Angstroms
        rjk = system.constants.uniqueAngles[it].rjk;
        //theta_ijk = system.constants.uniqueAngles[it].theta_ijk; // in rads
        C2 = system.constants.uniqueAngles[it].C2; // 1/rad^2
        C1 = system.constants.uniqueAngles[it].C1; // 1/rad
        C0 = system.constants.uniqueAngles[it].C0; // 1
        t1 = system.constants.uniqueAngles[it].t1;
        t2 = system.constants.uniqueAngles[it].t2;
        t3 = system.constants.uniqueAngles[it].t3;


        angle = get_angle(system, i, j, l, m);
        rik = get_rik(system, rij, rjk, t3, angle); // r_ik (A-C) is computed differently than r_ij (A-B) and r_jk (B-C)
        K_ijk = get_Kijk(system, rik, t1,t2,t3);
        if (K_ijk==0) continue; // skip 0-energy

        potential += K_ijk*(C0 + C1*cos(angle) + C2*cos(2.0*angle)); // in kcal/mol
        //printf("angle potential of %i %i %i = %f\n", j,l,m,K_ijk*(C0 + C1*cos(angle) + C2*cos(2.0*angle)));
    }
    potential /= system.constants.kbk; // kcal/mol -> K
    system.stats.Uangles.value = potential;
    return potential; // in kcal/mol
} // end angle bend energy


double get_dihedral_angle(System &system, unsigned int mol, unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
    // current dihedral angle for i--j--k--l atoms
    /*
     *          L    <- 0/180? degree dihedral (in plane of screen)..
     *         /
     *    J---K
     *   /
     *  I
     *
     * */
    // 4 components of a plane are A,B,C,D in Ax + By + Cz + D = 0
    // The D parameter is irrelevant here so we just need a vector holding A,B,C
    double Plane1[3]; // build from atoms i,j,k
    double Plane2[3]; // build from atoms j,k,l
    double v1a[3], v1b[3]; // for plane 1
    double v2a[3], v2b[3]; // for plane 2
    double norm1[3], norm2[3];

	/*
    double xi,xj,xk,xl;
    double yi,yj,yk,yl;
    double zi,zj,zk,zl;
	*/

    double* distancesJI = getDistanceXYZ(system,mol,j,mol,i);
    for (int n=0; n<3; n++)
        v1a[n] = distancesJI[n];

    double* distancesKI = getDistanceXYZ(system,mol,k,mol,i);
    for (int n=0;n<3;n++)
        v1b[n] = distancesKI[n];

    double* distancesKJ = getDistanceXYZ(system,mol,k,mol,j);
    for (int n=0;n<3;n++)
        v2a[n] = distancesKJ[n];

    double* distancesLJ = getDistanceXYZ(system,mol,l,mol,j);
    for (int n=0;n<3;n++)
        v2b[n] = distancesLJ[n];

/*
        v1a[n] = system.molecules[mol].atoms[j].pos[n] - system.molecules[mol].atoms[i].pos[n];
        v1b[n] = system.molecules[mol].atoms[k].pos[n] - system.molecules[mol].atoms[i].pos[n];

        v2a[n] = system.molecules[mol].atoms[k].pos[n] - system.molecules[mol].atoms[j].pos[n];
        v2b[n] = system.molecules[mol].atoms[l].pos[n] - system.molecules[mol].atoms[j].pos[n];
    */

    norm1[0] = v1a[1]*v1b[2] - v1a[2]*v1b[1];
    norm1[1] = v1a[2]*v1b[0] - v1a[0]*v1b[2];
    norm1[2] = v1a[0]*v1b[1] - v1a[1]*v1b[0];

    norm2[0] = v2a[1]*v2b[2] - v2a[2]*v2b[1];
    norm2[1] = v2a[2]*v2b[0] - v2a[0]*v2b[2];
    norm2[2] = v2a[0]*v2b[1] - v2a[1]*v2b[0];

    for (int n=0;n<3;n++) {
        Plane1[n] = norm1[n];
        Plane2[n] = norm2[n];
    }

    // both planes done; now get angle
    const double dotplanes = dddotprod(Plane1,Plane2);
    const double mag1 = sqrt(dddotprod(Plane1,Plane1));
    const double mag2 = sqrt(dddotprod(Plane2,Plane2));

    double arg = dotplanes/(mag1*mag2);
    if (arg > 1.0) arg = 1.0;
    else if (arg < -1.0) arg = -1.0;

    return acos(arg);
}

double * get_torsion_params(System &system, string a1, string a2) {
    static double o[3] = {0,0,0};
    // output 0 -- equilibrium angle
    // output 1 -- V_jk (energy param)
    // output 2 -- n (periodicity)

    // based on the shared bond of the torsion
    // sp3--sp3
    if (a1.find("_3") != std::string::npos && a2.find("_3") != std::string::npos) {
        if (a1.find("O") != std::string::npos && a2.find("O") != std::string::npos) {
            o[0] = 90;
            o[1] = 2.0; //sqrt(2.*2.);
            o[2] = 2;
        }
        else if (a1.find("O") != std::string::npos &&
            (a2.find("S_") != std::string::npos ||
             a2.find("Se") != std::string::npos ||
             a2.find("Te") != std::string::npos ||
             a2.find("Po") != std::string::npos)  )
        {
            o[0] = 90;
            o[1] = sqrt(2.0 * 6.8);
            o[2] = 2;

        }
        else if (a2.find("O") != std::string::npos &&
            (a1.find("S_") != std::string::npos ||
             a1.find("Se") != std::string::npos ||
             a1.find("Te") != std::string::npos ||
             a1.find("Po") != std::string::npos)  )
        {
            o[0] = 90;
            o[1] = 6.8; // sqrt(6.8*6.8)
            o[2] = 2;
        }
        else {
            o[0] = 60; // "or 180"
            o[1] = sqrt(system.constants.UFF_torsions[a1] * system.constants.UFF_torsions[a2]);
            o[2] = 3;
        }
    }
    // sp3--sp2
    else if ((a1.find("_3") != std::string::npos && a2.find("_2") != std::string::npos) ||
            (a1.find("_2") != std::string::npos && a2.find("_3") != std::string::npos) ||
            (a1.find("_R") != std::string::npos && a2.find("_3") != std::string::npos) ||
            (a1.find("_3") != std::string::npos && a2.find("_R") != std::string::npos)     ) {
        o[0]=0;
        o[1]=1.0;
        o[2] = 6;
    }
    // sp2--sp2
    else if ((a1.find("_2") != std::string::npos && a2.find("_2") != std::string::npos) ||
            (a1.find("_R") != std::string::npos && a2.find("_R") != std::string::npos) ||
            (a1.find("_2") != std::string::npos && a2.find("_R") != std::string::npos) ||
            (a1.find("_R") != std::string::npos && a2.find("_2") != std::string::npos) ) {
        o[0]=180; // "or 60"
        const double BO = 1.5; // assume resonance bond order..
        double Uj=1.25,Uk=1.25; // assume second period...
        o[1] = 5.0*sqrt(Uj*Uk)*(1.0 + 4.18*log(BO));
        o[2] = 2;
    }
    // sp--sp
    else if (a1.find("_1") != std::string::npos && a2.find("_1") != std::string::npos) {
        o[0]=180;
        o[1]=0;
        o[2] = 0;
    }
    else {
        // main group garbage
        o[0] = 0;
        o[1] = 0;
        o[2] = 0;
    }

    return o;
}

// get the total potential from torsions
// via simple Fourier small cosine expansion
double torsions_energy(System &system) {
    double potential=0;
    double vjk, n, dihedral, phi_ijkl; // vj,vk; // n is periodicity (integer quantity)
    //const double deg2rad = M_PI/180.0;
    unsigned int i,j,l,m,p; // molecule i, atoms (j,l,m and p)
    for (unsigned int it=0; it<system.constants.uniqueDihedrals.size(); it++) {

        vjk = system.constants.uniqueDihedrals[it].vjk;
        if (vjk==0) continue; // skip  0-energy
        i = system.constants.uniqueDihedrals[it].mol;
        j = system.constants.uniqueDihedrals[it].atom1;
        l = system.constants.uniqueDihedrals[it].atom2;
        m = system.constants.uniqueDihedrals[it].atom3;
        p = system.constants.uniqueDihedrals[it].atom4;


        phi_ijkl = system.constants.uniqueDihedrals[it].phi_ijkl;
        n = system.constants.uniqueDihedrals[it].n;

        dihedral = get_dihedral_angle(system, i, j,l,m,p);
        potential += 0.5*vjk*(1.0 - cos(n*phi_ijkl)*cos(n*dihedral));//0.5*vjk;
        //printf("dihedral %i %i %i %i = %f; phi_goal = %f; actual_phi = %f\n", j,l,m,p, 0.5*vjk*(1.0 - cos(n*phi_ijkl)*cos(n*dihedral)), phi_ijkl, dihedral);
    }

    potential /= system.constants.kbk; // kcal/mol -> K
    system.stats.Udihedrals.value = potential;
    return potential; // in kcal/mol
}


// Morse potential gradient for all bonded atoms, to minimize energy
double morse_gradient(System &system) {
    // x,y,z is the point of interest for this gradient
    // ij tells us whether it's the 1st atom or 2nd within definition of delta x (+1 or -1)
    double alpha,Dij,rij; //kij; // bond params
    unsigned int i,j,l; // atom indices
    double r; // actual, current distance for pair.
    /* ============================ */
    //double BO;
    /* ============================ */
    double prefactor,grad,delta;
        // typical gradient elements (e.g. dE/dx_i) are ~10^2 in these units.


            for (unsigned int it=0; it<system.constants.uniqueBonds.size(); it++) {

                i = system.constants.uniqueBonds[it].mol;
                j = system.constants.uniqueBonds[it].atom1;
                l = system.constants.uniqueBonds[it].atom2;

                rij = system.constants.uniqueBonds[it].rij; // in A
                //kij = system.constants.uniqueBonds[it].kij; // in kcal/molA^2
//          printf("rij = %f; kij = %f\n", rij, kij);

                //BO = system.constants.uniqueBonds[it].BO;
                Dij = system.constants.uniqueBonds[it].Dij;// in kcal/mol
                alpha = system.constants.uniqueBonds[it].alpha; // in 1/A


                double* distances = getDistanceXYZ(system, i,j,i,l);
                r = distances[3];

                // gradient for a single bond is 6D (3D on each atom, 1 for each D.O.F.)
                prefactor = 2*alpha*Dij*exp(alpha*(rij-r))/r;
                if (prefactor==0) continue; // skip 0-contributions
                for (int n=0;n<3;n++) {
                    delta = distances[n]; //system.molecules[i].atoms[j].pos[n] - system.molecules[i].atoms[l].pos[n];
                    grad = prefactor * delta;
                    grad *= (1 - exp(alpha*(rij-r)));
                        system.molecules[i].atoms[j].force[n] -= grad / system.constants.kbk;
                        system.molecules[i].atoms[l].force[n] += grad / system.constants.kbk;
                }
                // xj, yj, zj
                // since gradient of the other atom is just minus the other, we apply a Newton-pair style thing above
                // instead of recomputing.
                /*
                    for (int n=0;n<3;n++) {
                    delta = system.molecules[i].atoms[l].pos[n] - system.molecules[i].atoms[j].pos[n];
                    grad = prefactor * delta;
                    grad *= (1 - exp(alpha*(rij-r)));
                    printf("%f\n", grad);
                    // move the atom position element in direction of energy minimum
                    system.molecules[i].atoms[l].pos[n] += grad*move_factor;
                }
                */

            }

    return 0; //.5*potential; // in kcal/mol
}

double simple_r(double xi, double xj, double yi, double yj, double zi, double zj) {
    return sqrt((xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj));
}

// get the total potential from angle bends
// via simple Fourier small cosine expansion
double angle_bend_gradient(System &system) {
//    const double deg2rad = M_PI/180.0;
    unsigned int i,j,l,m;
    double rij, rjk, rik, K_ijk, C1, C2;// theta_ijk; // C0; angle-bend params
    double angle; // the actual angle IJK
    double grad;
    double xi, yi, zi, xj, yj, zj, xk, yk, zk;
    double t1,t2,t3;
        // cos-derivative terms in the gradient, + other terms.


    for (unsigned int it=0; it<system.constants.uniqueAngles.size(); it++) {
        i = system.constants.uniqueAngles[it].mol;
        j = system.constants.uniqueAngles[it].atom1;
        l = system.constants.uniqueAngles[it].atom2;
        m = system.constants.uniqueAngles[it].atom3;
        rij = system.constants.uniqueAngles[it].rij; // in Angstroms
        rjk = system.constants.uniqueAngles[it].rjk;
        //theta_ijk = system.constants.uniqueAngles[it].theta_ijk; // in rads
        C2 = system.constants.uniqueAngles[it].C2; // 1/rad^2
        C1 = system.constants.uniqueAngles[it].C1; // 1/rad
        //C0 = system.constants.uniqueAngles[it].C0; // 1
        t1 = system.constants.uniqueAngles[it].t1;
        t2 = system.constants.uniqueAngles[it].t2;
        t3 = system.constants.uniqueAngles[it].t3;


        angle = get_angle(system, i, j, l, m);
        //printf("Angle %i %i %i = %f; real angle = %f\n", j,l,m,theta_ijk/deg2rad, angle/deg2rad);
        rik = get_rik(system, rij, rjk, t3,angle); // r_ik (A-C) is computed differently than r_ij (A-B) and r_jk (B-C)
        K_ijk = get_Kijk(system, rik, t1,t2,t3);
        if (K_ijk==0) continue; // skip 0-contrib

        // compute the gradient for all angle components (xyz for 3 atoms = 9)
        // gradient is [dE/dx_i...] which ends up as a sum of two cosine derivs (d/dx C0 term -> 0)

        // based on atom i, we need to make periodic "ghosts" j and k.
        // this is way easier than trying to deal with periodic distances
        // inside the derivative..
        xi = system.molecules[i].atoms[j].pos[0];
        yi = system.molecules[i].atoms[j].pos[1];
        zi = system.molecules[i].atoms[j].pos[2];

        double* distances_ij = getDistanceXYZ(system,i,j,i,l);
        xj = xi - distances_ij[0];
        yj = yi - distances_ij[1];
        zj = zi - distances_ij[2];

        double* distances_ik = getDistanceXYZ(system,i,j,i,m);
        xk = xi - distances_ik[0];
        yk = yi - distances_ik[1];
        zk = zi - distances_ik[2];
        // done with ghosts

        /*
         * old manual crap
        B = (yi-yj)*(yk-yj) + (zi-zj)*(zk-zj);
        C = (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj);

        dij = simple_r(xi,xj,yi,yj,zi,zj);
        djk = simple_r(xj,xk,yj,yk,zj,zk);

        cos1 = ((B*(xj-xi) + C*(xk-xj))/(djk*dij*dij*dij));
        cos2 = 2*(((xk-xj)/(djk*dij)) - (( (xi-xj)*(B + (xk-xj)*(xi-xj)))/(djk*dij*dij*dij))) * sin(2*acos((B + (xk-xj)*(xi-xj))/(djk*dij))) / sqrt(1 - pow((B + (xk-xj)*(xi-xj)),2)/(djk*djk*dij*dij));

        grad = K_ijk*(C1*cos1 + C2*cos2);
  */
        // MATLAB GENERATED PARTIALS...
        // recall J,L,M are i,j,k
        double ij2sum = (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) + (zi-zj)*(zi-zj);
        double jk2sum = (xj-xk)*(xj-xk) + (yj-yk)*(yj-yk) + (zj-zk)*(zj-zk);
        double rootij = sqrt(ij2sum);
        double rootjk = sqrt(jk2sum);
        double mixer = ((xi-xj)*(xj-xk)+(yi-yj)*(yj-yk)+(zi-zj)*(zj-zk));
        double thing = -mixer*mixer/(ij2sum*jk2sum)+1.0;
        double arg = 1.0/rootij*1.0/rootjk*mixer;
        if (arg>1.0) arg=1.0;
        else if (arg<-1.0) arg=-1.0;

        if (thing > 0) {
        // x_i
        grad =
        K_ijk*(C2*sin(acos(arg)*2.0)*
        ((xj-xk)*1.0/rootij*1.0/rootjk-(xi*2.0-xj*2.0)*1.0/pow(ij2sum,3.0/2.0)*1.0/rootjk*mixer*0.5)*
        1.0/sqrt(thing)*2.0-
        C1*(xj-xk)*1.0/rootij*1.0/rootjk+
        C1*(xi*2.0-xj*2.0)*1.0/pow(ij2sum,3.0/2.0)*1.0/rootjk*mixer*0.5);
         system.molecules[i].atoms[j].force[0] -= grad / system.constants.kbk;

        // y_i
        grad =
        K_ijk*(C2*sin(acos(arg)*2.0)*
        ((yj-yk)*1.0/rootij*1.0/rootjk-(yi*2.0-yj*2.0)*1.0/pow(ij2sum,3.0/2.0)*1.0/rootjk*mixer*0.5)*
        1.0/sqrt(thing)*2.0-
        C1*(yj-yk)*1.0/rootij*1.0/rootjk+
        C1*(yi*2.0-yj*2.0)*1.0/pow(ij2sum,3.0/2.0)*1.0/rootjk*mixer*0.5);
         system.molecules[i].atoms[j].force[1] -= grad / system.constants.kbk;

        // z_i
        grad =
        K_ijk*(C2*sin(acos(arg)*2.0)*
        ((zj-zk)*1.0/rootij*1.0/rootjk-(zi*2.0-zj*2.0)*1.0/pow(ij2sum,3.0/2.0)*1.0/rootjk*mixer*0.5)*
        1.0/sqrt(thing)*2.0-
        C1*(zj-zk)*1.0/rootij*1.0/rootjk+
        C1*(zi*2.0-zj*2.0)*1.0/pow(ij2sum,3.0/2.0)*1.0/rootjk*mixer*0.5);
         system.molecules[i].atoms[j].force[2] -= grad / system.constants.kbk;

        // x_j
        grad =
        K_ijk*(C2*sin(acos(arg)*2.0)*
        1.0/sqrt(thing)*((xi-xj*2.0+xk)*1.0/rootij*
        1.0/rootjk+(xi*2.0-xj*2.0)*1.0/pow(ij2sum,3.0/2.0)*1.0/rootjk*mixer*0.5-(xj*2.0-xk*2.0)*
        1.0/rootij*1.0/pow(jk2sum,3.0/2.0)*mixer*0.5)*2.0-
        C1*(xi-xj*2.0+xk)*1.0/rootij*1.0/rootjk-
        C1*(xi*2.0-xj*2.0)*1.0/pow(ij2sum,3.0/2.0)*1.0/rootjk*mixer*0.5+
        C1*(xj*2.0-xk*2.0)*1.0/rootij*1.0/pow(jk2sum,3.0/2.0)*mixer*0.5);
         system.molecules[i].atoms[l].force[0] -= grad / system.constants.kbk;

        // y_j
        grad =
        K_ijk*(C2*sin(acos(arg)*2.0)*
        1.0/sqrt(thing)*((yi-yj*2.0+yk)*1.0/rootij*
        1.0/rootjk+(yi*2.0-yj*2.0)*1.0/pow(ij2sum,3.0/2.0)*1.0/rootjk*mixer*0.5-(yj*2.0-yk*2.0)*
        1.0/rootij*1.0/pow(jk2sum,3.0/2.0)*mixer*0.5)*2.0-
        C1*(yi-yj*2.0+yk)*1.0/rootij*1.0/rootjk-
        C1*(yi*2.0-yj*2.0)*1.0/pow(ij2sum,3.0/2.0)*1.0/rootjk*mixer*0.5+
        C1*(yj*2.0-yk*2.0)*1.0/rootij*1.0/pow(jk2sum,3.0/2.0)*mixer*0.5);
         system.molecules[i].atoms[l].force[1] -= grad / system.constants.kbk;

        // z_j
        grad =
        K_ijk*(C2*sin(acos(arg)*2.0)*
        1.0/sqrt(thing)*((zi-zj*2.0+zk)*1.0/rootij*
        1.0/rootjk+(zi*2.0-zj*2.0)*1.0/pow(ij2sum,3.0/2.0)*1.0/rootjk*mixer*0.5-(zj*2.0-zk*2.0)*
        1.0/rootij*1.0/pow(jk2sum,3.0/2.0)*mixer*0.5)*2.0-
        C1*(zi-zj*2.0+zk)*1.0/rootij*1.0/rootjk-
        C1*(zi*2.0-zj*2.0)*1.0/pow(ij2sum,3.0/2.0)*1.0/rootjk*mixer*0.5+
        C1*(zj*2.0-zk*2.0)*1.0/rootij*1.0/pow(jk2sum,3.0/2.0)*mixer*0.5);
         system.molecules[i].atoms[l].force[2] -= grad / system.constants.kbk;

        // x_k
        grad =
        -K_ijk*(C2*sin(acos(arg)*2.0)*
        ((xi-xj)*1.0/rootij*1.0/rootjk-(xj*2.0-xk*2.0)*1.0/rootij*
        1.0/pow(jk2sum,3.0/2.0)*mixer*0.5)*1.0/sqrt(thing)*2.0-
        C1*(xi-xj)*1.0/rootij*1.0/rootjk+
        C1*(xj*2.0-xk*2.0)*1.0/rootij*1.0/pow(jk2sum,3.0/2.0)*mixer*0.5);
         system.molecules[i].atoms[m].force[0] -= grad / system.constants.kbk;

        // y_k
        grad =
        -K_ijk*(C2*sin(acos(arg)*2.0)*
        ((yi-yj)*1.0/rootij*1.0/rootjk-(yj*2.0-yk*2.0)*1.0/rootij*
        1.0/pow(jk2sum,3.0/2.0)*mixer*0.5)*1.0/sqrt(thing)*2.0-
        C1*(yi-yj)*1.0/rootij*1.0/rootjk+
        C1*(yj*2.0-yk*2.0)*1.0/rootij*1.0/pow(jk2sum,3.0/2.0)*mixer*0.5);
         system.molecules[i].atoms[m].force[1] -= grad / system.constants.kbk;

        // z_k
        grad =
        -K_ijk*(C2*sin(acos(arg)*2.0)*
        ((zi-zj)*1.0/rootij*1.0/rootjk-(zj*2.0-zk*2.0)*1.0/rootij*
        1.0/pow(jk2sum,3.0/2.0)*mixer*0.5)*1.0/sqrt(thing)*2.0-
        C1*(zi-zj)*1.0/rootij*1.0/rootjk+
        C1*(zj*2.0-zk*2.0)*1.0/rootij*1.0/pow(jk2sum,3.0/2.0)*mixer*0.5);
         system.molecules[i].atoms[m].force[2] -= grad / system.constants.kbk;


       } // end if sqrt argument is > 0 (to avoid nan)
        //double POT=K_ijk*(C0 + C1*cos(angle) + C2*cos(2.0*angle)); // in kcal/mol

    }

    return 0.; // in kcal/mol
} // end angle bend gradient

double torsions_gradient(System &system) {
    double vjk, n, dihedral, phi_ijkl; // n is periodicity (integer quantity)
    //const double deg2rad = M_PI/180.0;
    double xi,xj,xk,xl, yi,yj,yk,yl, zi,zj,zk,zl;
    unsigned int i,j,l,m,p; // molecule i, atoms (j,l,m and p)
    double grad;
    double arg, dih_thing1, dih_thing2, dih_thing3, dih_thing4, dih_thing5;
    for (unsigned int it=0; it<system.constants.uniqueDihedrals.size(); it++) {
        vjk = system.constants.uniqueDihedrals[it].vjk;
        n = system.constants.uniqueDihedrals[it].n;
        if (vjk==0 || n==0) continue; // skip  0-gradients
        i = system.constants.uniqueDihedrals[it].mol;
        j = system.constants.uniqueDihedrals[it].atom1;
        l = system.constants.uniqueDihedrals[it].atom2;
        m = system.constants.uniqueDihedrals[it].atom3;
        p = system.constants.uniqueDihedrals[it].atom4;
        phi_ijkl = system.constants.uniqueDihedrals[it].phi_ijkl;
        dihedral = get_dihedral_angle(system, i, j,l,m,p);
        if (dihedral==0.0 || fabs(dihedral) == 180.0) continue; // this will avoid NAN/INF errors

        // there are 12 gradients for each dihedral
        // 3 for each of the 4 atoms
        // recall that i,j,k,l (conventional atoms) --> j,l,m,p here
        // based on atom i, we need to make periodic "ghosts" j, k, l.
        // this is way easier than trying to deal with periodic distances
        // inside the derivative..
        xi = system.molecules[i].atoms[j].pos[0];
        yi = system.molecules[i].atoms[j].pos[1];
        zi = system.molecules[i].atoms[j].pos[2];

        double* distances_ij = getDistanceXYZ(system,i,j,i,l);
        xj = xi - distances_ij[0];
        yj = yi - distances_ij[1];
        zj = zi - distances_ij[2];

        double* distances_ik = getDistanceXYZ(system,i,j,i,m);
        xk = xi - distances_ik[0];
        yk = yi - distances_ik[1];
        zk = zi - distances_ik[2];

        double* distances_il = getDistanceXYZ(system,i,j,i,p);
        xl = xi - distances_il[0];
        yl = yi - distances_il[1];
        zl = zi - distances_il[2];
        // done with ghosts

        // MATLAB GENERATED GRADIENTS

        dih_thing1 = -pow(((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))+((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj))*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk))+((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj))*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk)),2.0)/((pow((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj),2.0)+pow((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj),2.0)+pow((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj),2.0))*(pow((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk),2.0)+pow((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk),2.0)+pow((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk),2.0)))+1.0;

        if (dih_thing1 > 0.0) { // prevent NAN and INF
          dih_thing2 = sqrt(pow((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj),2.0)+pow((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj),2.0)+pow((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj),2.0));
          dih_thing3 = sqrt(pow((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk),2.0)+pow((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk),2.0)+pow((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk),2.0));
          dih_thing4 = pow(pow((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj),2.0)+pow((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj),2.0)+pow((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj),2.0),3.0/2.0);
          dih_thing5 = (((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))+((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj))*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk))+((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj))*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk)));

          arg = 1.0/dih_thing2*1.0/dih_thing3*dih_thing5;
          if (arg>1.0) arg=1.0;
          else if (arg<-1.0) arg=-1.0;

          // xi
          grad = n*vjk*cos(n*phi_ijkl)*sin(n*acos(arg))*
          1.0/sqrt(dih_thing1)*(((yj-yk)*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))+(zj-zk)*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk)))*
          1.0/dih_thing2*
          1.0/dih_thing3-((yj-yk)*((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))*2.0+(zj-zk)*((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj))*2.0)*1.0/dih_thing4*
          1.0/dih_thing3*dih_thing5*(1.0/2.0))*(-1.0/2.0);
           system.molecules[i].atoms[j].force[0] -= grad / system.constants.kbk;

          // yi
          grad = n*vjk*cos(n*phi_ijkl)*sin(n*acos(arg))*1.0/sqrt(dih_thing1)*(((xj-xk)*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))-(zj-zk)*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk)))*1.0/dih_thing2*1.0/dih_thing3-((xj-xk)*((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))*2.0-(zj-zk)*((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj))*2.0)*1.0/dih_thing4*1.0/dih_thing3*dih_thing5*(1.0/2.0))*(1.0/2.0);
           system.molecules[i].atoms[j].force[1] -= grad / system.constants.kbk;

          // zi
          grad = n*vjk*cos(n*phi_ijkl)*sin(n*acos(arg))*1.0/sqrt(dih_thing1)*(((xj-xk)*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk))+(yj-yk)*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk)))*1.0/dih_thing2*1.0/dih_thing3-((xj-xk)*((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj))*2.0+(yj-yk)*((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj))*2.0)*1.0/dih_thing4*1.0/dih_thing3*dih_thing5*(1.0/2.0))*(1.0/2.0);
           system.molecules[i].atoms[j].force[2] -= grad / system.constants.kbk;

          // xj
          grad = n*vjk*cos(n*phi_ijkl)*sin(n*acos(arg))*1.0/sqrt(dih_thing1)*(1.0/dih_thing2*1.0/dih_thing3*((yk-yl)*((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))-(yi-yk)*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))+(zk-zl)*((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj))-(zi-zk)*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk)))+((yi-yk)*((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))*2.0+(zi-zk)*((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj))*2.0)*1.0/dih_thing4*1.0/dih_thing3*dih_thing5*(1.0/2.0)-((yk-yl)*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))*2.0+(zk-zl)*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk))*2.0)*1.0/dih_thing2*1.0/pow(pow((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk),2.0)+pow((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk),2.0)+pow((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk),2.0),3.0/2.0)*dih_thing5*(1.0/2.0))*(-1.0/2.0);
           system.molecules[i].atoms[l].force[0] -= grad / system.constants.kbk;

          // yj
          grad = n*vjk*cos(n*phi_ijkl)*sin(n*acos(arg))*1.0/sqrt(dih_thing1)*(1.0/dih_thing2*1.0/dih_thing3*((xk-xl)*((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))-(xi-xk)*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))-(zk-zl)*((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj))+(zi-zk)*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk)))+((xi-xk)*((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))*2.0-(zi-zk)*((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj))*2.0)*1.0/dih_thing4*1.0/dih_thing3*dih_thing5*(1.0/2.0)-((xk-xl)*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))*2.0-(zk-zl)*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk))*2.0)*1.0/dih_thing2*1.0/pow(pow((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk),2.0)+pow((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk),2.0)+pow((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk),2.0),3.0/2.0)*dih_thing5*(1.0/2.0))*(1.0/2.0);
           system.molecules[i].atoms[l].force[1] -= grad / system.constants.kbk;

          // zj
          grad = n*vjk*cos(n*phi_ijkl)*sin(n*acos(arg))*1.0/sqrt(dih_thing1)*(1.0/dih_thing2*1.0/dih_thing3*((xk-xl)*((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj))-(xi-xk)*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk))+(yk-yl)*((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj))-(yi-yk)*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk)))+((xi-xk)*((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj))*2.0+(yi-yk)*((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj))*2.0)*1.0/dih_thing4*1.0/dih_thing3*dih_thing5*(1.0/2.0)-((xk-xl)*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk))*2.0+(yk-yl)*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk))*2.0)*1.0/dih_thing2*1.0/pow(pow((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk),2.0)+pow((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk),2.0)+pow((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk),2.0),3.0/2.0)*dih_thing5*(1.0/2.0))*(1.0/2.0);
           system.molecules[i].atoms[l].force[2] -= grad / system.constants.kbk;

          // xk
          grad = n*vjk*cos(n*phi_ijkl)*sin(n*acos(arg))*1.0/sqrt(dih_thing1)*(1.0/dih_thing2*1.0/dih_thing3*((yj-yl)*((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))-(yi-yj)*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))+(zj-zl)*((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj))-(zi-zj)*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk)))+((yi-yj)*((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))*2.0+(zi-zj)*((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj))*2.0)*1.0/dih_thing4*1.0/dih_thing3*dih_thing5*(1.0/2.0)-((yj-yl)*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))*2.0+(zj-zl)*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk))*2.0)*1.0/dih_thing2*1.0/pow(pow((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk),2.0)+pow((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk),2.0)+pow((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk),2.0),3.0/2.0)*dih_thing5*(1.0/2.0))*(1.0/2.0);
           system.molecules[i].atoms[m].force[0] -= grad / system.constants.kbk;

          // yk
          grad = n*vjk*cos(n*phi_ijkl)*sin(n*acos(arg))*1.0/sqrt(dih_thing1)*(1.0/dih_thing2*1.0/dih_thing3*((xj-xl)*((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))-(xi-xj)*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))-(zj-zl)*((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj))+(zi-zj)*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk)))+((xi-xj)*((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))*2.0-(zi-zj)*((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj))*2.0)*1.0/dih_thing4*1.0/dih_thing3*dih_thing5*(1.0/2.0)-((xj-xl)*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))*2.0-(zj-zl)*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk))*2.0)*1.0/dih_thing2*1.0/pow(pow((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk),2.0)+pow((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk),2.0)+pow((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk),2.0),3.0/2.0)*dih_thing5*(1.0/2.0))*(-1.0/2.0);
           system.molecules[i].atoms[m].force[1] -= grad / system.constants.kbk;

          // zk
          grad = n*vjk*cos(n*phi_ijkl)*sin(n*acos(arg))*1.0/sqrt(dih_thing1)*(1.0/dih_thing2*1.0/dih_thing3*((xj-xl)*((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj))-(xi-xj)*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk))+(yj-yl)*((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj))-(yi-yj)*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk)))+((xi-xj)*((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj))*2.0+(yi-yj)*((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj))*2.0)*1.0/dih_thing4*1.0/dih_thing3*dih_thing5*(1.0/2.0)-((xj-xl)*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk))*2.0+(yj-yl)*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk))*2.0)*1.0/dih_thing2*1.0/pow(pow((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk),2.0)+pow((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk),2.0)+pow((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk),2.0),3.0/2.0)*dih_thing5*(1.0/2.0))*(-1.0/2.0);
           system.molecules[i].atoms[m].force[2] -= grad / system.constants.kbk;

          // xl
          grad = n*vjk*cos(n*phi_ijkl)*sin(n*acos(arg))*1.0/sqrt(dih_thing1)*(((yj-yk)*((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))+(zj-zk)*((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj)))*1.0/dih_thing2*1.0/dih_thing3-((yj-yk)*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))*2.0+(zj-zk)*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk))*2.0)*1.0/dih_thing2*1.0/pow(pow((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk),2.0)+pow((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk),2.0)+pow((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk),2.0),3.0/2.0)*dih_thing5*(1.0/2.0))*(-1.0/2.0);
           system.molecules[i].atoms[p].force[0] -= grad / system.constants.kbk;

          // yl
          grad = n*vjk*cos(n*phi_ijkl)*sin(n*acos(arg))*1.0/sqrt(dih_thing1)*(((xj-xk)*((xi-xj)*(yi-yk)-(xi-xk)*(yi-yj))-(zj-zk)*((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj)))*1.0/dih_thing2*1.0/dih_thing3-((xj-xk)*((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk))*2.0-(zj-zk)*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk))*2.0)*1.0/dih_thing2*1.0/pow(pow((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk),2.0)+pow((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk),2.0)+pow((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk),2.0),3.0/2.0)*dih_thing5*(1.0/2.0))*(1.0/2.0);
           system.molecules[i].atoms[p].force[1] -= grad / system.constants.kbk;

          // zl
          grad = n*vjk*cos(n*phi_ijkl)*sin(n*acos(arg))*1.0/sqrt(dih_thing1)*(((xj-xk)*((xi-xj)*(zi-zk)-(xi-xk)*(zi-zj))+(yj-yk)*((yi-yj)*(zi-zk)-(yi-yk)*(zi-zj)))*1.0/dih_thing2*1.0/dih_thing3-((xj-xk)*((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk))*2.0+(yj-yk)*((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk))*2.0)*1.0/dih_thing2*1.0/pow(pow((xj-xk)*(yj-yl)-(xj-xl)*(yj-yk),2.0)+pow((xj-xk)*(zj-zl)-(xj-xl)*(zj-zk),2.0)+pow((yj-yk)*(zj-zl)-(yj-yl)*(zj-zk),2.0),3.0/2.0)*dih_thing5*(1.0/2.0))*(1.0/2.0);
           system.molecules[i].atoms[p].force[2] -= grad / system.constants.kbk;

          //potential += 0.5*vjk*(1.0 - cos(n*phi_ijkl)*cos(n*dihedral));//0.5*vjk;
      } // end if dih_thing1 > 0 (prevent NAN/INF)
    }
    return 0;
}


double LJ_intramolec_energy(System &system) {
    // this is a non-bonding potential, but I'm including it here as a separate function
    // for optimizations, with unique units (kcal/mol) from the MC/MD lj.cpp (Kelvin)
    // the latter of which excludes all intramolecular contributions.
    double potential=0;
    unsigned int mol,i,j;
    double eps, sig, r, sr6;
    for (unsigned int it=0; it<system.constants.uniqueLJNonBonds.size(); it++) {
        mol = system.constants.uniqueLJNonBonds[it].mol;
        i = system.constants.uniqueLJNonBonds[it].atom1;
        j = system.constants.uniqueLJNonBonds[it].atom2;
        eps = system.constants.uniqueLJNonBonds[it].eps;
        sig = system.constants.uniqueLJNonBonds[it].sig;
        if (eps==0 || sig==0) continue; // skip 0-energy

                    double* distances = getDistanceXYZ(system, mol,i,mol,j);
                    r = distances[3];
                    if (r > system.pbc.cutoff) continue;
                    //printf("pair %i %i LJ\n",i,j);
                    sr6 = sig/r;
                    sr6 *= sr6;
                    sr6 *= sr6*sr6;
                    potential += 4.0*eps*(sr6*sr6 - sr6);
                        //printf("LJ %i %i = %f\n", i,j, 4.0*eps*(sr6*sr6 - sr6));
    }

    system.stats.UintraLJ.value = potential;
    return potential;
} // LJ intramolecular potential function


double LJ_intramolec_gradient(System &system) {
    // this is a non-bonding potential, but I'm including it here as a separate function
    // for optimizations, with unique units (kcal/mol) from the MC/MD lj.cpp (Kelvin)
    // the latter of which excludes all intramolecular contributions.
    unsigned int mol,i,j;
    double eps, sig, r,rsq,r6,s6;
    double grad;
    for (unsigned int it=0; it<system.constants.uniqueLJNonBonds.size(); it++) {
        mol = system.constants.uniqueLJNonBonds[it].mol;
        i = system.constants.uniqueLJNonBonds[it].atom1;
        j = system.constants.uniqueLJNonBonds[it].atom2;
        eps = system.constants.uniqueLJNonBonds[it].eps;
        sig = system.constants.uniqueLJNonBonds[it].sig;
        if (eps==0 || sig==0) continue; // skip 0-contributions
                    double* distances = getDistanceXYZ(system, mol,i,mol,j);
                    r = distances[3];
                    if (r > system.pbc.cutoff) continue;
                    rsq= r*r;
                    r6 = rsq*rsq*rsq;
                    s6 = sig*sig;
                    s6 *= s6*s6;

                    // 6 gradients (xyz for 2 atoms)
                    for (int n=0;n<3;n++) {
                        grad = -24.0*distances[n]*eps*(2*(s6*s6)/(r6*r6*rsq) - s6/(r6*rsq));
                        system.molecules[mol].atoms[i].force[n] -= grad;
                        system.molecules[mol].atoms[j].force[n] += grad;


                    }
    }

    return 0;
} // LJ intramolecular gradient function

double ES_intramolec_energy(System &system) {
    // this is a non-bonding potential, but I'm including it here as a separate function
    double potential=0;
    unsigned int mol,i,j;
    double qq,r;
    for (unsigned int it=0; it<system.constants.uniqueChargeNonBonds.size(); it++) {
        mol = system.constants.uniqueChargeNonBonds[it].mol;
        i = system.constants.uniqueChargeNonBonds[it].atom1;
        j = system.constants.uniqueChargeNonBonds[it].atom2;
        qq = system.constants.uniqueChargeNonBonds[it].chargeprod;

        if (qq==0) continue; // skip 0-energy

                    double* distances = getDistanceXYZ(system, mol,i,mol,j);
                    r = distances[3];
                    if (r > system.pbc.cutoff) continue;
                        potential += qq/r;
    }

    system.stats.UintraES.value = potential;
    return potential;
} // LJ intramolecular potential function

double ES_intramolec_gradient(System &system) {
    // this is a non-bonding potential, but I'm including it here as a separate function
    unsigned int mol,i,j;
    double qq,r;
    for (unsigned int it=0; it<system.constants.uniqueChargeNonBonds.size(); it++) {
        mol = system.constants.uniqueChargeNonBonds[it].mol;
        i = system.constants.uniqueChargeNonBonds[it].atom1;
        j = system.constants.uniqueChargeNonBonds[it].atom2;
        qq = system.constants.uniqueChargeNonBonds[it].chargeprod;

        if (qq==0) continue; // skip 0-force

                    double* distances = getDistanceXYZ(system, mol,i,mol,j);
                    r = distances[3];
                    if (r > system.pbc.cutoff) continue;
                    for (int n=0;n<3;n++) {
                        system.molecules[mol].atoms[i].force[n] += distances[n]*qq/(r*r*r);
                        system.molecules[mol].atoms[j].force[n] -= distances[n]*qq/(r*r*r);
                    }
    }

    return 0; //potential*system.constants.kbk; // to kcal/mol
} // LJ intramolecular potential function


// function to find all bonds (and angles) for all atoms.
void findBonds(System &system) {

    unsigned int i,j,l,m,p; // i=mol, j,l,m,p are atoms (conventionally IJKL)
    double r, ra, rh; // bond r, angle bond ra, dihedral bond rh
    unsigned int local_bonds=0;
    unsigned int duplicateFlag=0; //bond dupes
    unsigned int duplicateAngleFlag=0; // angle dupes
    unsigned int duplicateDihFlag=0; // dihedral dupes
    unsigned int molecule_limit = 1;
    if (system.constants.md_mode == MD_FLEXIBLE) molecule_limit = system.molecules.size();

    for (i=0; i<molecule_limit; i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            local_bonds = 0;
            // for each atom, we find its bonded neighbors by a distance search
            // (by only searching atoms on this same molecule)
            for (l=0; l<system.molecules[i].atoms.size(); l++) {
               if (j==l) continue; // don't do self-atom
               double* distances = getDistanceXYZ(system, i,j,i,l);
               r = distances[3];

               if (qualify_bond(system, r, i, j, l)) {
                    local_bonds++;

                    system.molecules[i].atoms[j].bonds.push_back(l);

                    // check for duplicate bond.
                    duplicateFlag=0;
                    for (unsigned int n=0;n<system.constants.uniqueBonds.size();n++) {
                        if (system.constants.uniqueBonds[n].mol == i &&
                            system.constants.uniqueBonds[n].atom1 == l &&
                            system.constants.uniqueBonds[n].atom2 == j) {

                            duplicateFlag=1; break;
                        }
                    }
                    if (!duplicateFlag) {
                        // this bond is unique
                        Constants::UniqueBond tmp; tmp.mol=i; tmp.atom1=j; tmp.atom2=l;
                        tmp.value = r;
                        system.constants.uniqueBonds.push_back(tmp);
                    }

                    // now loop for angles
                    for (m=0; m<system.molecules[i].atoms.size(); m++) {
                        if (j==m) continue; // don't do ABA angles (0 degrees)
                        if (l==m) continue; // don't do ABB angles (0 degrees)
                        double* distancesa = getDistanceXYZ(system, i,l,i,m);
                        ra = distancesa[3];

                        if (qualify_bond(system, ra, i, l, m)) {

                            // check for duplicate angles
                            duplicateAngleFlag = 0;
                            for (unsigned int n=0;n<system.constants.uniqueAngles.size();n++) {
                                if (system.constants.uniqueAngles[n].mol == i &&
                                    system.constants.uniqueAngles[n].atom1 == m &&
                                    system.constants.uniqueAngles[n].atom2 == l &&
                                    system.constants.uniqueAngles[n].atom3 == j) {

                                    duplicateAngleFlag=1; break;

                                }
                            }
                            if (!duplicateAngleFlag) {
                                // this angle is unique
                                Constants::UniqueAngle tmp; tmp.mol=i; tmp.atom1=j; tmp.atom2=l; tmp.atom3=m;
                                tmp.value = get_angle(system, i, j,l,m);
                                system.constants.uniqueAngles.push_back(tmp);
                            }

                            // now loop for dihedrals
                            for (p=0; p<system.molecules[i].atoms.size(); p++) {
                                double * distancesh = getDistanceXYZ(system,i,m,i,p);
                                rh = distancesh[3];
                                if (j==l) continue; // don't do AABC
                                if (j==p) continue; // don't do ABCA
                                if (l==m) continue; // don't do ABBC
                                if (m==p) continue; // don't do ABCC
                                if (j==m) continue; // don't do ABAC
                                if (l==p) continue; // don't do ABCB

                                if (qualify_bond(system, rh, i, m, p)) {
                                    // check duplicate dihedral
                                    duplicateDihFlag = 0;
                                    for (unsigned int n=0;n<system.constants.uniqueDihedrals.size();n++) {
                                        if ((system.constants.uniqueDihedrals[n].mol==i &&
                                            system.constants.uniqueDihedrals[n].atom1==p &&
                                            system.constants.uniqueDihedrals[n].atom2==m &&
                                            system.constants.uniqueDihedrals[n].atom3==l &&
                                            system.constants.uniqueDihedrals[n].atom4==j) ||
                                            (system.constants.uniqueDihedrals[n].mol==i &&
                                            system.constants.uniqueDihedrals[n].atom1==p &&
                                            system.constants.uniqueDihedrals[n].atom2==m &&
                                            system.constants.uniqueDihedrals[n].atom3==l &&
                                            system.constants.uniqueDihedrals[n].atom4==j )) {

                                                duplicateDihFlag = 1;break;
                                        }
                                    }
                                    // this dihedral is unique
                                    if (!duplicateDihFlag) {
                                        Constants::UniqueDihedral tmp; tmp.mol=i; tmp.atom1=j;
                                        tmp.atom2=l; tmp.atom3=m; tmp.atom4=p;
                                        tmp.value = get_dihedral_angle(system, i, j,l,m,p);
                                        system.constants.uniqueDihedrals.push_back(tmp);
                                    }
                                }

                            }
                        }
                    }

               } // end if r < bond-length
            } // end pair (i,j) -- (i,l)

        } // end j
    } // end i


    // get UFF atom labels for all atoms
    for (i=0;i<molecule_limit;i++) {
        for (j=0;j<system.molecules[i].atoms.size();j++) {
            // based on the total number of bonds to this atom,
            // determine the atom-type from UFF.
            system.molecules[i].atoms[j].UFFlabel = getUFFlabel(system, system.molecules[i].atoms[j].name, system.molecules[i].atoms[j].bonds.size(), i,j);
        }
    }


    // get unique qualified LJ/ES non-bond pairs (beyond 1,3)
    unsigned int mol,qualified, y,z;
    double rlj;
    const double r_c = (system.pbc.cutoff==0) ? 12.0 : (system.pbc.cutoff); // default 12A if non-periodic
    for (mol=0; mol<molecule_limit; mol++) {
        // all pairs inside the molecule
        for (i=0; i<system.molecules[mol].atoms.size(); i++) {
            for (j=i+1; j<system.molecules[mol].atoms.size(); j++) {
                // need to check if beyond 2 bonds -- i.e. no 1-2 or 1-3 interactions.
                qualified = 1;
                for (y=0; y<system.constants.uniqueBonds.size();y++) {
                    if (system.constants.uniqueBonds[y].mol==mol &&
                        ((
                        system.constants.uniqueBonds[y].atom1==i &&
                        system.constants.uniqueBonds[y].atom2==j
                        ) ||
                        (
                        system.constants.uniqueBonds[y].atom1==j &&
                        system.constants.uniqueBonds[y].atom2==i
                        ))) {
                        qualified=0;  break;
                    } // end if bonded therefore unqualified
                } // end bonds loop
                for (z=0; z<system.constants.uniqueAngles.size();z++) {
                    if (system.constants.uniqueAngles[z].mol==mol &&
                       ((
                        system.constants.uniqueAngles[z].atom1==i &&
                        system.constants.uniqueAngles[z].atom3==j
                       ) ||
                       (
                        system.constants.uniqueAngles[z].atom3==i &&
                        system.constants.uniqueAngles[z].atom1==j
                       ))) {
                       qualified=0; break;
                    } // end if 1--3 therefore unqualified
                }

                double* distanceslj = getDistanceXYZ(system,mol,i,mol,j);
                rlj = distanceslj[3];

                // apply cutoff..
                if (rlj > r_c) qualified=0;

                if (qualified) {
                    Constants::UniqueLJNonBond tmp;
                    tmp.mol = mol; tmp.atom1=i; tmp.atom2=j;
                    tmp.sig = 0.5*(system.molecules[mol].atoms[i].sig + system.molecules[mol].atoms[j].sig);
                    tmp.eps = sqrt(system.molecules[mol].atoms[i].eps * system.molecules[mol].atoms[j].eps);
                    system.constants.uniqueLJNonBonds.push_back(tmp);

                    //also coulombic pairs
                    Constants::UniqueChargeNonBond tmp2;
                    tmp2.mol = mol; tmp2.atom1 = i; tmp2.atom2 = j;
                    tmp2.chargeprod = system.molecules[mol].atoms[i].C * system.molecules[mol].atoms[j].C;
                    system.constants.uniqueChargeNonBonds.push_back(tmp2);
                }
            } // end pair-atom j
        } // end atom loop i
    } // end molecule loop mol

}


void setBondingParameters(System &system) {
    // save all bond/angle/torsion/non-bond parameters to memory
    // (before running optimization)
    unsigned int it, mol,atom1,atom2,atom3;//,atom4;

    // 1) bonds
    for (it=0; it<system.constants.uniqueBonds.size(); it++) {
        mol = system.constants.uniqueBonds[it].mol;
        atom1 = system.constants.uniqueBonds[it].atom1;
        atom2 = system.constants.uniqueBonds[it].atom2;

        system.constants.uniqueBonds[it].BO = get_BO(system.molecules[mol].atoms[atom1].UFFlabel, system.molecules[mol].atoms[atom2].UFFlabel);
        system.constants.uniqueBonds[it].rij = get_rij(system,mol,atom1,mol,atom2);
        system.constants.uniqueBonds[it].kij = get_kij(system,mol,atom1,mol,atom2, system.constants.uniqueBonds[it].rij);
        system.constants.uniqueBonds[it].Dij = system.constants.uniqueBonds[it].BO*70.0; // kcal/mol
        system.constants.uniqueBonds[it].alpha = sqrt(0.5*system.constants.uniqueBonds[it].kij/system.constants.uniqueBonds[it].Dij);

    }

    // 2) angles
    for (it=0; it<system.constants.uniqueAngles.size(); it++) {
        mol = system.constants.uniqueAngles[it].mol;
        atom1 = system.constants.uniqueAngles[it].atom1;
        atom2 = system.constants.uniqueAngles[it].atom2;
        atom3 = system.constants.uniqueAngles[it].atom3;

        system.constants.uniqueAngles[it].rij = get_rij(system,mol,atom1,mol,atom2);
        double rij = system.constants.uniqueAngles[it].rij;
        system.constants.uniqueAngles[it].rjk = get_rij(system,mol,atom2,mol,atom3);
        double rjk = system.constants.uniqueAngles[it].rjk;
        system.constants.uniqueAngles[it].theta_ijk = M_PI/180.0*system.constants.UFF_angles[system.molecules[mol].atoms[atom2].UFFlabel.c_str()];
        double theta_ijk = system.constants.uniqueAngles[it].theta_ijk;
        system.constants.uniqueAngles[it].C2 = 1.0/(4.0*sin(theta_ijk)*sin(theta_ijk));
        double C2 = system.constants.uniqueAngles[it].C2;
        system.constants.uniqueAngles[it].C1 = -4.0*C2*cos(theta_ijk);
        system.constants.uniqueAngles[it].C0 = C2*(2.0*cos(theta_ijk)*cos(theta_ijk) + 1.0);

        double Zi = system.constants.UFF_Z[system.molecules[mol].atoms[atom1].UFFlabel.c_str()];
        double Zk = system.constants.UFF_Z[system.molecules[mol].atoms[atom3].UFFlabel.c_str()];

        // stuff for computing other terms (which have dihedral dependence)
        system.constants.uniqueAngles[it].t1 = 664.12/(rij*rjk)*Zi*Zk*rij*rjk;
        system.constants.uniqueAngles[it].t2 = cos(theta_ijk);
        system.constants.uniqueAngles[it].t3 = rij*rjk;
    }

    // 3) dihedrals
    for (it=0; it<system.constants.uniqueDihedrals.size(); it++) {
        mol = system.constants.uniqueDihedrals[it].mol;
        atom1 = system.constants.uniqueDihedrals[it].atom1;
        atom2 = system.constants.uniqueDihedrals[it].atom2;
        atom3 = system.constants.uniqueDihedrals[it].atom3;
        //atom4 = system.constants.uniqueDihedrals[it].atom4;

        double * params = get_torsion_params(system, system.molecules[mol].atoms[atom2].UFFlabel, system.molecules[mol].atoms[atom3].UFFlabel);
        system.constants.uniqueDihedrals[it].phi_ijkl = params[0]*M_PI/180.0;
        system.constants.uniqueDihedrals[it].vjk = params[1];
        system.constants.uniqueDihedrals[it].n = params[2];
    }


}



double totalBondedEnergy(System &system) {
    // each function here saves the component energies to Stats class. (system.stats)
    system.stats.Ubonded_tot.value = 0;
    if (system.constants.opt_bonds)
        system.stats.Ubonded_tot.value += stretch_energy(system);
    if (system.constants.opt_angles)
        system.stats.Ubonded_tot.value += angle_bend_energy(system);
    if (system.constants.opt_dihedrals)
        system.stats.Ubonded_tot.value += torsions_energy(system);
    if (system.constants.mode == "opt") {
        if (system.constants.opt_LJ)
            system.stats.Ubonded_tot.value += LJ_intramolec_energy(system);
        if (system.constants.opt_ES)
            system.stats.Ubonded_tot.value += ES_intramolec_energy(system);
    }
    return system.stats.Ubonded_tot.value;
}
