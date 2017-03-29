#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>

// by giving molecule/atom IDs
double * getDistanceXYZ(System &system, int i, int j, int k, int l) {
    if (system.constants.mc_pbc == "on") {
   // calculate distance between atoms
        double rimg;
        double d[3],di[3],img[3],dimg[3];
        int p,q;
        double r,r2,ri,ri2;

        // get dx dy dz
        for (int n=0; n<3; n++) d[n] = system.molecules[i].atoms[j].pos[n] - system.molecules[k].atoms[l].pos[n];

        // images from reciprocal basis.
        for (p=0; p<3; p++) {
            img[p] = 0;
            for (q=0; q<3; q++) {
                img[p] += system.pbc.reciprocal_basis[q][p]*d[q];
            }
            img[p] = rint(img[p]);
        }

        // get d_image
        for (p=0; p<3; p++) {
            di[p]=0;
            for (q=0; q<3; q++) {
                di[p] += system.pbc.basis[q][p]*img[q];
            }
        }
    
        // correct displacement
        for (p=0; p<3; p++)
            di[p] = d[p] - di[p];

        // pythagorean terms
        r2=0; ri2=0;
        for (p=0; p<3; p++) {
            r2 += d[p]*d[p];
            ri2 += di[p]*di[p];
        }
        r = sqrt(r2);
        ri = sqrt(ri2);

        if (isnan(ri) != 0 ) {
            rimg = r;
            for (p=0; p<3; p++)
                dimg[p] = d[p];
        } else {
            rimg = ri;
            for (p=0; p<3; p++)
                dimg[p] = di[p];
        }

        static double output[4];
        for (p=0;p<3;p++) output[p] = dimg[p];
        output[3] = rimg;
        return output;
    }
    else
    {
        //printf("no pbc r\n");
        // no PBC r
        double r, d[3];
        for (int n=0; n<3; n++) d[n] = system.molecules[i].atoms[j].pos[n] - system.molecules[k].atoms[l].pos[n];
        static double output[4];
        for (int p=0; p<3; p++) output[p] = d[p];
        output[3] = sqrt(dddotprod(d, d));
        return output;
    }
}

// by giving two r vectors.
double * getR(System &system, double * com1, double * com2) {
        double rimg;
        double d[3],di[3],img[3],dimg[3];
        int p,q;
        double r,r2,ri,ri2;

        // get dx dy dz
        for (int n=0; n<3; n++) d[n] = com1[n] - com2[n];

        // images from reciprocal basis.
        for (p=0; p<3; p++) {
            img[p] = 0;
            for (q=0; q<3; q++) {
                img[p] += system.pbc.reciprocal_basis[q][p]*d[q];
            }
            img[p] = rint(img[p]);
        }

        // get d_image
        for (p=0; p<3; p++) {
            di[p]=0;
            for (q=0; q<3; q++) {
                di[p] += system.pbc.basis[q][p]*img[q];
            }
        }

        // correct displacement
        for (p=0; p<3; p++)
            di[p] = d[p] - di[p];

        // pythagorean terms
        r2=0; ri2=0;
        for (p=0; p<3; p++) {
            r2 += d[p]*d[p];
            ri2 += di[p]*di[p];
        }
        r = sqrt(r2);
        ri = sqrt(ri2);

        if (isnan(ri) != 0 ) {
            rimg = r;
            for (p=0; p<3; p++)
                dimg[p] = d[p];
        } else {
            rimg = ri;
            for (p=0; p<3; p++)
                dimg[p] = di[p];
        }

        
        static double output[4];
        for (p=0; p<3; p++) output[p] = dimg[p];
        output[3] = rimg;
        return output;
    
}


void computeDistances(System &system) {
    int i,j,k,l,p,q,r;
    int molsize = (int)system.molecules.size();
    int frozenatoms = system.stats.count_frozens;
    
    double whatever[molsize][frozenatoms][molsize][frozenatoms];

    // gets all the pairwise distances for the entire system in one shot.
    for (i=0; i<molsize; i++) {
        for (j=0; j<system.molecules[i].atoms.size(); j++) {
            for (k=0; k<molsize; k++) {
                for (l=0; l<system.molecules[k].atoms.size(); l++) {
                    if (i==k && j==l) continue; // always skip self-atom distance
                    double *distances = getDistanceXYZ(system, i,j,k,l);
                    whatever[i][j][k][l] = distances[3];
                    printf("whatever[%i][%i][%i][%i] = %f\n", i,j,k,l, whatever[i][j][k][l]);
                }
            }
        }
    }

}

