#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>

// by giving molecule/atom IDs
double * getDistanceXYZ(System &system, int i, int j, int k, int l) {
    if (system.constants.all_pbc) {
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
                //if (i==0 && j==0 && k==1 && l==0) printf("img[%i] = reciprocal_basis[%i][%i]*d[%i] = %f\n",p,q,p,q,system.pbc.reciprocal_basis[q][p]*d[q]);
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

        if (ri != ri) {
        //if (isnan(ri) != 0 ) {
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
    else // NO PBC calculation.
    {
        // no PBC r
        double d[3];
        for (int n=0; n<3; n++) d[n] = system.molecules[i].atoms[j].pos[n] - system.molecules[k].atoms[l].pos[n];
        static double output[4];
        for (int p=0; p<3; p++) output[p] = d[p];
        output[3] = sqrt(dddotprod(d, d));
        return output;
    }
}

// by giving two r vectors.
double * getR(System &system, double * com1, double * com2, int pbcflag) {
    if (pbcflag) {
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

        if (ri != ri) {
        //if (isnan(ri) != 0 ) {
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
        
    } else {
        // no pbc
        double d[3];
        for (int n=0; n<3; n++) d[n] = com1[n] - com2[n];
        static double output[4];
        for (int p=0; p<3; p++) output[p] = d[p];
        output[3] = sqrt(dddotprod(d, d));
        return output;
    }


}



