#include <iostream>
#include <string>
#include <strings.h>
#include <algorithm>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <map>
#include <string>
#include <stdlib.h>
#include <ostream>
#include <sstream>
#include <fstream>

/* quick function for getting volume of a sphere (for normalization of distribution) */
double sphere_volume(double r) {
    return (4.0/3.0)*M_PI*(r*r*r);
}

/* THIS FUNCTION SETS UP THE INITIAL VECTOR FOR RADIAL DISTRIBUTION */
void setupRadialDist(System &system) {
    
    int num_bins = ceil(system.stats.radial_max_dist / system.stats.radial_bin_size);
    printf("The number of radial bins to create is %i\n", num_bins);
    printf("The number of different g(r) pair calculations is %i\n", (int)system.stats.radial_centroid.size());

    vector<long unsigned int> dummy; // = vector<long unsigned int>(num_bins); // dummy to push into the vector of g(r)'s

    // make vectors to hold each g(r), and fill each with zero in all bins
    for (int z=0; z<system.stats.radial_centroid.size();z++) {
        for (int i=0; i<num_bins; i++) {
            dummy.push_back(0); // fill each g(r) bin with zero
        }
        system.stats.radial_bins.push_back(dummy);
        dummy.clear();
    }
    // so if max = 10 and size = 0.2, 50 bins are created with index 0->49.
    // if dist between 0.0 and 0.2, index 0++, etc.
    // ... if dist between 9.8 and 10.0, index 49++.

  //  printf("Testing radial_bins[48] = %i\n",system.stats.radial_bins[48]);
    return;   
}


/* THIS FUNCTION WILL BE CALLED EVERY CORRTIME AND WILL ADD TO BINS AS NEEDED */ 
/* every step is excessive and increases step runtime by ~x15        */
void radialDist(System &system) {
    const double bin_size = system.stats.radial_bin_size;
    const double max_dist = system.stats.radial_max_dist;
    
    for (int y=0; y<(int)system.stats.radial_bins.size(); y++) {
    
    string centroid = system.stats.radial_centroid[y];
    string counterpart = system.stats.radial_counterpart[y];

    //system.checkpoint("starting loop");
    // loop through all the atom pairs. Doing intramolecular too b/c MD needs it sometimes.
    for (int i=0; i<system.molecules.size(); i++) {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
            for (int k=0; k<system.molecules.size(); k++) { // k=i+1 would skip intramolec
                for (int l=0; l<system.molecules[k].atoms.size(); l++) {
                    // and only if its a centroid/counterpart pair
                    if (
                        !(i==k && j==l) // don't do self interaction (r=0)
                     && ((system.molecules[i].atoms[j].name == centroid && system.molecules[k].atoms[l].name == counterpart)
                     || (system.molecules[i].atoms[j].name == counterpart && system.molecules[k].atoms[l].name == centroid))) 
                    {
                        double* distances = getDistanceXYZ(system, i, j, k, l);
                        double r = distances[3];     
                        //printf("distance = %f\n",dist);      
                        //system.checkpoint("getting index"); 
                        if (r < max_dist) {
                            // determine index of radial_bins
                            int indexr = floor(r / bin_size);  // so 0.02/0.2 -> index 0; 0.25/0.2 -> index 1..
                            //system.checkpoint("adding to bin");
                            system.stats.radial_bins[y][indexr]++;
                        } // end dist<max_dist

                    } // end if proper pair.
                } // end atoms-in-k loop l
            } // end molecules loop k
        } // end atoms-in-i loop j
    } // end molecules loop i
    //system.checkpoint("ending radialDist");
    
    } // end bins loop for different g(r)'s
    return;
}

/* THIS FUNCTION WILL BE CALLED AFTER EACH CORRTIME AND WRITE THE RADIAL DISTRIBUTION DATA TO FILE */
void writeRadialDist(System &system) {
        
    for (int y=0; y<(int)system.stats.radial_bins.size(); y++) {
        string suffix = to_string(y);
        string radfilename = system.stats.radial_file;
        radfilename = radfilename + suffix;
        remove(radfilename.c_str()); // JIC
    

    ofstream radfile;
    radfile.open (radfilename, ios_base::app);
    radfile << "#r_(" << system.stats.radial_centroid[y].c_str() << "--" << system.stats.radial_counterpart[y].c_str() << ")   #count(normalized%)\n";
    radfile << "0       0\n";

    double spherev = 0.0;
    double prevspherev = 0.0;
    //double sum = 0.0;

    //loop to generate sum
    /*
    for (int i=0; i<system.stats.radial_bins[y].size(); i++) {
        spherev = sphere_volume((i+1)*system.stats.radial_bin_size);
        sum += system.stats.radial_bins[y][i]/(spherev - prevspherev);
        prevspherev = spherev;
    }
    */
    // reset prevspherev
    //prevspherev=0.0;
    
    // loop to write normalized counts
    for (int i=0; i<system.stats.radial_bins[y].size(); i++) {
        spherev = sphere_volume((i+1)*system.stats.radial_bin_size);
        radfile << ((double)(i+1) * system.stats.radial_bin_size);
        radfile << "        ";
            // normalize as density of sorbates in selected r-region (N/V) as a percent
            // i.e. the integral of g(r) from 0 -> maximum r = 100
            // NEVERMIND. sum is unconventional and I'm turning it off.
            radfile << system.stats.radial_bins[y][i]/(spherev - prevspherev);// /sum *100;         
            radfile << "\n";
        prevspherev = spherev;
    }      
    radfile.close();

    } // end loop y for different g(r)'s
    return;
}

void countAtomInRadius(System &system, string atomname, double radius) {
    int count=0;
    double r=0;
    // loop through all atoms and count the number that are within radius
    for (int i=0; i<system.molecules.size(); i++) {
        for (int j=0; j<system.molecules[i].atoms.size(); j++) {
            if (system.molecules[i].atoms[j].name == atomname) {
                for (int n=0; n<3; n++) 
                    r += system.molecules[i].atoms[j].pos[n] * system.molecules[i].atoms[j].pos[n];

                r = sqrt(r);
                if (r <= radius) count++;
            }
        }
    }

    system.stats.dist_within.value = count;
}

