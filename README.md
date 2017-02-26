# mcmd
This is a Monte Carlo / Molecular Dynamics Simulation software.

--> MC Currently fully supporting NPT,NVT,NVE,uVT ensembles.
--> MD supporting NVT, basically.
--> POLARIZATION IS NOT WORKING so avoid using
    potential_form   ljespolar

PRE-COMPILED EXECUTABLE WORKS WITH THE FOLLOWING COMPILERS:
    -> gcc compiler 6.2.0 (circe)
    -> gcc compiler 4.9.3 (stampede)

To compile:
g++ main.cpp -lm -o my_executable -I. -std=c++11

To run
./my_executable

------------------------------------------

TODO

-> make averages moving
-> write .pdb traj (or just restart.pdb) with BOX included
-> include more-than-static polarization energy
	-> right now I just use V = -0.5 sum{u.E}
-> add rotation to MD displacements
    --> Added but a bit dysfunctional. 
        Mols tend to rotate into oblivion (lose rigidity) and 
        explode so I have a cap of 0.1 rad/fs now.
        They still disform. Maybe shouldn't do com() after xyz rot.
--> Use GPU for MD force calculations? (add option)
--> Phast2 needs definite fix. Units probably wrong.
