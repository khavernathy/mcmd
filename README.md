# mcmd
This is a Monte Carlo / Molecular Dynamics Simulation software.

--> MC Currently fully supporting NPT,NVT,NVE,uVT ensembles.  
--> MD supporting NVT, basically.  
--> POLARIZATION IS NOT WORKING so avoid using  
    potential_form   ljespolar  

PRE-COMPILED EXECUTABLE WORKS WITH THE FOLLOWING COMPILERS:  
    -> gcc compiler 6.2.0 (circe)  
    -> gcc compiler 4.9.3 (stampede)  

To compile:  (in src dir)
g++ main.cpp -lm -o ../t -I. -std=c++11  

To run  
./my_executable myinput.inp  
  
------------------------------------------
  
TODO

-> make universal unit-cell  
    -> Need to make PDB restart output with the box, first.
-> make correct pressure calculator (Frenkel method, check MPMC to compare)  
-> make correct thermostat for NVT dynamics  
-> write .pdb traj (or just restart.pdb) with BOX included  
-> include more-than-static polarization energy  
	-> right now I just use V = -0.5 sum{u.E}  
-> add rotation to MD displacements  
    -> Added but a bit dysfunctional.   
    -> maybe torque is calc'd wrong. should be local force on atoms or something  
-> Use GPU for MD force calculations? (add option)  
    -> Use GPU for polarization routine (later)  
-> Implement Phast2 model?  
