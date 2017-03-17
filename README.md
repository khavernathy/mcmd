# mcmd
This is a Monte Carlo / Molecular Dynamics Simulation software.<br />

--> MC Currently fully supporting NPT,NVT,NVE,uVT ensembles.  <br />
--> MD supporting NVT, basically.  <br />
--> POLARIZATION IS NOT WORKING so avoid using  <br />
    `potential_form   ljespolar`<br /><br />

PRE-COMPILED EXECUTABLE WORKS WITH THE FOLLOWING COMPILERS:  <br />
    -> gcc compiler 6.2.0 (circe)  <br />
    -> gcc compiler 4.9.3 (stampede)  <br />

To compile:  (in src dir)<br />
g++ main.cpp -lm -o ../t -I. -std=c++11  <br />

To run  <br />
./t myinput.inp<br /><br />  
  
<hr />
  
TODO<br /><br />

-> for some reason my RD energy is 0.001% off from MPMC every time. Maybe self energy is diff?<br />
-> allow 0 sorbate molecules in uvt<br />
-> make pair lists for MC to run faster. <br />
-> add input for easy model selection e.g. co2*, co2, h2_bssp, etc<br />
-> make averages update every step instead of every corrtime..<br />
    -> fix SD's maybe?<br />
-> make universal unit-cell<br />  
    -> NEED TO CHECK IF 4TH POINT OF PLANE MAKES THE SAME PLANE AS PREVIOUS 3 BEFORE MOVING ON<br />
-> REDUCE CHARGE CALCULATIONS BY CONVERTING TO AU BEFOREHAND AND CONVERTING BACK IN OUTPUT<br />
-> make correct pressure calculator (Frenkel method, check MPMC to compare)<br />  
-> make correct thermostat for NVT dynamics  <br />
-> include more-than-static polarization energy  <br />
	-> right now I just use V = -0.5 sum{u.E}  <br />
-> add rotation to MD displacements  <br />
    -> Added but a bit dysfunctional.   <br />
    -> maybe torque is calc'd wrong. should be local force on atoms or something  <br />
-> Use GPU for MD force calculations? (add option)  <br />
    -> Use GPU for polarization routine (later)  <br />
-> Implement Phast2 model?  <br />
