# mcmd
This is a Monte Carlo / Molecular Dynamics Simulation software.<br />

--> MC Currently fully supporting NPT,NVT,NVE,uVT ensembles.  <br />
--> MD supporting NVT.  <br />
--> POLARIZATION IS NOT WORKING so avoid using  <br />
    `potential_form   ljespolar`<br /><br />

PRE-COMPILED EXECUTABLE WORKS WITH THE FOLLOWING COMPILERS:  <br />
    -> gcc compiler 6.2.0 (circe)  <br />
    -> gcc compiler 4.9.3 (stampede)  <br />

To download: <br />
git clone https://github.com/khavernathy/mcmd<br />

To compile:  (in src dir "src")<br />
`g++ main.cpp -lm -o ../t -I. -std=c++11`  <br />

To run (in base dir "mcmd") <br />
`./t mcmd.inp`<br /><br />  
  
<hr />
  
TODO<br /><br />

-> MD: rotation cap is on theta, so rotational velocity can still get big. This throws off rotational KE<br />
-> for some reason my RD energy is 0.001% off from MPMC every time. Maybe self energy is diff?<br />
-> allow 0 sorbate molecules in uvt<br />
-> make pair lists for MC to run faster AND to do stuff like pressure calculation in MD. <br />
-> fix MC SD's maybe? Seem to be wrong<br />
-> NEED TO CHECK IF 4TH POINT OF PLANE MAKES THE SAME PLANE AS PREVIOUS 3 BEFORE MOVING ON<br />
    -> although I think it's always the case with crystalline systems.
-> REDUCE CHARGE CALCULATIONS BY CONVERTING TO AU BEFOREHAND AND CONVERTING BACK IN OUTPUT<br />
-> make correct pressure calculator (Frenkel method, check MPMC to compare?)<br />  
-> include more-than-static polarization energy  <br />
    -> right now I just use V = -0.5 sum{u.E}  <br />
-> Use GPU for MD force calculations? (add option)  <br />
-> Use GPU for polarization routine (later)  <br />
-> Implement Phast2 model?  <br />
