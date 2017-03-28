# mcmd
This is a Monte Carlo / Molecular Dynamics Simulation software.<br />

--> MC Currently fully supporting NPT,NVT,NVE,uVT ensembles.  <br />
--> MD supporting NVT.  <br />

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

-> MC: got polarization working but:<br />
    -> 4x slower than mpmc. (def. need to make pair lists)<br />
-> MC: add S.A. linear? (exponential already there)<br />
-> MC: fix atom-stuck-at-origin issue in uVT<br />
-> MC: for some reason my RD energy is 0.001% off from MPMC every time. Maybe self energy is diff?<br />
-> MC: allow 0 sorbate molecules in uvt<br />
-> MC: make pair lists for MC to run faster AND to do stuff like pressure calculation in MD. <br />
-> MC: fix MC SD's maybe? Seem to be wrong<br />
-> both: NEED TO CHECK IF 4TH POINT OF PLANE MAKES THE SAME PLANE AS PREVIOUS 3 BEFORE MOVING ON<br />
    -> although I think it's always the case with crystalline systems.
-> both: make correct pressure calculator (Frenkel method, check MPMC to compare?)<br />  
-> MD: Use GPU for MD force calculations? (add option)  <br />
-> MC: Use GPU for polarization routine (later)  <br />
-> MC: Implement Phast2 model?  <br />
