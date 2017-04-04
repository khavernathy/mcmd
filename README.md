# mcmd
This is a Monte Carlo / Molecular Dynamics Simulation software.<br />
It works on Linux (tested on Ubuntu) and Mac (tested on OS X El Capitan v10.11.6).<br />

--> MC supports NPT, NVT, NVE, uVT ensembles.  <br />
--> MD supports NVT.  <br />
--> Any periodic cell is supported for both MC and MD; non-periodic systems also supported.<br />
--> Force-fields available are Lennard-Jones (12-6), (Ewald) electrostatics, and Thole-Applequist polarization.<br />
--> Sample inputs are included. The program takes just one argument: the input file.<br />

PRE-COMPILED EXECUTABLE WORKS WITH THE FOLLOWING COMPILERS:  <br />
&emsp;-> gcc compiler 6.2.0 (circe)  <br />
&emsp;-> gcc compiler 4.9.3 (stampede)  <br />

Quick start:<br />
1) Download: <br />
`git clone https://github.com/khavernathy/mcmd` <br />
`cd mcmd` <br />

2) Compile: <br />
`cd src` <br />
`g++ main.cpp -lm -o ../t -I. -std=c++11`, or simply  <br />
`bash compile.sh` <br />
`cd ..` <br />

3) Run: <br />
`./t mcmd.inp`<br /><br />  
  
<hr />
  
TODO<br /><br />

-> MC: add dipole output and histogram output functions<br />
-> MC: add multi-sorb Qst calculator<br />
-> MC: add LJ Feynmann-Hibbs corrections for hydrogen<br />
-> MC: got polarization working but:<br />
&emsp;-> 4x slower than mpmc. (def. need to make pair lists)<br />
&emsp;-> Aperiodic energy is super low. (high magnitude, negative). Thole field for no PBC may be issue.<br />
-> MC: make pair lists for MC to run faster AND to do stuff like pressure calculation in MD. <br />
&emsp;-> although I think it's always the case with crystalline systems.<br />
-> MD: Emergent pressure seems large.. (Frenkel method p.84)<br />
-> MD: Use GPU for MD force calculations? (add option)<br />
-> MC: Use GPU for polarization routine (later)<br />
-> MC: Implement Phast2 model?<br />
