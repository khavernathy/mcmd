# mcmd
This is a Monte Carlo / Molecular Dynamics Simulation software used primarily for gas sorption in crystalline materials.<br />
It works on Linux (tested on Ubuntu 16.04 and ) and Mac (tested on OS X El Capitan v10.11.6).<br />

--> MC supports NPT, NVT, NVE, uVT ensembles.  <br />
--> MD supports NVT.  <br />
--> Any periodic cell is supported for both MC and MD; non-periodic systems also supported.<br />
--> Force-fields available are Lennard-Jones (12-6), (Ewald) electrostatics, and Thole-Applequist polarization.<br />
--> Sample inputs are included. The program takes just one argument: the input file.<br />

On computer clusters, you may need to load a compiler module first:  <br />
&emsp;-> gcc compiler 6.2.0 (circe)  <br />
&emsp;-> gcc compiler 4.9.3 (stampede)  <br />

Quick start:<br />
1) Download: <br />
`git clone https://github.com/khavernathy/mcmd` <br />

2) Compile: <br />
`cd mcmd` <br />
`cd src` <br />
`g++ main.cpp -lm -o ../mcmd -I. -std=c++11`, or simply `bash compile.sh` <br />
`cd ..` <br />

3) Run: <br />
`./mcmd mcmd.inp`<br /><br />  
  
<hr />
  
TODO<br /><br />
-> MC: reduce energy calculations to adjust by the 1 particle that is changed.<br /> 
-> MC: add multi-sorb Qst calculator<br />
-> MC: add LJ Feynmann-Hibbs corrections for hydrogen<br />
-> MD: Use GPU for MD force calculations? (add option)<br />
-> MC: Use GPU for polarization routine (later)<br />
-> MC: Implement Phast2 model?<br />
