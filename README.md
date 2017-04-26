# mcmd
This is a Monte Carlo / Molecular Dynamics Simulation software used primarily for gas sorption in crystalline materials. It has several features that make simulation easier for the user:<br />
&emsp;->Automatic radial distribution calculator<br />
&emsp;->Pre-programmed sorbate models for easy input, including multi-sorbate support<br />
&emsp;->Easy system basis parametrization via a,b,c,alpha,beta,gamma crystal values<br />
&emsp;->Quick routines for energy computation<br />

The software outputs several averaged quantities of interest:<br />
&emsp;->Uptake in wt%, reduced wt%, and mmol/g (MC)<br />
&emsp;->Excess adsorption ratio in mg/g (MC)<br />
&emsp;->Qst (heat of sorption) for sorbate (MC)<br />
&emsp;->Sorbate occupation about some site/atom (g(r)) (MC)<br />
&emsp;->Selectivity for multi-sorbate simulations (MC)<br />
&emsp;->Compressibility factor Z for homogenous gases as indicator of ideality (MC)<br />
&emsp;->Monte Carlo statistics (accepted moves, Boltzmann Factor and Acceptance Ratios) (MC)<br />
&emsp;->Diffusion coefficient and specific heat (MD)<br />
&emsp;->Emergent pressure approximation (MD)<br />

It works on Linux (tested on Ubuntu 16.04) and Mac (tested on OS X El Capitan v10.11.6).<br />

--> MC supports NPT, NVT, NVE, uVT ensembles.  <br />
--> MD supports NVT.  <br />
--> Any periodic cell is supported for both MC and MD; non-periodic systems also supported.<br />
--> Force-fields available are Lennard-Jones (12-6), (Ewald) electrostatics, and Thole-Applequist polarization.<br />
--> Sample inputs are included. The program takes just one argument: the input file (which itself usually points to a file containing starting atoms).<br />

On HPC clusters, you may need to load a compiler module first:  <br />
&emsp;-> `module load compilers/gcc/6.2.0` (circe)  <br />
&emsp;-> `module load gcc/6.3.0` (bridges) <br />

Quick start:<br />
1) Download: <br />
`git clone https://github.com/khavernathy/mcmd` <br />

2) Compile: <br />
`cd mcmd` <br />
`cd src` <br />
`bash compile.sh` <br />
`cd ..` <br />

3) Run: <br />
`./mcmd mcmd.inp`<br /><br />  
  
<hr />
  
TODO<br /><br />
-> MC: add "desired N" option to fill a system until a desired N is reached, then do NVT<br />
-> MC: speed up by adjusting energy by the 1 particle that is changed.<br /> 
-> MC: add multi-sorb Qst calculator<br />
-> MD: speed up via GPU for MD force calculations<br />
-> MC: speed up via GPU for polarization routine <br />
-> MC: add Phast2 model<br />
-> MD: Make correct(ed) force (most important), pressure, specific heat, temperature calculators<br />
&emsp;-> Corrected periodic (e.g. ewald) force would probably solve several MD bugs at once<br />
