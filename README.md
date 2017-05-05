# MCMD
This is a Monte Carlo / Molecular Dynamics Simulation software used primarily for gas sorption in crystalline materials. 

# Quick start:
On HPC clusters, you may need to load a compiler module first:  <br />
&emsp;-> `module load compilers/gcc/6.2.0` (circe)  <br />
&emsp;-> `module load gcc/6.3.0` (bridges) <br />
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

# Features
MCMD has several features that make simulation easier for the user:<br />
&emsp;-> MC supports NPT, NVT, NVE, uVT ensembles.  <br />
&emsp;-> MD supports NVT.  <br />
&emsp;->Automatic radial distribution calculator<br />
&emsp;->Pre-programmed sorbate models for easy input, including multi-sorbate support<br />
&emsp;->Easy system basis parametrization via a,b,c,alpha,beta,gamma crystal values<br />
&emsp;->Quick routines for energy computation<br />
&emsp;->Faster than most serial MC/MD programs to-date<br />
&emsp;-> Simulated annealing
&emsp;-> Any periodic cell is supported for both MC and MD; non-periodic systems also supported.<br />
&emsp;-> Force-fields available are Lennard-Jones (12-6), Ewald electrostatics, and Thole-Applequist polarization.<br />
&emsp;-> Sample inputs are included. The program takes just one argument: the input file (which itself usually points to a file containing starting atoms).<br />

# What can be obtained from this software
The program outputs several quantities of interest:<br />
&emsp;->Uptake of sorbates in wt%, reduced wt%, cm^3/g and mmol/g (MC)<br />
&emsp;->Excess adsorption ratio in mg/g (MC)<br />
&emsp;->Selectivities for multi-sorbate simulation (MC)<br />
&emsp;->Qst (heat of sorption) for sorbate (MC)<br />
&emsp;->Sorbate occupation about some site/atom (g(r)) (MC)<br />
&emsp;->Compressibility factor Z for homogenous gases as indicator of ideality (MC)<br />
&emsp;->Monte Carlo statistics (accepted moves, Boltzmann Factor and Acceptance Ratios) (MC)<br />
&emsp;->Diffusion coefficient and specific heat (MD)<br />
&emsp;->Emergent pressure approximation (MD)<br />
&emsp;->Trajectory and restart files to easily pickup a halted job and visualize simulation<br />
&emsp;->3D histogram data for visualization of sorbate occupation in a material (density visualization).<br />
&emsp;->Induced dipole strengths for polarization simulations<br />

# Operating System requirements
MCMD works out-of-the-box on Linux (tested on Ubuntu 16.04) and Mac (tested on OS X El Capitan v10.11.6).<br />
We recommend Visual Molecular Dynamics for data visualization, but the output is compatible with most other software<br />

# TODO
-> MC: add "desired N" option to fill a system until a desired N is reached (accelerate uptake), then do NVT<br />
-> MC: speed up by adjusting energy by the 1 particle that is changed.<br /> 
-> MC: add multi-sorb Qst calculator<br />
-> MD: speed up via GPU for MD force calculations<br />
-> MC: speed up via GPU for polarization routine <br />
-> MC: add Phast2 model<br />
-> MD: Make correct(ed) force (most important), pressure, specific heat, temperature calculators<br />
-> MD: add polarizable force calculator <br />
-> MC: add more sorbate models from literature<br />
-> MC/MD: Flexible materials (bonding potential)<br />
-> MC: Use openMP or fork() for energy tasks<br />
-> change timing/string stuff for icpc compatibility<br />
-> try out float instead of double optimization<br />

<br />
TESTING TODO:<br /><br />
-> MC: replay (frame-by-frame) energy calculations compared to MPMC<br />
-> MC: benchmarks: LJ, LJ+ES, LJ+ES+POLAR+FH of MOF-5+H2 vs. MPMC<br />
-> MC: non orthorhombic benchmarks: MPM-1-Br, Mn-Formate<br />
-> MC: ensemble tests: NPT, uVT, NVT, NVE<br />
-> MD: NVT test<br />
-> Mac/Linux tests. HPC platform tests<br />
