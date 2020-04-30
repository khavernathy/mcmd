# MCMD
This is a Monte Carlo and Molecular Dynamics Simulation software used primarily for gas sorption in crystalline materials. It is a project that began as a re-write and expansion of <a href="https://github.com/mpmccode/mpmc">Massively Parallel Monte Carlo (MPMC)</a>, another code developed and maintained by our laboratory, led by <a href="http://chemistry.usf.edu/faculty/space/">Brian Space</a> at the University of South Florida, <a href="http://chemistry.usf.edu/">Dept. of Chemistry</a>, <a href="http://chemistry.usf.edu/smmartt/">Smart Metal-organic Materials Advanced Research and Technology Transfer (SMMARTT)</a>.

![MCMD simulation screenshot](https://github.com/khavernathy/mcmd/blob/master/gui/images/mof%2Bco2.png)

# Quick start:
<!--
On HPC clusters, you may need to load a compiler module first:  <br />
&emsp;-> `module load compilers/gcc/6.2.0` (circe)  <br />
&emsp;-> `module load gcc/6.3.0` (bridges) <br />
-->
Using a terminal,<br />
&nbsp;&nbsp;&nbsp;&nbsp;0. (On Windows only) Get the <a href="https://docs.microsoft.com/en-us/windows/wsl/install-win10">Linux Subsystem</a> (<a href="https://www.howtogeek.com/249966/how-to-install-and-use-the-linux-bash-shell-on-windows-10/">easier instructions for beginners</a>) (or equivalent software, e.g. CygWin): 
1. Download: <br />
`git clone https://github.com/khavernathy/mcmd` or <a href="https://github.com/khavernathy/mcmd/archive/master.zip">download .zip file</a><br />

2. Compile: <br />
`cd mcmd` <br />
`cd src` <br />
`bash compile.sh   [ options ]` <br />
`cd ..` <br />
`export PATH=$PATH:/path/to/mcmd/`<br />

3. Run: <br />
`mcmd mcmd.inp`<br /><br />  
  
# Update
`cd mcmd` <br />
`git pull` <br />
`cd src` <br />
`bash compile.sh   [ options ]` <br />

# Advanced compilation
Take a look at mcmd/src/compile.sh for different options in compilation (OS-specific, CUDA implementation, OpenMP, optimization on different HPC systems, etc.)

<hr />

# Docs
You can find details on available options, built-in potentials, etc. on the wiki page: https://github.com/khavernathy/mcmd/wiki

# Contact
Douglas Franz: dfranz@mail.usf.edu
University of South Florida
Dept. of Chemistry

[![MCMD](https://img.youtube.com/vi/rSj4Q_VtO-Y/0.jpg)](https://www.youtube.com/watch?v=rSj4Q_VtO-Y)

# Features
&emsp;-> Monte Carlo simulation in NPT, NVT, NVE, and &mu;VT ensembles.  <br />
&emsp;-> Molecular Dynamics simulation in NVT, NVE, and &mu;VT ensembles.  <br />
&emsp;-> A crystal builder to create fully parameterized supercells from unit cells. <br />
&emsp;-> A fragment creator based around uniquely named atoms. <br /> 
&emsp;-> A LAMMPS input file exporter. <br />
&emsp;-> Trajectories and restart files in various formats. <br />
&emsp;-> Automatic radial distribution calculator<br />
&emsp;-> Hard-coded molecular models for easy input, including multi-molecule support<br />
&emsp;-> Easy system basis parametrization via a, b, c, &alpha;, &beta;, &gamma; crystal parameters, or basis vectors<br />
&emsp;-> Quick routines for energy/force computation<br />
&emsp;-> Simulated annealing<br />
&emsp;-> Any periodic cell is supported for both MC and MD; non-periodic systems also supported.<br />
&emsp;-> Potentials available are Lennard-Jones (12-6), Tang-Toennies (6-8-10), Ewald electrostatics, and Thole-Applequist polarization.<br />
&emsp;-> Built-in force fields from UFF, OPLS, and other sources<br />
&emsp;-> Sample inputs are included. The program takes just one argument: the input file (which itself usually points to a file containing starting atoms).<br />

# What can be obtained from this software
The program outputs several quantities of interest:<br />
&emsp;->Uptake of sorbates in wt%, reduced wt%, cm^3/g, mmol/g and mg/g<br />
&emsp;->Excess adsorption ratio in mg/g<br />
&emsp;->Selectivities for multi-sorbate simulation<br />
&emsp;->Qst (heat of sorption) for sorbate<br />
&emsp;->Sorbate occupation about some site/atom (g(r))<br />
&emsp;->Diffusion coefficient and specific heat<br />
&emsp;->Trajectory and restart files to easily pickup a halted job and visualize simulation<br />
&emsp;->3D histogram data for visualization of sorbate occupation in a material (density visualization).<br />
&emsp;->Induced dipole strengths for polarization simulations<br />

# Operating System requirements
MCMD works out-of-the-box on <br />
&emsp;-> Linux (tested on Ubuntu 16.04)<br />
&emsp;-> Mac (tested on OS X El Capitan v10.11.6)<br />
&emsp;-> Windows (tested using Cygwin and Windows 7 with gcc 5.4.0 installed)<br />
&emsp;-> Raspberry Pi (3, using Raspian OS).<br /><br />

# Visualization
We recommend Visual Molecular Dynamics (VMD) for data visualization, but the output is compatible with most other software, e.g. Avogadro, Molden, Ovito, etc.<br />

# Cite MCMD

<a href="https://onlinelibrary.wiley.com/doi/full/10.1002/adts.201900113">Franz, D. M. et al. MPMC and MCMD: Free High‐Performance Simulation Software for Atomistic Systems. Adv. Theory Sim., 2019. DOI: 10.1002/adts.201900113</a>

# Selected Publications
Below is a list of scientific publications/presentations that have been facilitated by this software.

1. Mukherjee, S. et al. Trace CO2 Capture by an Ultramicroporous Physisorbent with Low Water Affinity. Science Advances, 2019. DOI:10.1126/sciadv.aax9171

2. Franz, D. M., Forrest, K. A., Pham, T., & Space, B. (2016). Accurate H2 Sorption Modeling in the rht-MOF NOTT-112 Using Explicit Polarization. Crystal Growth & Design, 16(10), 6024-6032. DOI:10.1021/acs.cgd.6b01058

3. Mukherjee et al. Halogen‐C2H2 Binding in Ultramicroporous MOFs for Benchmark C2H2/CO2 Separation Selectivity. Chem. Eur. J. 2020. DOI:10.1002/chem.202000008 

4. Pham, T., Forrest, K. A., Franz, D. M., Guo, Z., Chen, B., & Space, B. (2017). Predictive models of gas sorption in a metal–organic framework with open-metal sites and small pore sizes. Physical Chemistry Chemical Physics, 19(28), 18587-18602. DOI:10.1039/C7CP02767B

5. Pham, T., Forrest, K. A., Franz, D. M., & Space, B. (2017). Experimental and theoretical investigations of the gas adsorption sites in rht-metal-organic frameworks. CrystEngComm, 19, 4646-4665. DOI:10.1039/C7CE01032J

6. Franz, D. M.; Dyott, Z.; Forrest, K. A.; Hogan, A.; Pham, T.; Space, B. Simulations of hydrogen, carbon dioxide, and small hydrocarbon sorption in a nitrogen-rich rht-metal–organic framework. Phys. Chem. Chem. Phys. 2017. DOI:10.1039/C7CP06885A

7. Franz, D. M.; Djulbegovic, M.; Pham, T.; Space, B. Theoretical study of the effect of halogen substitution in molecular porous materials for CO2 and C2H2 sorption. AIMS Materials Science. 2017. DOI:10.3934/matersci.2018.2.226

8. Forrest, K. A.; Franz, D. M.; Pham, T. ; Space, B. Investigating C2H2 sorption in a-[M3(O2CH)6] (M = Mg, Mn) Through Theoretical Studies. Cryst. Growth Des. 2018. DOI:10.1021/acs.cgd.8b00770

9. Yu et al. Enhanced Gas Uptake in a Microporous Metal–Organic Framework via a Sorbate Induced-Fit Mechanism. J. Amer. Chem. Soc. 2019. DOI:10.1021/jacs.9b07807

10. Pal et al. A Microporous Co-MOF for Highly Selective CO2 Sorption in High Loadings Involving Aryl C–H...O=C=O Interactions: Combined Simulation and Breakthrough Studies" ACS Inorganic Chem., 2019. DOI:10.1021/acs.inorgchem.9b01402



<!--
pending...
MPM-1-TIFSIX paper
Mn-Formate paper
Zn-datzbdc paper
-->

<!--
# TODO
-->
