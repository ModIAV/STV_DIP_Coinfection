# STV_DIP_Coinfection
This code implements a multiscale model of influenza A virus and defective interfering particle co-infection developed at the MPI Magdeburg. The current model version is documented in [1]. A detailed description of the intra- and extracellular model is provided in [2] and [3], respectively.

The code provided here allows to simulate the model with different MOI and MODIP conditions and to conduct the parameter estimation performed in [1]. 

## References
1. Rüdiger D, Pelz L, Mein MD, Kupke SY, Reichl U. MMultiscale model of defective interfering particle replication for influenza A virus infection in animal cell culture. Unpublished.
2. Laske T, Heldt FS, Hoffmann H, Frensing T, Reichl U. Modeling the intracellular replication of influenza A virus in the presence of defective interfering RNAs. Virus Research. 2016;213:90-99. https://www.doi.org/10.1016/j.virusres.2015.11.016
3. Rüdiger D, Kupke SY, Laske T, Zmora P, Reichl U. Multiscale modeling of influenza A virus replication in cell cultures predicts infection dynamics for highly different infection conditions. PLOS Computational Biology. 2019;15(2):e1006819. https://www.doi.org/10.1371/journal.pcbi.1006819

## Requirements
- MATLAB (MathWorks, Inc.)

- IQM Toolbox for MATLAB by Schmidt and Jirstrand (Bioinformatics, 2006), available at https://iqmtools.intiquan.com/main.html

## Optional programs (for faster simulation and optimization)
- C/C++ compiler: Creates MEX-files for a faster simulation with the SB Toolbox (e.g. MinGW 6.3 C/C++ for Windows or GCC for Linux)

- CVODE solver from SUNDIALS: Simulates MEX-files. Cohen and Hindmarsh (Computers in Physics, 1996), available at https://computation.llnl.gov/projects/sundials/sundials-software

## Running the code and main options
The function `STV_DIP_CoinfectionModel_Main.m` is used for model simulation and parameter estimation. A model simulation can be done by running this script. The following main options are available:
-	Define if a model simulation, parameter estimation or just the plotting part should be conducted. Set the variable `p.Task` to either 'simulate', 'optimize' or 'plot'. 

-	Define the MOI and MODIP conditions that should be used in the variables `p.SimCond` or `p.OptCond`.

-	Choose the optimization algorithm. Set `p.Solver` to 'fminsearch', 'fmincon', or 'CMAES'.

A simulation of the model for conditions 2 to 13 reproduces the main results provided in [1].

## Contributors
The code base was written by Stefan Heldt. Code and model extension was performed by Daniel Rüdiger. 

## Citation
If you use this code, we ask you to cite the appropriate papers in your publication.
