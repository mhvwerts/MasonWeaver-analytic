# MasonWeaver-analytic: Numerical evaluation of analytic solutions of the Mason-Weaver equation

by Lancelot Barthe and Martinus Werts, 2022.

*ENS Rennes, CNRS.*

[![DOI](https://zenodo.org/badge/569429237.svg)](https://zenodo.org/badge/latestdoi/569429237)

The Mason-Weaver equation (MWE) is a partial differential equation that describes the sedimentation of small particles in a fluid with the particles being subject to Brownian motion.[1] Understanding the evolution of the vertical concentration profile of an initially homogeneous solution of nanoparticles undergoing sedimentation requires solving the MWE. Previously, we used a numerical finite-difference scheme to find concentration profiles obeying the Mason-Weaver equation.[2][3]

We have now developed Python code that evaluates directly the analytic solutions given  in the original publication by Mason and Weaver.[1] The task of numerically evaluating the analytic expressions was not as straight-forward as initially expected, but we have finally arrived at an efficient and robust program for obtaining sedimentation concentration profiles from the analytic solutions of the Mason-Weaver equation. The details of the computations are described in a technical note.[4]


[1] Mason, M.; Weaver, W. "The Settling of Small Particles in a Fluid".
    Phys. Rev. 1924, 23, 412.
    [doi:10.1103/PhysRev.23.412](https://doi.org/10.1103/PhysRev.23.412)

[2] Midelet, J.; El-Sagheer, A. H.; Brown, T.; Kanaras, A. G.;
    Werts, M. H. V. "The Sedimentation of Colloidal Nanoparticles in 
    Solution and Its Study Using Quantitative Digital Photography".
    Part. & Part. Syst. Charact. 2017, 34, 1700095. 
    [doi:10.1002/ppsc.201700095](https://doi.org/10.1002/ppsc.201700095)

[3] [`MasonWeaver-finite-diff`](https://github.com/mhvwerts/MasonWeaver-finite-diff)

[4] Barthe, L.; Werts, M. H. V. "Sedimentation of colloidal nanoparticles
	in fluids: efficient and robust numerical evaluation of analytic solutions
	of the Mason-Weaver equation". ChemRXiv 2022.
	[doi:10.26434/chemrxiv-2022-91vrq](https://doi.org/10.26434/chemrxiv-2022-91vrq)

## Installation and usage

No specific installation is necessary. The package runs with Python 3, and depends only on `numpy` and `scipy`. To use the package, just copy the `masonweaver_analytic.py` file to the working directory for the project at hand and import the required functions from this module into your Python calculation code.

The main functionality in the `masonweaver_analytic` module is provided through two functions:

- `MW_c_profile` calculates the concentration profile for a sedimenting solution at a given time, based on dimensional parameters for the system ('real-world' calculation)
- `MW_adim` is the underlying dimensionless calculation

The script `example.py` gives an example of how to use `MW_c_profile`. This script is an adaptation of the example script demonstrating our finite-difference Mason-Weaver solver, [`MasonWeaver-finite-diff`](https://github.com/mhvwerts/MasonWeaver-finite-diff). The new script yields very similar curves and sedimentation profiles, but these are now calculated with the exact analytic expressions instead of the approximate numerical finite-difference scheme. Mass conservation is perfect for the analytic solution, and the sedimentation profile can be calculated directly for any instant in time, without iterating over many time-steps.

The folder `Figures` contains scripts that do calculations using the `masonweaver_analytic` module and then plot the results in figures for the scientific document which is in preparation. These scripts also provide a reasonably thorough test of the code base.


## To do

- Extend this README, in particular by adding more detailed descriptions for the `MW_c_profile` and `MW_adim` functions.
- Further improve the example script.
- Finish preparing the accompanying scientific document.

