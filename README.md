# LBM-Thomas  
Source code repository for the Bachelor's Thesis of Thomas Pereira do Carmo.

## Description  
This project implements the Lattice Boltzmann Method (D2Q9 model, BGK operator) for two-dimensional flow simulation in a channel with a square obstacle. The code allows Reynolds number variation via a parabolic inlet profile and variations of the outlet boundary condition.  

## Repository structure  
- `LBMOpen/` — Main code. Contain the files to compile and run the simulation.  
- `ReadFiles/` — Codes to read diferent aspects of the outputs files from LBMOpen.  
- `SimCode/` — main C++ implementation (without segmentation of functions in subfiles): initialization, collision, streaming, boundary conditions (bounce-back, prescribed velocity, convective outlet), statistics, and drag computation.  
- `README.md` — this file.  

## Dependencies and compilation  
- C++ compiler supporting C++14 or higher.  
- Makefile available for compilation (`make all`).  
- Only standard libraries required.  
- The Python post-processing module (≈ 12.8%) requires:  
  - Python 3.x  
  - `numpy`, `vtk` (or `pyvista`), and `matplotlib`.  

<!-- ===========================
## How to run  
1. Adjust simulation parameters in the configuration file (Reynolds number, domain size, obstacle location, etc.).  
2. Compile the code with `make all`.  
3. Run the resulting executable; it will produce VTK output files containing the 2D velocity field (u, v).  
4. Use the Python post-processing scripts to visualize flow fields, plot velocity profiles, or generate animations.  
5. Check statistical output files for drag coefficient and other quantities to validate results.  

## Reynolds ramp-up  
The code supports Reynolds number ramp-up starting from Re = 70 and increasing linearly through 80 → 90 → 100 → 150 → 200 → 250 → 300 over a fixed number of timesteps (e.g. 20 000), followed by 100 000 timesteps at constant velocity. The parameter α = 0.007698 and characteristic length D = 40 are used to compute viscosity.  
=========================== -->

## Boundary conditions  
- Top, bottom walls and Square obstacle: **Half-way bounce-back**.  
- Inlet (left): **prescribed velocity** (parabolic profile).  
<!-- - Outlet (right): **convective boundary condition** (Neumann approximation).  -->

## Academic goal  
The objective of this work is to analyze vortex shedding behind a square obstacle and to investigate why the conventional LBM model may fail to capture secondary frequencies at critical Reynolds numbers.  

## License  
This project is licensed under the MIT License — see the `LICENSE` file for details.  

![Re150_Mrk6](https://github.com/user-attachments/assets/1902d5e0-5d7a-43f3-ac63-485ebfb58922)
