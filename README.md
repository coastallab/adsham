# [ADSHAM_v1.0] A robust HLL-based shceme for capturing contact-discontinuity in scalar transport by shallow water flow
## Description
This archive contains the source code for a robust HLL-based shceme for capturing contact-discontinuity in scalar transport by shallow water flow.
The source code is developed by Sooncheol Hwang. 
The entire project was advised by Sangyoung Son, PhD.

## Authors
- Sooncheol Hwang [leneords@hanmail.net](mailto:leneords@hanmail.net)
- Sangyoung Son [sson@korea.ac.kr](mailto:sson@korea.ac.kr)

## Mathematical Model
This code solves the nonlinear shallow water equations (NSWEs) and the scalar transport equation (STE).
Governing equations are discretized using a hybrid finite volume-finite difference shceme (Kurganov and Petrova, 2007).
An anti-diffusion function is implemented in HLL Riemann solver for scalar concentration.
Diffusion terms in the STE are discretized using the finite difference method.

## References
- Kurganov, A., & Petrova, G. (2007). A second-order well-balanced positivity preserving central-upwind scheme for the Saint-Venant system. Communications in Mathematical Sciences, 5(1), 133-160.

## Code Architecture
The source code is written by MATLAB.
It consists 10 MATLAB scripts, and the main script is Engine.m.
The description of each file is as follows:
- Engine.m : Main engine for the numerical simulation. It contains every computation such as governing equation and time integration. Also it involves the following 9 subscripts.
- setParams.m : Takes physical and numerical parameters from input file.
- reconstB.m : Automatically reconstruct bathymetry matrix with the given grid resolution (dx, dy) from original bathymetry file.
- BoundaryCondition.m : Calculate the variables at the boundary cells with the given boundary conditions. Note that fully reflective wall boundary condition (B.C.) and periodic B.C. are available now.
- CalcUVC.m : Calculates variables (h, u, v, c) from the conservative variables (w, hu, hv, hc) at the cell interfaces.
- CAntiDissipation.m : Calculates the anti-diffusion function for minimizing the numerical diffusion of the scalar concentration.
- reconstruction.m : Calculates variables (w, hu, hv, hc) at the cell interfaces from the cell-averaged conservative variables (w, hu, hv, hc) at the cell center. It contains minmod.m which can be replaced by another limiter.
- minmod.m : Calculate a generalized minmod limiter.
- correction.m : Corrects basic piecewise linear reconstruction to guarantee the non-negativity of the computed water depth. ![Fig1](https://user-images.githubusercontent.com/44635414/161505692-548fc66b-aa9c-45be-be17-9ea05e050a5b.png)
- HLLsolver.m : Calculate a HLL Riemann solver to estimate the numerical flux at the cell interface.
