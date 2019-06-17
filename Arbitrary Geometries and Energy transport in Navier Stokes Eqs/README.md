Arbitrary Geometries and Energy Transport for the Navier-Stokes Equations

Code asks user to choose 1 out of 6 geometries given in Worksheet 2. Can also work for arbitrary geometries as long as following conditions are satisfied:
1. No-slip condition on inner obstacles
2. No obstacle cell can share 2 opposite, 3 or 4 edges with a fluid cell, i.e., only single-edge sharing and corner obstacle cells are allowed.


For other arbitrary geometries, 
please provide a .pgm and a .dat file in the same folder. 
Mention name of geometry in .dat file (first line) and while asked in input as well. 
keeping in mind:
*For temperature independent problems, please put beta=0, Pr=1, TI=Te=Tw=Tn=Ts=0 in .dat file. 
*For Temperature dependent problems, if wall is adiabatic, put T=1000 for that wall. If wall is kept at a constant temperature, insert that value in Te, Tn, Tw or Ts for directional sense instead of T_h and T_c.
*Inflow velocity are recorded as U_in and V_in. 
*For inital conditions, set UI, VI, PI, TI.

The binary coding of flags for each cell is mentioned in cfd_ws2_flags.pdf. This helps us distinguish between boundary, obstacle and fluid cells and specify the geometry and boundary conditions. Matrix to store: flag[i][j]. Option to print the matrix is commented out in the main.c file. Similarly option to print iteration step, dt and loop step as well.



Major Function definitions follow:
- calculate_dt() Determine the maximal time step size for stable convergence in an explicit scheme.
- To Set the boundary values for the next time step.
  for temeprature independent problem: boundaryvalues(),
  for temeprature dependent problem: boundaryvalues_t(),

- To Set the obstacle boundary values for the next time step.
  for temeprature independent problem: spec_boundary_val(),
  for temeprature dependent problem: spec_boundary_val_t(),

- calculate_temp- updates temperature value after one time step using forward euler on Dimensionless Energy equation obtained from Boussinesq Approximation. The updates happen not immediately but only after the entire matrix is traversed in the loop. (Gauss-Siedel style)


  - calculate_fg() Determine the values of F and G (diffusion and convection). Remember to put beta=0 and Pr=1 for temeprature dependent problems.
 
-calcualte_rs: This is the right hand side of the pressure equation and used later on for the time step transition.

- sor(): - Iterate the pressure poisson equation until the residual becomes smaller than eps or the maximal number of iterations is performed. Within the iteration loop the operation sor() is used. Whenever itermax is reached, warning message is displayed to show no convergence on that time step
 
- calculate_uv() Calculate the velocity at the next time step.

- Screenshots of worksheet situations saved in tar file name :  screenshots_results.tar.xz 

 