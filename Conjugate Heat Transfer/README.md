# CFDLabGroupJ_WS3

To run this code:
1. Download the entire master branch.
2. Unzip openfoam_files.zip files and cut-paste all contents (3 folders and 3 bash files) in the main folder.
3. Alternatively, you can get these 3 folders and 3 bash files from the original resources folder uploaded on moodle
In that case, In your original bash file run_solid_plate, please change line 8 to:
ln -s -f precice-configs/precice_config_plate_parallel.xml precice-config.xml
to:
ln -s -f precice-configs/precice_config_plate_explicit.xml precice-config.xml
 because this is sequential explicit coupling method.


4. We have already added precice_config_convection_implicit.xml file in Resources/precice-configs folder. Just a note, nothing for you to do here.
5. While solving problems, first run solid solver (OpenFOAM) in one terminal and then run Fluid Solver in another terminal.
6. In example 3, run F1, F2 and Solid Plate Exchanger problems in three different terminals simultaneously.

 
 
 Code Explanation:
 The code takes in input the choice to select which example you want to run. 
 You can run your own example as well (similar natured) but you need to make necessary .pgm, .dat and .xml files. You can also choose between implicit and explicit coupling options.
 For arbitrary geometries these conditions need to be satisfied:
 1. No-slip condition on inner obstacles.
2. No obstacle cell can share 2 opposite, 3 or 4 edges with a fluid cell, i.e., only single-edge sharing and corner obstacle cells are allowed.
3. Temperature conditions need to be specified. (Code is only for temperature dependent problems)
*if wall is adiabatic, put T=1000 for that wall. If wall is kept at a constant temperature, insert that value in Te, Tn, Tw or Ts for directional sense instead of T_h and T_c.
*Inflow velocity are recorded as U_in and V_in.
*For inital conditions, set UI, VI, PI, TI.


The binary coding of flags for each cell is mentioned in cfd_ws3_flags.pdf. This helps us distinguish between boundary, obstacle, coupled and fluid cells and specify the geometry and boundary conditions. Matrix to store: flag[i][j]. Option to print the matrix is commented out in the main.c file. Similarly option to print iteration step, dt and loop step as well.
 
Screenshots of worksheet situations saved in tar file name :  screenshots_results.tar.xz
 
Source Code for any precice functions help: https://github.com/precice/precice/blob/develop/src/precice/SolverInterface.hpp

