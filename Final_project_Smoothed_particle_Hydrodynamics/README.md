#to see results
1. Check the report


#to run heated plate problem:
1.	Make clean 
2.	Make 
3.	Press “1”
4.	Simulation will run 

To see results in paraview: 
1.	Load VTKs
2.  Press apply.
3.	Selects Glyph (Sphere), with scaling factor = 0.5 and Glyph mode =”All points”
4.  Press apply.
5.	Select “Temperature” from the drop down list 
6.	Select Background “Black” to see the Temperature profile
7.	Animations can be seen using “Animation tool bar ” using “real time” mode using 10 as a Duration

#to run Cavity problem with 9 particles visulization

1.	Make clean 
2.	Make 
3.	Press “2”
4.	Simulation will run 

To see results in paraview: 
1.	load boundaries.vtk available in repository
2.	Press apply
3.	Load all VTKs
4.  Press apply
5.	Selects Glyph (arrow), with scaling factor = 0.1*(10^-4) and Glyph mode =”All points”. (to see spehere: select "Sphere" glyph)
6.	Press apply 
7.	Select Background “Black” to see nice results
8.  Select "glyph vector" to see velocity profile and "Pressure" to see pressure values.
9.	Animations can be seen using “Animation tool bar ” using “real time” mode using 10 as a Duration

#to run Cavity problem with all particles visulization

1.	Make clean 
2.	Make 
3.	Press “3”
4.	Simulation will run 

To see results in paraview: 
1.	load boundaries.vtk available in repository
2.	Press apply
3.	Load all VTKs
4.  Press apply
5.	Selects Glyph (arrow), with scaling factor = 9*(10^-5) and Glyph mode =”All points”. (to see spehere: select "Sphere" glyph)
6.	Press apply 
7.	Select Background “Black” to see nice results
8.  Select "glyph vector" to see velocity profile and "Pressure" to see pressure values.
9.	Animations can be seen using “Animation tool bar ” using “real time” mode using 10 as a Duration


#relevant links
1. https://github.com/DualSPHysics/DualSPHysics/wiki/3.-SPH-formulation#38-boundary-conditions
2. https://github.com/bigthetaio/mueller-sph
3. https://github.com/phanix5/SPH-solver/blob/master/SPH/SPH/kernel.cu
4. Lid driven :  http://eprints.utem.edu.my/4245/1/SPH_Paper.pdf
5. for boundary: https://github.com/momohuang/SPH-Fluid-Simulation/blob/master/SPHfluid.cpp, https://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf
6. https://cg.informatik.uni-freiburg.de/course_notes/sim_10_sph.pdf


