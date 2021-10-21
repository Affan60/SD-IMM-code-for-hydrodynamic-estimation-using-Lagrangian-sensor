The sub-files details are as follows:
1- IMM_Exp.m : it contains matlab code for SD-IMM with Kalman filter to estimate water level and water velocites in irrigation channel 
2- matrix2.m: it contains Matlab code for 1D saint venant model. 
3- steady_values.m : it contains Matlab code for linearization of 1D saint venant model. it is used in matrix2.m file. 
4- robertsidae.m: it has some funtions to solve ODEs required for linearization. It is used in steady_values.m file.
5- boundary_cond.m : it has Matlab code to compute the upstream and downstream boundary conditions. It uses water levels and gate equations to compute water velocities. It considers undershot gate at upstream and overshot gate at downstream end.   
6- interpolation.m : it has Matlab code for Lagrange interpolation, used in boundary_cond.m file. 
7- MBL_rangedata.txt : This text file contains Eulerian sensor values of water levels at upstream end of undershoot gate at upstream end of irrigation channel. This is only used in boundary_cond.m file. 
8- 0KM_rangedata.txt : This text file contains Eulerian sensor values of water levels at downstream side of udnershot gate at upstream end of irrigation channel. This file is also used as upstream water level measurements in estimation. 
9- 1KM_rangedata.txt : This text file contains Eulerian sensor values placed at 1 km distance from upstream end of irrigation channel. 
10- 2KM_rangedata.txt : This text file contains Eulerian sensor values placed at 2 km distance from upstream end of irrigation channel. 
11- 3KM_rangedata.txt : This text file contains Eulerian sensor values placed at downstream end (at 3 km distance from upstream end) of irrigation channel. 
12- position2.txt : This text file contains Lagrangian sensor GPS position data. 