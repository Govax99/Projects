%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%              SPACECRAFT ATTITUDE DYNAMICS PROJECT                       %
%                    Academic year 2021/2022                              %
%                                                                         %
%              Description and usage of project functions                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DESCRIPTION OF THE PROJECT STRUCTURE:
The main contents of this project are distributed as:

1)final_sim: folder containing the simulation of detumbling/slew/Earth_pointing/Fixed_pointing/desaturation

2)LIBRARIES: folder containing all the basic blocks used in the simulation



HOW TO USE THE FILES:

1)run the file called RUN_FIRST.m to add all the files and the functions to the path

2)run final_sim\variables.m to initialize all the data for the simulation

3)run final_sim\final_simulation.m

CHOOSING SIMULATION MODE:
To test the various modes do the following (the switches are found in the Controller subsystem) and then RERUN variables.m:

1)Slew or inertial pointing
In rand_q_generator switch the manual switch from q0 to the other option (random quaternion) or directly insert the wanted q in q_ref.
Set omega0 to [0,0,0]' and omegar0 to [0,0,0]'. Manual switch of slew ON. All other controller switches OFF.
Recommended simulation time: 2 or 3 minutes for slews, whatever you want for intertial pointing.

2)Detumbling
Set omega0 to [6,5,4]' (or any value of your choice) and omegar0 to [0,0,0]'. Manual switch of detumbling ON. All other controller switches OFF.
Recommended simulation time: 1 hour.

3)Desaturation
Set omega0 to [0,0,0]' and omegar to [+-20,+-20,+-20] or any value of your choice between 20 and -20. Manual switch of desaturation ON. All other controller switches OFF.
Recommended simulation time: 1 hour.

4)Earth pointing
uncomment lines 110 and 111 of the variable file (to set the right q0). Set omega0 to [0,0,0]' and omegar0 to [0,0,0]'. Manual switch of earth_pointing ON. All other controller switches OFF.
Recommended simulation time: variable "T" = approx a little more than 2 hours for a single orbit or whatever you want.

5)Free motion
Set omega0 to whatever you'd like and omegar0 to [0,0,0]'. Set all of the manual switches OFF.
Recommended simulation time: any.