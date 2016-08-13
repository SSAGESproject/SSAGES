#RUN PARALLEL TEMPERING
1) lammps input file
replace temperature values with keyword ptemp in Butane_SSAGES.in

2) temeplate json file
set the parameters in the Template_Input.json
MDSteps = total timesteps to run
frequency = timestep frequency for swaping configurations
CVS = it is a dummpy CVs and does have any influence on simulation

3) generate json file
set the parameters in the ParallelTemp_Input_Generator.py
nwalker = number of walkers
min_temp = minimum temperature
max_temp = maximum temperature

run
python ParallelTemp_Input_Generator.py

5) run SSAGES
mpirun -np nwalker ../../../build/ssages ParallelTemp.json

