How to run SSAGES on LAMMPS
---------------------------

SSAGES can be compiled with LAMMPS using LAMMPS by pointing to a LAMMPS source via

	cmake -DLAMMPS_SRC=/path/to/LAMMPS/src
	
or

	cmake -DLAMMPS=ON

With the first option, the LAMMPS source directory is specified. In the second option 
the LAMMPS source file is downloaded automatically.  SSAGES is then compiled using 
	make

Using SSAGES with LAMMPS is as simple as running LAMMPS by itself. SSAGES requires the same 
files and commands that a stand-alone LAMMPS simulation requires, with one exception. 
LAMMPS input files usually contain a ‘run’ command. Instead of a run command, SSAGES requires an 
extra ‘fix’ in the LAMMPS input file:
	fix ssages all ssages

Additionally, SSAGES requires a .json file which specifies the method that will be used, 
collective variables, etc. The MD Engine should also be specified

	"type": "LAMMPS"


The number of MD steps of the simulation are specified in "MDSteps" of the .json file. LAMMPS log files for each 
simulation are specified in the "logfile" of the .json file. If they are specified as "none", no LAMMPS log
files will be generated. 

All values in the .json file will be in the units specified in the LAMMPS input file. 




