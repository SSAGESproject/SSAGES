This example calculates the MFEP of the isomerization of alanine dipeptide in water, using Gromacs.  
First, you must copy the .tpr file obtained from `gmx grompp`.  The included Python script can be run as:

`python copytpr.py`

To do this for you.

To run the example, execute:

```
mpirun -np 22 ./ssages Swarm.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the Gromacs input files can be modified to change any underlying properties of the simulation.  
After changing any of these parameters, the .tpr file must be regenerated using `gmx grompp`, as follow:

```
gmx_mpi grompp -f nvt.mdp -c adp_water.gro -p topol_water.top -o adp_H2O.tpr
```

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf
```
```
