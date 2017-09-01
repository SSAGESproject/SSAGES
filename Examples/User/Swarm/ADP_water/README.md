This example calculates the MFEP of the isomerization of alanine dipeptide in water, using LAMMPS.  
To run the example, execute:

```
mpirun -np 22 ./ssages Swarm.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the LAMMPS input files can be directly modified to change any underlying properties of the simulation.

An image of the string can be overlaid on a free energy surface (data contained in F_out) by executing:
```
python plotter.py 22 abf
```
