This example calculates the MFEP of a particle moving on a 2D energy surface with two wells (at [-1, -1] and [1, 1]) and a barrier at [0,0].
By default, this example runs with 16 images.  

```
mpirun -np 16 ./ssages FTS.json
```

Using the file Input_Generator.py, you can modify how many images the example will use and the initial string images.  

Parameters for the method can be modified in Template_Input.json to experiment with different spring constants, etc.

Lastly, the LAMMPS input files can be directly modified to change any underlying properties of the simulation.
