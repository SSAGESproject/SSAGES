This example demonstrates the CFF method by calculating the potential of mean force
alanine dipeptide in water, using Gromacs.

Inside the directory `1walker`, run CFF using a single walker.

To run this example, execute:

```
{SSAGES_path}/build/hooks/gromacs/gromacs/bin/gmx_mpi grompp -f npt.mdp -p topol.top -c npt.gro -o adp.tpr
mpirun -np ${nprocessors} {SSAGES_path}/build/ssages CFF.json
```

This run will produce the main output file called `CFF.out`, similar to the file
`CFF.out1000` located in this folder.

Please see the SSAGES documentation under `Methods/CFF` for more information on definition
of input variables.

**Additional notes**

The Gromacs input files can be modified to change any underlying properties of the
simulation. After changing any of these parameters, the `.tpr` file must be regenerated
using `gmx grompp`, as follows:

```
{SSAGES_path}/build/hooks/gromacs/gromacs/bin/gmx_mpi grompp -f npt.mdp -p topol.top -c npt.gro -o adp.tpr
```

Free energy surface or the potential of mean force can be plotted from one of the CFF
output data files, `CFF.out`. First two columns in `CFF.out` are the peptide torsional
angles, phi and psi (radians), and the last column is the potential of mean force (kJ/mol)
