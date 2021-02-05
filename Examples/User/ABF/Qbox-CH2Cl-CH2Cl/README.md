This example demonstrates the free energy surface (FES) of the Cl-C-C-Cl
torsional angle in 1,2-dichloroethane (CH2Cl-CH2Cl) in vacuum, using Qbox as
the engine and ABF as the method.

This run will produce three output files: F_out, Fworld_cv0, Nworld.
The FES can be created by running the script from the Tools folder:

    python /path/to/ssages/Tools/ABF_integrator.py -i F_out -o G -p True

**Additional notes**

In the example input.json, `"md_iterations"` is set to be 80000, which can
create a rough estimate of the FES with respect to the torsional angle.
To get a more accurate estimate of the FES, increase `"md_iterations"` to
improve sampling.



