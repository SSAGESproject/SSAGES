Umbrella Sampling
============
This is an example of how to run **Umbrella Sampling** using **SSAGES** and **JSON**.

The template file illustrates how you would set up a single umbrella for a single run.

<a id="Generating Input From Templates"></a>
##Generating Input From Templates
The Umbrella_Input_Generator.py is a python script that will take the template file in as a JSON file, change it into a python dictionary. 
 Python allows a high level of flexibility for the user to modify their python scripts. 
 In the end the dictionary is converted back into a JSON file

 For example Umbrella_Input_Generator.py takes Template_Input.json and creates Umbrella.json which can be used with SSAGES.

<a id="Torsional Angle Butane"></a>
##Torsional Angle Butane
 An example LAMMPS file (Butane_SSAGES.in) is provided, which reads the initial structure from Butane_VMD.data.

<a id="Running SSAGES"></a>
##Running SSAGES
 mpiexec -np x ./ssages Umbrella.json

 where x is the number of processors defined in the input file. This example creates 12 different drivers each using one umbrella with one processor, so the total number of processors is 12.

 mpiexec -np 12 ./ssages Umbrella.json


<a id="Results"></a>
##Results
Each driver will print out the results for each step:
iteration umbrella_center cv_value