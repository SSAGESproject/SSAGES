set xc PBE
set cell 40 0 0 0 40 0 0 0 40 
species sodium    Na_ONCV_PBE-1.0.xml
species chlorine  Cl_ONCV_PBE-1.0.xml
move Na to 0. 0. 0.
move Cl to 2.0 2.0 2.0 
set wf_dyn PSDA
set scf_tol 1.e-5
set atoms_dyn MD
set dt 10
set thermostat BDP
set th_temp 300
randomize_v 300
run 1 20 1  
