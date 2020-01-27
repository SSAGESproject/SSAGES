set cell 20 0 0 0 20 0 0 0 20 
species sodium    Na_ONCV_PBE-1.0.xml
species chlorine  Cl_ONCV_PBE-1.0.xml
atom Na sodium 0.0 0.0 0.0
atom Cl chlorine 2.0 2.0 2.0
set ecut 35
set xc PBE
set wf_dyn PSDA
set scf_tol 1.e-5
# ground state calculation
run -atomic_density 0 200 10
# optimize geometry
set atoms_dyn CG
run 20 10 10
# prepare for MD driven by SSAGES
set atoms_dyn MD
set dt 10
randomize_v 300
# optional equilibration in Qbox prior to SSAGES
#run 50 10 0
