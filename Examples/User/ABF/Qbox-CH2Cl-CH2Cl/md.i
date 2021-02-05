set cell 20 0 0 0 20 0 0 0 20 
species chlorine Cl_ONCV_PBE-1.2.xml
species carbon   C_ONCV_PBE-1.2.xml
species hydrogen H_ONCV_PBE-1.2.xml
atom  Cl1  chlorine   6.818780755130579 12.505473751842473 9.408518840470162
atom  Cl2  chlorine   13.185234513776031 7.5052847802260105 10.591481159529838
atom  C1   carbon     8.59889338221399 9.714362976680903 10.024566310140218         
atom  C2   carbon     11.40512188669262 10.29450583922295 9.975433689859782
atom  H1   hydrogen   8.009301938848786 8.996270834120715 11.870819002985751
atom  H2   hydrogen   8.090559733927964 8.323531879511698 8.580823160361314
atom  H3   hydrogen   11.913455534978645 11.687226652556786 11.419176839638686
atom  H4   hydrogen   11.996603046222457 11.014487697947768 8.129180997014249

set ecut 35
set xc PBE
set wf_dyn PSDA
set scf_tol 1.e-6

# ground state calculation
run -atomic_density 0 100 10
torsion Cl1 C1 C2 Cl2

# optimize geometry
set atoms_dyn CG
run 100 10
torsion Cl1 C1 C2 Cl2
save gs.xml

# prepare for MD driven by SSAGES
set atoms_dyn MD
set scf_tol 1.e-6
set dt 20
rseed
randomize_v 300
reset_rotation
set thermostat BDP
set th_temp 150
# equilibration in Qbox prior to SSAGES
run 500 20 5
