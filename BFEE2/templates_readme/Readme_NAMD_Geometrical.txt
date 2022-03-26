To calculate the binding free energy:

1.1   run 000_eq/000.1_eq.conf
1.2.  (optionally) run 000_eq/000.2_updateCenters.py
      (If this is done, then there is usually no need to change restraining centers
       after steps 3.x.)
2.    run 001_RMSDBound/001_abf_1.conf
3.1.  run 002_EulerTheta/002_abf_1.conf
3.2.  change restraining center of "eulerTheta" in
      003_EulerPhi/colvars_1.in to
      the CV value corresponding to dG=0 in 
      002_EulerTheta/output/abf_1.czar.pmf,
3.3.  run 003_EulerPhi/003_abf_1.conf
3.4.  similarly, change "eulerTheta" and "eulerPhi" in 
      004_EulerPsi/colvars_1.in
3.5.  run 004_EulerPsi/004_abf_1.conf
3.6.  change 005_PolarTheta/colvars_1.in correspondingly
3.7.  run 005_PolarTheta/005_abf_1.conf
3.8.  change 006_PolarPhi/colvars_1.in correspondingly
3.9.  run 006_PolarPhi/006_abf_1.conf
3.10. change 007_r/colvars_eq.in and 007_r/colvars_1.in correspondingly
3.11. (if you didn't link BFEE with VMD before generating inputs)
      run 007_r/007.0_solvate.tcl using VMD
3.12. run 007_r/007.1_eq.conf
3.13. run 007_r/007.2_abf_1.conf
4.1.  CHARMM user:
        (if you didn't link BFEE with VMD before generating inputs)
        run 008_RMSDUnbound/008.0_removeProtein.tcl using VMD
      Amber user:
        run 008.0_removeProtein.cpptraj using cpptraj
4.2.  run 008_RMSDUnbound/008.1_eq.conf
4.3.  run 008_RMSDUnbound/008.2_abf_1.conf
(steps 2-4 can be done in parallel)
5.    do post-treatment using BFEE