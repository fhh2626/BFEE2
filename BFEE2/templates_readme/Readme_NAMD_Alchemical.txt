To calculate the binding free energy:

1.1.  run 000_eq/000.1_eq.conf
1.2.  (optionally) run 000_eq/000.3_updateCenters.py
2.1.  run 001_MoleculeBound/001.1_fep_backward.conf
2.2.  run 001_MoleculeBound/001.2_fep_forward.conf
3.1.  run 002_RestraintBound/002.1_ti_backward.conf
3.2.  run 002_RestraintBound/002.2_ti_forward.conf
4.1.  CHARMM user:
        (if you didn't link BFEE with VMD before generating inputs)
        run 002.3_removeProtein.tcl using VMD
      Amber user:
        run 002.3_removeProtein.cpptraj using cpptraj
4.2   run 000_eq/000.2_eq_ligandOnly.conf
4.3.  run 003_MoleculeUnbound/003.1_fep_backward.conf
4.4.  run 003_MoleculeUnbound/003.2_fep_forward.conf
5.1.  run 004_RestraintUnbound/004.1_ti_backward.conf
5.2.  run 004_RestraintUnbound/004.2_ti_forward.conf
(steps 2-5 can be done in parallel)
6.    do post-treatment using BFEE