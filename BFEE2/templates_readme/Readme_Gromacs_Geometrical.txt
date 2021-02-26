To calculate the binding free energy:

1.1.  run 000_eq/000_generate_tpr.sh
1.2.  run 000_eq/000_eq.tpr
2.1.  run 001_RMSD_bound/001_generate_tpr.sh
2.2.  run 001_RMSD_bound/001_RMSD_bound.tpr
3.1.  run 002_euler_theta/002_generate_tpr.sh
3.2.  run 002_euler_theta/002_euler_theta.tpr
4.1.  run 003_euler_phi/003_generate_tpr.sh
4.2.  run 003_euler_phi/003_euler_phi.tpr
5.1.  run 004_euler_psi/004_generate_tpr.sh
5.2.  run 004_euler_psi/004_euler_psi.tpr
6.1.  run 005_polar_theta/005_generate_tpr.sh
6.2.  run 005_polar_theta/005_polar_theta.tpr
7.1.  run 006_polar_phi/006_generate_tpr.sh
7.2.  run 006_polar_phi/006_polar_phi.tpr
8.1.  run 007_r/007.1_generate_eq_tpr.sh
8.2.  run 007_r/007_r.eq.tpr
8.3.  run 007_r/007.2_generate_tpr.sh
8.4   run 007_r/007_r.tpr
9.1.  run 008_RMSD_unbound/008.1_generate_eq_tpr.sh
9.2.  run 008_RMSD_unbound/008_RMSD_unbound.eq.tpr
9.3.  run 008_RMSD_unbound/008.2_generate_tpr.sh
9.4.  run 008_RMSD_unbound/008_RMSD_unbound.tpr
10.   do post-treatment using BFEE