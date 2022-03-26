import sys

sys.path.append("..")
import os, shutil

import BFEE2.inputGenerator as inputGenerator

import commonTools
import userDefinedPath

targetFiles = [
    r'000_eq/000.1_eq.conf',
    r'000_eq/000.2_eq_ligandOnly.conf',
    r'000_eq/colvars.in',
    r'000_eq/colvars_ligandOnly.in',
    r'001_MoleculeBound/001.1_fep_backward.conf',
    r'001_MoleculeBound/001.2_fep_forward.conf',
    r'001_MoleculeBound/colvars.in',
    r'002_RestraintBound/002.1_ti_backward.conf',
    r'002_RestraintBound/002.2_ti_forward.conf',
    r'002_RestraintBound/colvars_backward.in',
    r'002_RestraintBound/colvars_forward.in',
    r'003_MoleculeUnbound/003.1_fep_backward.conf',
    r'003_MoleculeUnbound/003.2_fep_forward.conf',
    r'003_MoleculeUnbound/colvars.in',
    r'004_RestraintUnbound/004.1_ti_backward.conf',
    r'004_RestraintUnbound/004.2_ti_forward.conf',
    r'004_RestraintUnbound/colvars_backward.in',
    r'004_RestraintUnbound/colvars_forward.in',
    r'002.3.1_removeProtein.tcl',
    r'002.3.2_neutrilize.tcl',
    r'complex.ndx',
    r'complex.pdb',
    r'complex.psf',
    r'complex.xyz',
    r'fep.pdb',
    r'fep.tcl',
    r'fep_ligandOnly.pdb',
    r'ligandOnly.ndx',
    r'ligandOnly.pdb',
    r'ligandOnly.psf',
    r'ligandOnly.xyz',
    r'par_all36m_prot.prm',
    r'toppar_water_ions.prm'
]

refPath = r'exampleOutputs/002_abl_sh3_p41_alchemical'
targetPath = r'test_output/002_abl_sh3_p41_alchemical'

if os.path.exists(targetPath):
    shutil.rmtree(targetPath)
os.makedirs(targetPath)

iGenerator = inputGenerator.inputGenerator()
iGenerator.generateNAMDAlchemicalFiles(
    userDefinedPath.UNITTESTPATH + '/' + targetPath,
    userDefinedPath.UNITTESTPATH + '/' + r'exampleInputs/abl_sh3_p41/complex.psf',
    userDefinedPath.UNITTESTPATH + '/' + r'exampleInputs/abl_sh3_p41/complex.pdb',
    'charmm',
    [
        userDefinedPath.UNITTESTPATH + '/' + r'exampleInputs/abl_sh3_p41/par_all36m_prot.prm',
        userDefinedPath.UNITTESTPATH + '/' + r'exampleInputs/abl_sh3_p41/toppar_water_ions.prm'
    ],
    300.0,
    'segid SH3D',
    'segid PPRO',
    [100,100,100,50],
    vmdPath = userDefinedPath.VMDPATH
)

commonTools.checkAllFiles(refPath, targetPath, targetFiles)
