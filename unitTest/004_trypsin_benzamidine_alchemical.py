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
    r'002.3.1_removeProtein.cpptraj',
    r'complex.ndx',
    r'complex.pdb',
    r'complex.parm7',
    r'complex.xyz',
    r'fep.pdb',
    r'fep.tcl',
    r'fep_ligandOnly.pdb',
    r'ligandOnly.ndx',
    r'ligandOnly.pdb',
    r'ligandOnly.xyz',
]

refPath = r'exampleOutputs/004_trypsin_benzamidine_alchemical'
targetPath = r'test_output/004_trypsin_benzamidine_alchemical'

if os.path.exists(targetPath):
    shutil.rmtree(targetPath)
os.makedirs(targetPath)

iGenerator = inputGenerator.inputGenerator()
iGenerator.generateNAMDAlchemicalFiles(
    userDefinedPath.UNITTESTPATH + '/' + targetPath,
    userDefinedPath.UNITTESTPATH + '/' + r'exampleInputs/trypsin_benzamidine/small_box.parm7',
    userDefinedPath.UNITTESTPATH + '/' + r'exampleInputs/trypsin_benzamidine/small_box.rst7',
    'amber',
    [],
    300.0,
    'protein',
    'resname BAMI',
    [100,100,50,50],
    pinDownPro = False,
    useOldCv = False,
)

commonTools.checkAllFiles(refPath, targetPath, targetFiles)
