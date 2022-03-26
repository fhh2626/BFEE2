import sys

sys.path.append("..")
import os, shutil

import BFEE2.inputGenerator as inputGenerator

import commonTools
import userDefinedPath

targetFiles = [
    r'000_eq/000.1_eq.conf',
    r'000_eq/000.2_updateCenters.py',
    r'000_eq/colvars.in',
    r'001_RMSDBound/001_abf_1.conf',
    r'001_RMSDBound/001_abf_1.extend.conf',
    r'001_RMSDBound/001_abf_2.conf',
    r'001_RMSDBound/001_abf_2.extend.conf',
    r'001_RMSDBound/colvars_1.in',
    r'001_RMSDBound/colvars_2.in',
    r'002_EulerTheta/002_abf_1.conf',
    r'002_EulerTheta/002_abf_1.extend.conf',
    r'002_EulerTheta/002_abf_2.conf',
    r'002_EulerTheta/002_abf_2.extend.conf',
    r'002_EulerTheta/colvars_1.in',
    r'002_EulerTheta/colvars_2.in',
    r'003_EulerPhi/003_abf_1.conf',
    r'003_EulerPhi/003_abf_1.extend.conf',
    r'003_EulerPhi/003_abf_2.conf',
    r'003_EulerPhi/003_abf_2.extend.conf',
    r'003_EulerPhi/colvars_1.in',
    r'003_EulerPhi/colvars_2.in',
    r'004_EulerPsi/004_abf_1.conf',
    r'004_EulerPsi/004_abf_1.extend.conf',
    r'004_EulerPsi/004_abf_2.conf',
    r'004_EulerPsi/004_abf_2.extend.conf',
    r'004_EulerPsi/colvars_1.in',
    r'004_EulerPsi/colvars_2.in',
    r'005_PolarTheta/005_abf_1.conf',
    r'005_PolarTheta/005_abf_1.extend.conf',
    r'005_PolarTheta/005_abf_2.conf',
    r'005_PolarTheta/005_abf_2.extend.conf',
    r'005_PolarTheta/colvars_1.in',
    r'005_PolarTheta/colvars_2.in',
    r'006_PolarPhi/006_abf_1.conf',
    r'006_PolarPhi/006_abf_1.extend.conf',
    r'006_PolarPhi/006_abf_2.conf',
    r'006_PolarPhi/006_abf_2.extend.conf',
    r'006_PolarPhi/colvars_1.in',
    r'006_PolarPhi/colvars_2.in',
    r'007_r/007.1_eq.conf',
    r'007_r/007.2_abf_1.conf',
    r'007_r/007.2_abf_1.extend.conf',
    r'007_r/007.2_abf_2.conf',
    r'007_r/007.2_abf_2.extend.conf',
    r'007_r/colvars_1.in',
    r'007_r/colvars_2.in',
    r'007_r/colvars_eq.in',
    r'007_r/complex_largeBox.pdb',
    r'007_r/complex_largeBox.parm7',
    r'007_r/complex_largeBox.xyz',
    r'008_RMSDUnbound/008.0.1_removeProtein.cpptraj',
    r'008_RMSDUnbound/008.1_eq.conf',
    r'008_RMSDUnbound/008.2_abf_1.conf',
    r'008_RMSDUnbound/008.2_abf_1.extend.conf',
    r'008_RMSDUnbound/008.2_abf_2.conf',
    r'008_RMSDUnbound/008.2_abf_2.extend.conf',
    r'008_RMSDUnbound/colvars_1.in',
    r'008_RMSDUnbound/colvars_2.in',
    r'008_RMSDUnbound/ligandOnly.ndx',
    r'008_RMSDUnbound/ligandOnly.pdb',
    r'008_RMSDUnbound/ligandOnly.xyz',
    r'complex.ndx',
    r'complex.pdb',
    r'complex.parm7',
    r'complex.xyz',
]

refPath = r'exampleOutputs/003_trypsin_benzamidine_geometrical'
targetPath = r'test_output/003_trypsin_benzamidine_geometrical'

if os.path.exists(targetPath):
    shutil.rmtree(targetPath)
os.makedirs(targetPath)

iGenerator = inputGenerator.inputGenerator()
iGenerator.generateNAMDGeometricFiles(
    userDefinedPath.UNITTESTPATH + '/' + targetPath,
    userDefinedPath.UNITTESTPATH + '/' + r'exampleInputs/trypsin_benzamidine/small_box.parm7',
    userDefinedPath.UNITTESTPATH + '/' + r'exampleInputs/trypsin_benzamidine/small_box.rst7',
    'amber',
    [],
    300.0,
    'protein',
    'resname BAMI',
    'protein and (resid 116:118 or resid 139:142 or resid 163:165)',
    userDefinedPath.UNITTESTPATH + '/' + r'exampleInputs/trypsin_benzamidine/large_box.parm7',
    userDefinedPath.UNITTESTPATH + '/' + r'exampleInputs/trypsin_benzamidine/large_box.rst7',
    [2,2,2,2,2,2,2,2],
    pinDownPro = False,
    useOldCv = False,
)

commonTools.checkAllFiles(refPath, targetPath, targetFiles)
