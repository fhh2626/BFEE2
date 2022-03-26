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
    r'001_RMSDBound/001_abf_3.conf',
    r'001_RMSDBound/001_abf_3.extend.conf',
    r'001_RMSDBound/colvars_1.in',
    r'001_RMSDBound/colvars_2.in',
    r'001_RMSDBound/colvars_3.in',
    r'002_EulerTheta/002_abf_1.conf',
    r'002_EulerTheta/002_abf_1.extend.conf',
    r'002_EulerTheta/colvars_1.in',
    r'003_EulerPhi/003_abf_1.conf',
    r'003_EulerPhi/003_abf_1.extend.conf',
    r'003_EulerPhi/colvars_1.in',
    r'004_EulerPsi/004_abf_1.conf',
    r'004_EulerPsi/004_abf_1.extend.conf',
    r'004_EulerPsi/colvars_1.in',
    r'005_PolarTheta/005_abf_1.conf',
    r'005_PolarTheta/005_abf_1.extend.conf',
    r'005_PolarTheta/colvars_1.in',
    r'006_PolarPhi/006_abf_1.conf',
    r'006_PolarPhi/006_abf_1.extend.conf',
    r'006_PolarPhi/colvars_1.in',
    r'007_r/007.0_solvate.tcl',
    r'007_r/007.1_eq.conf',
    r'007_r/007.2_abf_1.conf',
    r'007_r/007.2_abf_1.extend.conf',
    r'007_r/007.2_abf_2.conf',
    r'007_r/007.2_abf_2.extend.conf',
    r'007_r/007.2_abf_3.conf',
    r'007_r/007.2_abf_3.extend.conf',
    r'007_r/007.2_abf_4.conf',
    r'007_r/007.2_abf_4.extend.conf',
    r'007_r/007.2_abf_5.conf',
    r'007_r/007.2_abf_5.extend.conf',
    r'007_r/colvars_1.in',
    r'007_r/colvars_2.in',
    r'007_r/colvars_3.in',
    r'007_r/colvars_4.in',
    r'007_r/colvars_5.in',
    r'007_r/colvars_eq.in',
    r'007_r/complex_largeBox.pdb',
    r'007_r/complex_largeBox.psf',
    r'007_r/complex_largeBox.xyz',
    r'008_RMSDUnbound/008.0.1_removeProtein.tcl',
    r'008_RMSDUnbound/008.0.2_neutrilize.tcl',
    r'008_RMSDUnbound/008.1_eq.conf',
    r'008_RMSDUnbound/008.2_abf_1.conf',
    r'008_RMSDUnbound/008.2_abf_1.extend.conf',
    r'008_RMSDUnbound/008.2_abf_2.conf',
    r'008_RMSDUnbound/008.2_abf_2.extend.conf',
    r'008_RMSDUnbound/008.2_abf_3.conf',
    r'008_RMSDUnbound/008.2_abf_3.extend.conf',
    r'008_RMSDUnbound/colvars_1.in',
    r'008_RMSDUnbound/colvars_2.in',
    r'008_RMSDUnbound/colvars_3.in',
    r'008_RMSDUnbound/ligandOnly.ndx',
    r'008_RMSDUnbound/ligandOnly.pdb',
    r'008_RMSDUnbound/ligandOnly.psf',
    r'008_RMSDUnbound/ligandOnly.xyz',
    r'complex.ndx',
    r'complex.pdb',
    r'complex.psf',
    r'complex.xyz',
    r'par_all36m_prot.prm',
    r'toppar_water_ions.prm',
]

refPath = r'exampleOutputs/001_abl_sh3_p41_geometrical'
targetPath = r'test_output/001_abl_sh3_p41_geometrical'

if os.path.exists(targetPath):
    shutil.rmtree(targetPath)
os.makedirs(targetPath)

iGenerator = inputGenerator.inputGenerator()
iGenerator.generateNAMDGeometricFiles(
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
    '',
    '',
    '',
    [3,1,1,1,1,1,5,3],
    vmdPath = userDefinedPath.VMDPATH
)

commonTools.checkAllFiles(refPath, targetPath, targetFiles)
