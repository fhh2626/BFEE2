# generate all the inputs and define corresponding slots

import os
import shutil
import subprocess
import sys

import numpy as np

from BFEE2.commonTools import fileParser
from BFEE2.templates_gromacs.BFEEGromacs import BFEEGromacs
from BFEE2.templates_namd import configTemplate, scriptTemplate

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from BFEE2 import templates_namd, templates_readme


# an runtime error
# directory already exists
class DirectoryExistError(RuntimeError):
    def __init__(self, arg):
        self.args = arg

# an runtime error
# file type is unknown
class FileTypeUnknownError(RuntimeError):
    def __init__(self, arg):
        self.args = arg

class inputGenerator():
    """generate all the inputs and define corresponding slots
    """

    def __init__(self):
        """simply initialize the configTemplate objects
        """

        self.cTemplate = configTemplate.configTemplate()

    def generateNAMDAlchemicalFiles(
        self,
        path,
        topFile,
        coorFile,
        forceFieldType,
        forceFieldFiles,
        temperature,
        selectionPro,
        selectionLig,
        stratification = [1,1,1,1],
        doubleWide = False,
        minBeforeSample = False,
        membraneProtein = False,
        neutralizeLigOnly = 'NaCl',
        pinDownPro = True,
        useOldCv = True,
        vmdPath = '',
        OPLSMixingRule = False
    ):
        """generate all the input files for NAMD alchemical simulation

        Args:
            path (str): the directory for generation of all the files
            topFile (str): the path of the topology (psf, parm) file
            coorFile (str): the path of the coordinate (pdb, rst7) file
            forceFieldType (str): 'charmm' or 'amber'
            forceFieldFiles (list of str): list of CHARMM force field files
            temperature (float): temperature of the simulation
            selectionPro (str): MDAnalysis-style selection of the protein
            selectionLig (str): MDAnalysis-style selection of the ligand
            stratification (list of int, 8, optional): number of windows for each simulation.
                                                       Defaults to [1,1,1,1].
            doubleWide (bool, optional): whether double-wide simulations are carried out. 
                                         Defaults to False.
            minBeforeSample (bool, optional): minimization before sampling in each FEP window.
                                              Defaults to False.
            membraneProtein (bool, optional): whether simulating a membrane protein.
                                              Defaults to False.
            neutralizeLigOnly (str or None, optional): 'NaCl', 'KCl', 'CaCl2' or None.
                                                       neutralize the lig-only system using the salt.
                                                       Defaluts to NaCl.
            pinDownPro (bool, optional): whether pinning down the protien. Defaults to True.
            useOldCv (bool, optional): whether used old, custom-function-based cv. Defaults to True.
            vmdPath (str, optional): path to vmd. Defaults to ''.
            OPLSMixingRule (bool, optional): whether use the OPLS mixing rules. Defaults to False.
        """

        assert(len(stratification) == 4)
        assert(forceFieldType == 'charmm' or forceFieldType == 'amber')

        # determine the type of input top and coor file
        topType, coorType = self._determineFileType(topFile, coorFile)

        self._makeDirectories(path, 'alchemical')
        # after copying files
        # the coordinate file is converted to pdb
        self._copyFiles(
            path, topFile, topType, coorFile, coorType, forceFieldType, forceFieldFiles,
            selectionPro, selectionLig, selectionPro, '', '',
            'alchemical', membraneProtein, neutralizeLigOnly, vmdPath)

        # get relative force field path
        relativeFFPath = []
        for item in forceFieldFiles:
            _, name = os.path.split(item)
            relativeFFPath.append(f'../{name}')
        
        self._generateAlchemicalNAMDConfig(
            path, forceFieldType, relativeFFPath, temperature, stratification, doubleWide, minBeforeSample,
            membraneProtein, OPLSMixingRule=OPLSMixingRule
        )
        self._generateAlchemicalColvarsConfig(
            path, topType, 'pdb', selectionPro, selectionLig, selectionPro, stratification, pinDownPro, useOldCv
        )

    def generateNAMDGeometricFiles(
        self,
        path,
        topFile,
        coorFile,
        forceFieldType,
        forceFieldFiles,
        temperature,
        selectionPro,
        selectionLig,
        selectionRef = '',
        userProvidedPullingTop = '',
        userProvidedPullingCoor = '',
        stratification = [1,1,1,1,1,1,1,1],
        membraneProtein = False,
        neutralizeLigOnly = 'NaCl',
        pinDownPro = True,
        useOldCv = True,
        parallelRuns = 1,
        vmdPath = '',
        reflectionBoundary = False,
        MDEngine = 'namd',
        OPLSMixingRule = False,
        considerRMSDCV = True,
        GaWTM = False
    ):
        """generate all the input files for NAMD geometric simulation

        Args:
            path (str): the directory for generation of all the files
            topFile (str): the path of the topology (psf, parm) file
            coorFile (str): the path of the coordinate (pdb, rst7) file
            forceFieldType (str): 'charmm' or 'amber'
            forceFieldFiles (list of str): list of CHARMM force field files
            temperature (float): temperature of the simulation
            selectionPro (str): MDAnalysis-style selection of the protein
            selectionLig (str): MDAnalysis-style selection of the ligand
            selectionRef (str, optional): MDAnalysis-style selection of the reference group for pulling simulation,
                                          by default, this is the protein. 
                                          Defaults to ''.
            userProvidedPullingTop (str, optional): user-provided large solvation box for pulling simulation.
                                                    Defaults to ''.
            userProvidedPullingCoor (str, optional): user-provided large solvation box for pulling simulation. 
                                                     Defaults to ''.
            stratification (list, optional): number of windows for each simulation. 
                                             Defaults to [1,1,1,1,1,1,1,1].
            membraneProtein (bool, optional): whether simulation a membrane protein. Defaults to False.
            neutralizeLigOnly (str or None, optional): 'NaCl', 'KCl', 'CaCl2' or None.
                                                       neutralize the lig-only system using the salt.
                                                       Defaluts to NaCl.
            pinDownPro (bool, optional): whether pinning down the protien. Defaults to True.
            useOldCv (bool, optional): whether used old, custom-function-based cv. Defaults to True.
            parallelRuns (int, optional): generate files for duplicate runs. Defaults to 1.
            vmdPath (str, optional): path to vmd. Defaults to ''.
            reflectingBoundary (bool, optional): Whether use reflecting boundaries, requires setBoundary on. Default to False.
            MDEngine (str, optional): namd or gromacs. Default to namd.
            OPLSMixingRule (bool, optional): whether use the OPLS mixing rules. Defaults to False.
            considerRMSDCV (bool, optional): Whethre consider the RMSD CV. Default to True.
            GaWTM (bool, optional): Whether performing GaWTM-eABF simulations. Default to False
        """ 

        assert(len(stratification) == 8)
        assert(forceFieldType == 'charmm' or forceFieldType == 'amber')

        if selectionRef == '':
            selectionRef = selectionPro

        # determine the type of input top and coor file
        # However, whatever coorType is, NAMD only read pdb as coordinate
        topType, coorType = self._determineFileType(topFile, coorFile)

        self._makeDirectories(path, 'geometric', considerRMSDCV=considerRMSDCV)
        self._copyFiles(
            path, topFile, topType, coorFile, coorType, forceFieldType, forceFieldFiles,
            selectionPro, selectionLig, selectionRef, userProvidedPullingTop, userProvidedPullingCoor,
            'geometric', membraneProtein, neutralizeLigOnly, vmdPath, considerRMSDCV=considerRMSDCV)
        
        if MDEngine == 'gromacs':
            self._makeGromacsTopGro(path, forceFieldType, forceFieldFiles, 'geometric', vmdPath, considerRMSDCV=considerRMSDCV)

        # get relative force field path
        relativeFFPath = []
        for item in forceFieldFiles:
            _, name = os.path.split(item)
            relativeFFPath.append(f'../{name}')
        
        if MDEngine == 'namd':
            self._generateGeometricNAMDConfig(
                path, forceFieldType, relativeFFPath, temperature, stratification, membraneProtein,
                OPLSMixingRule=OPLSMixingRule, considerRMSDCV=considerRMSDCV, GaWTM=GaWTM
            )
        elif MDEngine == 'gromacs':
            self._generateGeometricGromacsConfig(
                path, forceFieldType, temperature, membraneProtein, considerRMSDCV=considerRMSDCV
            )
            
        self._generateGeometricColvarsConfig(
            path, topType, coorType, selectionPro, selectionLig, selectionRef, 
            stratification, pinDownPro, useOldCv, reflectionBoundary, MDEngine,
            considerRMSDCV=considerRMSDCV, reweightAMD=False
        )
        # for GaWTM-eABF simulations, two colvars config files are needed. One for pre-equilibration
        # and the other for production
        if GaWTM:
            self._generateGeometricColvarsConfig(
                path, topType, coorType, selectionPro, selectionLig, selectionRef, 
                stratification, pinDownPro, useOldCv, reflectionBoundary, MDEngine,
                considerRMSDCV=considerRMSDCV, reweightAMD=True
            )

        self._duplicateFileFolder(path, parallelRuns)

    def generateGromacsGeometricFiles(
        self, 
        path, 
        topFile, 
        pdbFile, 
        pdbFileFormat,
        ligandOnlyTopFile, 
        ligandOnlyPdbFile,
        ligandOnlyPdbFileFormat,
        selectionPro,
        selectionLig,
        selectionSol='resname TIP3* or resname SPC* or resname HOH or resname WAT',
        temperature=300.0
    ):
        """generate all the input files for Gromacs Geometric simulation
           This function is based on BFEEGromacs.py
           contributed by Haochuan Chen (yjcoshc_at_mail.nankai.edu.cn)
           
        Args:
            path (str): the directory for generation of all the files
            topFile (str): the path of the topology (top) file
            pdbFile (str): the path of the coordinate (pdb) file
            ligandOnlyTopFile (str): the path of the ligand-only topology (top) file
            ligandOnlyPdbFile (str): the path of the ligand-only coordinate (pdb) file
            selectionPro (str): MDAnalysis-style selection of the protein
            selectionLig (str): MDAnalysis-style selection of the ligand
            selectionSol (str, optional): MDAnalysis-style selection of the solvent. 
                                          Defaults to 'resname TIP3* or resname SPC*'.
            temperature (float, optional): Temperature of the simulation. Defaults to 300.0.

        Raises:
            DirectoryExistError: when {path}/BFEE exists
        """        
        
        # if the folder already exists, raise an error
        if os.path.exists(f'{path}/BFEE'):
            raise DirectoryExistError('Directory exists')
        os.mkdir(f'{path}/BFEE')

        # readme file 
        with pkg_resources.path(templates_readme, 'Readme_Gromacs_Geometrical.txt') as p:
            shutil.copyfile(p, f'{path}/BFEE/Readme.txt')

        bfee = BFEEGromacs(
            structureFile=pdbFile,
            topologyFile=topFile,
            ligandOnlyStructureFile=ligandOnlyPdbFile,
            ligandOnlyTopologyFile=ligandOnlyTopFile,
            baseDirectory=f'{path}/BFEE',
            structureFormat=pdbFileFormat,
            ligandOnlyStructureFileFormat=ligandOnlyPdbFileFormat
        )
        bfee.setProteinHeavyAtomsGroup(f'{selectionPro} and not (name H*)')
        bfee.setLigandHeavyAtomsGroup(f'{selectionLig} and not (name H*)')
        bfee.setSolventAtomsGroup(selectionSol)
        bfee.setTemperature(temperature)
        bfee.generate000()
        bfee.generate001()
        bfee.generate002()
        bfee.generate003()
        bfee.generate004()
        bfee.generate005()
        bfee.generate006()
        bfee.generate007()
        bfee.generate008()


    def _determineFileType(self, topFile, coorFile):
        """determine the file type of topology and coordinate files

        Args:
            topFile (str): the path of the topology (psf, parm) file
            coorFile (str): the path of the coordinate (pdb, rst7) file

        Raises:
            FileTypeUnknownError: if the postfix of topFile is not 
                                  psf, parm, prm7, parm7, prmtop
                                  or the postfix of coorFile is not
                                  pdb, coor, rst7, rst, inpcrd

        Returns:
            tuple of str: (topType, coorType), e.g., (psf, pdb)
        """

        topPostfix = os.path.splitext(topFile)[-1]
        coorPostfix = os.path.splitext(coorFile)[-1]

        topType = ''
        coorType = ''

        if topPostfix == '.psf':
            topType = 'psf'
        elif topPostfix == '.parm' or topPostfix == '.prm7' or topPostfix == '.parm7' or topPostfix == '.prmtop':
            topType = 'parm7'

        if coorPostfix == '.pdb':
            coorType = 'pdb'
        elif coorPostfix == '.coor':
            coorType = 'namdbin'
        elif coorPostfix == '.rst7' or coorPostfix == '.rst' or coorPostfix == '.inpcrd':
            coorType = 'inpcrd'

        if topType == '' or coorType == '':
            raise FileTypeUnknownError('File type unknown')

        return topType, coorType

    def _makeDirectories(self, path, jobType='geometric', considerRMSDCV=True):
        """make directories for BFEE calculation

        Args:
            path (str): the path for putting BFEE input files into
            jobType (str, optional): geometric or alchemical. Defaults to 'geometric'.
            considerRMSDCV (bool, optional): Whethre consider the RMSD CV. Default to True.

        Raises:
            DirectoryExistError: if {path}/BFEE exists
        """        

        # if the folder already exists, raise an error
        if os.path.exists(f'{path}/BFEE'):
            raise DirectoryExistError('Directory exists')
        os.mkdir(f'{path}/BFEE')
        os.mkdir(f'{path}/BFEE/000_eq')
        os.mkdir(f'{path}/BFEE/000_eq/output')
        if jobType == 'geometric':
            if considerRMSDCV:
                os.mkdir(f'{path}/BFEE/001_RMSDBound')
                os.mkdir(f'{path}/BFEE/001_RMSDBound/output')
                os.mkdir(f'{path}/BFEE/008_RMSDUnbound')
                os.mkdir(f'{path}/BFEE/008_RMSDUnbound/output')
            
            os.mkdir(f'{path}/BFEE/002_EulerTheta')
            os.mkdir(f'{path}/BFEE/003_EulerPhi')
            os.mkdir(f'{path}/BFEE/004_EulerPsi')
            os.mkdir(f'{path}/BFEE/005_PolarTheta')
            os.mkdir(f'{path}/BFEE/006_PolarPhi')
            os.mkdir(f'{path}/BFEE/007_r')

            os.mkdir(f'{path}/BFEE/002_EulerTheta/output')
            os.mkdir(f'{path}/BFEE/003_EulerPhi/output')
            os.mkdir(f'{path}/BFEE/004_EulerPsi/output')
            os.mkdir(f'{path}/BFEE/005_PolarTheta/output')
            os.mkdir(f'{path}/BFEE/006_PolarPhi/output')
            os.mkdir(f'{path}/BFEE/007_r/output')
            

        if jobType == 'alchemical':
            os.mkdir(f'{path}/BFEE/001_MoleculeBound')
            os.mkdir(f'{path}/BFEE/002_RestraintBound')
            os.mkdir(f'{path}/BFEE/003_MoleculeUnbound')
            os.mkdir(f'{path}/BFEE/004_RestraintUnbound')
            os.mkdir(f'{path}/BFEE/001_MoleculeBound/output')
            os.mkdir(f'{path}/BFEE/002_RestraintBound/output')
            os.mkdir(f'{path}/BFEE/003_MoleculeUnbound/output')
            os.mkdir(f'{path}/BFEE/004_RestraintUnbound/output')


    def _copyFiles(
        self, 
        path, 
        topFile, 
        topType, 
        coorFile, 
        coorType, 
        forceFieldType, 
        forceFieldFiles,
        selectionPro,
        selectionLig,
        selectionRef,
        userProvidedPullingTop = '',
        userProvidedPullingCoor = '',
        jobType = 'geometric',
        membraneProtein = False,
        neutralizeLigOnly = 'NaCl',
        vmdPath = '',
        considerRMSDCV = True
    ):
        """copy original and generate necessary topology/structure files

        Args:
            path (str): the directory for generation of all the files
            topFile (str): the path of the topology (psf, parm) file
            topType (str): the type (psf, parm) of the topology file
            coorFile (str): the path of the coordinate (pdb, rst7) file
            coorType (str): the type (pdb, rst) of the coordinate file
            forceFieldType (str): 'charmm' or 'amber'
            forceFieldFiles (list of str): list of CHARMM force field files
            selectionPro (str): MDAnalysis-style selection of the protein
            selectionLig (str): MDAnalysis-style selection of the ligand
            selectionRef (str): MDAnalysis-style selection of the reference group for pulling simulation
            userProvidedPullingTop (str, optional): user-provided large solvation box for pulling simulation. 
                                                    Defaults to ''.
            userProvidedPullingCoor (str, optional): user-provided large solvation box for pulling simulation. 
                                                     Defaults to ''.
            jobType (str, optional): 'geometric' or 'alchemical'. Defaults to 'geometric'.
            membraneProtein (bool, optional): whether simulating a membrane protein. Defaults to False.
            neutralizeLigOnly (str or None, optional): 'NaCl', 'KCl', 'CaCl2' or None.
                                                       neutralize the lig-only system using the salt.
                                                       Defaluts to NaCl.
            vmdPath (str, optional): path to vmd, space is forbidden. Defaults to ''.
            considerRMSDCV (bool, optional): Whethre consider the RMSD CV. Default to True.
        """
        
        # copy force fields
        if forceFieldType == 'charmm':
            for item in forceFieldFiles:
                _, fileName = os.path.split(item)
                shutil.copyfile(item, f'{path}/BFEE/{fileName}')

        # read the original topology and coordinate file
        fParser = fileParser.fileParser(topFile, coorFile)
            
        # check polar angles
        # if theta < 30 or > 150, then rotate 90 degrees to avoid polar angle invarience
        if (fParser.measurePolarAngles(selectionRef, selectionLig)[0] > 150 \
            or fParser.measurePolarAngles(selectionRef, selectionLig)[0] < 30):
            fParser.rotateSystem('x', 90)

        # normalize/centerize coordinate
        fParser.centerSystem()
        
        # write the new topology and structure file
        fParser.saveFile(
            'all', f'{path}/BFEE/complex.pdb', 'pdb', True,
            f'{path}/BFEE/complex.{topType}'
        )
        # xyz file
        fParser.saveFile(
            'all', f'{path}/BFEE/complex.xyz', 'xyz'
        )
        # ndx file
        fParser.saveNDX(
            [selectionPro, selectionLig, selectionRef], ['protein', 'ligand', 'reference'], 
            f'{path}/BFEE/complex.ndx', True
        )
        # cv definition file
        #shutil.copyfile(f'{sys.path[0]}/scripts/CVs.tcl', f'{path}/BFEE/CVs.tcl')
        
        # determine cation and anion
        if neutralizeLigOnly == 'NaCl':
            cation = 'SOD'
            anion = 'CLA'
        elif neutralizeLigOnly == 'KCl':
            cation = 'POT'
            anion = 'CLA'
        elif neutralizeLigOnly == 'CaCl2':
            cation = 'CAL'
            anion = 'CLA'
        else:
            neutralizeLigOnly = None

        if jobType == 'geometric':

            # readme file 
            with pkg_resources.path(templates_readme, 'Readme_NAMD_Geometrical.txt') as p:
                shutil.copyfile(p, f'{path}/BFEE/Readme.txt')
                
            # script to update the center after equilibration
            with pkg_resources.path(templates_namd, 'updateCenters.py') as p:
                shutil.copyfile(p, f'{path}/BFEE/000_eq/000.2_updateCenters.py')

            # remove protein for step 8
            # this cannot be done in pure python
            if considerRMSDCV:
                if not membraneProtein:
                    fParser.saveFile(
                        f'not {selectionPro}', f'{path}/BFEE/008_RMSDUnbound/ligandOnly.pdb', 'pdb'
                    )
                else:
                    # membrane protein
                    fParser.saveFile(
                        f'{selectionLig}', f'{path}/BFEE/008_RMSDUnbound/ligandOnly.pdb', 'pdb'
                    )
                # connect to VMD
                if forceFieldType == 'charmm':
                    if not membraneProtein:
                        with open( f'{path}/BFEE/008_RMSDUnbound/008.0.1_removeProtein.tcl', 'w') as rScript:
                            rScript.write(
                                scriptTemplate.removeProteinTemplate.substitute(
                                    path='../complex', selectionPro=f'{selectionPro}'.replace('segid', 'segname'),
                                    outputPath=f'./ligandOnly'
                                )
                            )
                    else:
                        # membrane protein
                        with open( f'{path}/BFEE/008_RMSDUnbound/008.0.1_removeProtein.tcl', 'w') as rScript:
                            rScript.write(
                                scriptTemplate.removeMemProteinTemplate.substitute(
                                    path='../complex', selectionLig=f'{selectionLig}'.replace('segid', 'segname'),
                                    outputPath=f'./ligandOnly'
                                )
                            )
                    
                    # neutralization
                    if neutralizeLigOnly is not None:
                        with open(f'{path}/BFEE/008_RMSDUnbound/008.0.2_neutrilize.tcl', 'w') as rScript:
                            rScript.write(
                                scriptTemplate.neutralizeSystempTemplate.substitute(
                                    path='./ligandOnly', cationName=cation, anionName=anion,
                                    extraCommand=''
                                )
                            )
                            
                    # if vmd path is defined
                    # then execute vmd automatically
                    if vmdPath != '':
                        subprocess.run(
                            [vmdPath, '-dispdev', 'text', '-e', f'{path}/BFEE/008_RMSDUnbound/008.0.1_removeProtein.tcl'],
                            cwd=f'{path}/BFEE/008_RMSDUnbound'
                        )
                        if neutralizeLigOnly is not None:
                            subprocess.run(
                                [vmdPath, '-dispdev', 'text', '-e', f'{path}/BFEE/008_RMSDUnbound/008.0.2_neutrilize.tcl'],
                                cwd=f'{path}/BFEE/008_RMSDUnbound'
                            )
                        
                elif forceFieldType == 'amber':
                    with open( f'{path}/BFEE/008_RMSDUnbound/008.0.1_removeProtein.cpptraj', 'w') as rScript:
                        rScript.write(
                            scriptTemplate.removeProteinAmberTemplate.substitute(
                                path='../complex', 
                                residueNum=fParser.getResid(selectionPro),
                                outputPath=f'./ligandOnly'
                            )
                        )
                # xyz and ndx for step 8
                fParserStep8 = fileParser.fileParser(f'{path}/BFEE/008_RMSDUnbound/ligandOnly.pdb')
                if not membraneProtein:
                    # otherwise the xyz file will be generated by vmd
                    if (vmdPath == '') or (forceFieldType == 'amber'):
                        fParserStep8.saveFile('all', f'{path}/BFEE/008_RMSDUnbound/ligandOnly.xyz', 'xyz')
                fParserStep8.saveNDX(
                    [selectionLig], ['ligand'], f'{path}/BFEE/008_RMSDUnbound/ligandOnly.ndx', True
                )
            
            # add water for step 7
            # this cannot be done in pure python
            # connect to VMD
            if userProvidedPullingTop == '' and userProvidedPullingCoor == '':
                if forceFieldType == 'charmm':
                    if not membraneProtein:
                        with pkg_resources.path(templates_namd, 'solvate.tcl') as p:
                            shutil.copyfile(p, f'{path}/BFEE/007_r/007.0_solvate.tcl')
                    else:
                        # membrane protein
                        with pkg_resources.path(templates_namd, 'solvate_mem.tcl') as p:
                            shutil.copyfile(p, f'{path}/BFEE/007_r/007.0_solvate.tcl')
                    # if vmd path is defined
                    # then execute vmd automatically
                    if vmdPath != '':
                        subprocess.run(
                            [vmdPath, '-dispdev', 'text', '-e', f'{path}/BFEE/007_r/007.0_solvate.tcl'],
                            cwd=f'{path}/BFEE/007_r'
                        )
            else:
                # copy the user-provided large box
                fParserStep7 = fileParser.fileParser(userProvidedPullingTop, userProvidedPullingCoor)

                # check polar angles
                # if theta < 30 or > 150, then rotate 90 degrees to avoid polar angle invarience
                if (fParserStep7.measurePolarAngles(selectionRef, selectionLig)[0] > 150 \
                    or fParserStep7.measurePolarAngles(selectionRef, selectionLig)[0] < 30):
                    fParserStep7.rotateSystem('x', 90)
                fParserStep7.centerSystem()

                fParserStep7.saveFile(
                    'all', f'{path}/BFEE/007_r/complex_largeBox.pdb', 'pdb',
                    True, f'{path}/BFEE/007_r/complex_largeBox.{topType}'
                )

                fParserStep7.saveFile(
                    'all', f'{path}/BFEE/007_r/complex_largeBox.xyz', 'xyz'
                )

        if jobType == 'alchemical':

            # readme file 
            with pkg_resources.path(templates_readme, 'Readme_NAMD_Alchemical.txt') as p:
                shutil.copyfile(p, f'{path}/BFEE/Readme.txt')
                
            # script to update the center after equilibration
            with pkg_resources.path(templates_namd, 'updateCenters.py') as p:
                shutil.copyfile(p, f'{path}/BFEE/000_eq/000.3_updateCenters.py')

            # fep file
            fParser.setBeta('all', 0)
            fParser.setBeta(selectionLig, 1)
            fParser.saveFile(
                'all', f'{path}/BFEE/fep.pdb', 'pdb'
            )
            with pkg_resources.path(templates_namd, 'fep.tcl') as p:
                shutil.copyfile(p, f'{path}/BFEE/fep.tcl')
                
            if not membraneProtein:
                # otherwise the these files will be generated by vmd
                fParser.saveFile(
                    f'not {selectionPro}', f'{path}/BFEE/ligandOnly.pdb', 'pdb'
                )
                fParser.saveFile(
                    f'not {selectionPro}', f'{path}/BFEE/fep_ligandOnly.pdb', 'pdb'
                )
            # remove protein for the unbound state
            if forceFieldType == 'charmm':
                if not membraneProtein:
                    with open(f'{path}/BFEE/002.3.1_removeProtein.tcl', 'w') as rScript:
                        rScript.write(
                            scriptTemplate.removeProteinTemplate.substitute(
                                path='./complex', selectionPro=f'{selectionPro}'.replace('segid', 'segname'),
                                outputPath=f'./ligandOnly'
                            )
                        )
                else:
                    # membrane protein
                    # fake ligandOnly.pdb, only used to generate ndx
                    # putting this here to guarantee ligandOnly.pdb is currected overwritten
                    fParser.saveFile(
                        f'{selectionLig}', f'{path}/BFEE/ligandOnly.pdb', 'pdb'
                    )
                    with open( f'{path}/BFEE/002.3.1_removeProtein.tcl', 'w') as rScript:
                        rScript.write(
                            scriptTemplate.removeMemProteinFepTemplate.substitute(
                                path='./complex', selectionLig=f'{selectionLig}'.replace('segid', 'segname'),
                                outputPath=f'./ligandOnly', outputFepPath=f'./fep_ligandOnly'
                            )
                        )
                        
                # neutralization
                if neutralizeLigOnly is not None:
                    with open(f'{path}/BFEE/002.3.2_neutrilize.tcl', 'w') as rScript:
                        rScript.write(
                            scriptTemplate.neutralizeSystempTemplate.substitute(
                                path='./ligandOnly', cationName=cation, anionName=anion,
                                extraCommand='$aa writepdb fep_ligandOnly.pdb'
                            )
                        )
                
                # if vmd path is defined
                # then execute vmd automatically
                if vmdPath != '':
                    subprocess.run(
                        [vmdPath, '-dispdev', 'text', '-e', f'{path}/BFEE/002.3.1_removeProtein.tcl'],
                        cwd=f'{path}/BFEE'
                    )
                    if neutralizeLigOnly is not None:
                        subprocess.run(
                            [vmdPath, '-dispdev', 'text', '-e', f'{path}/BFEE/002.3.2_neutrilize.tcl'],
                            cwd=f'{path}/BFEE'
                        )
            elif forceFieldType == 'amber':
                with open( f'{path}/BFEE/002.3.1_removeProtein.cpptraj', 'w') as rScript:
                    rScript.write(
                        scriptTemplate.removeProteinAmberTemplate.substitute(
                            path='./complex', 
                            residueNum=fParser.getResid(selectionPro),
                            outputPath=f'./ligandOnly'
                        )
                    )
            
            # xyz and ndx
            fParserLigandOnly = fileParser.fileParser( f'{path}/BFEE/ligandOnly.pdb')
            if not membraneProtein:
                # otherwise the xyz file will be generated by vmd
                if (vmdPath == '') or (forceFieldType == 'amber'):
                    fParserLigandOnly.saveFile('all', f'{path}/BFEE/ligandOnly.xyz', 'xyz')
            fParserLigandOnly.saveNDX(
                [selectionLig], ['ligand'], f'{path}/BFEE/ligandOnly.ndx', True
            )
            
    def _makeGromacsTopGro(
        self,
        path,
        forceFieldType,
        forceFieldFiles,
        jobType = 'geometric',
        vmdPath = '',
        considerRMSDCV = True
    ):
        """generate topology and coordinate files from CHARMM/Amber files for gromacs simulation

        Args:
            path (str): the root directory of simulation files
            forceFieldType (str): 'charmm' or 'amber'
            forceFieldFiles (list of str): list of CHARMM force field files
            jobType (str, optional): 'geometric' or 'alchemical'. Defaults to 'geometric'.
            vmdPath (str, optional): path to vmd, space is forbidden. Defaults to ''.
            considerRMSDCV (bool, optional): Whethre consider the RMSD CV. Default to True.
        """
        
        if forceFieldType == 'charmm':
            topType = 'psf'
        elif forceFieldType == 'amber':
            topType = 'parm7'
        coorType = 'pdb'

        # read the original topology and coordinate file
        fParser = fileParser.fileParser(f'{path}/BFEE/complex.{topType}', f'{path}/BFEE/complex.{coorType}')
        # pbc
        pbc = fParser.measurePBC()
        
        # gromacs assume the center of system is at (x/2, y/2, z/2)
        fParser.moveSystem(-pbc[1] + pbc[0] / 2)
        fParser.saveFile('all', f'{path}/BFEE/complex.{coorType}', coorType)
        
        if forceFieldType == 'charmm':
            fileParser.charmmToGromacs(
                f'{path}/BFEE/complex.psf', 
                f'{path}/BFEE/complex.pdb',
                forceFieldFiles,
                pbc[0],
                f'{path}/BFEE/complex_gmx'
            )
        elif forceFieldType == 'amber':
            fileParser.amberToGromacs(
                f'{path}/BFEE/complex.parm7', f'{path}/BFEE/complex.pdb', pbc[0], f'{path}/BFEE/complex_gmx'
            )
            
        
        if forceFieldType == 'charmm':
            if vmdPath != '':
                fParser_7 = fileParser.fileParser(
                    f'{path}/BFEE/007_r/complex_largeBox.psf', 
                    f'{path}/BFEE/007_r/complex_largeBox.pdb'
                )
                pbc_7 = fParser_7.measurePBC()
                fParser_7.moveSystem(-pbc_7[1] + pbc_7[0] / 2)
                fParser_7.saveFile('all', f'{path}/BFEE/007_r/complex_largeBox.pdb', 'pdb')

                if considerRMSDCV:
                    fParser_8 = fileParser.fileParser(
                        f'{path}/BFEE/008_RMSDUnbound/ligandOnly.psf', 
                        f'{path}/BFEE/008_RMSDUnbound/ligandOnly.pdb'
                    )
                    pbc_8 = fParser_8.measurePBC()
                    fParser_8.moveSystem(-pbc_8[1] + pbc_8[0] / 2)
                    fParser_8.saveFile('all', f'{path}/BFEE/008_RMSDUnbound/ligandOnly.pdb', 'pdb')
            
                fileParser.charmmToGromacs(
                    f'{path}/BFEE/007_r/complex_largeBox.psf', 
                    f'{path}/BFEE/007_r/complex_largeBox.pdb',
                    forceFieldFiles,
                    pbc_7[0],
                    f'{path}/BFEE/007_r/complex_largeBox_gmx'
                )
                
                if considerRMSDCV:
                    fileParser.charmmToGromacs(
                        f'{path}/BFEE/008_RMSDUnbound/ligandOnly.psf', 
                        f'{path}/BFEE/008_RMSDUnbound/ligandOnly.pdb',
                        forceFieldFiles,
                        pbc_8[0],
                        f'{path}/BFEE/008_RMSDUnbound/ligandOnly_gmx'
                    )
            else:
                with open( f'{path}/BFEE/007_r/007.0.2_genGromacsTop.py', 'w') as rScript:
                    rScript.write(
                        scriptTemplate.charmmToGromacsTemplate(
                            inputPrefix='./complex_largeBox',
                            forceFieldList=forceFieldFiles,
                        )
                    )
                
                if considerRMSDCV:
                    with open( f'{path}/BFEE/008_RMSDUnbound/008.0.3_genGromacsTop.py', 'w') as rScript:
                        rScript.write(
                            scriptTemplate.charmmToGromacsTemplate(
                                inputPrefix='./ligandOnly',
                                forceFieldList=forceFieldFiles,
                            )
                        )
                
        if forceFieldType == 'amber':
            
            fParser_7 = fileParser.fileParser(
                f'{path}/BFEE/007_r/complex_largeBox.parm7', 
                f'{path}/BFEE/007_r/complex_largeBox.pdb'
            )
            pbc_7 = fParser_7.measurePBC()
            fParser_7.moveSystem(-pbc_7[1] + pbc_7[0] / 2)
            fParser_7.saveFile('all', f'{path}/BFEE/007_r/complex_largeBox.pdb', 'pdb')
            
            fileParser.amberToGromacs(
                f'{path}/BFEE/007_r/complex_largeBox.parm7', 
                f'{path}/BFEE/007_r/complex_largeBox.pdb',
                pbc_7[0],
                f'{path}/BFEE/007_r/complex_largeBox_gmx'
            )
            
            if considerRMSDCV:
                with open( f'{path}/BFEE/008_RMSDUnbound/008.0.2_genGromacsTop.py', 'w') as rScript:
                    rScript.write(
                        scriptTemplate.amberToGromacsTemplate.substitute(
                            inputPrefix='./ligandOnly'
                        )
                    )   
            

    def _generateAlchemicalNAMDConfig(
        self,
        path,
        forceFieldType,
        forceFields,
        temperature,
        stratification = [1,1,1,1],
        doubleWide = False,
        minBeforeSample = False,
        membraneProtein = False,
        OPLSMixingRule = False
    ):
        """generate NAMD config fils for the alchemical route

        Args:
            path (str): the directory for generation of all the files
            forceFieldType (str): 'charmm' or 'amber'
            forceFields (list of str): list of CHARMM force field files
            temperature (float): temperature of the simulation
            stratification (list of int, 4, optional): number of windows for each simulation. 
                                                       Defaults to [1,1,1,1].
            doubleWide (bool, optional): whether double-wide simulations are carried out. 
                                         Defaults to False.
            minBeforeSample (bool, optional): minimization before sampling in each FEP window. 
                                              Defaults to False.
            membraneProtein (bool, optional): whether simulating a membrane protein. 
                                              Defaults to False.
            OPLSMixingRule (bool, optional): whether use the OPLS mixing rules. Defaults to False.
        """

        if forceFieldType == 'charmm':
            topType = 'psf'
        elif forceFieldType == 'amber':
            topType = 'parm7'
        coorType = 'pdb'

        # read the original topology and coordinate file
        fParser = fileParser.fileParser(f'{path}/BFEE/complex.{topType}', f'{path}/BFEE/complex.{coorType}')
        # pbc
        pbc = fParser.measurePBC()

        # read the ligandOnly topology and coordinate file
        if os.path.exists(f'{path}/BFEE/ligandOnly.{topType}'):
            fParserLig = fileParser.fileParser(f'{path}/BFEE/ligandOnly.{topType}', f'{path}/BFEE/ligandOnly.{coorType}')
            # pbc
            pbcLig = fParserLig.measurePBC()
        else:
            pbcLig = pbc


        # 000_eq
        with open(f'{path}/BFEE/000_eq/000.1_eq.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    '', '', '', pbc,
                    'output/eq', temperature, 5000000, 'colvars.in', '', membraneProtein=membraneProtein,
                    OPLSMixingRule=OPLSMixingRule
                )
            )
        with open(f'{path}/BFEE/000_eq/000.2_eq_ligandOnly.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../ligandOnly.{topType}', f'../ligandOnly.pdb',
                    '', '', '', pbcLig,
                    'output/eq_ligandOnly', temperature, 1000000, 'colvars_ligandOnly.in',
                    '', OPLSMixingRule=OPLSMixingRule
                )
            )

        # 001_MoleculeBound
        with open(f'{path}/BFEE/001_MoleculeBound/001.2_fep_forward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'output/fep_backward.coor', f'output/fep_backward.vel', f'output/fep_backward.xsc', '',
                    'output/fep_forward', temperature, 0, 'colvars.in', '', '', '../fep.pdb', 
                    stratification[0], True, False, minBeforeSample, membraneProtein=membraneProtein,
                    OPLSMixingRule=OPLSMixingRule
                )
            )
        with open(f'{path}/BFEE/001_MoleculeBound/001.1_fep_backward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc', '',
                    'output/fep_backward', temperature, 0, 'colvars.in', '', '', '../fep.pdb', 
                    stratification[0], False, False, minBeforeSample, membraneProtein=membraneProtein,
                    OPLSMixingRule=OPLSMixingRule
                )
            )
        
        if doubleWide:
            with open(f'{path}/BFEE/001_MoleculeBound/001_fep_doubleWide.conf', 'w') as namdConfig:
                namdConfig.write(
                    self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc', '',
                        'output/fep_doubleWide', temperature, 0, 'colvars.in', '', '', '../fep.pdb', 
                        stratification[0], False, True, membraneProtein=membraneProtein,
                        OPLSMixingRule=OPLSMixingRule
                    )
                )

        # 002_RestraintBound
        with open(f'{path}/BFEE/002_RestraintBound/002.2_ti_forward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'output/ti_backward.coor', f'output/ti_backward.vel', f'output/ti_backward.xsc', '',
                    'output/ti_forward', temperature, f'{500000*(stratification[1]+1)}', 'colvars_forward.in', 
                    '', membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule
                )
            )
        with open(f'{path}/BFEE/002_RestraintBound/002.1_ti_backward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc', '',
                    'output/ti_backward', temperature, f'{500000*(stratification[1]+1)}', 'colvars_backward.in', 
                    '', membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule
                )
            )

        # 003_MoleculeUnbound
        with open(f'{path}/BFEE/003_MoleculeUnbound/003.2_fep_forward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../ligandOnly.{topType}', f'../ligandOnly.pdb',
                    f'output/fep_backward.coor', f'output/fep_backward.vel', 
                    f'output/fep_backward.xsc', '',
                    'output/fep_forward', temperature, 0, 'colvars.in', '', '', '../fep_ligandOnly.pdb', 
                    stratification[2], True, False, minBeforeSample, OPLSMixingRule=OPLSMixingRule
                )
            )
        with open(f'{path}/BFEE/003_MoleculeUnbound/003.1_fep_backward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../ligandOnly.{topType}', f'../ligandOnly.pdb',
                    f'../000_eq/output/eq_ligandOnly.coor', f'../000_eq/output/eq_ligandOnly.vel', 
                    f'../000_eq/output/eq_ligandOnly.xsc', '',
                    'output/fep_backward', temperature, 0, 'colvars.in', '', '', '../fep_ligandOnly.pdb', 
                    stratification[2], False, False, minBeforeSample, OPLSMixingRule=OPLSMixingRule
                )
            )

        if doubleWide:
            with open(f'{path}/BFEE/003_MoleculeUnbound/003_fep_doubleWide.conf', 'w') as namdConfig:
                namdConfig.write(
                    self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../ligandOnly.{topType}', f'../ligandOnly.pdb',
                        f'../000_eq/output/eq_ligandOnly.coor', f'../000_eq/output/eq_ligandOnly.vel', 
                        f'../000_eq/output/eq_ligandOnly.xsc', '',
                        'output/fep_doubleWide', temperature, 0, 'colvars.in', '', '', '../fep_ligandOnly.pdb', 
                        stratification[2], False, True, OPLSMixingRule=OPLSMixingRule
                    )
                )

        # 004_RestraintUnbound
        with open(f'{path}/BFEE/004_RestraintUnbound/004.2_ti_forward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../ligandOnly.{topType}', f'../ligandOnly.pdb',
                    f'output/ti_backward.coor', f'output/ti_backward.vel', 
                    f'output/ti_backward.xsc', '',
                    'output/ti_forward', temperature, f'{500000*(stratification[3]+1)}', 'colvars_forward.in', 
                    '', OPLSMixingRule=OPLSMixingRule
                )
            )
        with open(f'{path}/BFEE/004_RestraintUnbound/004.1_ti_backward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../ligandOnly.{topType}', f'../ligandOnly.pdb',
                    f'../000_eq/output/eq_ligandOnly.coor', f'../000_eq/output/eq_ligandOnly.vel', 
                    f'../000_eq/output/eq_ligandOnly.xsc', '',
                    'output/ti_backward', temperature, f'{500000*(stratification[3]+1)}', 'colvars_backward.in', 
                    '', OPLSMixingRule=OPLSMixingRule
                )
            )

    def _generateAlchemicalColvarsConfig(
        self, path, topType, coorType, selectionPro, selectionLig, selectionRef, 
        stratification=[1,1,1,1], pinDownPro=True, useOldCv=True
    ):
        """generate Colvars config fils for geometric route

        Args:
            path (str): the directory for generation of all the files
            topType (str): the type (psf, parm) of the topology file
            coorType (str): the type (pdb, rst) of the topology file
            selectionPro (str): MDAnalysis-style selection of the protein
            selectionLig (str): MDAnalysis-style selection of the ligand
            selectionRef (str): MDAnalysis-style selection of the reference group for pulling simulation
            stratification (list of int, 4, optional): number of windows for each simulation. 
                                                       Defaults to [1,1,1,1].
            pinDownPro (bool, optinal): Whether pinning down the protein. Defaults to True.
            useOldCv (bool, optional): whether used old, custom-function-based cv. Defaults to True.
        """

        assert(len(stratification) == 4)

        # read the original topology and coordinate file
        fParser = fileParser.fileParser(f'{path}/BFEE/complex.{topType}', f'{path}/BFEE/complex.{coorType}')
        center = fParser.measureCenter(selectionPro)
        polarAngles = fParser.measurePolarAngles(selectionRef, selectionLig)
        distance = fParser.measureDistance(selectionRef, selectionLig)

        # 000_eq
        with open(f'{path}/BFEE/000_eq/colvars.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../complex.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(
                    True, 0, 10, '../complex.xyz',
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    True, -90, 90, 'eulerTheta', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    True, -180, 180, 'eulerPhi', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    True, -180, 180, 'eulerPsi', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    True, 0, 180, 'polarTheta', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    True, -180, 180, 'polarPhi', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvRTemplate(
                    True, (distance - 5 if distance - 5 > 0 else 0), 
                    distance + 22, extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvHistogramTemplate('RMSD')
            )
            colvarsConfig.write(
                self.cTemplate.cvHistogramTemplate('eulerTheta')
            )
            colvarsConfig.write(
                self.cTemplate.cvHistogramTemplate('eulerPhi')
            )
            colvarsConfig.write(
                self.cTemplate.cvHistogramTemplate('eulerPsi')
            )
            colvarsConfig.write(
                self.cTemplate.cvHistogramTemplate('polarTheta')
            )
            colvarsConfig.write(
                self.cTemplate.cvHistogramTemplate('polarPhi')
            )
            colvarsConfig.write(
                self.cTemplate.cvHistogramTemplate('r')
            )
            colvarsConfig.write(
                self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
            )

        with open(f'{path}/BFEE/000_eq/colvars_ligandOnly.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../ligandOnly.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(
                    False, '', '', '../ligandOnly.xyz',
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('RMSD', 10, 0)
            )

        # 001_MoleculeBound
        with open(f'{path}/BFEE/001_MoleculeBound/colvars.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../complex.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(
                    False, '', '', '../complex.xyz',
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'eulerTheta', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'eulerPhi', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'eulerPsi', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'polarTheta', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'polarPhi', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvRTemplate(
                    False, 0, 0, extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('RMSD', 10, 0)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('eulerTheta', 0.1, 0)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('eulerPhi', 0.1, 0)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('eulerPsi', 0.1, 0)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('polarTheta', 0.1, polarAngles[0])
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('polarPhi', 0.1, polarAngles[1])
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('r', 10, distance)
            )
            if pinDownPro:
                colvarsConfig.write(
                    self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
                )

        # 002_RestraintBound
        with open(f'{path}/BFEE/002_RestraintBound/colvars_forward.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../complex.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(
                    False, '', '', '../complex.xyz',
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'eulerTheta', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'eulerPhi', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'eulerPsi', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'polarTheta', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'polarPhi', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvRTemplate(
                    False, 0, 0, extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('RMSD', 0, 0, stratification[1], True, 10)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('eulerTheta', 0, 0, stratification[1], True, 0.1)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('eulerPhi', 0, 0, stratification[1], True, 0.1)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('eulerPsi', 0, 0, stratification[1], True, 0.1)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('polarTheta', 0, polarAngles[0], stratification[1], True, 0.1)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('polarPhi', 0, polarAngles[1], stratification[1], True, 0.1)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('r', 0, distance, stratification[1], True, 10)
            )
            if pinDownPro:
                colvarsConfig.write(
                    self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
                )
        with open(f'{path}/BFEE/002_RestraintBound/colvars_backward.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../complex.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(
                    False, '', '', '../complex.xyz', extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'eulerTheta', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'eulerPhi', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'eulerPsi', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'polarTheta', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvAngleTemplate(
                    False, 0, 0, 'polarPhi', '../complex.xyz', useOldCv,
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvRTemplate(
                    False, 0, 0, extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('RMSD', 0, 0, stratification[1], False, 10)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('eulerTheta', 0, 0, stratification[1], False, 0.1)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('eulerPhi', 0, 0, stratification[1], False, 0.1)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('eulerPsi', 0, 0, stratification[1], False, 0.1)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('polarTheta', 0, polarAngles[0], stratification[1], False, 0.1)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('polarPhi', 0, polarAngles[1], stratification[1], False, 0.1)
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('r', 0, distance, stratification[1], False, 10)
            )
            if pinDownPro:
                colvarsConfig.write(
                    self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
                )

        # 003_MoleculeUnbound
        with open(f'{path}/BFEE/003_MoleculeUnbound/colvars.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../ligandOnly.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(
                    False, '', '', '../ligandOnly.xyz',
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('RMSD', 10, 0)
            )

        # 004_RestraintUnbound
        with open(f'{path}/BFEE/004_RestraintUnbound/colvars_forward.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../ligandOnly.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(
                    False, '', '', '../ligandOnly.xyz',
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('RMSD', 0, 0, stratification[3], True, 10)
            )
        with open(f'{path}/BFEE/004_RestraintUnbound/colvars_backward.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../ligandOnly.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(
                    False, '', '', '../ligandOnly.xyz',
                    extendedLagrangian = False
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('RMSD', 0, 0, stratification[3], False, 10)
            )

    def _generateGeometricNAMDConfig(
        self,
        path,
        forceFieldType,
        forceFields,
        temperature,
        stratification = [1,1,1,1,1,1,1,1],
        membraneProtein = False,
        OPLSMixingRule = False,
        considerRMSDCV = True,
        GaWTM = False
    ):
        """generate NAMD config fils for the geometric route

        Args:
            path (str): the directory for generation of all the files
            forceFieldType (str): 'charmm' or 'amber'
            forceFieldFiles (list of str): list of CHARMM force field files
            temperature (float): temperature of the simulation
            stratification (list of int, 8): number of windows for each simulation. 
                                             Defaults to [1,1,1,1,1,1,1,1].
            membraneProtein (bool, optional): whether simulating a membrane protein. Defaults to False.
            OPLSMixingRule (bool, optional): whether use the OPLS mixing rules. Defaults to False.
            considerRMSDCV (bool, optional): Whethre consider the RMSD CV. Default to True.
            GaWTM (bool, optional): Whether this is an GaWTM-eABF simulation. Default to False
        """

        if forceFieldType == 'charmm':
            topType = 'psf'
        elif forceFieldType == 'amber':
            topType = 'parm7'
        coorType = 'pdb'

        # read the original topology and coordinate file
        fParser = fileParser.fileParser(f'{path}/BFEE/complex.{topType}', f'{path}/BFEE/complex.{coorType}')

        # pbc
        pbc = fParser.measurePBC()
        
        # the extended water box
        # read the ligandOnly topology and coordinate file
        if os.path.exists(f'{path}/BFEE/007_r/complex_largeBox.{topType}'):
            fParserLargeBox = fileParser.fileParser(
                f'{path}/BFEE/007_r/complex_largeBox.{topType}', 
                f'{path}/BFEE/007_r/complex_largeBox.{coorType}'
            )
            # pbc
            pbcStep7 = fParserLargeBox.measurePBC()
        else:
            if not membraneProtein:
                pbcStep7 = pbc + np.array([[24,24,24],[0,0,0]])
            else:
                pbcStep7 = pbc + np.array([[0,0,32],[0,0,0]])

        # read the ligandOnly topology and coordinate file
        if os.path.exists(f'{path}/BFEE/008_RMSDUnbound/ligandOnly.{topType}'):
            fParserLig = fileParser.fileParser(
                f'{path}/BFEE/008_RMSDUnbound/ligandOnly.{topType}', 
                f'{path}/BFEE/008_RMSDUnbound/ligandOnly.{coorType}'
            )
            # pbc
            pbcLig = fParserLig.measurePBC()
        else:
            pbcLig = pbc
        

        # 000_eq
        with open(f'{path}/BFEE/000_eq/000.1_eq.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    '', '', '', pbc,
                    'output/eq', temperature, 5000000, 'colvars.in',
                    membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                    GaWTM=False
                )
            )

        # RMSD bound
        if considerRMSDCV:
            with open(f'{path}/BFEE/001_RMSDBound/001_abf_1.conf', 'w') as namdConfig:
                namdConfig.write(
                    self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc',
                        '', 'output/abf_1', temperature, 5000000, 'colvars_1.in',
                        membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                        GaWTM=GaWTM
                    )
                )
            with open(f'{path}/BFEE/001_RMSDBound/001_abf_1.extend.conf', 'w') as namdConfig:
                namdConfig.write(
                    self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_1.restart.coor', f'output/abf_1.restart.vel', f'output/abf_1.restart.xsc',
                        '', 'output/abf_1.extend', temperature, 5000000, 'colvars_1.in', 
                        CVRestartFile='output/abf_1.restart', membraneProtein=membraneProtein,
                        OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                    )
                )

            # stratification
            if stratification[0] > 1:
                for i in range(1, stratification[0]):
                    with open(f'{path}/BFEE/001_RMSDBound/001_abf_{i+1}.conf', 'w') as namdConfig:
                        namdConfig.write(
                            self.cTemplate.namdConfigTemplate(
                            forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                            f'output/abf_{i}.restart.coor', f'output/abf_{i}.restart.vel', 
                            f'output/abf_{i}.restart.xsc',
                            '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in',
                            membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                            GaWTM=GaWTM
                        )
                    )
                    with open(f'{path}/BFEE/001_RMSDBound/001_abf_{i+1}.extend.conf', 'w') as namdConfig:
                        namdConfig.write(
                            self.cTemplate.namdConfigTemplate(
                            forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                            f'output/abf_{i+1}.restart.coor', f'output/abf_{i+1}.restart.vel', 
                            f'output/abf_{i+1}.restart.xsc',
                            '', f'output/abf_{i+1}.extend', temperature, 5000000, f'colvars_{i+1}.in', 
                            CVRestartFile=f'output/abf_{i+1}.restart', membraneProtein=membraneProtein, 
                            OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                        )
                    )

        # Theta
        with open(f'{path}/BFEE/002_EulerTheta/002_abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc',
                    '', 'output/abf_1', temperature, 5000000, 'colvars_1.in', '',
                    membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                    GaWTM=GaWTM
                )
            )
        with open(f'{path}/BFEE/002_EulerTheta/002_abf_1.extend.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'output/abf_1.restart.coor', f'output/abf_1.restart.vel', f'output/abf_1.restart.xsc',
                    '', 'output/abf_1.extend', temperature, 5000000, 'colvars_1.in', '',
                    CVRestartFile='output/abf_1.restart', membraneProtein=membraneProtein,
                    OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                )
            )

        # stratification
        if stratification[1] > 1:
            for i in range(1, stratification[1]):
                with open(f'{path}/BFEE/002_EulerTheta/002_abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i}.restart.coor', f'output/abf_{i}.restart.vel', 
                        f'output/abf_{i}.restart.xsc',
                        '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in', '',
                        membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                        GaWTM=GaWTM
                    )
                )
                with open(f'{path}/BFEE/002_EulerTheta/002_abf_{i+1}.extend.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i+1}.restart.coor', f'output/abf_{i+1}.restart.vel', 
                        f'output/abf_{i+1}.restart.xsc',
                        '', f'output/abf_{i+1}.extend', temperature, 5000000, f'colvars_{i+1}.in', '',
                        CVRestartFile=f'output/abf_{i+1}.restart', membraneProtein=membraneProtein,
                        OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                    )
                )

        # Phi
        with open(f'{path}/BFEE/003_EulerPhi/003_abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc',
                    '', 'output/abf_1', temperature, 5000000, 'colvars_1.in', '',
                    membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                    GaWTM=GaWTM
                )
            )
        with open(f'{path}/BFEE/003_EulerPhi/003_abf_1.extend.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'output/abf_1.restart.coor', f'output/abf_1.restart.vel', f'output/abf_1.restart.xsc',
                    '', 'output/abf_1.extend', temperature, 5000000, 'colvars_1.in', '',
                    CVRestartFile='output/abf_1.restart', membraneProtein=membraneProtein,
                    OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                )
            )

        # stratification
        if stratification[2] > 1:
            for i in range(1, stratification[2]):
                with open(f'{path}/BFEE/003_EulerPhi/003_abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i}.restart.coor', f'output/abf_{i}.restart.vel', 
                        f'output/abf_{i}.restart.xsc',
                        '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in', '',
                        membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                        GaWTM=GaWTM
                    )
                )
                with open(f'{path}/BFEE/003_EulerPhi/003_abf_{i+1}.extend.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i+1}.restart.coor', f'output/abf_{i+1}.restart.vel', 
                        f'output/abf_{i+1}.restart.xsc',
                        '', f'output/abf_{i+1}.extend', temperature, 5000000, f'colvars_{i+1}.in', '',
                        CVRestartFile=f'output/abf_{i+1}.restart', membraneProtein=membraneProtein,
                        OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                    )
                )

        # Psi
        with open(f'{path}/BFEE/004_EulerPsi/004_abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc',
                    '', 'output/abf_1', temperature, 5000000, 'colvars_1.in', '', 
                    membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                    GaWTM=GaWTM
                )
            )
        with open(f'{path}/BFEE/004_EulerPsi/004_abf_1.extend.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'output/abf_1.restart.coor', f'output/abf_1.restart.vel', f'output/abf_1.restart.xsc',
                    '', 'output/abf_1.extend', temperature, 5000000, 'colvars_1.in', '',
                    CVRestartFile='output/abf_1.restart', membraneProtein=membraneProtein,
                    OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                )
            )

        # stratification
        if stratification[3] > 1:
            for i in range(1, stratification[3]):
                with open(f'{path}/BFEE/004_EulerPsi/004_abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i}.restart.coor', f'output/abf_{i}.restart.vel', 
                        f'output/abf_{i}.restart.xsc',
                        '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in', '',
                        membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                        GaWTM=GaWTM
                    )
                )
                with open(f'{path}/BFEE/004_EulerPsi/004_abf_{i+1}.extend.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i+1}.restart.coor', f'output/abf_{i+1}.restart.vel', 
                        f'output/abf_{i+1}.restart.xsc',
                        '', f'output/abf_{i+1}.extend', temperature, 5000000, f'colvars_{i+1}.in', '',
                        CVRestartFile=f'output/abf_{i+1}.restart', membraneProtein=membraneProtein,
                        OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                    )
                )

        # theta
        with open(f'{path}/BFEE/005_PolarTheta/005_abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc',
                    '', 'output/abf_1', temperature, 5000000, 'colvars_1.in', '',
                    membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                    GaWTM=GaWTM
                )
            )
        with open(f'{path}/BFEE/005_PolarTheta/005_abf_1.extend.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'output/abf_1.restart.coor', f'output/abf_1.restart.vel', f'output/abf_1.restart.xsc',
                    '', 'output/abf_1.extend', temperature, 5000000, 'colvars_1.in', '',
                    CVRestartFile='output/abf_1.restart', membraneProtein=membraneProtein,
                    OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                )
            )


        # stratification
        if stratification[4] > 1:
            for i in range(1, stratification[4]):
                with open(f'{path}/BFEE/005_PolarTheta/005_abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i}.restart.coor', f'output/abf_{i}.restart.vel', 
                        f'output/abf_{i}.restart.xsc',
                        '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in', '',
                        membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                        GaWTM=GaWTM
                    )
                )
                with open(f'{path}/BFEE/005_PolarTheta/005_abf_{i+1}.extend.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i+1}.restart.coor', f'output/abf_{i+1}.restart.vel', 
                        f'output/abf_{i+1}.restart.xsc',
                        '', f'output/abf_{i+1}.extend', temperature, 5000000, f'colvars_{i+1}.in', '',
                        CVRestartFile=f'output/abf_{i+1}.restart', membraneProtein=membraneProtein,
                        OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                    )
                )

        # phi
        with open(f'{path}/BFEE/006_PolarPhi/006_abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc',
                    '', 'output/abf_1', temperature, 5000000, 'colvars_1.in', '',
                    membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                    GaWTM=GaWTM
                )
            )
        with open(f'{path}/BFEE/006_PolarPhi/006_abf_1.extend.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'output/abf_1.restart.coor', f'output/abf_1.restart.vel', f'output/abf_1.restart.xsc',
                    '', 'output/abf_1.extend', temperature, 5000000, 'colvars_1.in', '',
                    CVRestartFile='output/abf_1.restart', membraneProtein=membraneProtein,
                    OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                )
            )

        # stratification
        if stratification[5] > 1:
            for i in range(1, stratification[5]):
                with open(f'{path}/BFEE/006_PolarPhi/006_abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i}.restart.coor', f'output/abf_{i}.restart.vel', 
                        f'output/abf_{i}.restart.xsc',
                        '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in', '',
                        membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                        GaWTM=GaWTM
                    )
                )
                with open(f'{path}/BFEE/006_PolarPhi/006_abf_{i+1}.extend.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i+1}.restart.coor', f'output/abf_{i+1}.restart.vel', 
                        f'output/abf_{i+1}.restart.xsc',
                        '', f'output/abf_{i+1}.extend', temperature, 5000000, f'colvars_{i+1}.in', '',
                        CVRestartFile=f'output/abf_{i+1}.restart', membraneProtein=membraneProtein,
                        OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                    )
                )

        # r
        # eq
        with open(f'{path}/BFEE/007_r/007.1_eq.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'./complex_largeBox.{topType}', f'./complex_largeBox.pdb',
                    '', '', '', 
                    pbcStep7,
                    'output/eq', temperature, 5000000, 'colvars_eq.in', '',
                    membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                    GaWTM=False
                )
            )
        # abf
        with open(f'{path}/BFEE/007_r/007.2_abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'./complex_largeBox.{topType}', f'./complex_largeBox.pdb',
                    'output/eq.coor', 'output/eq.vel', 'output/eq.xsc', '',
                    'output/abf_1', temperature, 20000000, 'colvars_1.in', '',
                    membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                    GaWTM=GaWTM
                )
            )
        with open(f'{path}/BFEE/007_r/007.2_abf_1.extend.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'./complex_largeBox.{topType}', f'./complex_largeBox.pdb',
                    'output/abf_1.restart.coor', 'output/abf_1.restart.vel', 'output/abf_1.restart.xsc', '',
                    'output/abf_1.extend', temperature, 20000000, 'colvars_1.in', '',
                    CVRestartFile=f'output/abf_1.restart', membraneProtein=membraneProtein,
                    OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                )
            )

        # stratification
        if stratification[6] > 1:
            for i in range(1, stratification[6]):
                with open(f'{path}/BFEE/007_r/007.2_abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'./complex_largeBox.{topType}', f'./complex_largeBox.pdb',
                        f'output/abf_{i}.restart.coor', f'output/abf_{i}.restart.vel', 
                        f'output/abf_{i}.restart.xsc',
                        '', f'output/abf_{i+1}', temperature, 20000000, f'colvars_{i+1}.in', '',
                        membraneProtein=membraneProtein, OPLSMixingRule=OPLSMixingRule,
                        GaWTM=GaWTM
                    )
                )
                with open(f'{path}/BFEE/007_r/007.2_abf_{i+1}.extend.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'./complex_largeBox.{topType}', f'./complex_largeBox.pdb',
                        f'output/abf_{i+1}.restart.coor', f'output/abf_{i+1}.restart.vel', 
                        f'output/abf_{i+1}.restart.xsc',
                        '', f'output/abf_{i+1}.extend', temperature, 20000000, f'colvars_{i+1}.in', '',
                        CVRestartFile=f'output/abf_{i+1}.restart', membraneProtein=membraneProtein,
                        OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                    )
                )

        # RMSD unbound
        if considerRMSDCV:
            # eq
            with open(f'{path}/BFEE/008_RMSDUnbound/008.1_eq.conf', 'w') as namdConfig:
                namdConfig.write(
                    self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'./ligandOnly.{topType}', f'./ligandOnly.pdb',
                        '', '', '', 
                        pbcLig,
                        'output/eq', temperature, 1000000, OPLSMixingRule=OPLSMixingRule,
                        GaWTM=False
                    )
                )
            # abf
            with open(f'{path}/BFEE/008_RMSDUnbound/008.2_abf_1.conf', 'w') as namdConfig:
                namdConfig.write(
                    self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'./ligandOnly.{topType}', f'./ligandOnly.pdb',
                        'output/eq.coor', 'output/eq.vel', 'output/eq.xsc', '',
                        'output/abf_1', temperature, 5000000, 'colvars_1.in',
                        OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                    )
                )
            with open(f'{path}/BFEE/008_RMSDUnbound/008.2_abf_1.extend.conf', 'w') as namdConfig:
                namdConfig.write(
                    self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'./ligandOnly.{topType}', f'./ligandOnly.pdb',
                        'output/abf_1.restart.coor', 'output/abf_1.restart.vel', 'output/abf_1.restart.xsc', '',
                        'output/abf_1.extend', temperature, 5000000, 'colvars_1.in',
                        CVRestartFile=f'output/abf_1.restart', OPLSMixingRule=OPLSMixingRule,
                        GaWTM=GaWTM
                    )
                )

            # stratification
            if stratification[7] > 1:
                for i in range(1, stratification[7]):
                    with open(f'{path}/BFEE/008_RMSDUnbound/008.2_abf_{i+1}.conf', 'w') as namdConfig:
                        namdConfig.write(
                            self.cTemplate.namdConfigTemplate(
                            forceFieldType, forceFields, f'./ligandOnly.{topType}', f'./ligandOnly.pdb',
                            f'output/abf_{i}.restart.coor', f'output/abf_{i}.restart.vel', 
                            f'output/abf_{i}.restart.xsc',
                            '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in',
                            OPLSMixingRule=OPLSMixingRule, GaWTM=GaWTM
                        )
                    )
                    with open(f'{path}/BFEE/008_RMSDUnbound/008.2_abf_{i+1}.extend.conf', 'w') as namdConfig:
                        namdConfig.write(
                            self.cTemplate.namdConfigTemplate(
                            forceFieldType, forceFields, f'./ligandOnly.{topType}', f'./ligandOnly.pdb',
                            f'output/abf_{i+1}.restart.coor', f'output/abf_{i+1}.restart.vel', 
                            f'output/abf_{i+1}.restart.xsc',
                            '', f'output/abf_{i+1}.extend', temperature, 5000000, f'colvars_{i+1}.in',
                            CVRestartFile=f'output/abf_{i+1}.restart', OPLSMixingRule=OPLSMixingRule,
                            GaWTM=GaWTM
                        )
                    )
                    
    def _generateGeometricGromacsConfig(
        self,
        path,
        forceFieldType,
        temperature,
        membraneProtein = False,
        considerRMSDCV = True
    ):
        """generate NAMD config fils for the geometric route

        Args:
            path (str): the directory for generation of all the files
            forceFieldType (str): 'charmm' or 'amber'
            temperature (float): temperature of the simulation
            membraneProtein (bool, optional): whether simulating a membrane protein. Defaults to False.
            considerRMSDCV (bool, optional): Whethre consider the RMSD CV. Default to True.
        """
        
        # 000_eq
        with open(f'{path}/BFEE/000_eq/000.1_min.mdp', 'w') as gromacsConfig:
            gromacsConfig.write(
                self.cTemplate.gromacsMinimizeConfigTemplate()
            )
        with open(f'{path}/BFEE/000_eq/000.2_eq.mdp', 'w') as gromacsConfig:
            gromacsConfig.write(
                self.cTemplate.gromacsConfigTemplate(
                    forceFieldType, temperature, 5000000,
                    membraneProtein=membraneProtein,
                    generateVelocities=True
                )
            )

        # RMSD bound
        if considerRMSDCV:
            with open(f'{path}/BFEE/001_RMSDBound/001_abf.mdp', 'w') as gromacsConfig:
                gromacsConfig.write(
                    self.cTemplate.gromacsConfigTemplate(
                        forceFieldType, temperature, 5000000,
                        membraneProtein=membraneProtein,
                        generateVelocities=False
                    )
                )

        # Theta
        with open(f'{path}/BFEE/002_EulerTheta/002_abf.mdp', 'w') as gromacsConfig:
            gromacsConfig.write(
                self.cTemplate.gromacsConfigTemplate(
                    forceFieldType, temperature, 5000000,
                    membraneProtein=membraneProtein,
                    generateVelocities=False
                )
            )

        # Phi
        with open(f'{path}/BFEE/003_EulerPhi/003_abf.mdp', 'w') as gromacsConfig:
            gromacsConfig.write(
                self.cTemplate.gromacsConfigTemplate(
                    forceFieldType, temperature, 5000000,
                    membraneProtein=membraneProtein,
                    generateVelocities=False
                )
            )

        # Psi
        with open(f'{path}/BFEE/004_EulerPsi/004_abf.mdp', 'w') as gromacsConfig:
            gromacsConfig.write(
                self.cTemplate.gromacsConfigTemplate(
                    forceFieldType, temperature, 5000000,
                    membraneProtein=membraneProtein,
                    generateVelocities=False
                )
            )

        # theta
        with open(f'{path}/BFEE/005_PolarTheta/005_abf.mdp', 'w') as gromacsConfig:
            gromacsConfig.write(
                self.cTemplate.gromacsConfigTemplate(
                    forceFieldType, temperature, 5000000,
                    membraneProtein=membraneProtein,
                    generateVelocities=False
                )
            )

        # phi
        with open(f'{path}/BFEE/006_PolarPhi/006_abf.mdp', 'w') as gromacsConfig:
            gromacsConfig.write(
                self.cTemplate.gromacsConfigTemplate(
                    forceFieldType, temperature, 5000000,
                    membraneProtein=membraneProtein,
                    generateVelocities=False
                )
            )

        # r
        # eq
        with open(f'{path}/BFEE/007_r/007.1_min.mdp', 'w') as gromacsConfig:
            gromacsConfig.write(
                self.cTemplate.gromacsMinimizeConfigTemplate()
            )
        with open(f'{path}/BFEE/007_r/007.2_eq.mdp', 'w') as gromacsConfig:
            gromacsConfig.write(
                self.cTemplate.gromacsConfigTemplate(
                    forceFieldType, temperature, 5000000,
                    membraneProtein=membraneProtein,
                    generateVelocities=True
                )
            )
        # abf
        with open(f'{path}/BFEE/007_r/007.3_abf.mdp', 'w') as gromacsConfig:
            gromacsConfig.write(
                self.cTemplate.gromacsConfigTemplate(
                    forceFieldType, temperature, 20000000,
                    membraneProtein=membraneProtein,
                    generateVelocities=False
                )
            )

        # RMSD unbound
        if considerRMSDCV:
            # eq
            with open(f'{path}/BFEE/008_RMSDUnbound/008.1_min.mdp', 'w') as gromacsConfig:
                gromacsConfig.write(
                    self.cTemplate.gromacsMinimizeConfigTemplate()
                )
            with open(f'{path}/BFEE/008_RMSDUnbound/008.2_eq.mdp', 'w') as gromacsConfig:
                gromacsConfig.write(
                    self.cTemplate.gromacsConfigTemplate(
                        forceFieldType, temperature, 1000000,
                        membraneProtein=membraneProtein,
                        generateVelocities=True
                    )
                )
            # abf
            with open(f'{path}/BFEE/008_RMSDUnbound/008.3_abf.mdp', 'w') as gromacsConfig:
                gromacsConfig.write(
                    self.cTemplate.gromacsConfigTemplate(
                        forceFieldType, temperature, 5000000,
                        membraneProtein=membraneProtein,
                        generateVelocities=False
                    )
                )

    def _generateGeometricColvarsConfig(
        self, 
        path, 
        topType, 
        coorType, 
        selectionPro, 
        selectionLig, 
        selectionRef, 
        stratification=[1,1,1,1,1,1,1,1], 
        pinDownPro=True,
        useOldCv=True,
        reflectingBoundary=False,
        unit='namd',
        considerRMSDCV=True,
        reweightAMD=False
    ):
        """generate Colvars config fils for geometric route

        Args:
            path (str): the directory for generation of all the files
            topType (str): the type (psf, parm) of the topology file
            coorType (str): the type (pdb, rst) of the topology file
            selectionPro (str): MDAnalysis-style selection of the protein
            selectionLig (str): MDAnalysis-style selection of the ligand
            selectionRef (str): MDAnalysis-style selection of the reference group for pulling simulation
            stratification (list of int, 8, optional): number of windows for each simulation. 
                                                       Defaults to [1,1,1,1,1,1,1,1].
            pinDownPro (bool, optional): Whether pinning down the protein. Defaults to True.
            useOldCv (bool, optional): whether used old, custom-function-based cv. Defaults to True.
            reflectingBoundary (bool, optional): Whether use reflecting boundaries, requires setBoundary on. Default to False
            unit (str, optional): unit, 'namd' or 'gromacs'. Default to namd.
            considerRMSDCV (bool, optional): Whether consider the RMSD CV. Default to True.
            reweightAMD (bool, optional): Whether add reweightAMD action. If this option is on, colvars config files for
                                          equilibration will not be generated. Default to False
        """    

        assert(len(stratification) == 8)
        
        if unit == 'namd':
            if not reweightAMD:
                colvarPostfix = 'in'
            else:
                colvarPostfix = 'in.amd'
        elif unit == 'gromacs':
            colvarPostfix = 'dat'
            assert(not reweightAMD)

        # read the original topology and coordinate file
        fParser = fileParser.fileParser(f'{path}/BFEE/complex.{topType}', f'{path}/BFEE/complex.pdb')
        center = fParser.measureCenter(selectionPro)
        polarAngles = fParser.measurePolarAngles(selectionRef, selectionLig)
        distance = fParser.measureDistance(selectionRef, selectionLig)
        
        if os.path.exists(f'{path}/BFEE/007_r/complex_largeBox.{topType}'):
            fParserLargeBox = fileParser.fileParser(f'{path}/BFEE/007_r/complex_largeBox.pdb')
            centerLargeBox = fParserLargeBox.measureCenter(selectionPro)
        else:
            centerLargeBox = center

        # 000_eq
        if not reweightAMD:
            with open(f'{path}/BFEE/000_eq/colvars.{colvarPostfix}', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx')
                )
                # monitor CVs during equilibration
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(
                        True, 0, 10, '../complex.xyz',
                        extendedLagrangian = False,
                        unit = unit
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        True, -90, 90, 'eulerTheta', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        True, -180, 180, 'eulerPhi', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        True, -180, 180, 'eulerPsi', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        True, 0, 180, 'polarTheta', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        True, -180, 180, 'polarPhi', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvRTemplate(
                        True, (distance - 5 if distance - 5 > 0 else 0), 
                        distance + 22, extendedLagrangian = False,
                        unit = unit
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvHistogramTemplate('RMSD')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHistogramTemplate('eulerTheta')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHistogramTemplate('eulerPhi')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHistogramTemplate('eulerPsi')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHistogramTemplate('polarTheta')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHistogramTemplate('polarPhi')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHistogramTemplate('r')
                )
                colvarsConfig.write(
                    self.cTemplate.cvProteinTemplate(
                        center, '../complex.xyz', unit = unit
                    )
                )

        # 001_RMSDBound
        if considerRMSDCV:
            for i in range(stratification[0]):
                with open(f'{path}/BFEE/001_RMSDBound/colvars_{i+1}.{colvarPostfix}', 'w') as colvarsConfig:
                    colvarsConfig.write(
                        self.cTemplate.cvHeadTemplate('../complex.ndx', reweightAMD=reweightAMD)
                    )
                    colvarsConfig.write(
                        self.cTemplate.cvRMSDTemplate(
                            True, float(i)/stratification[0] * 3.0, float(i+1)/stratification[0] * 3.0, '../complex.xyz',
                            extendedLagrangian = True, reflectingBoundary = reflectingBoundary,
                            unit = unit
                        )
                    )
                    colvarsConfig.write(
                        self.cTemplate.cvABFTemplate('RMSD', unit = unit)
                    )
                    if reweightAMD:
                        colvarsConfig.write(
                            self.cTemplate.cvReweightAMDTemplate('RMSD')
                        )
                    if not reflectingBoundary:
                        colvarsConfig.write(
                            self.cTemplate.cvHarmonicWallsTemplate(
                                'RMSD', float(i)/stratification[0] * 3.0, float(i+1)/stratification[0] * 3.0,
                                unit = unit
                            )
                        )
                    if pinDownPro:
                        colvarsConfig.write(
                            self.cTemplate.cvProteinTemplate(center, '../complex.xyz', unit = unit)
                        )

        # 002_Theta
        for i in range(stratification[1]):
            with open(f'{path}/BFEE/002_EulerTheta/colvars_{i+1}.{colvarPostfix}', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx', reweightAMD=reweightAMD)
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(
                        False, '', '', '../complex.xyz',
                        extendedLagrangian = False,
                        unit = unit
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        True, float(i)/stratification[1] * 20 - 10, float(i+1)/stratification[1] * 20 - 10, 
                        'eulerTheta', '../complex.xyz', useOldCv, extendedLagrangian = True,
                        reflectingBoundary = reflectingBoundary
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('eulerTheta', unit = unit)
                )
                if reweightAMD:
                    colvarsConfig.write(
                        self.cTemplate.cvReweightAMDTemplate('eulerTheta')
                    )
                if not reflectingBoundary:
                    colvarsConfig.write(
                        self.cTemplate.cvHarmonicWallsTemplate(
                            'eulerTheta', 
                            float(i)/stratification[1] * 20 - 10, 
                            float(i+1)/stratification[1] * 20 - 10,
                            unit = unit,
                        )
                    )
                if considerRMSDCV:
                    colvarsConfig.write(
                        self.cTemplate.cvHarmonicTemplate('RMSD', 10, 0, unit = unit)
                    )
                if pinDownPro:
                    colvarsConfig.write(
                        self.cTemplate.cvProteinTemplate(center, '../complex.xyz', unit = unit)
                    )

        # 003_Phi
        for i in range(stratification[2]):
            with open(f'{path}/BFEE/003_EulerPhi/colvars_{i+1}.{colvarPostfix}', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx', reweightAMD=reweightAMD)
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(
                        False, '', '', '../complex.xyz',
                        extendedLagrangian = False,
                        unit = unit
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerTheta', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        True, float(i)/stratification[2] * 20 - 10, float(i+1)/stratification[2] * 20 - 10, 
                        'eulerPhi', '../complex.xyz', useOldCv, extendedLagrangian = True,
                         reflectingBoundary = reflectingBoundary
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('eulerPhi', unit = unit)
                )
                if reweightAMD:
                    colvarsConfig.write(
                        self.cTemplate.cvReweightAMDTemplate('eulerPhi')
                    )
                if not reflectingBoundary:
                    colvarsConfig.write(
                        self.cTemplate.cvHarmonicWallsTemplate(
                            'eulerPhi', 
                            float(i)/stratification[2] * 20 - 10, 
                            float(i+1)/stratification[2] * 20 - 10,
                            unit = unit,
                        )
                    )
                if considerRMSDCV:
                    colvarsConfig.write(
                        self.cTemplate.cvHarmonicTemplate('RMSD', 10, 0, unit = unit)
                    )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerTheta', 0.1, 0, unit = unit)
                )
                if pinDownPro:
                    colvarsConfig.write(
                        self.cTemplate.cvProteinTemplate(center, '../complex.xyz', unit = unit)
                    )

        # 004_Psi
        for i in range(stratification[3]):
            with open(f'{path}/BFEE/004_EulerPsi/colvars_{i+1}.{colvarPostfix}', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx', reweightAMD=reweightAMD)
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(
                        False, '', '', '../complex.xyz',
                        extendedLagrangian = False,
                        unit = unit
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerTheta', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerPhi', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        True, float(i)/stratification[3] * 20 - 10, float(i+1)/stratification[3] * 20 - 10, 
                        'eulerPsi', '../complex.xyz', useOldCv, extendedLagrangian = True,
                        reflectingBoundary = reflectingBoundary
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('eulerPsi', unit = unit)
                )
                if reweightAMD:
                    colvarsConfig.write(
                        self.cTemplate.cvReweightAMDTemplate('eulerPsi')
                    )
                if not reflectingBoundary:
                    colvarsConfig.write(
                        self.cTemplate.cvHarmonicWallsTemplate(
                            'eulerPsi', 
                            float(i)/stratification[3] * 20 - 10, 
                            float(i+1)/stratification[3] * 20 - 10,
                            unit = unit
                        )
                    )
                if considerRMSDCV:
                    colvarsConfig.write(
                        self.cTemplate.cvHarmonicTemplate('RMSD', 10, 0, unit = unit)
                    )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerTheta', 0.1, 0, unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerPhi', 0.1, 0, unit = unit)
                )
                if pinDownPro:
                    colvarsConfig.write(
                        self.cTemplate.cvProteinTemplate(center, '../complex.xyz', unit = unit)
                    )

        # 005_polarTheta
        for i in range(stratification[4]):
            with open(f'{path}/BFEE/005_PolarTheta/colvars_{i+1}.{colvarPostfix}', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx', reweightAMD=reweightAMD)
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(
                        False, '', '', '../complex.xyz',
                        extendedLagrangian = False,
                        unit = unit
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerTheta', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerPhi', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerPsi', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        True, float(i)/stratification[4] * 20 - 10 + polarAngles[0], 
                        float(i+1)/stratification[4] * 20 - 10 + polarAngles[0], 'polarTheta', 
                        '../complex.xyz', useOldCv, extendedLagrangian = True,
                        reflectingBoundary = reflectingBoundary
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('polarTheta', unit = unit)
                )
                if reweightAMD:
                    colvarsConfig.write(
                        self.cTemplate.cvReweightAMDTemplate('polarTheta')
                    )
                if not reflectingBoundary:
                    colvarsConfig.write(
                        self.cTemplate.cvHarmonicWallsTemplate(
                            'polarTheta', float(i)/stratification[4] * 20 - 10 + polarAngles[0],
                            float(i+1)/stratification[4] * 20 - 10 + polarAngles[0],
                            unit = unit
                        )
                    )
                if considerRMSDCV:
                    colvarsConfig.write(
                        self.cTemplate.cvHarmonicTemplate('RMSD', 10, 0, unit = unit)
                    )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerTheta', 0.1, 0, unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerPhi', 0.1, 0, unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerPsi', 0.1, 0, unit = unit)
                )
                if pinDownPro:
                    colvarsConfig.write(
                        self.cTemplate.cvProteinTemplate(center, '../complex.xyz', unit = unit)
                    )

        # 006_polarPhi
        for i in range(stratification[5]):
            with open(f'{path}/BFEE/006_PolarPhi/colvars_{i+1}.{colvarPostfix}', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx', reweightAMD=reweightAMD)
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(
                        False, '', '', '../complex.xyz',
                        extendedLagrangian = False,
                        unit = unit
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerTheta', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerPhi', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerPsi', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'polarTheta', '../complex.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        True, float(i)/stratification[5] * 20 - 10 + polarAngles[1], 
                        float(i+1)/stratification[5] * 20 - 10 + polarAngles[1], 'polarPhi', 
                        '../complex.xyz', useOldCv, extendedLagrangian = True,
                        reflectingBoundary = reflectingBoundary
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('polarPhi', unit = unit)
                )
                if reweightAMD:
                    colvarsConfig.write(
                        self.cTemplate.cvReweightAMDTemplate('polarPhi')
                    )
                if not reflectingBoundary:
                    colvarsConfig.write(
                        self.cTemplate.cvHarmonicWallsTemplate(
                            'polarPhi', float(i)/stratification[5] * 20 - 10 + polarAngles[1],
                            float(i+1)/stratification[5] * 20 - 10 + polarAngles[1],
                            unit = unit
                        )
                    )
                if considerRMSDCV:
                    colvarsConfig.write(
                        self.cTemplate.cvHarmonicTemplate('RMSD', 10, 0, unit = unit)
                    )   
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerTheta', 0.1, 0, unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerPhi', 0.1, 0, unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerPsi', 0.1, 0, unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('polarTheta', 0.1, polarAngles[0], unit = unit)
                )
                if pinDownPro:
                    colvarsConfig.write(
                        self.cTemplate.cvProteinTemplate(center, '../complex.xyz', unit = unit)
                    )

        # 007_r
        if not reweightAMD:
            # eq
            with open(f'{path}/BFEE/007_r/colvars_eq.{colvarPostfix}', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx')
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(
                        False, '', '', './complex_largeBox.xyz',
                        extendedLagrangian = False,
                        unit = unit
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerTheta', './complex_largeBox.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerPhi', './complex_largeBox.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerPsi', './complex_largeBox.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'polarTheta', './complex_largeBox.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'polarPhi', './complex_largeBox.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvRTemplate(
                        False, 0, 0, extendedLagrangian = False, unit = unit
                    )
                )
                if considerRMSDCV:
                    colvarsConfig.write(
                        self.cTemplate.cvHarmonicTemplate('RMSD', 10, 0, unit = unit)
                    )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerTheta', 0.1, 0, unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerPhi', 0.1, 0, unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerPsi', 0.1, 0, unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('polarTheta', 0.1, polarAngles[0], unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('polarPhi', 0.1, polarAngles[1], unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvProteinTemplate(centerLargeBox, './complex_largeBox.xyz', unit = unit)
                )

        # abf
        for i in range(stratification[6]):
            with open(f'{path}/BFEE/007_r/colvars_{i+1}.{colvarPostfix}', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx', reweightAMD=reweightAMD)
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(
                        False, '', '', './complex_largeBox.xyz',
                        extendedLagrangian = False,
                        unit = unit
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerTheta', './complex_largeBox.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerPhi', './complex_largeBox.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'eulerPsi', './complex_largeBox.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'polarTheta', './complex_largeBox.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvAngleTemplate(
                        False, 0, 0, 'polarPhi', './complex_largeBox.xyz', useOldCv,
                        extendedLagrangian = False
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvRTemplate(
                        True, float(i)/stratification[6] * 24 - 2 + distance, 
                        float(i+1)/stratification[6] * 24 - 2 + distance,
                        extendedLagrangian = True, reflectingBoundary = reflectingBoundary,
                        unit = unit
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('r', unit = unit)
                )
                if reweightAMD:
                    colvarsConfig.write(
                        self.cTemplate.cvReweightAMDTemplate('r')
                    )
                if not reflectingBoundary:
                    colvarsConfig.write(
                        self.cTemplate.cvHarmonicWallsTemplate(
                            'r', float(i)/stratification[6] * 24 - 2 + distance, 
                            float(i+1)/stratification[6] * 24 - 2 + distance,
                            unit = unit
                        )
                    )
                if considerRMSDCV:
                    colvarsConfig.write(
                        self.cTemplate.cvHarmonicTemplate('RMSD', 10, 0, unit = unit)
                    )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerTheta', 0.1, 0, unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerPhi', 0.1, 0, unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerPsi', 0.1, 0, unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('polarTheta', 0.1, polarAngles[0], unit = unit)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('polarPhi', 0.1, polarAngles[1], unit = unit)
                )
                if pinDownPro:
                    colvarsConfig.write(
                        self.cTemplate.cvProteinTemplate(centerLargeBox, './complex_largeBox.xyz', unit = unit)
                    )

        # 008_RMSDUnbound
        if considerRMSDCV:
            for i in range(stratification[7]):
                with open(f'{path}/BFEE/008_RMSDUnbound/colvars_{i+1}.{colvarPostfix}', 'w') as colvarsConfig:
                    colvarsConfig.write(
                        self.cTemplate.cvHeadTemplate('./ligandOnly.ndx', reweightAMD=reweightAMD)
                    )
                    colvarsConfig.write(
                        self.cTemplate.cvRMSDTemplate(
                            True, float(i)/stratification[7] * 3.0, float(i+1)/stratification[7] * 3.0, './ligandOnly.xyz',
                            extendedLagrangian = True, reflectingBoundary = reflectingBoundary,
                            unit = unit
                        )
                    )
                    colvarsConfig.write(
                        self.cTemplate.cvABFTemplate('RMSD', unit = unit)
                    )
                    if reweightAMD:
                        colvarsConfig.write(
                            self.cTemplate.cvReweightAMDTemplate('RMSD')
                        )
                    if not reflectingBoundary:
                        colvarsConfig.write(
                            self.cTemplate.cvHarmonicWallsTemplate(
                                'RMSD', float(i)/stratification[7] * 3.0, float(i+1)/stratification[7] * 3.0,
                                unit = unit
                            )
                        )

    def _duplicateFileFolder(self, path, number):
        """duplicate the ./BFEE folder

        Args:
            path (string): the directory for generation of all the files
            number (int): the number of copies

        Raises:
            DirectoryExistError: if {path}/BFEE_{i} exists
        """        
        for i in range(number - 1):
            if os.path.exists(f'{path}/BFEE_{i}'):
                raise DirectoryExistError('Directory exists')
            shutil.copytree(f'{path}/BFEE', f'{path}/BFEE_{i}')
