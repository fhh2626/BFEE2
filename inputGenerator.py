# generate all the inputs and define corresponding slots

import os, sys, shutil, subprocess
import numpy as np
from commonTools import fileParser
from templates import configTemplate, scriptTemplate
from templates.gromacs.BFEEGromacs import BFEEGromacs

class inputGenerator():
    ''' generate all the inputs and define corresponding slots '''

    def __init__(self):
        ''' simply initialize the configTemplate objects '''

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
        vmdPath = ''
    ):
        ''' generate all the input files for NAMD alchemical simulation
            Inputs:
                path (string): the directory for generation of all the files
                topFile (string): the path of the topology (psf, parm) file
                coorFile (string): the path of the coordinate (pdb, rst7) file
                forceFieldType (string): 'charmm' or 'amber'
                forceFieldFiles (list of strings): list of CHARMM force field files
                temperature (float): temperature of the simulation
                selectionPro (string): MDAnalysis-style selection of the protein
                selectionLig (string): MDAnalysis-style selection of the ligand
                stratification (list of int, 8): number of windows for each simulation
                vmdPath (string): path to vmd '''

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
            'alchemical', vmdPath)

        # get relative force field path
        relativeFFPath = []
        for item in forceFieldFiles:
            _, name = os.path.split(item)
            relativeFFPath.append(f'../{name}')
        
        self._generateAlchemicalNAMDConfig(
            path, forceFieldType, relativeFFPath, temperature, stratification
        )
        self._generateAlchemicalColvarsConfig(
            path, topType, 'pdb', selectionPro, selectionLig, selectionPro, stratification
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
        vmdPath = ''
    ):
        ''' generate all the input files for NAMD Geometric simulation
            Inputs:
                path (string): the directory for generation of all the files
                topFile (string): the path of the topology (psf, parm) file
                coorFile (string): the path of the coordinate (pdb, rst7) file
                forceFieldType (string): 'charmm' or 'amber'
                forceFieldFiles (list of strings): list of CHARMM force field files
                temperature (float): temperature of the simulation
                selectionPro (string): MDAnalysis-style selection of the protein
                selectionLig (string): MDAnalysis-style selection of the ligand
                selectionRef (string): MDAnalysis-style selection of the reference group for pulling simulation,
                                       by default, this is the protein
                userProvidedPullingTop (string): user-provided large solvation box for pulling simulation
                userProvidedPullingCoor (string): user-provided large solvation box for pulling simulation
                stratification (list of int, 8): number of windows for each simulation
                vmdPath (string): path to vmd '''

        assert(len(stratification) == 8)
        assert(forceFieldType == 'charmm' or forceFieldType == 'amber')

        if selectionRef == '':
            selectionRef = selectionPro

        # determine the type of input top and coor file
        # However, whatever coorType is, NAMD only read pdb as coordinate
        topType, coorType = self._determineFileType(topFile, coorFile)

        self._makeDirectories(path, 'geometric')
        self._copyFiles(
            path, topFile, topType, coorFile, coorType, forceFieldType, forceFieldFiles,
            selectionPro, selectionLig, selectionRef, userProvidedPullingTop, userProvidedPullingCoor,
            'geometric', vmdPath)

        # get relative force field path
        relativeFFPath = []
        for item in forceFieldFiles:
            _, name = os.path.split(item)
            relativeFFPath.append(f'../{name}')
        
        self._generateGeometricNAMDConfig(
            path, forceFieldType, relativeFFPath, temperature, stratification
        )
        self._generateGeometricColvarsConfig(
            path, topType, coorType, selectionPro, selectionLig, selectionRef, stratification
        )

    def generateGromacsGeometricFiles(
        self, 
        path, 
        topFile, 
        pdbFile, 
        ligandOnlyTopFile, 
        ligandOnlyPdbFile,
        selectionPro,
        selectionLig,
        selectionSol='resname TIP3* or resname SPC*',
    ):
        ''' generate all the input files for Gromacs Geometric simulation
            This function is based on BFEEGromacs.py
            contributed by Haochuan Chen (yjcoshc_at_mail.nankai.edu.cn)
            Inputs:
                path (string): the directory for generation of all the files
                topFile (string): the path of the topology (top) file
                pdbFile (string): the path of the coordinate (pdb) file
                ligandOnlyTopFile (string): the path of the ligand-only topology (top) file
                ligandOnlyPdbFile (string): the path of the ligand-only coordinate (pdb) file
                selectionPro (string): MDAnalysis-style selection of the protein
                selectionLig (string): MDAnalysis-style selection of the ligand
                selectionSol (string): MDAnalysis-style selection of the solvent '''
        
        # if the folder already exists, simply remove it
        if os.path.exists(f'{path}/BFEE'):
            shutil.rmtree(f'{path}/BFEE')
        os.mkdir(f'{path}/BFEE')

        bfee = BFEEGromacs(
            pdbFile, 
            topFile, 
            ligandOnlyPdbFile, 
            ligandOnlyTopFile, 
            f'{path}/BFEE'
        )
        bfee.setProteinHeavyAtomsGroup(f'{selectionPro} and not (name H*)')
        bfee.setLigandHeavyAtomsGroup(f'{selectionLig} and not (name H*)')
        bfee.setSolventAtomsGroup(selectionSol)
        bfee.generate001()
        bfee.generate002()
        bfee.generate003()
        bfee.generate004()
        bfee.generate005()
        bfee.generate006()
        bfee.generate007()
        bfee.generate008()


    def _determineFileType(self, topFile, coorFile):
        ''' determine the file type of topology and coordinate files
            Inputs:
                topFile (string): the path of the topology (psf, parm) file
                coorFile (string): the path of the coordinate (pdb, rst7) file
            Return:
                tuple of string, 2: (topType, coorType), e.g., (psf, pdb) '''

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
        elif coorPostfix == '.rst7' or coorPostfix == '.rst' or coorPostfix == '.inpcrd':
            coorType = 'inpcrd'

        assert(topType != '' and coorType != '')

        return topType, coorType

    def _makeDirectories(self, path, jobType='geometric'):
        ''' make directories for BFEE calculation
            Inputs:
                path (string): the path for putting BFEE input files into
                jobType (string): geometric or alchemical '''

        # if the folder already exists, simply remove it
        if os.path.exists(f'{path}/BFEE'):
            shutil.rmtree(f'{path}/BFEE')
        os.mkdir(f'{path}/BFEE')
        os.mkdir(f'{path}/BFEE/000_eq')
        os.mkdir(f'{path}/BFEE/000_eq/output')
        if jobType == 'geometric':
            os.mkdir(f'{path}/BFEE/001_RMSDBound')
            os.mkdir(f'{path}/BFEE/002_EulerTheta')
            os.mkdir(f'{path}/BFEE/003_EulerPhi')
            os.mkdir(f'{path}/BFEE/004_EulerPsi')
            os.mkdir(f'{path}/BFEE/005_PolarTheta')
            os.mkdir(f'{path}/BFEE/006_PolarPhi')
            os.mkdir(f'{path}/BFEE/007_r')
            os.mkdir(f'{path}/BFEE/008_RMSDUnbound')
            os.mkdir(f'{path}/BFEE/001_RMSDBound/output')
            os.mkdir(f'{path}/BFEE/002_EulerTheta/output')
            os.mkdir(f'{path}/BFEE/003_EulerPhi/output')
            os.mkdir(f'{path}/BFEE/004_EulerPsi/output')
            os.mkdir(f'{path}/BFEE/005_PolarTheta/output')
            os.mkdir(f'{path}/BFEE/006_PolarPhi/output')
            os.mkdir(f'{path}/BFEE/007_r/output')
            os.mkdir(f'{path}/BFEE/008_RMSDUnbound/output')

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
        jobType='geometric',
        vmdPath = ''
    ):
        ''' copy original and generate necessary topology/structure files
            Inputs:
                path (string): the directory for generation of all the files
                topFile (string): the path of the topology (psf, parm) file
                topType (string): the type (psf, parm) of the topology file
                coorFile (string): the path of the coordinate (pdb, rst7) file
                coorType (string): the type (pdb, rst) of the coordinate file
                forceFieldType (string): 'charmm' or 'amber'
                forceFieldFiles (list of strings): list of CHARMM force field files
                selectionPro (string): MDAnalysis-style selection of the protein
                selectionLig (string): MDAnalysis-style selection of the ligand
                selectionRef (string): MDAnalysis-style selection of the reference group for pulling simulation
                userProvidedPullingTop (string): user-provided large solvation box for pulling simulation
                userProvidedPullingCoor (string): user-provided large solvation box for pulling simulation
                jobType (string): 'geometric' or 'alchemical'
                vmdPath (string): path to vmd, space is forbidden '''
        
        # copy force fields
        if forceFieldType == 'charmm':
            for item in forceFieldFiles:
                _, fileName = os.path.split(item)
                shutil.copyfile(item, f'{path}/BFEE/{fileName}')

        # read the original topology and coordinate file
        fParser = fileParser.fileParser(topFile, coorFile)
            
        # check polar angles
        # if theta < 30 or > 150, then rotate 90 degrees to avoid polar angle invarience
        if (fParser.measurePolarAngles(selectionPro, selectionLig)[0] > 150 \
            or fParser.measurePolarAngles(selectionPro, selectionLig)[0] < 30):
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


        if jobType == 'geometric':

            # remove protein for step 8
            # this cannot be done in pure python
            fParser.saveFile(
                f'not {selectionPro}', f'{path}/BFEE/008_RMSDUnbound/ligandOnly.pdb', 'pdb'
            )
            # connect to VMD
            if forceFieldType == 'charmm':
                with open( f'{path}/BFEE/008_RMSDUnbound/000_removeProtein.tcl', 'w') as rScript:
                    rScript.write(
                        scriptTemplate.removeProteinTemplate.substitute(
                            path='../complex', selectionPro=f'{selectionPro}'.replace('segid', 'segname'),
                            outputPath=f'./ligandOnly'
                        )
                    )
                # if vmd path is defined
                # then execute vmd automatically
                if vmdPath != '':
                    subprocess.run(
                        [vmdPath, '-dispdev', 'text', '-e', f'{path}/BFEE/008_RMSDUnbound/000_removeProtein.tcl'],
                        cwd=f'{path}/BFEE/008_RMSDUnbound'
                    )
            elif forceFieldType == 'amber':
                with open( f'{path}/BFEE/008_RMSDUnbound/000_removeProtein.cpptraj', 'w') as rScript:
                    rScript.write(
                        scriptTemplate.removeProteinAmberTemplate.substitute(
                            path='../complex', 
                            residueNum=fParser.getResid(selectionPro),
                            outputPath=f'./ligandOnly'
                        )
                    )
            # xyz and ndx for step 8
            fParserStep8 = fileParser.fileParser(f'{path}/BFEE/008_RMSDUnbound/ligandOnly.pdb')
            fParserStep8.saveFile('all', f'{path}/BFEE/008_RMSDUnbound/ligandOnly.xyz', 'xyz')
            fParserStep8.saveNDX(
                [selectionLig], ['ligand'], f'{path}/BFEE/008_RMSDUnbound/ligandOnly.ndx', True
            )
            
            # add water for step 7
            # this cannot be done in pure python
            # connect to VMD
            if userProvidedPullingTop == '' and userProvidedPullingCoor == '':
                if forceFieldType == 'charmm':
                    shutil.copyfile(f'{sys.path[0]}/templates/scripts/solvate.tcl', f'{path}/BFEE/007_r/000_solvate.tcl')
                    # if vmd path is defined
                    # then execute vmd automatically
                    if vmdPath != '':
                        subprocess.run(
                            [vmdPath, '-dispdev', 'text', '-e', f'{path}/BFEE/007_r/000_solvate.tcl'],
                            cwd=f'{path}/BFEE/007_r'
                        )
            else:
                # copy the user-provided large box
                fParserStep7 = fileParser.fileParser(userProvidedPullingTop, userProvidedPullingCoor)

                # check polar angles
                # if theta < 30 or > 150, then rotate 90 degrees to avoid polar angle invarience
                if (fParserStep7.measurePolarAngles(selectionPro, selectionLig)[0] > 150 \
                    or fParserStep7.measurePolarAngles(selectionPro, selectionLig)[0] < 30):
                    fParserStep7.rotateSystem('x', 90)

                fParserStep7.saveFile(
                    'all', f'{path}/BFEE/007_r/complexLargeBox.pdb', 'pdb'
                )

        if jobType == 'alchemical':
            # fep file
            fParser.setBeta('all', 0)
            fParser.setBeta(selectionLig, 1)
            fParser.saveFile(
                'all', f'{path}/BFEE/fep.pdb', 'pdb'
            )
            shutil.copyfile(f'{sys.path[0]}/templates/scripts/fep.tcl', f'{path}/BFEE/fep.tcl')

            # remove protein for the unbound state
            if forceFieldType == 'charmm':
                with open(f'{path}/BFEE/000_removeProtein.tcl', 'w') as rScript:
                    rScript.write(
                        scriptTemplate.removeProteinTemplate.substitute(
                            path='./complex', selectionPro=f'{selectionPro}'.replace('segid', 'segname'),
                            outputPath=f'./ligandOnly'
                        )
                    )
                # if vmd path is defined
                # then execute vmd automatically
                if vmdPath != '':
                    subprocess.run(
                        [vmdPath, '-dispdev', 'text', '-e', f'{path}/BFEE/000_removeProtein.tcl'],
                        cwd=f'{path}/BFEE'
                    )
            elif forceFieldType == 'amber':
                with open( f'{path}/BFEE/000_removeProtein.cpptraj', 'w') as rScript:
                    rScript.write(
                        scriptTemplate.removeProteinAmberTemplate.substitute(
                            path='./complex', 
                            residueNum=fParser.getResid(selectionPro),
                            outputPath=f'./ligandOnly'
                        )
                    )
            fParser.saveFile(
                f'not {selectionPro}', f'{path}/BFEE/ligandOnly.pdb', 'pdb'
            )
            fParser.saveFile(
                f'not {selectionPro}', f'{path}/BFEE/fep_ligandOnly.pdb', 'pdb'
            )
            # xyz and ndx
            fParserLigandOnly = fileParser.fileParser( f'{path}/BFEE/ligandOnly.pdb')
            fParserLigandOnly.saveFile('all', f'{path}/BFEE/ligandOnly.xyz', 'xyz')
            fParserLigandOnly.saveNDX(
                [selectionLig], ['ligand'], f'{path}/BFEE/ligandOnly.ndx', True
            )

    def _generateAlchemicalNAMDConfig(
        self,
        path,
        forceFieldType,
        forceFields,
        temperature,
        stratification = [1,1,1,1]
    ):
        ''' generate NAMD config fils for the alchemical route
            Inputs:
                path (string): the directory for generation of all the files
                forceFieldType (string): 'charmm' or 'amber'
                forceFieldFiles (list of strings): list of CHARMM force field files
                temperature (float): temperature of the simulation
                stratification (list of int, 4): number of windows for each simulation  '''

        if forceFieldType == 'charmm':
            topType = 'psf'
        elif forceFieldType == 'amber':
            topType = 'parm7'
        coorType = 'pdb'

        # read the original topology and coordinate file
        fParser = fileParser.fileParser(f'{path}/BFEE/complex.{topType}', f'{path}/BFEE/complex.{coorType}')

        # pbc
        pbc = fParser.measurePBC()

        # 000_eq
        with open(f'{path}/BFEE/000_eq/eq.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    '', '', '', pbc,
                    'output/eq', temperature, 5000000, 'colvars.in', ''
                )
            )
        with open(f'{path}/BFEE/000_eq/eq_ligandOnly.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../ligandOnly.{topType}', f'../ligandOnly.pdb',
                    '', '', '', pbc,
                    'output/eq_ligandOnly', temperature, 1000000, 'colvars_ligandOnly.in',
                    ''
                )
            )

        # 001_MoleculeBound
        with open(f'{path}/BFEE/001_MoleculeBound/fep_forward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc', '',
                    'output/fep_forward', temperature, 0, 'colvars.in', '', '', '../fep.pdb', 
                    stratification[0], True
                )
            )
        with open(f'{path}/BFEE/001_MoleculeBound/fep_backward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc', '',
                    'output/fep_backward', temperature, 0, 'colvars.in', '', '', '../fep.pdb', 
                    stratification[0], False
                )
            )

        # 002_RestraintBound
        with open(f'{path}/BFEE/002_RestraintBound/ti_forward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc', '',
                    'output/ti_forward', temperature, f'{500000*(stratification[1]+1)}', 'colvars_forward.in', 
                    ''
                )
            )
        with open(f'{path}/BFEE/002_RestraintBound/ti_backward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc', '',
                    'output/ti_backward', temperature, f'{500000*(stratification[1]+1)}', 'colvars_backward.in', 
                    ''
                )
            )

        # 003_MoleculeUnbound
        with open(f'{path}/BFEE/003_MoleculeUnbound/fep_forward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../ligandOnly.{topType}', f'../ligandOnly.pdb',
                    f'../000_eq/output/eq_ligandOnly.coor', f'../000_eq/output/eq_ligandOnly.vel', 
                    f'../000_eq/output/eq_ligandOnly.xsc', '',
                    'output/fep_forward', temperature, 0, 'colvars.in', '', '', '../fep_ligandOnly.pdb', 
                    stratification[0], True
                )
            )
        with open(f'{path}/BFEE/003_MoleculeUnbound/fep_backward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../ligandOnly.{topType}', f'../ligandOnly.pdb',
                    f'../000_eq/output/eq_ligandOnly.coor', f'../000_eq/output/eq_ligandOnly.vel', 
                    f'../000_eq/output/eq_ligandOnly.xsc', '',
                    'output/fep_backward', temperature, 0, 'colvars.in', '', '', '../fep_ligandOnly.pdb', 
                    stratification[0], False
                )
            )

        # 004_RestraintUnbound
        with open(f'{path}/BFEE/004_RestraintUnbound/ti_forward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../ligandOnly.{topType}', f'../ligandOnly.pdb',
                    f'../000_eq/output/eq_ligandOnly.coor', f'../000_eq/output/eq_ligandOnly.vel', 
                    f'../000_eq/output/eq_ligandOnly.xsc', '',
                    'output/ti_forward', temperature, f'{500000*(stratification[3]+1)}', 'colvars_forward.in', 
                    ''
                )
            )
        with open(f'{path}/BFEE/004_RestraintUnbound/ti_backward.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../ligandOnly.{topType}', f'../ligandOnly.pdb',
                    f'../000_eq/output/eq_ligandOnly.coor', f'../000_eq/output/eq_ligandOnly.vel', 
                    f'../000_eq/output/eq_ligandOnly.xsc', '',
                    'output/ti_backward', temperature, f'{500000*(stratification[3]+1)}', 'colvars_backward.in', 
                    ''
                )
            )

    def _generateAlchemicalColvarsConfig(
        self, path, topType, coorType, selectionPro, selectionLig, selectionRef, stratification=[1,1,1,1]
    ):
        ''' generate Colvars config fils for geometric route
            Inputs:
                path (string): the directory for generation of all the files
                topType (string): the type (psf, parm) of the topology file
                coorType (string): the type (pdb, rst) of the topology file
                selectionPro (string): MDAnalysis-style selection of the protein
                selectionLig (string): MDAnalysis-style selection of the ligand
                selectionRef (string): MDAnalysis-style selection of the reference group for pulling simulation
                stratification (list of int, 4): number of windows for each simulation ''' 

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
                self.cTemplate.cvRMSDTemplate(False, '', '', '../complex.xyz')
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerTheta', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerPhi', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerPsi', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvPolarAngleTemplate(
                    False, 0, 0, 'polarTheta', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvPolarAngleTemplate(
                    False, 0, 0, 'polarPhi', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvRTemplate(
                    False, 0, 0
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
            colvarsConfig.write(
                self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
            )

        with open(f'{path}/BFEE/000_eq/colvars_ligandOnly.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../ligandOnly.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(False, '', '', '../ligandOnly.xyz')
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
                self.cTemplate.cvRMSDTemplate(False, '', '', '../complex.xyz')
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerTheta', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerPhi', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerPsi', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvPolarAngleTemplate(
                    False, 0, 0, 'polarTheta', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvPolarAngleTemplate(
                    False, 0, 0, 'polarPhi', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvRTemplate(
                    False, 0, 0
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
            colvarsConfig.write(
                self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
            )

        # 002_RestraintBound
        with open(f'{path}/BFEE/002_RestraintBound/colvars_forward.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../complex.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(False, '', '', '../complex.xyz')
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerTheta', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerPhi', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerPsi', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvPolarAngleTemplate(
                    False, 0, 0, 'polarTheta', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvPolarAngleTemplate(
                    False, 0, 0, 'polarPhi', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvRTemplate(
                    False, 0, 0
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
            colvarsConfig.write(
                self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
            )
        with open(f'{path}/BFEE/002_RestraintBound/colvars_backward.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../complex.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(False, '', '', '../complex.xyz')
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerTheta', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerPhi', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerPsi', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvPolarAngleTemplate(
                    False, 0, 0, 'polarTheta', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvPolarAngleTemplate(
                    False, 0, 0, 'polarPhi', '../complex.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvRTemplate(
                    False, 0, 0
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
            colvarsConfig.write(
                self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
            )

        # 003_MoleculeUnbound
        with open(f'{path}/BFEE/003_MoleculeUnbound/colvars.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../ligandOnly.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(False, '', '', '../ligandOnly.xyz')
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
                self.cTemplate.cvRMSDTemplate(False, '', '', '../ligandOnly.xyz')
            )
            colvarsConfig.write(
                self.cTemplate.cvHarmonicTemplate('RMSD', 0, 0, stratification[3], True, 10)
            )
        with open(f'{path}/BFEE/004_RestraintUnbound/colvars_backward.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../ligandOnly.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(False, '', '', '../ligandOnly.xyz')
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
        stratification = [1,1,1,1,1,1,1,1]
    ):
        ''' generate NAMD config fils for the geometric route
            Inputs:
                path (string): the directory for generation of all the files
                forceFieldType (string): 'charmm' or 'amber'
                forceFieldFiles (list of strings): list of CHARMM force field files
                temperature (float): temperature of the simulation
                stratification (list of int, 8): number of windows for each simulation  '''

        if forceFieldType == 'charmm':
            topType = 'psf'
        elif forceFieldType == 'amber':
            topType = 'parm7'
        coorType = 'pdb'

        # read the original topology and coordinate file
        fParser = fileParser.fileParser(f'{path}/BFEE/complex.{topType}', f'{path}/BFEE/complex.{coorType}')

        # pbc
        pbc = fParser.measurePBC()

        # 000_eq
        with open(f'{path}/BFEE/000_eq/eq.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    '', '', '', pbc,
                    'output/eq', temperature, 5000000, 'colvars.in'
                )
            )

        # RMSD bound
        with open(f'{path}/BFEE/001_RMSDBound/abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc',
                    '', 'output/abf_1', temperature, 5000000, 'colvars_1.in'
                )
            )

        # stratification
        if stratification[0] > 1:
            for i in range(1, stratification[0]):
                with open(f'{path}/BFEE/001_RMSDBound/abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i}.coor', f'output/abf_{i}.vel', 
                        f'output/abf_{i}.xsc',
                        '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in'
                    )
                )

        # Theta
        with open(f'{path}/BFEE/002_EulerTheta/abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc',
                    '', 'output/abf_1', temperature, 5000000, 'colvars_1.in', ''
                )
            )

        # stratification
        if stratification[1] > 1:
            for i in range(1, stratification[1]):
                with open(f'{path}/BFEE/002_EulerTheta/abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i}.coor', f'output/abf_{i}.vel', 
                        f'output/abf_{i}.xsc',
                        '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in', ''
                    )
                )

        # Phi
        with open(f'{path}/BFEE/003_EulerPhi/abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc',
                    '', 'output/abf_1', temperature, 5000000, 'colvars_1.in', ''
                )
            )

        # stratification
        if stratification[2] > 1:
            for i in range(1, stratification[2]):
                with open(f'{path}/BFEE/003_EulerPhi/abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i}.coor', f'output/abf_{i}.vel', 
                        f'output/abf_{i}.xsc',
                        '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in', ''
                    )
                )

        # Psi
        with open(f'{path}/BFEE/004_EulerPsi/abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc',
                    '', 'output/abf_1', temperature, 5000000, 'colvars_1.in', ''
                )
            )

        # stratification
        if stratification[3] > 1:
            for i in range(1, stratification[3]):
                with open(f'{path}/BFEE/004_EulerPsi/abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i}.coor', f'output/abf_{i}.vel', 
                        f'output/abf_{i}.xsc',
                        '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in', ''
                    )
                )

        # theta
        with open(f'{path}/BFEE/005_PolarTheta/abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc',
                    '', 'output/abf_1', temperature, 5000000, 'colvars_1.in', ''
                )
            )

        # stratification
        if stratification[4] > 1:
            for i in range(1, stratification[4]):
                with open(f'{path}/BFEE/005_PolarTheta/abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i}.coor', f'output/abf_{i}.vel', 
                        f'output/abf_{i}.xsc',
                        '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in', ''
                    )
                )

        # phi
        with open(f'{path}/BFEE/006_PolarPhi/abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                    f'../000_eq/output/eq.coor', f'../000_eq/output/eq.vel', f'../000_eq/output/eq.xsc',
                    '', 'output/abf_1', temperature, 5000000, 'colvars_1.in', ''
                )
            )

        # stratification
        if stratification[5] > 1:
            for i in range(1, stratification[5]):
                with open(f'{path}/BFEE/006_PolarPhi/abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'../complex.{topType}', f'../complex.pdb',
                        f'output/abf_{i}.coor', f'output/abf_{i}.vel', 
                        f'output/abf_{i}.xsc',
                        '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in', ''
                    )
                )

        # r
        # eq
        with open(f'{path}/BFEE/007_r/eq.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'./complex_largeBox.{topType}', f'./complex_largeBox.pdb',
                    '', '', '', 
                    pbc + np.array([[22,22,22],[0,0,0]]),
                    'output/eq', temperature, 5000000, 'colvars_eq.in', ''
                )
            )
        # abf
        with open(f'{path}/BFEE/007_r/abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'./complex_largeBox.{topType}', f'./complex_largeBox.pdb',
                    'output/eq.coor', 'output/eq.vel', 'output/eq.xsc', '',
                    'output/abf_1', temperature, 5000000, 'colvars_1.in', ''
                )
            )

        # stratification
        if stratification[6] > 1:
            for i in range(1, stratification[6]):
                with open(f'{path}/BFEE/007_r/abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'./complex_largeBox.{topType}', f'./complex_largeBox.pdb',
                        f'output/abf_{i}.coor', f'output/abf_{i}.vel', 
                        f'output/abf_{i}.xsc',
                        '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in', ''
                    )
                )

        # RMSD unbound
        # eq
        with open(f'{path}/BFEE/008_RMSDUnbound/eq.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'./ligandOnly.{topType}', f'./ligandOnly.pdb',
                    '', '', '', 
                    pbc,
                    'output/eq', temperature, 1000000
                )
            )
        # abf
        with open(f'{path}/BFEE/008_RMSDUnbound/abf_1.conf', 'w') as namdConfig:
            namdConfig.write(
                self.cTemplate.namdConfigTemplate(
                    forceFieldType, forceFields, f'./ligandOnly.{topType}', f'./ligandOnly.pdb',
                    'output/eq.coor', 'output/eq.vel', 'output/eq.xsc', '',
                    'output/abf_1', temperature, 5000000, 'colvars_1.in'
                )
            )

        # stratification
        if stratification[7] > 1:
            for i in range(1, stratification[7]):
                with open(f'{path}/BFEE/008_RMSDUnbound/abf_{i+1}.conf', 'w') as namdConfig:
                    namdConfig.write(
                        self.cTemplate.namdConfigTemplate(
                        forceFieldType, forceFields, f'./ligandOnly.{topType}', f'./ligandOnly.pdb',
                        f'output/abf_{i}.coor', f'output/abf_{i}.vel', 
                        f'output/abf_{i}.xsc',
                        '', f'output/abf_{i+1}', temperature, 5000000, f'colvars_{i+1}.in'
                    )
                )

    def _generateGeometricColvarsConfig(
        self, path, topType, coorType, selectionPro, selectionLig, selectionRef, stratification=[1,1,1,1,1,1,1,1]
    ):
        ''' generate Colvars config fils for geometric route
            Inputs:
                path (string): the directory for generation of all the files
                topType (string): the type (psf, parm) of the topology file
                coorType (string): the type (pdb, rst) of the topology file
                selectionPro (string): MDAnalysis-style selection of the protein
                selectionLig (string): MDAnalysis-style selection of the ligand
                selectionRef (string): MDAnalysis-style selection of the reference group for pulling simulation
                stratification (list of int, 8): number of windows for each simulation ''' 

        assert(len(stratification) == 8)

        # read the original topology and coordinate file
        fParser = fileParser.fileParser(f'{path}/BFEE/complex.{topType}', f'{path}/BFEE/complex.pdb')
        center = fParser.measureCenter(selectionPro)
        polarAngles = fParser.measurePolarAngles(selectionRef, selectionLig)
        distance = fParser.measureDistance(selectionRef, selectionLig)

        # 000_eq
        with open(f'{path}/BFEE/000_eq/colvars.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../complex.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
            )

        # 001_RMSDBound
        for i in range(stratification[0]):
            with open(f'{path}/BFEE/001_RMSDBound/colvars_{i+1}.in', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx')
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(
                        True, float(i)/stratification[0] * 3.0, float(i+1)/stratification[0] * 3.0, '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('RMSD')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicWallsTemplate(
                        'RMSD', float(i)/stratification[0] * 3.0, float(i+1)/stratification[0] * 3.0
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
                )

        # 002_Theta
        for i in range(stratification[1]):
            with open(f'{path}/BFEE/002_EulerTheta/colvars_{i+1}.in', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx')
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(False, '', '', '../complex.xyz')
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        True, float(i)/stratification[1] * 20 - 10, float(i+1)/stratification[1] * 20 - 10, 'eulerTheta', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('eulerTheta')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicWallsTemplate(
                        'eulerTheta', float(i)/stratification[1] * 20 - 10, float(i+1)/stratification[1] * 20 - 10
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('RMSD', 10, 0)
                )
                colvarsConfig.write(
                    self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
                )

        # 003_Phi
        for i in range(stratification[2]):
            with open(f'{path}/BFEE/003_EulerPhi/colvars_{i+1}.in', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx')
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(False, '', '', '../complex.xyz')
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        False, 0, 0, 'eulerTheta', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        True, float(i)/stratification[2] * 20 - 10, float(i+1)/stratification[2] * 20 - 10, 'eulerPhi', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('eulerPhi')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicWallsTemplate(
                        'eulerPhi', float(i)/stratification[2] * 20 - 10, float(i+1)/stratification[2] * 20 - 10
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('RMSD', 10, 0)
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicTemplate('eulerTheta', 0.1, 0)
                )
                colvarsConfig.write(
                    self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
                )

        # 004_Psi
        for i in range(stratification[3]):
            with open(f'{path}/BFEE/004_EulerPsi/colvars_{i+1}.in', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx')
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(False, '', '', '../complex.xyz')
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        False, 0, 0, 'eulerTheta', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        False, 0, 0, 'eulerPhi', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        True, float(i)/stratification[3] * 20 - 10, float(i+1)/stratification[3] * 20 - 10, 'eulerPsi', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('eulerPsi')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicWallsTemplate(
                        'eulerPsi', float(i)/stratification[3] * 20 - 10, float(i+1)/stratification[3] * 20 - 10
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
                    self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
                )

        # 005_polarTheta
        for i in range(stratification[4]):
            with open(f'{path}/BFEE/005_PolarTheta/colvars_{i+1}.in', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx')
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(False, '', '', '../complex.xyz')
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        False, 0, 0, 'eulerTheta', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        False, 0, 0, 'eulerPhi', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        False, 0, 0, 'eulerPsi', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvPolarAngleTemplate(
                        True, float(i)/stratification[4] * 20 - 10 + polarAngles[0], 
                        float(i+1)/stratification[4] * 20 - 10 + polarAngles[0], 'polarTheta', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('polarTheta')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicWallsTemplate(
                        'polarTheta', float(i)/stratification[4] * 20 - 10 + polarAngles[0],
                        float(i+1)/stratification[4] * 20 - 10 + polarAngles[0]
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
                    self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
                )

        # 006_polarPhi
        for i in range(stratification[5]):
            with open(f'{path}/BFEE/006_PolarPhi/colvars_{i+1}.in', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx')
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(False, '', '', '../complex.xyz')
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        False, 0, 0, 'eulerTheta', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        False, 0, 0, 'eulerPhi', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        False, 0, 0, 'eulerPsi', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvPolarAngleTemplate(
                        False, 0, 0, 'polarTheta', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvPolarAngleTemplate(
                        True, float(i)/stratification[5] * 20 - 10 + polarAngles[1], 
                        float(i+1)/stratification[5] * 20 - 10 + polarAngles[1], 'polarPhi', '../complex.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('polarPhi')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicWallsTemplate(
                        'polarPhi', float(i)/stratification[5] * 20 - 10 + polarAngles[1],
                        float(i+1)/stratification[5] * 20 - 10 + polarAngles[1]
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
                    self.cTemplate.cvProteinTemplate(center, '../complex.xyz')
                )

        # 007_r
        # eq
        with open(f'{path}/BFEE/007_r/colvars_eq.in', 'w') as colvarsConfig:
            colvarsConfig.write(
                self.cTemplate.cvHeadTemplate('../complex.ndx')
            )
            colvarsConfig.write(
                self.cTemplate.cvRMSDTemplate(False, '', '', './complex_largeBox.xyz')
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerTheta', './complex_largeBox.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerPhi', './complex_largeBox.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvEulerAngleTemplate(
                    False, 0, 0, 'eulerPsi', './complex_largeBox.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvPolarAngleTemplate(
                    False, 0, 0, 'polarTheta', './complex_largeBox.xyz'
                )
            )
            colvarsConfig.write(
                self.cTemplate.cvPolarAngleTemplate(
                    False, 0, 0, 'polarPhi', './complex_largeBox.xyz'
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
                self.cTemplate.cvProteinTemplate(center, './complex_largeBox.xyz')
            )

        # abf
        for i in range(stratification[6]):
            with open(f'{path}/BFEE/007_r/colvars_{i+1}.in', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('../complex.ndx')
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(False, '', '', './complex_largeBox.xyz')
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        False, 0, 0, 'eulerTheta', './complex_largeBox.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        False, 0, 0, 'eulerPhi', './complex_largeBox.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvEulerAngleTemplate(
                        False, 0, 0, 'eulerPsi', './complex_largeBox.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvPolarAngleTemplate(
                        False, 0, 0, 'polarTheta', './complex_largeBox.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvPolarAngleTemplate(
                        False, 0, 0, 'polarPhi', './complex_largeBox.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvRTemplate(
                        True, float(i)/stratification[6] * 24 - 2 + distance, 
                        float(i+1)/stratification[6] * 24 - 2 + distance
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('r')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicWallsTemplate(
                        'r', float(i)/stratification[6] * 24 - 2 + distance, 
                        float(i+1)/stratification[6] * 24 - 2 + distance
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
                    self.cTemplate.cvProteinTemplate(center, './complex_largeBox.xyz')
                )

        # 008_RMSDUnbound
        for i in range(stratification[7]):
            with open(f'{path}/BFEE/008_RMSDUnbound/colvars_{i+1}.in', 'w') as colvarsConfig:
                colvarsConfig.write(
                    self.cTemplate.cvHeadTemplate('./ligandOnly.ndx')
                )
                colvarsConfig.write(
                    self.cTemplate.cvRMSDTemplate(
                        True, float(i)/stratification[7] * 3.0, float(i+1)/stratification[7] * 3.0, './ligandOnly.xyz'
                    )
                )
                colvarsConfig.write(
                    self.cTemplate.cvABFTemplate('RMSD')
                )
                colvarsConfig.write(
                    self.cTemplate.cvHarmonicWallsTemplate(
                        'RMSD', float(i)/stratification[7] * 3.0, float(i+1)/stratification[7] * 3.0
                    )
                )