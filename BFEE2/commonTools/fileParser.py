# file parser used in bfee

import MDAnalysis
from MDAnalysis import transformations
import numpy as np
import os, sys, math, shutil

# an runtime error
# selection corresponding to nothing
class SelectionError(RuntimeError):
    def __init__(self, arg):
        self.args = arg

class fileParser:
    ''' topology and coordinate parser,
        this class implements basic method for the topology and
        coordinate file, used by BFEE '''

    def __init__(self, topFile, coorFile=''):
        ''' initialize fileParser, on coorFile is not provided, then 
            topFile is considered as a top-coor file, such as mol2
            inputs: 
                topFile (string): path of the topology (psf, parm) file
                coorFile (string): path of the coordinate (pdb, rst) file '''

        coorPostfix = os.path.splitext(coorFile)[-1]
        if coorPostfix == '.rst7' or coorPostfix == '.rst' or coorPostfix == '.inpcrd':
            coorType = 'inpcrd'

        if coorFile == '':
            self.uObject = MDAnalysis.Universe(topFile)
        else:
            # MDAnalysis does not well recognize amber file type by postfix
            coorPostfix = os.path.splitext(coorFile)[-1]
            if coorPostfix == '.rst7' or coorPostfix == '.rst' or coorPostfix == '.inpcrd':
                coorType = 'INPCRD'
                self.uObject = MDAnalysis.Universe(topFile, coorFile, format=coorType)
            else:
                self.uObject = MDAnalysis.Universe(topFile, coorFile)
        self.uObject.add_TopologyAttr('tempfactors')
        self.topPath = topFile

    def saveFile(self, selection, targetPath, targetType, saveTop=False, topPath=''):
        ''' save the coordinate file to the target type 
            this function cannot really 'save' the topology file,
            it can only copy the original topology to the target path
            inputs:
                selection (string): selection of atoms to save
                targetPath (string): path for the coor file to be saved
                targetType (string): type of the coor file (pdb, xyz)
                saveTop (bool): whether the topology file will be saved
                topPath (string): path for the topology file to be saved '''

        atoms = self.uObject.select_atoms(selection)

        if len(atoms) == 0:
            raise SelectionError('Empty selection!')

        atoms.write(targetPath, targetType)
        if saveTop:
            assert(selection == 'all')
            shutil.copyfile(self.topPath, topPath)

    def saveNDX(self, selections, names, targetPath, nonHydrogen=True):
        ''' save an ndx file, including the selections
            inputs:
                selections (list of strings): the selections for the ndx file
                names (list of strings): the name in ndx of each selection
                targetPath (string): path for the file to be saved'''

        assert(len(selections) == len(names))

        if nonHydrogen == True:
            HString = 'and not (name H*)'
        else:
            HString = ''

        for selection, name in zip(selections, names):
            atoms = self.uObject.select_atoms(f'{selection} {HString}')
            
            if len(atoms) == 0:
                raise SelectionError('Empty selection!')
        
            atoms.write(targetPath, 'ndx', name=name, mode='a')

    def getResid(self, selection):
        ''' return a string listing the resid of the selection
            may be used to generate amber masks
            Inputs:
                selection (string): MDAnalysis selection
            Return:
                string: a list of resid, e.g. (4,5,6,7,8,9,10) '''
        atoms = self.uObject.select_atoms(selection)
        
        if len(atoms) == 0:
            raise SelectionError('Empty selection!')
        
        return ','.join([str(num+1) for num in atoms.residues.ix])

    def measureMinmax(self, selection):
        ''' mimic VMD measure minmax
            inputs:
                selection (string): selection of atoms
            return:
                np.array (2*3, float): ((minX, minY, minZ), (maxX, maxY, maxZ)) '''

        atoms = self.uObject.select_atoms(selection)
        
        if len(atoms) == 0:
            raise SelectionError('Empty selection!')
        
        atomPositions = atoms.positions
        xyz_array = np.transpose(atomPositions)
        min_x = np.min(xyz_array[0])
        max_x = np.max(xyz_array[0])
        min_y = np.min(xyz_array[1])
        max_y = np.max(xyz_array[1])
        min_z = np.min(xyz_array[2])
        max_z = np.max(xyz_array[2])
        
        return np.array([[min_x, min_y, min_z],[max_x, max_y, max_z]])

    def measureCenter(self, selection):
        ''' mimic vmd measure center
            inputs:
                selection (string): selection of atoms
            return:
                np.array (3, float): (x, y, z) '''

        atoms = self.uObject.select_atoms(selection)
        
        if len(atoms) == 0:
            raise SelectionError('Empty selection!')
        
        atomPositions = atoms.positions
        xyz_array = np.transpose(atomPositions)
        center_x = np.average(xyz_array[0])
        center_y = np.average(xyz_array[1])
        center_z = np.average(xyz_array[2])

        return np.array([center_x, center_y, center_z])

    def measureDistance(self, selection1, selection2):
        ''' measure the distance between the center of mass
            of two atom groups
            Inputs:
                selection1 (string): selection of atom group 1
                selection2 (string): selection of atom group 2
            Return:
                float: distance between the center of mass of the 
                       two atom groups '''

        center1 = self.measureCenter(selection1)
        center2 = self.measureCenter(selection2)
        return round(np.linalg.norm(center2 - center1), 1)

    def measurePBC(self):
        ''' measure periodic boundary condition of the file 
            return:
                np.array (2*3, float): ((lengthX, lengthY, lengthZ), (centerX, centerY, centerZ)) '''

        minmaxArray = self.measureMinmax('all')
        vec = minmaxArray[1] - minmaxArray[0]
        center = self.measureCenter('all')

        return np.array((vec, center))

    def measurePolarAngles(self, selectionPro, selectionLig):
        ''' calculation the polar angles based on selectionPro and selectionLig
            inputs:
                selectionPro (string): selection of the host molecule
                selectionLig (string): selection of the ligand molecule
            return:
                np.array (2, float): (theta, phi) in degrees'''
        
        vector = self.measureCenter(selectionLig) - self.measureCenter(selectionPro)
        vector /= np.linalg.norm(vector)

        return (float(int(math.degrees(np.arccos(-vector[1])))), float(int(math.degrees(np.arctan2(vector[2], vector[0])))))

    def setBeta(self, selection, beta):
        ''' set beta for the selected atoms
            Inputs:
                selection (string): selection of atoms to change beta
                beta (int): beta value '''

        atoms = self.uObject.select_atoms(selection)
        
        if len(atoms) == 0:
            raise SelectionError('Empty selection!')
        
        atoms.tempfactors = beta

    def moveSystem(self, moveVector):
        ''' move all the atoms in the loaded file
            inputs:
                moveVector (np.array, 3, float): the vector of moving '''

        allAtoms = self.uObject.select_atoms('all')
        transformations.translate(moveVector)(allAtoms)

    def rotateSystem(self, axis, degrees):
        ''' rotate all the atoms in the loaded file
            inputs:
                axis (string): 'x' or 'y' or 'z'
                degrees (float): degrees to move by '''

        assert(axis == 'x' or axis == 'y' or axis == 'z')

        allAtoms = self.uObject.select_atoms('all')
        if axis == 'x':
            axisVector = (1,0,0)
        elif axis == 'y':
            axisVector = (0,1,0)
        else:
            axisVector = (0,0,1)
        transformations.rotate.rotateby(degrees, axisVector, ag=allAtoms)(allAtoms)

    def centerSystem(self):
        ''' move all the atoms in the loaded file, 
            such that the center of system be (0,0,0) '''

        vec = self.measurePBC()[1] * -1.0
        self.moveSystem(vec)