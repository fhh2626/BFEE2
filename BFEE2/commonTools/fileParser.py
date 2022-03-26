# file parser used in bfee

import math
import os
import shutil
import sys

import MDAnalysis
import numpy as np
import parmed
from MDAnalysis import transformations


# an runtime error
# selection corresponding to nothing
class SelectionError(RuntimeError):
    def __init__(self, arg):
        self.args = arg

class fileParser:
    """topology and coordinate parser,
       this class implements basic method for the topology and
       coordinate file, used by BFEE
    """    

    def __init__(self, topFile, coorFile=''):
        """initialize fileParser, on coorFile is not provided, then 
           topFile is considered as a top-coor file, such as mol2

        Args:
            topFile (str): path of the topology (psf, parm) file
            coorFile (str): path of the coordinate (pdb, rst) file. Defaults to ''.
        """        

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
        """save the coordinate file to the target type 
           this function cannot really 'save' the topology file,
           it can only copy the original topology to the target path

        Args:
            selection (str): selection of atoms to save
            targetPath (str): path for the coor file to be saved
            targetType (str): type of the coor file (pdb, xyz)
            saveTop (bool): whether the topology file will be saved. Defaults to False.
            topPath (str, optional): path for the topology file to be saved. Defaults to ''.

        Raises:
            SelectionError: if the selection corresponds to nothing
        """ 

        atoms = self.uObject.select_atoms(selection)

        if len(atoms) == 0:
            raise SelectionError('Empty selection!')

        atoms.write(targetPath, targetType, bonds=None)
        if saveTop:
            assert(selection == 'all')
            shutil.copyfile(self.topPath, topPath)

    def saveNDX(self, selections, names, targetPath, nonHydrogen=True):
        """save an ndx file, including the selections

        Args:
            selections (list of str): the selections for the ndx file
            names (list of str): the name in ndx of each selection
            targetPath (str): path for the file to be saved
            nonHydrogen (bool, optional): whether only non-H atoms are considered. Defaults to True.

        Raises:
            SelectionError: if the selection corresponds to nothing
        """

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
        """return a string listing the resid of the selection
           may be used to generate amber masks

        Args:
            selection (str): MDAnalysis selection

        Raises:
            SelectionError: if the selection corresponds to nothing

        Returns:
            str: a list of resid, e.g. (4,5,6,7,8,9,10)
        """        

        atoms = self.uObject.select_atoms(selection)
        
        if len(atoms) == 0:
            raise SelectionError('Empty selection!')
        
        return ','.join([str(num+1) for num in atoms.residues.ix])

    def measureMinmax(self, selection):
        """mimic VMD measure minmax

        Args:
            selection (str): selection of atoms

        Raises:
            SelectionError: if the selection corresponds to nothing
            
        Returns:
            np.array (2*3, float): ((minX, minY, minZ), (maxX, maxY, maxZ))
        """        

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
        """mimic vmd measure center

        Args:
            selection (str): selection of atoms
            
        Raises:
            SelectionError: if the selection corresponds to nothing

        Returns:
            np.array (3, float): (x, y, z)
        """        

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
        """measure the distance between the center of mass
           of two atom groups

        Args:
            selection1 (str): selection of atom group 1
            selection2 (str): selection of atom group 2

        Returns:
            float: distance between the center of mass of the two atom groups
        """        

        center1 = self.measureCenter(selection1)
        center2 = self.measureCenter(selection2)
        return round(np.linalg.norm(center2 - center1), 1)

    def measurePBC(self):
        """measure periodic boundary condition of the file

        Returns:
            np.array (2*3, float): ((lengthX, lengthY, lengthZ), (centerX, centerY, centerZ))
        """        

        minmaxArray = self.measureMinmax('all')
        vec = minmaxArray[1] - minmaxArray[0]
        center = self.measureCenter('all')

        return np.array((vec, center))

    def measurePolarAngles(self, selectionPro, selectionLig):
        """calculation the polar angles based on selectionPro and selectionLig

        Args:
            selectionPro (str): selection of the host molecule
            selectionLig (str): selection of the ligand molecule

        Returns:
            np.array (2, float): (theta, phi) in degrees
        """        
        
        vector = self.measureCenter(selectionLig) - self.measureCenter(selectionPro)
        vector /= np.linalg.norm(vector)

        return (float(int(math.degrees(np.arccos(-vector[1])))), float(int(math.degrees(np.arctan2(vector[2], vector[0])))))

    def setBeta(self, selection, beta):
        """set beta for the selected atoms

        Args:
            selection (str): selection of atoms to change beta
            beta (int): beta value

        Raises:
            SelectionError: if the selection corresponds to nothing
        """        

        atoms = self.uObject.select_atoms(selection)
        
        if len(atoms) == 0:
            raise SelectionError('Empty selection!')
        
        atoms.tempfactors = beta

    def moveSystem(self, moveVector):
        """move all the atoms in the loaded file

        Args:
            moveVector (np.array, 3, float): the vector of moving
        """        

        allAtoms = self.uObject.select_atoms('all')
        transformations.translate(moveVector)(allAtoms)

    def rotateSystem(self, axis, degrees):
        """rotate all the atoms in the loaded file

        Args:
            axis (str): 'x' or 'y' or 'z'
            degrees (float): degrees to move by
        """        

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
        """move all the atoms in the loaded file, 
           such that the center of system be (0,0,0)
        """

        vec = self.measurePBC()[1] * -1.0
        self.moveSystem(vec)
        
def charmmToGromacs(psfFile, pdbFile, prmFiles, PBC, outputPrefix):
    """convert a set of CHARMM files (psf + pdb + prm + pbc) into the Gromacs format

    Args:
        psfFile (str): path of the psf file
        pdbFile (str): path of the pdb file
        prmFiles (list of str): pathes of the prm files
        PBC (list of flost): pbc information
        outputPrefix (str): path + prefix of the output file
    """
    struct = parmed.load_file(psfFile)
    struct.load_parameters(
        parmed.charmm.CharmmParameterSet(*prmFiles)
    )
    struct.coordinates = parmed.load_file(pdbFile).coordinates
    struct.box = [PBC[0], PBC[1], PBC[2], 90, 90, 90]
    struct.save(f'{outputPrefix}.top', format='gromacs')
    struct.save(f'{outputPrefix}.gro')

def amberToGromacs(parmFile, rstFile, PBC, outputPrefix):
    """convert a set of Amber files (parm7 + rst7) into the Gromacs format

    Args:
        parmFile (str): path of the parm7 file
        rstFile (str): path of the rst7 file
        PBC (list of flost): pbc information
        outputPrefix (str): path + prefix of the output file
    """
    struct = parmed.load_file(parmFile, xyz=rstFile)
    struct.box = [PBC[0], PBC[1], PBC[2], 90, 90, 90]
    struct.save(f'{outputPrefix}.top', format='gromacs')
    struct.save(f'{outputPrefix}.gro')

def gromacsToAmber(topFile, groFile, PBC, outputPrefix):
    """convert a set of Gromacs files (top + gro/pdb) into the Amber format

    Args:
        topFile (str): path of the top file
        groFile (str): path of the gro file
        PBC (list of flost): pbc information
        outputPrefix (str): path + prefix of the output file
    """
    struct = parmed.load_file(topFile, xyz=groFile)
    struct.box = [PBC[0], PBC[1], PBC[2], 90, 90, 90]
    struct.save(f'{outputPrefix}.parm7', format='amber')
    struct.save(f'{outputPrefix}.rst7', format='rst7') 
