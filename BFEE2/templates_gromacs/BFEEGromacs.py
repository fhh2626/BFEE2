#!/usr/bin/env python3
import numpy as np
import os
import sys
import posixpath
import string
import logging
import shutil
from MDAnalysis import Universe
from MDAnalysis.units import convert
from MDAnalysis import transformations
from math import isclose

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from BFEE2 import templates_gromacs

# an runtime error
# selection corresponding to nothing
class SelectionError(RuntimeError):
    def __init__(self, arg):
        self.args = arg

def scanGromacsTopologyInclude(gmxTopFile, logger=None):
    """scan the files included by a gromacs topology file

    Args:
        gmxTopFile (str): filename of the gromacs topology file
        logger (Logging.logger, optional): logger for debugging. Defaults to None.

    Returns:
        tuple: 
            tuple:
                a (list, list) tuple.
                The first list contains the absolute pathnames of included files.
                The second list contains the strings in the quotation marks for handling
                relative paths.
            Logging.logger:
                logger for debugging
    """    

    topology_dirpath = os.path.dirname(os.path.abspath(gmxTopFile))
    include_files = []
    include_strings = []
    with open(gmxTopFile, 'r') as finput:
        for line in finput:
            line = line.strip()
            if line.startswith('#include'):
                # this is an include line
                # find the first quote
                first_quote = line.find('"')
                # find the second quote
                second_quote = line.find('"', first_quote+1)
                # save the string of include
                include_str = line[first_quote+1:second_quote]
                # find the absolute path and save it to the list
                include_filename = os.path.join(topology_dirpath, include_str).replace('\\', '/')
                if (posixpath.exists(include_filename)):
                    include_files.append(include_filename)
                    include_strings.append(include_str)
                else:
                    if logger is None:
                        print(f'Warning: {include_filename} does not exists.')
                    else:
                        logger.warning(f'{include_filename} does not exists.')
    return include_files, include_strings


def measure_minmax(atom_positions):
    """mimic the VMD command "measure minmax"

    Args:
        atom_positions (numpy.array): a numpy array containing the XYZ coordinates of N atoms. The shape 
                                      should be (N, 3).

    Returns:
        Numpy.array: a shape of (2, 3) array, where the first subarray has the minimum 
                     values in XYZ directions, and the second subarray has the maximum
                     values in XYZ directions.
    """    

    xyz_array = np.transpose(atom_positions)
    min_x = np.min(xyz_array[0])
    max_x = np.max(xyz_array[0])
    min_y = np.min(xyz_array[1])
    max_y = np.max(xyz_array[1])
    min_z = np.min(xyz_array[2])
    max_z = np.max(xyz_array[2])
    return np.array([[min_x, min_y, min_z],[max_x, max_y, max_z]])


def measure_center(atom_positions):
    """mimic the VMD command "measure center"

    Args:
        atom_positions (numpy.array): a numpy array containing the XYZ coordinates of N atoms. The shape should be (N, 3).

    Returns:
        Numpy.array: a shape of (3,) array contains the geometric center
    """    

    xyz_array = np.transpose(atom_positions)
    center_x = np.average(xyz_array[0])
    center_y = np.average(xyz_array[1])
    center_z = np.average(xyz_array[2])
    return np.array([center_x, center_y, center_z])


def get_cell(atom_positions):
    """mimic the VMD script get_cell.tcl to calculate the cell units

    Args:
        atom_positions (numpy.array): a numpy array containing the XYZ coordinates of N atoms. The shape should be (N, 3).

    Returns:
         Numpy.array: a shape of (3,3) array contains periodic cell information
    """    

    minmax_array = measure_minmax(atom_positions)
    vec = minmax_array[1] - minmax_array[0]
    cell_basis_vector1 = np.array([vec[0], 0, 0])
    cell_basis_vector2 = np.array([0, vec[1], 0])
    cell_basis_vector3 = np.array([0, 0, vec[2]])
    return np.array([cell_basis_vector1,
                     cell_basis_vector2,
                     cell_basis_vector3])


def generateMDP(MDPTemplate, outputPrefix, timeStep, numSteps, temperature, pressure, logger=None):
    """generate a GROMACS mdp file from a template

    Args:
        MDPTemplate (str): template MDP file with $dt and $nsteps
        outputPrefix (str): prefix (no .mdp extension) of the output MDP file
        timeStep (float): timestep for running the simulation
        numSteps (int): number of steps for running the simulation
        temperature (float): simulation temperature
        pressure (float): simulation pressure
        logger (Logging.logger, optional): logger for debugging. Defaults to None.
    """    

    #if logger is None:
        #print(f'generateMDP: Generating {outputPrefix + ".mdp"} from template {MDPTemplate}...')
        #print(f'Timestep (dt): {timeStep}')
        #print(f'Number of simulation steps (nsteps): {numSteps}')
        #print(f'Temperature: {temperature}')
        #print(f'Pressure: {pressure}')
    #else:
        #logger.info(f'generateMDP: Generating {outputPrefix + ".mdp"} from template {MDPTemplate}...')
        #logger.info(f'Timestep (dt): {timeStep}')
        #logger.info(f'Number of simulation steps (nsteps): {numSteps}')
        #logger.info(f'Temperature: {temperature}')
        #logger.info(f'Pressure: {pressure}')
    #with open(MDPTemplate, 'r', newline='\n') as finput:
    MDP_content = string.Template(MDPTemplate)
    MDP_content = MDP_content.safe_substitute(dt=timeStep,
                                              nsteps=numSteps,
                                              temperature=temperature,
                                              pressure=pressure)
    with open(outputPrefix + '.mdp', 'w', newline='\n') as foutput:
        foutput.write(MDP_content)

def generateColvars(colvarsTemplate, outputPrefix, logger=None, **kwargs):
    """generate a Colvars configuration file from a template suffixed with '.dat'

    Args:
        colvarsTemplate (str): path to a colvars template
        outputPrefix (str): (no .dat extension) of the output Colvars configuration file
        logger (Logging.logger, optional): logger for debugging. Defaults to None.
        **kwargs: additional arguments passed to safe_substitute()
    """    

    #if logger is None:
        #print(f'generateColvars: Generating {outputPrefix + ".dat"} from template {colvarsTemplate}...')
        #print('Colvars parameters:')
    #else:
        #logger.info(f'generateColvars: Generating {outputPrefix + ".dat"} from template {colvarsTemplate}...')
        #logger.info('Colvars parameters:')
    for key, val in kwargs.items():
        if logger is None:
            print(f'{key} = {val}')
        else:
            logger.info(f'{key} = {val}')
    #with open(colvarsTemplate, 'r', newline='\n') as finput:
    content = string.Template(colvarsTemplate)
    content = content.safe_substitute(**kwargs)
    with open(outputPrefix + '.dat', 'w', newline='\n') as foutput:
        foutput.write(content)

def generateShellScript(shellTemplate, outputPrefix, logger=None, **kwargs):
    """generate a shell script from a template

    Args:
        shellTemplate (str): path to a shell template
        outputPrefix (str): prefix (no .sh extension) of the output Colvars configuration file
        logger (Logging.logger, optional): logger for debugging. Defaults to None.
        **kwargs: additional arguments passed to safe_substitute()
    """    

    #if logger is None:
        #print(f'generateShellScript: Generating {outputPrefix + ".sh"} from template {shellTemplate}...')
    #else:
        #logger.info(f'generateShellScript: Generating {outputPrefix + ".sh"} from template {shellTemplate}...')
    #with open(shellTemplate, 'r', newline='\n') as finput:
    content = string.Template(shellTemplate)
    content = content.safe_substitute(**kwargs)
    with open(outputPrefix + '.sh', 'w', newline='\n') as foutput:
        foutput.write(content)


def mearsurePolarAngles(proteinCenter, ligandCenter):
    """measure the polar angles between the protein and the ligand

    Args:
        proteinCenter (numpy.array): center-of-mass of the protein
        ligandCenter (numpy.array): center-of-mass of the ligand

    Returns:
        tuple: a (theta, phi) tuple where theta and phi are measured in degrees
    """    

    vector = ligandCenter - proteinCenter
    vector /= np.linalg.norm(vector)
    return (np.degrees(np.arccos(vector[2])),
            np.degrees(np.arctan2(vector[1], vector[0])))


class BFEEGromacs:
    """The entry class for handling gromacs inputs in BFEE.

    Attributes:
        logger (logging.Logger): logger object for debugging
        handler (logging.StreamHandler): output stream of the debug output
        baseDirectory (str): output directory of the generated files
        structureFile (str): filename of the structure file (either in PDB or GRO format) of the
                             protein-ligand binding complex
        topologyFile (str): filename of the GROMACS topology file of the protein-ligand binding
                            complex
        ligandOnlyStructureFile (str): filename of the structure file (either in PDB or GRO format) of the
                                       ligand-only system
        system (MDAnalysis.core.universe): MDAnalysis universe of the protein-ligand binding system
        ligandOnlySystem (MDAnalysis.core.universe): MDAnalysis universe of the ligand-only system
        basenames (str): subdirectory names of all eight steps
        ligand (MDAnalysis.core.groups.AtomGroup): selected HEAVY ATOMS of the ligand in the protein-ligand binding
                                                   complex. This attribute does not exist until the call of
                                                   setLigandHeavyAtomsGroup.
        ligandOnly (MDAnalysis.core.groups.AtomGroup): selected HEAVY ATOMS of the ligand in the ligand-only system.
                                                       This attribute does not exist until the call of 
                                                       setLigandHeavyAtomsGroup.
        protein (MDAnalysis.core.groups.AtomGroup): selected HEAVY ATOMS of the protein in the protein-ligand binding
                                                    complex. This attribute does not exist until the call of
                                                    setProteinHeavyAtomsGroup.
        solvent (MDAnalysis.core.groups.AtomGroup): selected atoms of the solvents in the protein-ligand binding complex.
                                                    This attribute does not exist until the call of setSolventAtomsGroup.
        temperature (float): the temperature of simulations (default : 300.0)
    """    

    def __init__(self, structureFile, topologyFile, ligandOnlyStructureFile, ligandOnlyTopologyFile, baseDirectory=None, structureFormat='pdb', ligandOnlyStructureFileFormat='pdb'):
        # setup the logger
        self.logger = logging.getLogger()
        self.logger.handlers.clear()
        self.handler = logging.StreamHandler(sys.stdout)
        self.handler.setFormatter(logging.Formatter('%(asctime)s [BFEEGromacs][%(levelname)s]:%(message)s'))
        self.logger.addHandler(self.handler)
        self.logger.setLevel(logging.INFO)
        self.logger.info('Initializing BFEEGromacs...')
        # setup class attributes
        self.structureFile = structureFile
        self.topologyFile = topologyFile
        if baseDirectory is None:
            self.baseDirectory = os.getcwd()
        else:
            self.baseDirectory = baseDirectory
        self.ligandOnlyStructureFile = ligandOnlyStructureFile
        self.ligandOnlyTopologyFile = ligandOnlyTopologyFile

        # start to load data into MDAnalysis
        self.logger.info(f'Calling MDAnalysis to load structure {self.structureFile} (format {structureFormat}).')
        self.system = Universe(self.structureFile, format=structureFormat)
        self.ligandOnlySystem = Universe(self.ligandOnlyStructureFile, format=structureFormat)

        # some PDB files do not have cell info
        # so we reset the cell by an estimation
        all_atoms = self.system.select_atoms("all")
        all_atoms_ligandOnly = self.ligandOnlySystem.select_atoms("all")
        self.system.trajectory[0].triclinic_dimensions = get_cell(all_atoms.positions)
        self.ligandOnlySystem.trajectory[0].triclinic_dimensions = get_cell(all_atoms_ligandOnly.positions)
        dim = self.system.dimensions
        ligandOnly_dim = self.ligandOnlySystem.dimensions
        # measure the cell
        volume = dim[0] * dim[1] * dim[2]
        self.logger.info(f'The volume of the simulation box is {volume} Å^3.')

        # set default temperature to 300.0 K
        self.temperature = 300.0
        self.logger.info(f'You have specified a new base directory at {self.baseDirectory}')
        if not posixpath.exists(self.baseDirectory):
            os.makedirs(self.baseDirectory)
        if not posixpath.exists(posixpath.join(self.baseDirectory, 'Protein')):
            os.makedirs(posixpath.join(self.baseDirectory, 'Protein'))
        if not posixpath.exists(posixpath.join(self.baseDirectory, 'Ligand')):
            os.makedirs(posixpath.join(self.baseDirectory, 'Ligand'))

        # check if the topologies have other itp files included
        topologyIncludeFiles, topologyIncludeStrings = scanGromacsTopologyInclude(self.topologyFile)
        for includeFile, includeString in zip(topologyIncludeFiles, topologyIncludeStrings):
            # handle something like "#include "toppar/xxx.itp""
            # in this case we need to create the directory named "toppar" in the base directory
            dest_dirname = posixpath.dirname(includeString)
            if dest_dirname:
                # if dest_dirname is not empty
                if not posixpath.exists(posixpath.join(self.baseDirectory, 'Protein', dest_dirname)):
                # if the destination directory does not exist
                    os.makedirs(posixpath.join(self.baseDirectory, 'Protein', dest_dirname))
            shutil.copy(includeFile, posixpath.join(self.baseDirectory, 'Protein', dest_dirname))
        # do the same thing to the ligand topology
        topologyIncludeFiles, topologyIncludeStrings = scanGromacsTopologyInclude(self.ligandOnlyTopologyFile)
        for includeFile, includeString in zip(topologyIncludeFiles, topologyIncludeStrings):
            # handle something like "#include "toppar/xxx.itp""
            # in this case we need to create the directory named "toppar" in the base directory
            dest_dirname = posixpath.dirname(includeString)
            if dest_dirname:
                # if dest_dirname is not empty
                if not posixpath.exists(posixpath.join(self.baseDirectory, 'Ligand', dest_dirname)):
                # if the destination directory does not exist
                    os.makedirs(posixpath.join(self.baseDirectory, 'Ligand', dest_dirname))
            shutil.copy(includeFile, posixpath.join(self.baseDirectory, 'Ligand', dest_dirname))

        #self.structureFile = shutil.copy(self.structureFile, self.baseDirectory)
        self.topologyFile = shutil.copy(self.topologyFile, posixpath.join(self.baseDirectory, 'Protein'))
        self.ligandOnlyStructureFile = shutil.copy(self.ligandOnlyStructureFile, posixpath.join(self.baseDirectory, 'Ligand'))
        self.ligandOnlyTopologyFile = shutil.copy(self.ligandOnlyTopologyFile, posixpath.join(self.baseDirectory, 'Ligand'))

        # move the system, so that the complex is at the center of the simulation box
        all_center = measure_center(all_atoms.positions)
        moveVector = (-all_center[0], -all_center[1], -all_center[2])
        transformations.translate(moveVector)(all_atoms)
        all_center = measure_center(all_atoms.positions)
        moveVector = (dim[0]/2, dim[1]/2, dim[2]/2)
        transformations.translate(moveVector)(all_atoms)
        all_center = measure_center(all_atoms.positions)

        ligandOnly_center = measure_center(all_atoms_ligandOnly.positions)
        moveVector = (-ligandOnly_center[0], -ligandOnly_center[1], -ligandOnly_center[2])
        transformations.translate(moveVector)(all_atoms_ligandOnly)
        ligandOnly_center = measure_center(all_atoms_ligandOnly.positions)
        moveVector = (ligandOnly_dim[0]/2, ligandOnly_dim[1]/2, ligandOnly_dim[2]/2)
        transformations.translate(moveVector)(all_atoms_ligandOnly)
        ligandOnly_center = measure_center(all_atoms_ligandOnly.positions)

        _, fileName = os.path.split(self.structureFile)
        all_atoms.write(self.baseDirectory + '/Protein/' + fileName)

        _, ligandOnly_fileName = os.path.split(self.ligandOnlyStructureFile)
        all_atoms_ligandOnly.write(self.baseDirectory + '/' + ligandOnly_fileName)

        #if isclose(volume, 0.0):
            #self.logger.warning(f'The volume is too small. Maybe the structure file is a PDB file without the unit cell.')
            #all_atoms = self.system.select_atoms("all")
            #self.logger.warning(f'The unit cell has been reset to {dim[0]:12.5f} {dim[1]:12.5f} {dim[2]:12.5f} .')
        newBasename = posixpath.splitext(fileName)[0]
        self.structureFile = self.baseDirectory + '/Protein/' + fileName + '.new.gro'
            #self.saveStructure(self.structureFile)
        all_atoms.write(self.structureFile)
        # measure the cell of the ligand-only system
        #dim = self.ligandOnlySystem.dimensions
        #volume = dim[0] * dim[1] * dim[2]
        #self.logger.info(f'The volume of the simulation box (ligand-only system) is {volume} Å^3.')
        #if isclose(volume, 0.0):
            #self.logger.warning(f'The volume is too small. Maybe the structure file is a PDB file without the unit cell.')
            #all_atoms = self.ligandOnlySystem.select_atoms("all")
            #self.ligandOnlySystem.trajectory[0].triclinic_dimensions = get_cell(all_atoms.positions)
            #dim = self.ligandOnlySystem.dimensions
            #self.logger.warning(f'The unit cell has been reset to {dim[0]:12.5f} {dim[1]:12.5f} {dim[2]:12.5f} .')
        newBasename = posixpath.splitext(ligandOnly_fileName)[0]
        self.ligandOnlyStructureFile = self.baseDirectory + '/' + ligandOnly_fileName + '.new.gro'
            #self.saveStructure(self.ligandOnlyStructureFile)
        all_atoms_ligandOnly.write(self.ligandOnlyStructureFile)

        self.basenames = [
                          '000_eq',
                          '001_RMSD_bound',
                          '002_euler_theta',
                          '003_euler_phi',
                          '004_euler_psi',
                          '005_polar_theta',
                          '006_polar_phi',
                          '007_r',
                          '008_RMSD_unbound'
                        ]
        self.stepnames = self.basenames.copy()
        self.basenames = [posixpath.join(self.baseDirectory, basename) for basename in self.basenames]
        self.logger.info('Initialization done.')

    def saveStructure(self, outputFile, selection='all'):
        """a helper method for selecting a group of atoms and save it

        Args:
            outputFile (str): filename of the output file
            selection (str, optional): MDAnalysis atom selection string. Defaults to 'all'.

        Raises:
            SelectionError: if the selection corresponds to nothing
        """        

        self.logger.info(f'Saving a new structure file at {outputFile} with selection ({selection}).')
        selected_atoms = self.system.select_atoms(selection)
        if len(selected_atoms) == 0:
            raise SelectionError('Empty selection!')
        selected_atoms.write(outputFile)

    def setProteinHeavyAtomsGroup(self, selection):
        """select the heavy atoms of the protein

        Args:
            selection (str): MDAnalysis atom selection string

        Raises:
            SelectionError: if the selection corresponds to nothing
        """        

        self.logger.info(f'Setup the atoms group of the protein by selection: {selection}')
        self.protein = self.system.select_atoms(selection)
        if len(self.protein) == 0:
            raise SelectionError('Empty selection!')
    
    def setLigandHeavyAtomsGroup(self, selection):
        """select the heavy atoms of the ligand in both the protein-ligand complex
           and the ligand-only systems

        Args:
            selection (str): MDAnalysis atom selection string

        Raises:
            SelectionError: if the selection corresponds to nothing
        """        

        self.logger.info(f'Setup the atoms group of the ligand by selection: {selection}')
        self.ligand = self.system.select_atoms(selection)
        if len(self.ligand) == 0:
            raise SelectionError('Empty selection!')
        self.ligandOnly = self.ligandOnlySystem.select_atoms(selection)
        if len(self.ligandOnly) == 0:
            raise SelectionError('Empty selection!')

    def setSolventAtomsGroup(self, selection):
        """select the solvent atoms

        Args:
            selection (str): MDAnalysis atom selection string

        Raises:
            SelectionError: if the selection corresponds nothing
        """        
        
        self.logger.info(f'Setup the atoms group of the solvent molecule by selection: {selection}')
        self.solvent = self.system.select_atoms(selection)
        if len(self.solvent) == 0:
            raise SelectionError('Empty selection!')

    def setTemperature(self, newTemperature):
        """set the temperature

        Args:
            newTemperature (float): new value of the temperature
        """        
        self.temperature = newTemperature

    def generateGromacsIndex(self, outputFile):
        """generate a GROMACS index file for atom selection in Colvars
        """

        self.system.select_atoms('all').write(outputFile, name='BFEE_all')
        if hasattr(self, 'ligand'):
            self.ligand.write(outputFile, name='BFEE_Ligand', mode='a')
        if hasattr(self, 'protein'):
            self.protein.write(outputFile, name='BFEE_Protein', mode='a')
        if hasattr(self, 'solvent'):
            self.solvent.write(outputFile, name='BFEE_Solvent', mode='a')

    def generate000(self):
        """generate files for running an equilibrium simulation
        """

        self.handler.setFormatter(logging.Formatter('%(asctime)s [BFEEGromacs][000][%(levelname)s]:%(message)s'))
        generate_basename = self.basenames[0]
        self.logger.info('=' * 80)
        self.logger.info(f'Generating simulation files for {generate_basename}...')
        if not posixpath.exists(generate_basename):
            self.logger.info(f'Making directory {os.path.abspath(generate_basename)}...')
            os.makedirs(generate_basename)
        # generate the MDP file
        generateMDP(pkg_resources.read_text(templates_gromacs, '000.mdp.template'),
                    posixpath.join(generate_basename, '000_eq'),
                    logger=self.logger,
                    timeStep=0.002,
                    numSteps=4000000,
                    temperature=self.temperature,
                    pressure=1.01325)
        # check if the ligand and protein is selected
        if not hasattr(self, 'ligand'):
            raise RuntimeError('The atoms of the ligand have not been selected.')
        if not hasattr(self, 'protein'):
            raise RuntimeError('The atoms of the protein have not been selected.')
        # measure the COM of the protein
        protein_center = measure_center(self.protein.positions)
        # convert angstrom to nanometer and format the string
        protein_center = convert(protein_center, "angstrom", "nm")
        protein_center_str = f'({protein_center[0]}, {protein_center[1]}, {protein_center[2]})'
        self.logger.info('COM of the protein: ' + protein_center_str + '.')
        # generate the index file
        self.generateGromacsIndex(posixpath.join(generate_basename, 'colvars.ndx'))
        # generate the colvars configuration
        colvars_inputfile_basename = posixpath.join(generate_basename, '000_colvars')
        generateColvars(pkg_resources.read_text(templates_gromacs, '000.colvars.template'),
                        colvars_inputfile_basename,
                        protein_selection='BFEE_Protein',
                        protein_center=protein_center_str,
                        logger=self.logger)
        # generate the reference file
        self.system.select_atoms('all').write(posixpath.join(generate_basename, 'reference.xyz'))
        # generate the shell script for making the tpr file
        generateShellScript(pkg_resources.read_text(templates_gromacs, '000.generate_tpr_sh.template'),
                            posixpath.join(generate_basename, '000_generate_tpr'),
                            logger=self.logger,
                            MDP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                '000_eq.mdp')),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            GRO_FILE_TEMPLATE=os.path.relpath(os.path.abspath(self.structureFile),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            TOP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(self.topologyFile),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            COLVARS_INPUT_TEMPLATE=os.path.relpath(os.path.abspath(colvars_inputfile_basename + '.dat'),
                                                                    os.path.abspath(generate_basename)).replace('\\', '/'))
        if not posixpath.exists(posixpath.join(generate_basename, 'output')):
            os.makedirs(posixpath.join(generate_basename, 'output'))
        self.logger.info(f"Generation of {generate_basename} done.")
        self.logger.info('=' * 80)


    def generate001(self):
        """generate files for determining the PMF along the RMSD of the ligand  
           with respect to its bound state
        """

        self.handler.setFormatter(logging.Formatter('%(asctime)s [BFEEGromacs][001][%(levelname)s]:%(message)s'))
        generate_basename = self.basenames[1]
        self.logger.info('=' * 80)
        self.logger.info(f'Generating simulation files for {generate_basename}...')
        if not posixpath.exists(generate_basename):
            self.logger.info(f'Making directory {os.path.abspath(generate_basename)}...')
            os.makedirs(generate_basename)
        # generate the MDP file
        generateMDP(pkg_resources.read_text(templates_gromacs, '001.mdp.template'),
                    posixpath.join(generate_basename, '001_PMF'),
                    logger=self.logger,
                    timeStep=0.002,
                    numSteps=4000000,
                    temperature=self.temperature,
                    pressure=1.01325)
        # check if the ligand and protein is selected
        if not hasattr(self, 'ligand'):
            raise RuntimeError('The atoms of the ligand have not been selected.')
        if not hasattr(self, 'protein'):
            raise RuntimeError('The atoms of the protein have not been selected.')
        # measure the COM of the protein
        protein_center = measure_center(self.protein.positions)
        # convert angstrom to nanometer and format the string
        protein_center = convert(protein_center, "angstrom", "nm")
        protein_center_str = f'({protein_center[0]}, {protein_center[1]}, {protein_center[2]})'
        self.logger.info('COM of the protein: ' + protein_center_str + '.')
        # generate the index file
        self.generateGromacsIndex(posixpath.join(generate_basename, 'colvars.ndx'))
        # generate the colvars configuration
        colvars_inputfile_basename = posixpath.join(generate_basename, '001_colvars')
        generateColvars(pkg_resources.read_text(templates_gromacs, '001.colvars.template'),
                        colvars_inputfile_basename,
                        rmsd_bin_width=0.005,
                        rmsd_lower_boundary=0.0,
                        rmsd_upper_boundary=0.5,
                        rmsd_wall_constant=0.8368,
                        ligand_selection='BFEE_Ligand',
                        protein_selection='BFEE_Protein',
                        protein_center=protein_center_str,
                        logger=self.logger)
        # generate the reference file
        self.system.select_atoms('all').write(posixpath.join(generate_basename, 'reference.xyz'))
        # generate the shell script for making the tpr file
        generateShellScript(pkg_resources.read_text(templates_gromacs, '001.generate_tpr_sh.template'),
                            posixpath.join(generate_basename, '001_generate_tpr'),
                            logger=self.logger,
                            MDP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                 '001_PMF.mdp')),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            GRO_FILE_TEMPLATE=os.path.relpath(os.path.abspath(f'{self.baseDirectory}/000_eq/output/000_eq.out.gro'),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            TOP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(self.topologyFile),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            COLVARS_INPUT_TEMPLATE=os.path.relpath(os.path.abspath(colvars_inputfile_basename + '.dat'),
                                                                     os.path.abspath(generate_basename)).replace('\\', '/'))
        if not posixpath.exists(posixpath.join(generate_basename, 'output')):
            os.makedirs(posixpath.join(generate_basename, 'output'))
        self.logger.info(f"Generation of {generate_basename} done.")
        self.logger.info('=' * 80)
    
    def generate002(self):
        """generate files for determining the PMF along the pitch (theta) angle of
           the ligand
        """

        self.handler.setFormatter(logging.Formatter('%(asctime)s [BFEEGromacs][002][%(levelname)s]:%(message)s'))
        generate_basename = self.basenames[2]
        self.logger.info('=' * 80)
        self.logger.info(f'Generating simulation files for {generate_basename}...')
        if not posixpath.exists(generate_basename):
            self.logger.info(f'Making directory {os.path.abspath(generate_basename)}...')
            os.makedirs(generate_basename)
        # generate the MDP file
        generateMDP(pkg_resources.read_text(templates_gromacs, '002.mdp.template'),
                    posixpath.join(generate_basename, '002_PMF'),
                    timeStep=0.002,
                    numSteps=4000000,
                    temperature=self.temperature,
                    pressure=1.01325,
                    logger=self.logger)
        # check if the ligand and protein is selected
        if not hasattr(self, 'ligand'):
            raise RuntimeError('The atoms of the ligand have not been selected.')
        if not hasattr(self, 'protein'):
            raise RuntimeError('The atoms of the protein have not been selected.')
        # measure the COM of the protein
        protein_center = measure_center(self.protein.positions)
        # convert angstrom to nanometer and format the string
        protein_center = convert(protein_center, "angstrom", "nm")
        protein_center_str = f'({protein_center[0]}, {protein_center[1]}, {protein_center[2]})'
        self.logger.info('COM of the protein: ' + protein_center_str + '.')
        # generate the index file
        self.generateGromacsIndex(posixpath.join(generate_basename, 'colvars.ndx'))
        # generate the colvars configuration
        colvars_inputfile_basename = posixpath.join(generate_basename, '002_colvars')
        generateColvars(pkg_resources.read_text(templates_gromacs, '002.colvars.template'),
                        colvars_inputfile_basename,
                        logger=self.logger,
                        eulerTheta_width=1,
                        eulerTheta_lower_boundary=-10.0,
                        eulerTheta_upper_boundary=10.0,
                        eulerTheta_wall_constant=0.8368,
                        ligand_selection='BFEE_Ligand',
                        protein_selection='BFEE_Protein',
                        protein_center=protein_center_str)
        # generate the reference file
        self.system.select_atoms('all').write(posixpath.join(generate_basename, 'reference.xyz'))
        # generate the shell script for making the tpr file
        generateShellScript(pkg_resources.read_text(templates_gromacs, '002.generate_tpr_sh.template'),
                            posixpath.join(generate_basename, '002_generate_tpr'),
                            logger=self.logger,
                            MDP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                 '002_PMF.mdp')),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            GRO_FILE_TEMPLATE=os.path.relpath(os.path.abspath(f'{self.baseDirectory}/000_eq/output/000_eq.out.gro'),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            TOP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(self.topologyFile),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            COLVARS_INPUT_TEMPLATE=os.path.relpath(os.path.abspath(colvars_inputfile_basename + '.dat'),
                                                                     os.path.abspath(generate_basename)).replace('\\', '/'))
        # also copy the awk script to modify the colvars configuration according to the PMF minima in previous stages
        with pkg_resources.path(templates_gromacs, 'find_min_max.awk') as p:
            shutil.copyfile(p, posixpath.join(generate_basename, 'find_min_max.awk'))
        if not posixpath.exists(posixpath.join(generate_basename, 'output')):
            os.makedirs(posixpath.join(generate_basename, 'output'))
        self.logger.info(f"Generation of {generate_basename} done.")
        self.logger.info('=' * 80)

    def generate003(self):
        """generate files for determining the PMF along the roll (phi) angle of
           the ligand
        """

        self.handler.setFormatter(logging.Formatter('%(asctime)s [BFEEGromacs][003][%(levelname)s]:%(message)s'))
        generate_basename = self.basenames[3]
        self.logger.info('=' * 80)
        self.logger.info(f'Generating simulation files for {generate_basename}...')
        if not posixpath.exists(generate_basename):
            self.logger.info(f'Making directory {os.path.abspath(generate_basename)}...')
            os.makedirs(generate_basename)
        # generate the MDP file
        generateMDP(pkg_resources.read_text(templates_gromacs, '003.mdp.template'),
                    posixpath.join(generate_basename, '003_PMF'),
                    timeStep=0.002,
                    numSteps=4000000,
                    temperature=self.temperature,
                    pressure=1.01325,
                    logger=self.logger)
        # check if the ligand and protein is selected
        if not hasattr(self, 'ligand'):
            raise RuntimeError('The atoms of the ligand have not been selected.')
        if not hasattr(self, 'protein'):
            raise RuntimeError('The atoms of the protein have not been selected.')
        # measure the COM of the protein
        protein_center = measure_center(self.protein.positions)
        # convert angstrom to nanometer and format the string
        protein_center = convert(protein_center, "angstrom", "nm")
        protein_center_str = f'({protein_center[0]}, {protein_center[1]}, {protein_center[2]})'
        self.logger.info('COM of the protein: ' + protein_center_str + '.')
        # generate the index file
        self.generateGromacsIndex(posixpath.join(generate_basename, 'colvars.ndx'))
        # generate the colvars configuration
        colvars_inputfile_basename = posixpath.join(generate_basename, '003_colvars')
        generateColvars(pkg_resources.read_text(templates_gromacs, '003.colvars.template'),
                        colvars_inputfile_basename,
                        logger=self.logger,
                        eulerPhi_width=1,
                        eulerPhi_lower_boundary=-10.0,
                        eulerPhi_upper_boundary=10.0,
                        eulerPhi_wall_constant=0.8368,
                        ligand_selection='BFEE_Ligand',
                        protein_selection='BFEE_Protein',
                        protein_center=protein_center_str)
        # generate the reference file
        self.system.select_atoms('all').write(posixpath.join(generate_basename, 'reference.xyz'))
        # generate the shell script for making the tpr file
        generateShellScript(pkg_resources.read_text(templates_gromacs, '003.generate_tpr_sh.template'),
                            posixpath.join(generate_basename, '003_generate_tpr'),
                            logger=self.logger,
                            BASENAME_002=self.stepnames[2],
                            MDP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                 '003_PMF.mdp')),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            GRO_FILE_TEMPLATE=os.path.relpath(os.path.abspath(f'{self.baseDirectory}/000_eq/output/000_eq.out.gro'),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            TOP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(self.topologyFile),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            COLVARS_INPUT_TEMPLATE=os.path.relpath(os.path.abspath(colvars_inputfile_basename + '.dat'),
                                                                     os.path.abspath(generate_basename)).replace('\\', '/'))
        # also copy the awk script to modify the colvars configuration according to the PMF minima in previous stages
        with pkg_resources.path(templates_gromacs, 'find_min_max.awk') as p:
            shutil.copyfile(p, posixpath.join(generate_basename, 'find_min_max.awk'))
        if not posixpath.exists(posixpath.join(generate_basename, 'output')):
            os.makedirs(posixpath.join(generate_basename, 'output'))
        self.logger.info(f"Generation of {generate_basename} done.")
        self.logger.info('=' * 80)

    def generate004(self):
        """generate files for determining the PMF along the yaw (psi) angle of
           the ligand
        """

        self.handler.setFormatter(logging.Formatter('%(asctime)s [BFEEGromacs][004][%(levelname)s]:%(message)s'))
        generate_basename = self.basenames[4]
        self.logger.info('=' * 80)
        self.logger.info(f'Generating simulation files for {generate_basename}...')
        if not posixpath.exists(generate_basename):
            self.logger.info(f'Making directory {os.path.abspath(generate_basename)}...')
            os.makedirs(generate_basename)
        # generate the MDP file
        generateMDP(pkg_resources.read_text(templates_gromacs, '004.mdp.template'),
                    posixpath.join(generate_basename, '004_PMF'),
                    timeStep=0.002,
                    numSteps=4000000,
                    temperature=self.temperature,
                    pressure=1.01325,
                    logger=self.logger)
        # check if the ligand and protein is selected
        if not hasattr(self, 'ligand'):
            raise RuntimeError('The atoms of the ligand have not been selected.')
        if not hasattr(self, 'protein'):
            raise RuntimeError('The atoms of the protein have not been selected.')
        # measure the COM of the protein
        protein_center = measure_center(self.protein.positions)
        # convert angstrom to nanometer and format the string
        protein_center = convert(protein_center, "angstrom", "nm")
        protein_center_str = f'({protein_center[0]}, {protein_center[1]}, {protein_center[2]})'
        self.logger.info('COM of the protein: ' + protein_center_str + '.')
        # generate the index file
        self.generateGromacsIndex(posixpath.join(generate_basename, 'colvars.ndx'))
        # generate the colvars configuration
        colvars_inputfile_basename = posixpath.join(generate_basename, '004_colvars')
        generateColvars(pkg_resources.read_text(templates_gromacs, '004.colvars.template'),
                        colvars_inputfile_basename,
                        logger=self.logger,
                        eulerPsi_width=1,
                        eulerPsi_lower_boundary=-10.0,
                        eulerPsi_upper_boundary=10.0,
                        eulerPsi_wall_constant=0.8368,
                        ligand_selection='BFEE_Ligand',
                        protein_selection='BFEE_Protein',
                        protein_center=protein_center_str)
        # generate the reference file
        self.system.select_atoms('all').write(posixpath.join(generate_basename, 'reference.xyz'))
        # generate the shell script for making the tpr file
        generateShellScript(pkg_resources.read_text(templates_gromacs, '004.generate_tpr_sh.template'),
                            posixpath.join(generate_basename, '004_generate_tpr'),
                            logger=self.logger,
                            BASENAME_002=self.stepnames[2],
                            BASENAME_003=self.stepnames[3],
                            MDP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                 '004_PMF.mdp')),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            GRO_FILE_TEMPLATE=os.path.relpath(os.path.abspath(f'{self.baseDirectory}/000_eq/output/000_eq.out.gro'),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            TOP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(self.topologyFile),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            COLVARS_INPUT_TEMPLATE=os.path.relpath(os.path.abspath(colvars_inputfile_basename + '.dat'),
                                                                     os.path.abspath(generate_basename)).replace('\\', '/'))
        # also copy the awk script to modify the colvars configuration according to the PMF minima in previous stages
        with pkg_resources.path(templates_gromacs, 'find_min_max.awk') as p:
            shutil.copyfile(p, posixpath.join(generate_basename, 'find_min_max.awk'))
        if not posixpath.exists(posixpath.join(generate_basename, 'output')):
            os.makedirs(posixpath.join(generate_basename, 'output'))
        self.logger.info(f"Generation of {generate_basename} done.")
        self.logger.info('=' * 80)

    def generate005(self):
        """generate files for determining the PMF along the polar theta angle of
           the ligand relative to the protein
        """

        self.handler.setFormatter(logging.Formatter('%(asctime)s [BFEEGromacs][005][%(levelname)s]:%(message)s'))
        generate_basename = self.basenames[5]
        self.logger.info('=' * 80)
        self.logger.info(f'Generating simulation files for {generate_basename}...')
        if not posixpath.exists(generate_basename):
            self.logger.info(f'Making directory {os.path.abspath(generate_basename)}...')
            os.makedirs(generate_basename)
        # generate the MDP file
        generateMDP(pkg_resources.read_text(templates_gromacs, '005.mdp.template'),
                    posixpath.join(generate_basename, '005_PMF'),
                    timeStep=0.002,
                    numSteps=4000000,
                    temperature=self.temperature,
                    pressure=1.01325,
                    logger=self.logger)
        # check if the ligand and protein is selected
        if not hasattr(self, 'ligand'):
            raise RuntimeError('The atoms of the ligand have not been selected.')
        if not hasattr(self, 'protein'):
            raise RuntimeError('The atoms of the protein have not been selected.')
        # measure the COM of the protein
        protein_center = measure_center(self.protein.positions)
        # convert angstrom to nanometer and format the string
        protein_center = convert(protein_center, "angstrom", "nm")
        protein_center_str = f'({protein_center[0]}, {protein_center[1]}, {protein_center[2]})'
        self.logger.info('COM of the protein: ' + protein_center_str + '.')
        # generate the index file
        self.generateGromacsIndex(posixpath.join(generate_basename, 'colvars.ndx'))
        # generate the colvars configuration
        colvars_inputfile_basename = posixpath.join(generate_basename, '005_colvars')
        # measure the current polar theta angles
        ligand_center = measure_center(self.ligand.positions)
        ligand_center = convert(ligand_center, "angstrom", "nm")
        ligand_center_str = f'({ligand_center[0]}, {ligand_center[1]}, {ligand_center[2]})'
        polar_theta, polar_phi = mearsurePolarAngles(protein_center, ligand_center)
        polar_theta_center = np.around(polar_theta, 1)
        self.logger.info(f'Measured polar angles: theta = {polar_theta:12.5f} ; phi = {polar_phi:12.5f}')
        polar_theta_width = 1
        polar_theta_lower = polar_theta_center - polar_theta_width * np.ceil(10 / polar_theta_width)
        polar_theta_upper = polar_theta_center + polar_theta_width * np.ceil(10 / polar_theta_width)
        generateColvars(pkg_resources.read_text(templates_gromacs, '005.colvars.template'),
                        colvars_inputfile_basename,
                        logger=self.logger,
                        polarTheta_width=polar_theta_width,
                        polarTheta_lower_boundary=np.around(polar_theta_lower, 2),
                        polarTheta_upper_boundary=np.around(polar_theta_upper, 2),
                        polarTheta_wall_constant=0.8368,
                        ligand_selection='BFEE_Ligand',
                        protein_selection='BFEE_Protein',
                        protein_center=protein_center_str)
        # generate the reference file
        self.system.select_atoms('all').write(posixpath.join(generate_basename, 'reference.xyz'))
        # generate the shell script for making the tpr file
        generateShellScript(pkg_resources.read_text(templates_gromacs, '005.generate_tpr_sh.template'),
                            posixpath.join(generate_basename, '005_generate_tpr'),
                            logger=self.logger,
                            BASENAME_002=self.stepnames[2],
                            BASENAME_003=self.stepnames[3],
                            BASENAME_004=self.stepnames[4],
                            MDP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                 '005_PMF.mdp')),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            GRO_FILE_TEMPLATE=os.path.relpath(os.path.abspath(f'{self.baseDirectory}/000_eq/output/000_eq.out.gro'),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            TOP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(self.topologyFile),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            COLVARS_INPUT_TEMPLATE=os.path.relpath(os.path.abspath(colvars_inputfile_basename + '.dat'),
                                                                     os.path.abspath(generate_basename)).replace('\\', '/'))
        # also copy the awk script to modify the colvars configuration according to the PMF minima in previous stages
        with pkg_resources.path(templates_gromacs, 'find_min_max.awk') as p:
            shutil.copyfile(p, posixpath.join(generate_basename, 'find_min_max.awk'))
        if not posixpath.exists(posixpath.join(generate_basename, 'output')):
            os.makedirs(posixpath.join(generate_basename, 'output'))
        self.logger.info(f"Generation of {generate_basename} done.")
        self.logger.info('=' * 80)

    def generate006(self):
        """generate files for determining the PMF along the polar phi angle of
           the ligand relative to the protein
        """

        self.handler.setFormatter(logging.Formatter('%(asctime)s [BFEEGromacs][006][%(levelname)s]:%(message)s'))
        generate_basename = self.basenames[6]
        self.logger.info('=' * 80)
        self.logger.info(f'Generating simulation files for {generate_basename}...')
        if not posixpath.exists(generate_basename):
            self.logger.info(f'Making directory {os.path.abspath(generate_basename)}...')
            os.makedirs(generate_basename)
        # generate the MDP file
        generateMDP(pkg_resources.read_text(templates_gromacs, '006.mdp.template'),
                    posixpath.join(generate_basename, '006_PMF'),
                    timeStep=0.002,
                    numSteps=4000000,
                    temperature=self.temperature,
                    pressure=1.01325,
                    logger=self.logger)
        # check if the ligand and protein is selected
        if not hasattr(self, 'ligand'):
            raise RuntimeError('The atoms of the ligand have not been selected.')
        if not hasattr(self, 'protein'):
            raise RuntimeError('The atoms of the protein have not been selected.')
        # measure the COM of the protein
        protein_center = measure_center(self.protein.positions)
        # convert angstrom to nanometer and format the string
        protein_center = convert(protein_center, "angstrom", "nm")
        protein_center_str = f'({protein_center[0]}, {protein_center[1]}, {protein_center[2]})'
        self.logger.info('COM of the protein: ' + protein_center_str + '.')
        # generate the index file
        self.generateGromacsIndex(posixpath.join(generate_basename, 'colvars.ndx'))
        # generate the colvars configuration
        colvars_inputfile_basename = posixpath.join(generate_basename, '006_colvars')
        # measure the current polar theta angles
        ligand_center = measure_center(self.ligand.positions)
        ligand_center = convert(ligand_center, "angstrom", "nm")
        ligand_center_str = f'({ligand_center[0]}, {ligand_center[1]}, {ligand_center[2]})'
        polar_theta, polar_phi = mearsurePolarAngles(protein_center, ligand_center)
        polar_phi_center = np.around(polar_phi, 1)
        self.logger.info(f'Measured polar angles: theta = {polar_theta:12.5f} ; phi = {polar_phi:12.5f}')
        polar_phi_width = 1
        polar_phi_lower = polar_phi_center - polar_phi_width * np.ceil(10 / polar_phi_width)
        polar_phi_upper = polar_phi_center + polar_phi_width * np.ceil(10 / polar_phi_width)
        generateColvars(pkg_resources.read_text(templates_gromacs, '006.colvars.template'),
                        colvars_inputfile_basename,
                        logger=self.logger,
                        polarPhi_width=polar_phi_width,
                        polarPhi_lower_boundary=np.around(polar_phi_lower, 2),
                        polarPhi_upper_boundary=np.around(polar_phi_upper, 2),
                        polarPhi_wall_constant=0.8368,
                        ligand_selection='BFEE_Ligand',
                        protein_selection='BFEE_Protein',
                        protein_center=protein_center_str)
        # generate the reference file
        self.system.select_atoms('all').write(posixpath.join(generate_basename, 'reference.xyz'))
        # generate the shell script for making the tpr file
        generateShellScript(pkg_resources.read_text(templates_gromacs, '006.generate_tpr_sh.template'),
                            posixpath.join(generate_basename, '006_generate_tpr'),
                            logger=self.logger,
                            BASENAME_002=self.stepnames[2],
                            BASENAME_003=self.stepnames[3],
                            BASENAME_004=self.stepnames[4],
                            BASENAME_005=self.stepnames[5],
                            MDP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                 '006_PMF.mdp')),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            GRO_FILE_TEMPLATE=os.path.relpath(os.path.abspath(f'{self.baseDirectory}/000_eq/output/000_eq.out.gro'),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            TOP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(self.topologyFile),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            COLVARS_INPUT_TEMPLATE=os.path.relpath(os.path.abspath(colvars_inputfile_basename + '.dat'),
                                                                     os.path.abspath(generate_basename)).replace('\\', '/'))
        # also copy the awk script to modify the colvars configuration according to the PMF minima in previous stages
        with pkg_resources.path(templates_gromacs, 'find_min_max.awk') as p:
            shutil.copyfile(p, posixpath.join(generate_basename, 'find_min_max.awk'))
        if not posixpath.exists(posixpath.join(generate_basename, 'output')):
            os.makedirs(posixpath.join(generate_basename, 'output'))
        self.logger.info(f"Generation of {generate_basename} done.")
        self.logger.info('=' * 80)

    def generate007(self):
        """generate files for determining the PMF along the distance between the
           the ligand and the protein
        """

        self.handler.setFormatter(logging.Formatter('%(asctime)s [BFEEGromacs][007][%(levelname)s]:%(message)s'))
        generate_basename = self.basenames[7]
        self.logger.info('=' * 80)
        self.logger.info(f'Generating simulation files for {generate_basename}...')
        if not posixpath.exists(generate_basename):
            self.logger.info(f'Making directory {os.path.abspath(generate_basename)}...')
            os.makedirs(generate_basename)
        # check if the topologies have other itp files included
        topologyIncludeFiles, topologyIncludeStrings = scanGromacsTopologyInclude(self.topologyFile)
        for includeFile, includeString in zip(topologyIncludeFiles, topologyIncludeStrings):
            # handle something like "#include "toppar/xxx.itp""
            # in this case we need to create the directory named "toppar" in the base directory
            dest_dirname = posixpath.dirname(includeString)
            if dest_dirname:
                # if dest_dirname is not empty
                if not posixpath.exists(posixpath.join(self.baseDirectory, generate_basename, dest_dirname)):
                # if the destination directory does not exist
                    os.makedirs(posixpath.join(self.baseDirectory, generate_basename, dest_dirname))
            shutil.copy(includeFile, posixpath.join(self.baseDirectory, generate_basename, dest_dirname))
        # generate the MDP file
        generateMDP(pkg_resources.read_text(templates_gromacs, '007_min.mdp.template'),
                    posixpath.join(generate_basename, '007_Minimize'),
                    timeStep=0.002,
                    numSteps=5000,
                    temperature=self.temperature,
                    pressure=1.01325,
                    logger=self.logger)
        # equilibration
        generateMDP(pkg_resources.read_text(templates_gromacs, '007.mdp.template'),
                    posixpath.join(generate_basename, '007_Equilibration'),
                    timeStep=0.002,
                    numSteps=5000000,
                    temperature=self.temperature,
                    pressure=1.01325,
                    logger=self.logger)
        # free-energy calculation
        generateMDP(pkg_resources.read_text(templates_gromacs, '007.mdp.template'),
                    posixpath.join(generate_basename, '007_PMF'),
                    timeStep=0.002,
                    numSteps=80000000,
                    temperature=self.temperature,
                    pressure=1.01325,
                    logger=self.logger)
        # check if the ligand, protein and solvent is selected
        if not hasattr(self, 'ligand'):
            raise RuntimeError('The atoms of the ligand have not been selected.')
        if not hasattr(self, 'protein'):
            raise RuntimeError('The atoms of the protein have not been selected.')
        if not hasattr(self, 'solvent'):
            raise RuntimeError('The atoms of the solvent have not been selected.')
        # measure the COM of the protein
        protein_center = measure_center(self.protein.positions)
        # convert angstrom to nanometer and format the string
        protein_center = convert(protein_center, "angstrom", "nm")
        protein_center_str = f'({protein_center[0]}, {protein_center[1]}, {protein_center[2]})'
        self.logger.info('COM of the protein: ' + protein_center_str + '.')
        # generate the index file
        self.generateGromacsIndex(posixpath.join(generate_basename, 'colvars.ndx'))
        # generate the colvars configuration
        colvars_inputfile_basename_eq = posixpath.join(generate_basename, '007_eq_colvars')
        colvars_inputfile_basename = posixpath.join(generate_basename, '007_colvars')
        # measure the current COM distance from the ligand to protein
        ligand_center = measure_center(self.ligand.positions)
        ligand_center = convert(ligand_center, "angstrom", "nm")
        ligand_center_str = f'({ligand_center[0]}, {ligand_center[1]}, {ligand_center[2]})'
        self.logger.info('COM of the ligand: ' + ligand_center_str + '.')
        r_center = np.sqrt(np.dot(ligand_center - protein_center, ligand_center - protein_center))
        # round r_center
        r_center = np.around(r_center, 2)
        self.logger.info('Distance of protein and ligand: ' + str(r_center) + ' nm.')
        r_width = 0.01
        # r_lower_boundary = r_center - r_lower_shift
        # r_lower_shift is default to 0.2 nm
        r_lower_shift = 0.2
        r_lower_boundary = r_center - r_lower_shift
        if r_lower_boundary < 0:
            r_lower_boundary = 0.0
        # r_upper_boundary = r_center + r_upper_shift
        # r_upper_shift is default to 2.1 nm
        # also we will need r_upper_shift to enlarge the solvent box
        r_upper_shift = 2.1
        r_upper_boundary = r_center + r_upper_shift
        # colvars file for equilibration
        generateColvars(pkg_resources.read_text(templates_gromacs, '007_eq.colvars.template'),
                        colvars_inputfile_basename_eq,
                        logger=self.logger,
                        ligand_selection='BFEE_Ligand',
                        protein_selection='BFEE_Protein',
                        protein_center=protein_center_str)
        # colvars file for free-energy calculation
        generateColvars(pkg_resources.read_text(templates_gromacs, '007.colvars.template'),
                        colvars_inputfile_basename,
                        logger=self.logger,
                        r_width=r_width,
                        r_lower_boundary=r_lower_boundary,
                        r_upper_boundary=r_upper_boundary,
                        r_wall_constant=0.5*4.184,
                        ligand_selection='BFEE_Ligand',
                        protein_selection='BFEE_Protein',
                        protein_center=protein_center_str)
        # generate the reference file
        self.system.select_atoms('all').write(posixpath.join(generate_basename, 'reference.xyz'))
        # write the solvent molecules
        self.solvent.write(posixpath.join(generate_basename, 'solvent.gro'))
        # generate the shell script for making the tpr file
        # further enlarge the water box by 10% since the size of box may be compressed under NPT
        new_box_x = np.around(convert(self.system.dimensions[0], 'angstrom', 'nm'), 2) + r_upper_shift * 1.1
        new_box_y = np.around(convert(self.system.dimensions[1], 'angstrom', 'nm'), 2) + r_upper_shift * 1.1
        new_box_z = np.around(convert(self.system.dimensions[2], 'angstrom', 'nm'), 2) + r_upper_shift * 1.1
        # generate shell script for equlibration
        generateShellScript(pkg_resources.read_text(templates_gromacs, '007_eq.generate_tpr_sh.template'),
                            posixpath.join(generate_basename, '007.1_generate_eq_tpr'),
                            logger=self.logger,
                            BASENAME_002=self.stepnames[2],
                            BASENAME_003=self.stepnames[3],
                            BASENAME_004=self.stepnames[4],
                            BASENAME_005=self.stepnames[5],
                            BASENAME_006=self.stepnames[6],
                            BOX_MODIFIED_GRO_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                         'box_modified.gro')),
                                                                        os.path.abspath(generate_basename)).replace('\\', '/'),
                            MODIFIED_TOP_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                     'solvated.top')),
                                                                    os.path.abspath(generate_basename)).replace('\\', '/'),
                            MODIFIED_GRO_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                     'solvated.gro')),
                                                                    os.path.abspath(generate_basename)).replace('\\', '/'),
                            NEW_BOX_X_TEMPLATE=f'{new_box_x:.5f}',
                            NEW_BOX_Y_TEMPLATE=f'{new_box_y:.5f}',
                            NEW_BOX_Z_TEMPLATE=f'{new_box_z:.5f}',
                            MIN_MDP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                     '007_Minimize.mdp')),
                                                                    os.path.abspath(generate_basename)).replace('\\', '/'),
                            MDP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                 '007_Equilibration.mdp')),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            GRO_FILE_TEMPLATE=os.path.relpath(os.path.abspath(self.structureFile),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            TOP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(self.topologyFile),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            COLVARS_INPUT_TEMPLATE=os.path.relpath(os.path.abspath(colvars_inputfile_basename_eq + '.dat'),
                                                                     os.path.abspath(generate_basename)).replace('\\', '/'))
        # generate shell script for free-energy calculation
        generateShellScript(pkg_resources.read_text(templates_gromacs, '007.generate_tpr_sh.template'),
                            posixpath.join(generate_basename, '007.2_generate_tpr'),
                            logger=self.logger,
                            BASENAME_002=self.stepnames[2],
                            BASENAME_003=self.stepnames[3],
                            BASENAME_004=self.stepnames[4],
                            BASENAME_005=self.stepnames[5],
                            BASENAME_006=self.stepnames[6],
                            MDP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                 '007_PMF.mdp')),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            GRO_FILE_TEMPLATE=os.path.relpath(os.path.abspath(f'{self.baseDirectory}/007_r/output/007_r_eq.out.gro'),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            TOP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(f'{self.baseDirectory}/007_r/solvated.top'),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            COLVARS_INPUT_TEMPLATE=os.path.relpath(os.path.abspath(colvars_inputfile_basename + '.dat'),
                                                                     os.path.abspath(generate_basename)).replace('\\', '/'))
        # also copy the awk script to modify the colvars configuration according to the PMF minima in previous stages
        with pkg_resources.path(templates_gromacs, 'find_min_max.awk') as p:
            shutil.copyfile(p, posixpath.join(generate_basename, 'find_min_max.awk'))
        if not posixpath.exists(posixpath.join(generate_basename, 'output')):
            os.makedirs(posixpath.join(generate_basename, 'output'))
        self.logger.info(f"Generation of {generate_basename} done.")
        self.logger.info('=' * 80)

    def generate008(self):
        """generate files for determining the PMF along the RMSD of the ligand  
           with respect to its unbound state
        """

        self.handler.setFormatter(logging.Formatter('%(asctime)s [BFEEGromacs][008][%(levelname)s]:%(message)s'))
        generate_basename = self.basenames[8]
        self.logger.info('=' * 80)
        self.logger.info(f'Generating simulation files for {generate_basename}...')
        if not posixpath.exists(generate_basename):
            self.logger.info(f'Making directory {os.path.abspath(generate_basename)}...')
            os.makedirs(generate_basename)
        # # generate the MDP file for equlibration
        generateMDP(pkg_resources.read_text(templates_gromacs, '008.mdp.template'),
                    posixpath.join(generate_basename, '008_Equilibration'),
                    logger=self.logger,
                    timeStep=0.002,
                    numSteps=1000000,
                    temperature=self.temperature,
                    pressure=1.01325)
        # generate the MDP file
        generateMDP(pkg_resources.read_text(templates_gromacs, '008.mdp.template'),
                    posixpath.join(generate_basename, '008_PMF'),
                    logger=self.logger,
                    timeStep=0.002,
                    numSteps=4000000,
                    temperature=self.temperature,
                    pressure=1.01325)
        # generate the index file
        if hasattr(self, 'ligandOnly'):
            self.ligandOnly.write(posixpath.join(generate_basename, 'colvars_ligand_only.ndx'), name='BFEE_Ligand_Only')
        # generate the colvars configuration
        colvars_inputfile_basename = posixpath.join(generate_basename, '008_colvars')
        generateColvars(pkg_resources.read_text(templates_gromacs, '008.colvars.template'),
                        colvars_inputfile_basename,
                        rmsd_bin_width=0.005,
                        rmsd_lower_boundary=0.0,
                        rmsd_upper_boundary=0.5,
                        rmsd_wall_constant=0.8368,
                        ligand_selection='BFEE_Ligand_Only',
                        logger=self.logger)
        # generate the reference file for ligand only
        # extract the positions from the host-guest binding system
        ligand_position_in_system = self.ligand.positions
        # modify positions in the ligand-only system
        self.ligandOnly.positions = ligand_position_in_system
        # write out the whole ligand-only system as reference
        self.ligandOnlySystem.select_atoms('all').write(posixpath.join(generate_basename, 'reference_ligand_only.xyz'))
        # generate the shell script for making the tpr file for equilibration
        generateShellScript(pkg_resources.read_text(templates_gromacs, '008_eq.generate_tpr_sh.template'),
                            posixpath.join(generate_basename, '008.1_generate_eq_tpr'),
                            logger=self.logger,
                            MDP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                 '008_Equilibration.mdp')),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            GRO_FILE_TEMPLATE=os.path.relpath(os.path.abspath(self.ligandOnlyStructureFile),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            TOP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(self.ligandOnlyTopologyFile),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),)
        # generate the shell script for making the tpr file for free-energy calculation
        generateShellScript(pkg_resources.read_text(templates_gromacs, '008.generate_tpr_sh.template'),
                            posixpath.join(generate_basename, '008.2_generate_tpr'),
                            logger=self.logger,
                            MDP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(posixpath.join(generate_basename,
                                                                                                 '008_PMF.mdp')),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            GRO_FILE_TEMPLATE=os.path.relpath(os.path.abspath(f'{self.baseDirectory}/008_RMSD_unbound/output/008_RMSD_unbound_eq.out.gro'),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            TOP_FILE_TEMPLATE=os.path.relpath(os.path.abspath(self.ligandOnlyTopologyFile),
                                                                os.path.abspath(generate_basename)).replace('\\', '/'),
                            COLVARS_INPUT_TEMPLATE=os.path.relpath(os.path.abspath(colvars_inputfile_basename + '.dat'),
                                                                     os.path.abspath(generate_basename)).replace('\\', '/'))
        if not posixpath.exists(posixpath.join(generate_basename, 'output')):
            os.makedirs(posixpath.join(generate_basename, 'output'))
        self.logger.info(f"Generation of {generate_basename} done.")
        self.logger.info('=' * 80)

if __name__ == "__main__":
    bfee = BFEEGromacs('p41-abl.pdb', 'p41-abl.top', 'ligand-only.pdb', 'ligand-only.top', 'p41-abl-test/abc/def')
    bfee.setProteinHeavyAtomsGroup('segid SH3D and not (name H*)')
    bfee.setLigandHeavyAtomsGroup('segid PPRO and not (name H*)')
    bfee.setSolventAtomsGroup('resname TIP3*')
    bfee.setTemperature(350.0)
    bfee.generate000()
    bfee.generate001()
    bfee.generate002()
    bfee.generate003()
    bfee.generate004()
    bfee.generate005()
    bfee.generate006()
    bfee.generate007()
    bfee.generate008()
