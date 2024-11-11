import string

removeProteinTemplate = string.Template('''
mol new ${path}.psf
mol addfile ${path}.pdb
set aa [atomselect top "not $selectionPro"]
$$aa writepsf ${outputPath}.psf
$$aa writepdb ${outputPath}.pdb
exit
''')

removeMemProteinTemplate = string.Template('''
mol new ${path}.psf
mol addfile ${path}.pdb
set aa [atomselect top "$selectionLig"]
$$aa writepsf ${outputPath}_base.psf
$$aa writepdb ${outputPath}_base.pdb
mol delete top
package require solvate
solvate ${outputPath}_base.psf ${outputPath}_base.pdb -x 20 -y 20 -z 20 +x 20 +y 20 +z 20 -o ${outputPath} -s WT -b 2.2
mol delete top
mol new ${outputPath}.psf
mol addfile ${outputPath}.pdb
set aa [atomselect top all]
$$aa writexyz ${outputPath}.xyz
exit
''')

removeMemProteinFepTemplate = string.Template('''
package require autoionize
mol new ${path}.psf
mol addfile ${path}.pdb
set aa [atomselect top "$selectionLig"]
$$aa writepsf ${outputPath}_base.psf
$$aa writepdb ${outputPath}_base.pdb
mol delete top
package require solvate
solvate ${outputPath}_base.psf ${outputPath}_base.pdb -x 20 -y 20 -z 20 +x 20 +y 20 +z 20 -o ${outputPath} -s WT -b 2.2
mol delete top
mol new ${outputPath}.psf
mol addfile ${outputPath}.pdb
set aa [atomselect top all]
$$aa writexyz ${outputPath}.xyz
$$aa set beta 0
set solute [atomselect top "not water"]
$$solute set beta 1
$$aa writepdb ${outputFepPath}.pdb
exit
''')

removeProteinAmberTemplate = string.Template('''
parm ${path}.parm7
trajin ${path}.pdb
strip :$residueNum parmout ${outputPath}.parm7
trajout ${outputPath}.pdb
go
''')

neutralizeSystempTemplate = string.Template('''
package require autoionize
autoionize -psf ${path}.psf -pdb ${path}.pdb -neutralize -cation ${cationName} -anion ${anionName} -seg IO2 -o ${path}
mol new ${path}.psf
mol addfile ${path}.pdb
set aa [atomselect top all]
$$aa set beta 0
set solute [atomselect top "not water and not ion"]
$$solute set beta 1
$$aa writexyz ${path}.xyz
${extraCommand}
exit
''')

def charmmToGromacsTemplate(inputPrefix, forceFieldList):
    return f'''
import numpy as np
import parmed, MDAnalysis
from MDAnalysis import transformations
def measurePBC(uObject):
    atoms = uObject.select_atoms('all')
    atomPositions = atoms.positions
    xyz_array = np.transpose(atomPositions)
    min_x = np.min(xyz_array[0])
    max_x = np.max(xyz_array[0])
    min_y = np.min(xyz_array[1])
    max_y = np.max(xyz_array[1])
    min_z = np.min(xyz_array[2])
    max_z = np.max(xyz_array[2])
    return [
        np.array([max_x, max_y, max_z]) - np.array([min_x, min_y, min_z]),
        (np.array([min_x, min_y, min_z]) + np.array([max_x, max_y, max_z])) / 2
    ]
def charmmToGromacs(psfFile, pdbFile, prmFiles, PBC, outputPrefix):
    struct = parmed.load_file(psfFile)
    struct.load_parameters(
        parmed.charmm.CharmmParameterSet(*prmFiles)
    )
    struct.coordinates = parmed.load_file(pdbFile).coordinates
    struct.box = [PBC[0], PBC[1], PBC[2], 90, 90, 90]
    struct.save(f'{{outputPrefix}}.top', format='gromacs')
    struct.save(f'{{outputPrefix}}.gro')
uObject = MDAnalysis.Universe('{inputPrefix}.psf', '{inputPrefix}.pdb')
pbcVector = measurePBC(uObject)
allAtoms = uObject.select_atoms('all')
transformations.translate((-pbcVector[1] + pbcVector[0] / 2))(allAtoms)
allAtoms.write('{inputPrefix}.pdb', 'pdb', bonds=None)
charmmToGromacs(
    '{inputPrefix}.psf', 
    '{inputPrefix}.pdb',
    {forceFieldList},
    pbcVector[0],
    '{inputPrefix}_gmx'
)'''

amberToGromacsTemplate = string.Template('''
import numpy as np
import parmed, MDAnalysis
from MDAnalysis import transformations
def measurePBC(uObject):
    atoms = uObject.select_atoms('all')
    atomPositions = atoms.positions
    xyz_array = np.transpose(atomPositions)
    min_x = np.min(xyz_array[0])
    max_x = np.max(xyz_array[0])
    min_y = np.min(xyz_array[1])
    max_y = np.max(xyz_array[1])
    min_z = np.min(xyz_array[2])
    max_z = np.max(xyz_array[2])
    return [
        np.array([max_x, max_y, max_z]) - np.array([min_x, min_y, min_z]),
        (np.array([min_x, min_y, min_z]) + np.array([max_x, max_y, max_z])) / 2
    ]
def amberToGromacs(parmFile, rstFile, PBC, outputPrefix):
    struct = parmed.load_file(parmFile, xyz=rstFile)
    struct.coordinates = parmed.load_file(rstFile).coordinates
    struct.box = [PBC[0], PBC[1], PBC[2], 90, 90, 90]
    struct.save(f'{outputPrefix}.top', format='gromacs')
    struct.save(f'{outputPrefix}.gro')
uObject = MDAnalysis.Universe('${inputPrefix}.parm7', '${inputPrefix}.pdb')
pbcVector = measurePBC(uObject)
allAtoms = uObject.select_atoms('all')
transformations.translate((-pbcVector[1] + pbcVector[0] / 2))(allAtoms)
allAtoms.write('${inputPrefix}.pdb', 'pdb', bonds=None)
amberToGromacs(
    '${inputPrefix}.parm7',
    '${inputPrefix}.pdb',
    pbcVector[0],
    '${inputPrefix}_gmx',
)
''')

genarateLDDMFilesTemplate = string.Template('''
import numpy as np
import os, sys

TEMPERATURE = ${temperature}
WINDOWS = ${windows}
BOLTZMANN = 0.0019872041

def parseDat(filename, temperature=300.0):
    """Parse a dat (histogram) file, return the most probable CV value
       and the corresponding force constant for restraints

    Args:
        filename (str): the dat file to be parsed with
        temperature (float): the temperature of the simulation
        
    Returns:
        Tuple(float, float): the most probable CV value and the force constant
    """
    
    data = np.loadtxt(filename)
    CVs = data[:,0]
    counts = data[:,1]

    beta = 1 / (-BOLTZMANN * temperature)
    
    maxCV = -1
    maxCount = -1
    maxIndex = -1
    for index, (i, j) in enumerate(zip(CVs, counts)):
        if j > maxCount:
            maxCV = i
            maxCount = j
            maxIndex = index
    
    # fitting
    assert (maxIndex - 2 >=0 and maxIndex + 2 < len(counts)), f"maxIndex is: {maxIndex}, length is: {len(counts)}"
    X = [CVs[maxIndex - 2], CVs[maxIndex], CVs[maxIndex + 2]]
    y = [beta * np.log(counts[maxIndex - 2]), beta * np.log(counts[maxIndex]), beta * np.log(counts[maxIndex + 2])]

    forceConstant = np.polyfit(X, y, 2)

    return maxCV, forceConstant[0]

def showForceConstant(temperature=300):
    """ Print force constants obtained from histogram.dat files

    Args:
        temperature (int, optional): _description_. Defaults to 300.
    """
    #print(parseDat("eq.histogram1.dat", temperature))
    print(parseDat("../000_eq/output/eq.histogram2.dat", temperature))
    print(parseDat("../000_eq/output/eq.histogram3.dat", temperature))
    print(parseDat("../000_eq/output/eq.histogram4.dat", temperature))
    print(parseDat("../000_eq/output/eq.histogram5.dat", temperature))
    print(parseDat("../000_eq/output/eq.histogram6.dat", temperature))
    print(parseDat("../000_eq/output/eq.histogram7.dat", temperature))

def updateForceConstant(filename, CVs, newForceConstants):
    """ Update force constant of the CVs in a colvars file named filename. 
        Generates a new files named filename.tmp.

    Args:
        filename (str): path to the colvars files
        CVs (list[str]): list of CVs
        newForceConstants (List[float]): new force constants
    """
    # whether parsing a harmonic block
    # either None or the index of the corresponding CV
    inBlock = None
    with open(filename, 'r') as oldColvarsFile:
        # Python cannot modify text file in place
        with open(filename + '.tmp', 'w') as newColvarsFile:
            for line in oldColvarsFile.readlines():
                if line.strip().lower().startswith('colvarstrajfrequency'):
                    newColvarsFile.write(f'colvarsTrajFrequency      50\\n')
                    continue

                splitedLine = line.strip().split()
                
                if inBlock is None:
                    newColvarsFile.write(line)
                else:
                    if not line.strip().lower().startswith('forceconstant'):
                        newColvarsFile.write(line)
                    else:
                        # 'colvars' not in CVs
                        if inBlock == -1:
                            newColvarsFile.write(line)
                        else:
                            newColvarsFile.write(f'    ForceConstant    $$afc_{newForceConstants[inBlock]}\\n')
                    
                    if len(splitedLine) > 0:
                        if splitedLine[0] == '}':
                            inBlock = None
                        
                if line.strip().lower().startswith('harmonic'):
                    inBlock = -1
                        
                if line.strip().lower().startswith('colvars'):
                    for index, CV in enumerate(CVs):
                        if splitedLine[1].lower() == CV.lower():
                            inBlock = index

def updateAllForceConstants(temperature=300):
    """ update all the force constants of the colvars.in file,
        based on histogram.dat files

    Args:
        temperature (int, optional): temperature. Defaults to 300.
    """
    newForceConstants = []
    newForceConstants.append(parseDat("../000_eq/output/eq.histogram2.dat", temperature)[1])
    newForceConstants.append(parseDat("../000_eq/output/eq.histogram3.dat", temperature)[1])
    newForceConstants.append(parseDat("../000_eq/output/eq.histogram4.dat", temperature)[1])
    newForceConstants.append(parseDat("../000_eq/output/eq.histogram5.dat", temperature)[1])
    newForceConstants.append(parseDat("../000_eq/output/eq.histogram6.dat", temperature)[1])
    newForceConstants.append(parseDat("../000_eq/output/eq.histogram7.dat", temperature)[1])

    updateForceConstant("./colvars.in", 
                        ["eulerTheta", "eulerPhi", "eulerPsi", "polarTheta", "polarPhi", "r"],
                        newForceConstants)

def generateColvarsFiles(template_file, windows, prefix):

    if os.path.exists(f"./{prefix}"):
        print(f"Error, ./{prefix} exists!")
        sys.exit(1)
    
    os.mkdir(prefix)

    with open(template_file, "r") as inputFile:
        lines = inputFile.readlines()

    for i in range(windows):
        with open(f"./{prefix}/colvars_{i}.in", "w") as new_colvars_file:
            for j in range(len(lines)):
                splitedLine = lines[j].strip().split()
                if len(splitedLine) < 2 or (not splitedLine[1].startswith("$$afc_")):
                    new_colvars_file.write(lines[j])
                else:
                    new_colvars_file.write(f"    forceConstant  {float(splitedLine[1].split('_')[1]) * (float(i) / windows)}\\n")
    
    with open(f"./{prefix}/colvars_{windows}.in", "w") as new_colvars_file:
        for j in range(len(lines)):
            splitedLine = lines[j].strip().split()
            if len(splitedLine) < 2 or (not splitedLine[1].startswith("$$afc_")):
                new_colvars_file.write(lines[j])
            else:
                new_colvars_file.write(f"    forceConstant  {float(splitedLine[1].split('_')[1])}\\n")

updateAllForceConstants(TEMPERATURE)
generateColvarsFiles("colvars.in.tmp", WINDOWS, "colvars_files")
''')