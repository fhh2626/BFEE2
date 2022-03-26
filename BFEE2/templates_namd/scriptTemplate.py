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