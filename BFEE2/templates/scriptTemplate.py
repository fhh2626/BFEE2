import string

removeProteinTemplate = string.Template('''
mol addfile ${path}.psf
mol addfile ${path}.pdb
set aa [atomselect top "not $selectionPro"]
$$aa writepsf ${outputPath}.psf
$$aa writepdb ${outputPath}.pdb
exit
''')

removeMemProteinTemplate = string.Template('''
mol addfile ${path}.psf
mol addfile ${path}.pdb
set aa [atomselect top "$selectionLig"]
$$aa writepsf ${outputPath}_base.psf
$$aa writepdb ${outputPath}_base.pdb
mol delete top
package require solvate
solvate ${outputPath}_base.psf ${outputPath}_base.pdb -x 20 -y 20 -z 20 +x 20 +y 20 +z 20 -o ${outputPath} -s WT -b 2.2
mol delete top
mol addfile ${outputPath}.psf
mol addfile ${outputPath}.pdb
set aa [atomselect top all]
$$aa writexyz ${outputPath}.xyz
exit
''')

removeMemProteinFepTemplate = string.Template('''
mol addfile ${path}.psf
mol addfile ${path}.pdb
set aa [atomselect top "$selectionLig"]
$$aa writepsf ${outputPath}_base.psf
$$aa writepdb ${outputPath}_base.pdb
mol delete top
package require solvate
solvate ${outputPath}_base.psf ${outputPath}_base.pdb -x 20 -y 20 -z 20 +x 20 +y 20 +z 20 -o ${outputPath} -s WT -b 2.2
mol delete top
mol addfile ${outputPath}.psf
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