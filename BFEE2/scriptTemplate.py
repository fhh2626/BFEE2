import string

removeProteinTemplate = string.Template('''
mol addfile ${path}.psf
mol addfile ${path}.pdb
set aa [atomselect top "not $selectionPro"]
$$aa writepsf ${outputPath}.psf
$$aa writepdb ${outputPath}.pdb
exit
''')

removeProteinAmberTemplate = string.Template('''
parm ${path}.parm7
trajin ${path}.pdb
strip :$residueNum parmout ${outputPath}.parm7
trajout ${outputPath}.pdb
go
''')