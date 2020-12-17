package require solvate
solvate ../complex.psf ../complex.pdb -x 12 -y 12 -z 12 +x 12 +y 12 +z 12 -o complex_largeBox -s W2 -b 2.2
mol addfile complex_largeBox.psf
mol addfile complex_largeBox.pdb
set all [atomselect top "all"]
$all writexyz complex_largeBox.xyz
$all delete
mol delete top
exit