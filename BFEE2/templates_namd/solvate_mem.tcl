package require solvate
solvate ../complex.psf ../complex.pdb -x 0 -y 0 -z 16 +x 0 +y 0 +z 16 -o complex_largeBox -s W2 -b 2.2
mol addfile complex_largeBox.psf
mol addfile complex_largeBox.pdb
set all [atomselect top "all"]
$all writexyz complex_largeBox.xyz
$all delete
mol delete top
exit