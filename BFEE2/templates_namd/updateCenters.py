# Read eq.histogramX.dat and update Centers in *.in files
# This step is optional but may improve the convergence

import os, sys
import numpy as np

def isGeometric():
    """Check the route of the simulation, True for geometric,
       False for alchemical.
       
    Returns:
        bool: the route of the free-energy calculation. True for geometric,
              False for alchemical
    """
    
    if os.path.exists('../fep.tcl'):
        return False
    else:
        return True
    
    
def parseDat(filename):
    """Parse a dat (histogram) file and return the most probable CV value

    Args:
        filename (str): the dat file to be parsed with
        
    Returns:
        float: the most probable CV value
    """
    
    data = np.loadtxt(filename)
    CVs = data[:,0]
    counts = data[:,1]
    
    maxCV = -1
    maxCount = -1
    for i, j in zip(CVs, counts):
        if j > maxCount:
            maxCV = i
            maxCount = j
    return maxCV

def findOptimalCVs():
    """ read *.dat files generated through equilibration
        and find out the most probable values of CVs
        
    Returns:
        list[Float]: the most probable values of RMSD, eulerTheta, eulerPhi, eulerPsi, polarTheta, polarPhi, r
    """
    
    optimalCVs = []
    files = [
        'output/eq.histogram1.dat',
        'output/eq.histogram2.dat',
        'output/eq.histogram3.dat',
        'output/eq.histogram4.dat',
        'output/eq.histogram5.dat',
        'output/eq.histogram6.dat',
        'output/eq.histogram7.dat'
    ]
    
    for file in files:
        optimalCVs.append(parseDat(file))
    
    return optimalCVs
 
def readCenters(filename, CVs):
    """ find the value of centers of provided CVs in a *.in file

    Args:
        filename (str): a Colvars config file
        CVs (List[str]): a list, showing the CVs whose centers are to be get
        
    Returns:
        List[float]: centers of the provided CVs
    """
    
    # whether parsing a harmonic block
    # either None or the index of the corresponding CV
    inBlock = None
    
    # the CV values
    CVValues = [0 for i in range(len(CVs))]
    
    with open(filename, 'r') as colvarsFile:
        for line in colvarsFile.readlines():
            splitedLine = line.strip().split()
                
            if inBlock is not None:
                if line.strip().lower().startswith('centers'):
                    if inBlock != -1:
                        CVValues[inBlock] = float(splitedLine[1])
                    
                    if splitedLine[0] == '}':
                        inBlock = None
                        
            if line.strip().lower().startswith('harmonic'):
                inBlock = -1
                        
            if line.strip().lower().startswith('colvars'):
                for index, CV in enumerate(CVs):
                    if splitedLine[1].lower() == CV.lower():
                        inBlock = index
    
    return CVValues
 
def changeCenters(filename, CVs, newCenters):
    """ change Centers of harmonic restraints in a *.in file

    Args:
        filename (str): a Colvars config file
        CVs (List[str]): a list, showing the Centers of which harmonic blocks are to be replaced
        newCenters (List[float]): a list, containing the new values of Centers
    """
    
    # whether parsing a harmonic block
    # either None or the index of the corresponding CV
    inBlock = None
    with open(filename, 'r') as oldColvarsFile:
        # Python cannot modify text file in place
        with open(filename + '.tmp', 'w') as newColvarsFile:
            for line in oldColvarsFile.readlines():
                splitedLine = line.strip().split()
                
                if inBlock is None:
                    newColvarsFile.write(line)
                else:
                    if not line.strip().lower().startswith('centers'):
                        newColvarsFile.write(line)
                    else:
                        # 'colvars' not in CVs
                        if inBlock == -1:
                            newColvarsFile.write(line)
                        else:
                            newColvarsFile.write(f'    Centers    {newCenters[inBlock]}\n')
                    
                    if len(splitedLine) > 0:
                        if splitedLine[0] == '}':
                            inBlock = None
                        
                if line.strip().lower().startswith('harmonic'):
                    inBlock = -1
                        
                if line.strip().lower().startswith('colvars'):
                    for index, CV in enumerate(CVs):
                        if splitedLine[1].lower() == CV.lower():
                            inBlock = index
    
    os.remove(filename)
    os.rename(filename + '.tmp', filename)
    
def changeABFRange(filename, CVs, newCenters, originalCenters):
    """ change lowerboundary, upperboundary, lowerWalls, upperWalls of ABF in a *.in file

    Args:
        filename (str): a Colvars config file
        CVs (List[str]): a list, showing the corresponding options of which CV blocks are to be replaced
        newCenters (List[float]): a list, containing the new values of Centers of these CVs
        originalCenters (List[float]): a list, containing the original values of Centers of these CVs
    """
    
    # whether parsing a harmonic block
    # either None or the index of the corresponding CV
    inBlock = None
    with open(filename, 'r') as oldColvarsFile:
        # Python cannot modify text file in place
        with open(filename + '.tmp', 'w') as newColvarsFile:
            for line in oldColvarsFile.readlines():
                splitedLine = line.strip().split()
                
                if inBlock is None:
                    newColvarsFile.write(line)
                else:
                    if (not line.strip().lower().startswith('lowerboundary')) and \
                        (not line.strip().lower().startswith('upperboundary')) and \
                        (not line.strip().lower().startswith('lowerwalls')) and \
                        (not line.strip().lower().startswith('upperwalls')):
                        newColvarsFile.write(line)
                    else:
                        # 'name' not in CVs
                        if inBlock == -1:
                            newColvarsFile.write(line)
                        else:
                            if line.strip().lower().startswith('lowerboundary'):
                                newColvarsFile.write(
                                    f'    Lowerboundary    {float(splitedLine[1]) + newCenters[inBlock] - originalCenters[inBlock]}\n'
                                )
                            elif line.strip().lower().startswith('upperboundary'):
                                newColvarsFile.write(
                                    f'    Upperboundary    {float(splitedLine[1]) + newCenters[inBlock] - originalCenters[inBlock]}\n'
                                    )
                            elif line.strip().lower().startswith('lowerwalls'):
                                newColvarsFile.write(
                                    f'    LowerWalls    {float(splitedLine[1]) + newCenters[inBlock] - originalCenters[inBlock]}\n'
                                    )
                            elif line.strip().lower().startswith('upperwalls'):
                                newColvarsFile.write(
                                    f'    UpperWalls    {float(splitedLine[1]) + newCenters[inBlock] - originalCenters[inBlock]}\n'
                                    )
                            else:
                                print("Selected CV(s) not for free-energy calculation!")
                    
                    if len(splitedLine) > 0:
                        if splitedLine[0] == '}':
                            inBlock = None
                
                if len(splitedLine) > 0:
                    if splitedLine[0].lower() == 'colvar' or line.strip().lower().startswith('harmonicwalls'):
                        inBlock = -1
                    
                    if line.strip().lower().startswith('name') or splitedLine[0].lower() == 'colvars':
                        for index, CV in enumerate(CVs):
                            if splitedLine[1].lower() == CV.lower():
                                inBlock = index
    
    os.remove(filename)
    os.rename(filename + '.tmp', filename)

def updateAlchemicalFiles():
    """Update the Centers of *.in files based on equilibration
    """
    
    files = [
        '../001_MoleculeBound/colvars.in',
        '../002_RestraintBound/colvars_forward.in',
        '../002_RestraintBound/colvars_backward.in'
    ]
    
    for file in files:
        changeCenters(
            file, 
            [
                'eulerTheta', 'eulerPhi', 'eulerPsi', 'polarTheta', 'polarPhi', 'r'
            ], 
            findOptimalCVs()[1:]
        )
        print(f'File {file} updated!')
        
def updateGeometricFiles():
    """Update the Centers of *.in files based on equilibration
    """
    
    CVNames = ['eulerTheta', 'eulerPhi', 'eulerPsi', 'polarTheta', 'polarPhi']
    optimalCVValues = findOptimalCVs()[1:6]
    originalCVValues = readCenters('../007_r/colvars_eq.in', CVNames)
    
    # in geometrical route, there may be many windows
    folders = [
        '../002_EulerTheta',
        '../003_EulerPhi',
        '../004_EulerPsi',
        '../005_PolarTheta',
        '../006_PolarPhi',
        '../007_r'
    ]
    
    for fIndex, folder in enumerate(folders):
        # *.in in folder 
        files = []
        for file in os.listdir(folder):
            if file.endswith('.in'):
                files.append(os.path.join(folder, file))
                
        for inFile in files:
            changeCenters(
                inFile, 
                CVNames[:fIndex], 
                optimalCVValues[:fIndex]
            )
            
            changeABFRange(
                inFile,
                CVNames[:fIndex+1],
                optimalCVValues[:fIndex+1],
                originalCVValues[:fIndex+1]
            )
            print(f'File {inFile} updated!')
        
if __name__ == '__main__':
    
    if isGeometric():
        updateGeometricFiles()
        print('All the *.in files updated for the Geometrical route!')
    else:
        updateAlchemicalFiles()
        print('All the *.in files updated for the Alchemical route!')