# plot figures

import math
import os
import pathlib

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate


# an runtime error
# does not have corresponding correction for a pmf
class NoCorrectionFileError(RuntimeError):
    def __init__(self, arg):
        self.args = arg

def isGaWTM(pmfFiles):
    """determine whether the input PMFs indicate Gs-WTM simulations

    Args:
        pmfFiles (list[str]): path to a set of PMFs (and pmf corrections)
    
    Returns:
        bool: GaWTM simulation or not
    """

    for file in pmfFiles:
        if ''.join(pathlib.Path(file).suffixes) == '.reweightamd1.cumulant.pmf':
            return True
    return False

def getGaWTMBaseName(filePath):
    """Extract the base name from a GaWTM PMF file or its correction file.
    
    For example:
        'path/to/001.czar.pmf' -> '001'
        'path/to/001.reweightamd1.cumulant.pmf' -> '001'
        'path/to/step1.czar.pmf' -> 'step1'
    
    Args:
        filePath (str): path to the file
    
    Returns:
        str: base name of the file (without extensions)
    """
    fileName = pathlib.Path(filePath).name
    # Remove known GaWTM suffixes
    if fileName.endswith('.reweightamd1.cumulant.pmf'):
        return fileName[:-len('.reweightamd1.cumulant.pmf')]
    elif fileName.endswith('.czar.pmf'):
        return fileName[:-len('.czar.pmf')]
    else:
        # For other PMF files, just remove .pmf extension
        return pathlib.Path(fileName).stem

def pairGaWTMFiles(pmfFiles):
    """Pair GaWTM PMF files with their corresponding correction files.
    
    This function looks for files with matching base names:
    - xxx.czar.pmf paired with xxx.reweightamd1.cumulant.pmf
    
    Files without a corresponding correction file are returned as unpaired.
    Orphan correction files (without matching czar.pmf) are also returned separately.
    
    Args:
        pmfFiles (list[str]): list of all PMF file paths (including correction files)
    
    Returns:
        tuple: (paired_list, unpaired_czar_list, orphan_correction_list)
            paired_list: list of tuples (pmf_file_path, correction_file_path)
            unpaired_czar_list: list of czar.pmf file paths without correction
            orphan_correction_list: list of correction file paths without matching czar.pmf
    """
    # Separate czar.pmf files and correction files
    czar_files = {}  # baseName -> filePath
    correction_files = {}  # baseName -> filePath
    other_files = []  # files that are neither
    
    for filePath in pmfFiles:
        fileName = pathlib.Path(filePath).name
        if fileName.endswith('.reweightamd1.cumulant.pmf'):
            baseName = getGaWTMBaseName(filePath)
            correction_files[baseName] = filePath
        elif fileName.endswith('.czar.pmf'):
            baseName = getGaWTMBaseName(filePath)
            czar_files[baseName] = filePath
        else:
            other_files.append(filePath)
    
    # Pair files by base name
    paired = []
    unpaired_czar = []
    
    for baseName, czarPath in czar_files.items():
        if baseName in correction_files:
            paired.append((czarPath, correction_files[baseName]))
        else:
            unpaired_czar.append(czarPath)
    
    # Find orphan correction files (correction without czar.pmf)
    orphan_corrections = []
    for baseName, corrPath in correction_files.items():
        if baseName not in czar_files:
            orphan_corrections.append(corrPath)
    
    return paired, unpaired_czar, orphan_corrections, other_files

def correctGaWTM(pmfFile, correctionFile=None):
    """read a 1D namd PMF file and correct it using cumulant.pmf file

    Args:
        pmfFile (str): path to the pmf File
        correctionFile (str, optional): path to the correction file. 
            If None, will try to find it in the same directory (legacy behavior).
    
    Returns:
        np.array (N*2): 1D PMF
    """

    pmf = np.loadtxt(pmfFile)
    
    # If correction file is not provided, try to find it (legacy behavior)
    if correctionFile is None:
        correctionFile = pmfFile.replace('.czar.pmf', '') + '.reweightamd1.cumulant.pmf'

    if not os.path.exists(correctionFile):
        raise NoCorrectionFileError(f'{pmfFile} does not have a corresponding correction!')

    correction_data = np.loadtxt(correctionFile)
    correction_interpolate = interpolate.interp1d(correction_data[:,0], correction_data[:,1], fill_value="extrapolate")

    pmf[:,1] += correction_interpolate(pmf[:,0])

    return pmf

def readPMF(pmfFile):
    """read a 1D namd PMF file

    Args:
        pmfFile (str): path to the pmf File

    Returns:
        np.array (N*2): 1D PMF
    """

    return np.loadtxt(pmfFile)

def mergePMF(pmfFiles):
    """merge several PMF files

    Args:
        pmfFiles (list of np.arrays): list of 1D pmfs

    Returns:
        np.array (N*2): merged PMF if the PMFs overlap, pmfFiles[0] otherwise
    """    
    
    numPmfs = len(pmfFiles)
    assert(numPmfs > 0)
    
    # sort pmfs
    pmfSort = [i for i in range(numPmfs)]
    pmfSort.sort(key=lambda x: pmfFiles[x][0][0])
    
    finalPMF = pmfFiles[pmfSort[0]]
    
    if len(pmfFiles) > 1:
        for i in range(1, len(pmfFiles)):
            for j in range(len(finalPMF)):
                if finalPMF[j][0] == pmfFiles[pmfSort[i]][0][0]:
                    # overlapped region
                    avgDifference = np.average(finalPMF[j:,1:] - pmfFiles[pmfSort[i]][0:len(finalPMF)-j,1:])
                    pmfFiles[pmfSort[i]][:,1:] += avgDifference
                    finalPMF[j:,1:] = (finalPMF[j:,1:] + pmfFiles[pmfSort[i]][0:len(finalPMF)-j,1:]) / 2
                    # other region
                    finalPMF = np.append(finalPMF, pmfFiles[pmfSort[i]][len(finalPMF)-j:], axis=0)
                    break

    finalPMF[:,1] -= finalPMF[:,1].min()

    return finalPMF

def writePMF(pmfFile, pmf):
    """write a 1D namd PMF file

    Args:
        pmfFile (str): path to the pmf File
        pmf (np.array, N*2): pmf to be written
    """

    np.savetxt(pmfFile, pmf, fmt='%g')

def plotPMF(pmf):
    """plot a pmf

    Args:
        pmf (np.array, N*2): pmf to be plotted
    """
    
    plt.plot(pmf[:,0], pmf[:,1])
    plt.xlabel('Transition coordinate')
    plt.ylabel('ΔG (kcal/mol)')
    plt.show()

def plotHysteresis(forwardProfile, backwardProfile):
    """plot the profile describing the hysteresis between forward and backward
       simulations

    Args:
        forwardProfile (np.array, N*2): forward free-energy profile to be plotted
        backwardProfile (np.array, N*2): backward free-energy profile to be plotted
    """
    
    plt.plot(forwardProfile[:,0], forwardProfile[:,1], label='Forward')
    plt.plot(backwardProfile[:,0], backwardProfile[:,1], label='Backward')
    plt.xlabel('Lambda')
    plt.ylabel('ΔG (kcal/mol)')
    plt.legend()
    plt.show()

def saveHysteresis(forwardProfile, backwardProfile, filePath):
    """save the hysteresis data to a text file

    Args:
        forwardProfile (np.array, N*2): forward free-energy profile data
        backwardProfile (np.array, N*2): backward free-energy profile data
        filePath (str): path to save the data file
    """
    
    # Combine forward and backward data into a single array
    # Format: Lambda, Forward_dG, Backward_dG
    combined_data = np.column_stack([
        forwardProfile[:,0], 
        forwardProfile[:,1], 
        backwardProfile[:,1]
    ])
    header = 'Lambda\tForward_dG(kcal/mol)\tBackward_dG(kcal/mol)'
    np.savetxt(filePath, combined_data, fmt='%g', header=header, delimiter='\t')

def calcRMSD(inputArray):
    """calculate RMSD of a np.array with respect to (0,0,0,...0)

    Args:
        inputArray (1D np.array): the input array

    Returns:
        float: RMSD of a np.array with respect to (0,0,0,...0)
    """    

    sumG2 = sum(map(lambda x: x * x, inputArray))
    return math.sqrt(sumG2 / len(inputArray))

def readFrame(input):
    """read a frame of Colvars hist file and calculate its RMSD with respect to zero array

    Args:
        input (python file object): input object

    Returns:
        float: RMSD with respect to zero array
    """

    G = []
    while True:
        line = input.readline()
        
        # end of file
        if not line:
            return False
            
        splitedLine = line.strip().split()
        if splitedLine == []:
            if G == []:
                continue
            else:
                break
        if splitedLine[0].startswith('#'):
            continue

        G.append(float(splitedLine[1]))

    if G != []:
        return calcRMSD(G)
    else:
        return None

def parseHistFile(histPath):
    """parse a hist.czar.pmf file and return frame-RMSD list

    Args:
        histPath (str): path to a hist.czar.pmf file

    Returns:
        1D np.array: time evolution of RMSD with respect to zero array
    """    
    
    rmsd = []
    with open(histPath, 'r') as ifile:
        while True:
            rmsdPerFrame = readFrame(ifile)
            if rmsdPerFrame is False:
                break
            rmsd.append(rmsdPerFrame)
    return rmsd

def plotConvergence(rmsdList):
    """plot the time evolution of PMF rmsd

    Args:
        rmsdList (list or 1D np.array, float): time evolution of RMSD with respect to zero array
    """    

    plt.plot(range(1, len(rmsdList) + 1), rmsdList)
    plt.xlabel('Frame')
    plt.ylabel('RMSD (Colvars Unit)')
    plt.show()

def saveConvergence(rmsdList, filePath):
    """save the PMF RMSD convergence data to a text file

    Args:
        rmsdList (list or 1D np.array, float): time evolution of RMSD with respect to zero array
        filePath (str): path to save the data file
    """    

    frames = np.arange(1, len(rmsdList) + 1)
    data = np.column_stack([frames, rmsdList])
    header = 'Frame\tRMSD(Colvars_Unit)'
    np.savetxt(filePath, data, fmt='%g', header=header, delimiter='\t')