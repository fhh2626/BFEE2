# plot figures

import os, math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def readPMF(pmfFile):
    """read a 1D namd PMF file

    Args:
        pmfFile (str): path to the pmf File

    Returns:
        np.array (N * 2): 1D PMF
    """

    return np.loadtxt(pmfFile)

def mergePMF(pmfFiles):
    """merge several PMF files

    Args:
        pmfFiles (list of np.arrays): list of 1D pmfs

    Returns:
        np.array (N * 2): merged PMF if the PMFs overlap, pmfFiles[0] otherwise
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
        pmf np.array (N * 2): pmf to be written
    """

    np.savetxt(pmfFile, pmf, fmt='%g')

def plotPMF(pmf):
    """plot a pmf

    Args:
        pmf np.array (N * 2): pmf to be plotted
    """
    
    plt.plot(pmf[:,0], pmf[:,1])
    plt.xlabel('Transition coordinate')
    plt.ylabel('ΔG (kcal/mol)')
    plt.show()

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