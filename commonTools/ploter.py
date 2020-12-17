# plot figures

import os, math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def readPMF(pmfFile):
    ''' read a 1D namd PMF file
        Inputs:
            pmfFile (string): path to the pmf File
        Return:
            np.array (N * 2): 1D PMF '''

    return np.loadtxt(pmfFile)

def mergePMF(pmfFiles):
    ''' merge several PMF files
        Inputs:
            pmfFiles (list of np.arrays): list of 1D pmfs
        Return:
            np.array (N * 2): merged PMF if the PMFs overlap, pmfFiles[0] otherwise '''
    
    assert(len(pmfFiles) > 0)
    finalPMF = pmfFiles[0]
    
    if len(pmfFiles) > 1:
        for i in range(1, len(pmfFiles)):
            for j in range(len(finalPMF)):
                if finalPMF[j][0] == pmfFiles[i][0][0]:
                    # overlapped region
                    avgDifference = np.average(finalPMF[j:,1:] - pmfFiles[i][0:len(finalPMF)-j,1:])
                    pmfFiles[i][:,1:] += avgDifference
                    finalPMF[j:,1:] = (finalPMF[j:,1:] + pmfFiles[i][0:len(finalPMF)-j,1:]) / 2
                    # other region
                    finalPMF = np.append(finalPMF, pmfFiles[i][len(finalPMF)-j:], axis=0)
                    break

    finalPMF[:,1] -= finalPMF[:,1].min()

    return finalPMF

def writePMF(pmfFile, pmf):
    ''' write a 1D namd PMF file
        Inputs:
            pmfFile (string): path to the pmf File
            pmf np.array (N * 2): pmf to be written '''

    return np.savetxt(pmfFile, pmf, fmt='%g')

def plotPMF(pmf):
    ''' plot a pmf
        Inputs:
            pmf np.array (N * 2): pmf to be plotted '''
    
    plt.plot(pmf[:,0], pmf[:,1])
    plt.xlabel('Transition coordinate')
    plt.ylabel('ΔG (kcal/mol)')
    plt.show()

def calcRMSD(inputArray):
    ''' calculate RMSD of a np.array with respect to (0,0,0,...0)
        Input:
            inputArray (1D np.array): the input array '''
    sumG2 = sum(map(lambda x: x * x, inputArray))
    return math.sqrt(sumG2 / len(inputArray))

def readFrame(input):
    ''' read a frame of Colvars hist file and calculate its RMSD with respect to zero array
        Input:
            input (python file object): input object
        Return:
            float: RMSD with respect to zero array '''

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
    ''' parse a hist.czar.pmf file and return frame-RMSD list
        Input:
            histPath (string): path to a hist.czar.pmf file
        Return:
            1D np.array: time evolution of RMSD with respect to zero array '''
    rmsd = []
    with open(histPath, 'r') as ifile:
        while True:
            rmsdPerFrame = readFrame(ifile)
            if rmsdPerFrame is False:
                break
            rmsd.append(rmsdPerFrame)
    return rmsd

def plotConvergence(rmsdList):
    ''' plot the time evolution of PMF rmsd
        Input:
            rmsdList (list or 1D np.array, float): 
                time evolution of RMSD with respect to zero array  '''
    plt.plot(range(1, len(rmsdList) + 1), rmsdList)
    plt.xlabel('Frame')
    plt.ylabel('RMSD (Colvars Unit)')
    plt.show()