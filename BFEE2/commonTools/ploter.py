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
    """determine whether the input PMFs indicate Ga-WTM simulations

    Args:
        pmfFiles (list[str]): path to a set of PMFs (and pmf corrections)
    
    Returns:
        bool: GaWTM simulation or not
    """

    for file in pmfFiles:
        fileName = pathlib.Path(file).name
        if fileName.endswith('.reweightamd1.cumulant.pmf') or \
           fileName.endswith('.reweightamd1.reweight.pmf'):
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

def getGaWTMBaseNames(filePath):
    """Extract all possible base names from a GaWTM PMF file for flexible matching.
    
    This function returns a list of candidate base names to support flexible pairing.
    For example:
        'path/to/abf_1.abf1.czar.pmf' -> ['abf_1.abf1', 'abf_1']
        'path/to/001.czar.pmf' -> ['001']
        'path/to/step1.abf2.czar.pmf' -> ['step1.abf2', 'step1']
    
    Args:
        filePath (str): path to the file
    
    Returns:
        list[str]: list of possible base names (primary first, then alternatives)
    """
    primaryBaseName = getGaWTMBaseName(filePath)
    candidates = [primaryBaseName]
    
    # Check if the base name contains patterns like .abf1, .abf2, etc.
    # Find the last dot and check if it's followed by 'abf' + digits
    dotIndex = primaryBaseName.rfind('.')
    if dotIndex > 0:
        suffix = primaryBaseName[dotIndex + 1:]
        if suffix.startswith('abf') and suffix[3:].isdigit():
            alternativeBaseName = primaryBaseName[:dotIndex]
            if alternativeBaseName not in candidates:
                candidates.append(alternativeBaseName)
    
    return candidates

def pairGaWTMFiles(pmfFiles):
    """Pair GaWTM PMF files with their corresponding correction files.
    
    This function looks for files with matching base names:
    - xxx.czar.pmf paired with xxx.reweightamd1.cumulant.pmf
    
    Files without a corresponding correction file are returned as unpaired.
    Orphan correction files (without matching czar.pmf) are also returned separately.
    Wrong correction files (.reweightamd1.reweight.pmf instead of .cumulant.pmf) are also detected.
    
    Args:
        pmfFiles (list[str]): list of all PMF file paths (including correction files)
    
    Returns:
        tuple: (paired_list, unpaired_czar_list, orphan_correction_list, other_files, wrong_correction_files)
            paired_list: list of tuples (pmf_file_path, correction_file_path)
            unpaired_czar_list: list of czar.pmf file paths without correction
            orphan_correction_list: list of correction file paths without matching czar.pmf
            other_files: list of other PMF files
            wrong_correction_files: list of .reweightamd1.reweight.pmf files (wrong type)
    """
    # Separate czar.pmf files and correction files
    czar_files = {}  # baseName -> filePath
    czar_file_candidates = {}  # baseName -> list of candidate base names for matching
    correction_files = {}  # baseName -> filePath
    other_files = []  # files that are neither
    wrong_correction_files = []  # .reweightamd1.reweight.pmf files (wrong type)
    
    for filePath in pmfFiles:
        fileName = pathlib.Path(filePath).name
        if fileName.endswith('.reweightamd1.cumulant.pmf'):
            baseName = getGaWTMBaseName(filePath)
            correction_files[baseName] = filePath
        elif fileName.endswith('.reweightamd1.reweight.pmf'):
            # Wrong correction file type
            wrong_correction_files.append(filePath)
        elif fileName.endswith('.czar.pmf'):
            baseName = getGaWTMBaseName(filePath)
            czar_files[baseName] = filePath
            czar_file_candidates[baseName] = getGaWTMBaseNames(filePath)
        else:
            other_files.append(filePath)
    
    # Pair files by base name (with flexible matching)
    paired = []
    unpaired_czar = []
    matched_corrections = set()  # Track which correction files have been matched
    
    for baseName, czarPath in czar_files.items():
        # Try all candidate base names for this czar file
        candidates = czar_file_candidates.get(baseName, [baseName])
        matched = False
        for candidate in candidates:
            if candidate in correction_files:
                paired.append((czarPath, correction_files[candidate]))
                matched_corrections.add(candidate)
                matched = True
                break
        if not matched:
            unpaired_czar.append(czarPath)
    
    # Find orphan correction files (correction without czar.pmf)
    orphan_corrections = []
    for baseName, corrPath in correction_files.items():
        if baseName not in matched_corrections:
            orphan_corrections.append(corrPath)
    
    return paired, unpaired_czar, orphan_corrections, other_files, wrong_correction_files

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

# ============== History PMF functions ==============

def isGaWTMHist(pmfFiles):
    """determine whether the input History PMFs indicate Ga-WTM simulations

    Args:
        pmfFiles (list[str]): path to a set of History PMFs (and pmf corrections)
    
    Returns:
        bool: GaWTM simulation or not
    """
    for file in pmfFiles:
        fileName = pathlib.Path(file).name
        # Check for correction files (correct or wrong type)
        if fileName.endswith('.reweightamd1.cumulant.hist.pmf') or \
           fileName.endswith('.reweightamd1.cumulant.pmf') or \
           fileName.endswith('.reweightamd1.reweight.hist.pmf') or \
           fileName.endswith('.reweightamd1.reweight.pmf'):
            return True
    return False

def getGaWTMHistBaseName(filePath):
    """Extract the base name from a GaWTM History PMF file or its correction file.
    
    For example:
        'path/to/001.hist.czar.pmf' -> '001'
        'path/to/001.reweightamd1.cumulant.hist.pmf' -> '001'
        'path/to/001.reweightamd1.cumulant.pmf' -> '001'
        'path/to/step1.hist.czar.pmf' -> 'step1'
    
    Args:
        filePath (str): path to the file
    
    Returns:
        str: base name of the file (without extensions)
    """
    fileName = pathlib.Path(filePath).name
    # Remove known GaWTM history suffixes (order matters - check longer first)
    if fileName.endswith('.reweightamd1.cumulant.hist.pmf'):
        return fileName[:-len('.reweightamd1.cumulant.hist.pmf')]
    elif fileName.endswith('.reweightamd1.cumulant.pmf'):
        return fileName[:-len('.reweightamd1.cumulant.pmf')]
    elif fileName.endswith('.hist.czar.pmf'):
        return fileName[:-len('.hist.czar.pmf')]
    else:
        # For other PMF files, try to remove common extensions
        if fileName.endswith('.hist.pmf'):
            return fileName[:-len('.hist.pmf')]
        return pathlib.Path(fileName).stem

def getGaWTMHistBaseNames(filePath):
    """Extract all possible base names from a GaWTM History PMF file for flexible matching.
    
    This function returns a list of candidate base names to support flexible pairing.
    For example:
        'path/to/abf_1.abf1.hist.czar.pmf' -> ['abf_1.abf1', 'abf_1']
        'path/to/001.hist.czar.pmf' -> ['001']
    
    Args:
        filePath (str): path to the file
    
    Returns:
        list[str]: list of possible base names (primary first, then alternatives)
    """
    primaryBaseName = getGaWTMHistBaseName(filePath)
    candidates = [primaryBaseName]
    
    # Check if the base name contains patterns like .abf1, .abf2, etc.
    # Find the last dot and check if it's followed by 'abf' + digits
    dotIndex = primaryBaseName.rfind('.')
    if dotIndex > 0:
        suffix = primaryBaseName[dotIndex + 1:]
        if suffix.startswith('abf') and suffix[3:].isdigit():
            alternativeBaseName = primaryBaseName[:dotIndex]
            if alternativeBaseName not in candidates:
                candidates.append(alternativeBaseName)
    
    return candidates

def pairGaWTMHistFiles(pmfFiles):
    """Pair GaWTM History PMF files with their corresponding correction files.
    
    This function looks for files with matching base names:
    - xxx.hist.czar.pmf paired with xxx.reweightamd1.cumulant.hist.pmf or xxx.reweightamd1.cumulant.pmf
    
    Files without a corresponding correction file are returned as unpaired.
    Orphan correction files (without matching hist.czar.pmf) are also returned separately.
    Wrong correction files (.reweightamd1.reweight.pmf instead of .cumulant.pmf) are also detected.
    
    Args:
        pmfFiles (list[str]): list of all History PMF file paths (including correction files)
    
    Returns:
        tuple: (paired_list, unpaired_czar_list, orphan_correction_list, other_files, wrong_correction_files)
            paired_list: list of tuples (pmf_file_path, correction_file_path, is_hist_correction)
                is_hist_correction: True if correction is .hist.pmf, False if single-frame .pmf
            unpaired_czar_list: list of hist.czar.pmf file paths without correction
            orphan_correction_list: list of correction file paths without matching hist.czar.pmf
            other_files: list of other PMF files
            wrong_correction_files: list of .reweightamd1.reweight.pmf files (wrong type)
    """
    czar_files = {}  # baseName -> filePath
    czar_file_candidates = {}  # baseName -> list of candidate base names for matching
    correction_files = {}  # baseName -> (filePath, is_hist_correction)
    other_files = []
    wrong_correction_files = []  # .reweightamd1.reweight.pmf files (wrong type)
    
    for filePath in pmfFiles:
        fileName = pathlib.Path(filePath).name
        if fileName.endswith('.reweightamd1.cumulant.hist.pmf'):
            baseName = getGaWTMHistBaseName(filePath)
            correction_files[baseName] = (filePath, True)
        elif fileName.endswith('.reweightamd1.cumulant.pmf'):
            baseName = getGaWTMHistBaseName(filePath)
            # Only add if not already have a hist correction (prefer hist over single-frame)
            if baseName not in correction_files:
                correction_files[baseName] = (filePath, False)
        elif fileName.endswith('.reweightamd1.reweight.hist.pmf') or \
             fileName.endswith('.reweightamd1.reweight.pmf'):
            # Wrong correction file type
            wrong_correction_files.append(filePath)
        elif fileName.endswith('.hist.czar.pmf'):
            baseName = getGaWTMHistBaseName(filePath)
            czar_files[baseName] = filePath
            czar_file_candidates[baseName] = getGaWTMHistBaseNames(filePath)
        else:
            other_files.append(filePath)
    
    # Pair files by base name (with flexible matching)
    paired = []
    unpaired_czar = []
    matched_corrections = set()  # Track which correction files have been matched
    
    for baseName, czarPath in czar_files.items():
        # Try all candidate base names for this czar file
        candidates = czar_file_candidates.get(baseName, [baseName])
        matched = False
        for candidate in candidates:
            if candidate in correction_files:
                corrPath, is_hist = correction_files[candidate]
                paired.append((czarPath, corrPath, is_hist))
                matched_corrections.add(candidate)
                matched = True
                break
        if not matched:
            unpaired_czar.append(czarPath)
    
    orphan_corrections = []
    for baseName, (corrPath, _) in correction_files.items():
        if baseName not in matched_corrections:
            orphan_corrections.append(corrPath)
    
    return paired, unpaired_czar, orphan_corrections, other_files, wrong_correction_files

def readHistPMF(histPmfFile):
    """Read a History PMF file and return a list of PMF frames.
    
    Each frame is a 2D numpy array with shape (N, 2) where N is the number of points.
    
    Args:
        histPmfFile (str): path to the history PMF file
    
    Returns:
        list[np.array]: list of PMF frames, each as (N, 2) array
    """
    frames = []
    current_frame = []
    
    with open(histPmfFile, 'r') as f:
        for line in f:
            stripped = line.strip()
            
            # Skip empty lines
            if not stripped:
                # If we have accumulated data, save the frame
                if current_frame:
                    frames.append(np.array(current_frame))
                    current_frame = []
                continue
            
            # Skip comment/header lines
            if stripped.startswith('#'):
                continue
            
            # Parse data line
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    current_frame.append([x, y])
                except ValueError:
                    continue
    
    # Don't forget the last frame if file doesn't end with empty line
    if current_frame:
        frames.append(np.array(current_frame))
    
    return frames

def correctGaWTMHist(histPmfFile, correctionFile, is_hist_correction=True):
    """Apply GaWTM correction to a History PMF file.
    
    Args:
        histPmfFile (str): path to the history PMF file (.hist.czar.pmf)
        correctionFile (str): path to the correction file
        is_hist_correction (bool): True if correction file is history format (.hist.pmf),
            False if single-frame format (.pmf)
    
    Returns:
        list[np.array]: list of corrected PMF frames
    """
    pmf_frames = readHistPMF(histPmfFile)
    
    if not os.path.exists(correctionFile):
        raise NoCorrectionFileError(f'{histPmfFile} does not have a corresponding correction!')
    
    if is_hist_correction:
        # History correction file - apply frame by frame
        correction_frames = readHistPMF(correctionFile)
        
        # Interpolate correction frames if needed to match PMF frames
        if len(correction_frames) != len(pmf_frames):
            correction_frames = interpolateHistPMFFrames([correction_frames], len(pmf_frames))[0]
        
        corrected_frames = []
        for i, pmf in enumerate(pmf_frames):
            correction = correction_frames[i]
            correction_interp = interpolate.interp1d(
                correction[:, 0], correction[:, 1], fill_value="extrapolate"
            )
            corrected = pmf.copy()
            corrected[:, 1] += correction_interp(pmf[:, 0])
            corrected_frames.append(corrected)
        
        return corrected_frames
    else:
        # Single-frame correction file - apply same correction to all frames
        correction_data = np.loadtxt(correctionFile)
        correction_interp = interpolate.interp1d(
            correction_data[:, 0], correction_data[:, 1], fill_value="extrapolate"
        )
        
        corrected_frames = []
        for pmf in pmf_frames:
            corrected = pmf.copy()
            corrected[:, 1] += correction_interp(pmf[:, 0])
            corrected_frames.append(corrected)
        
        return corrected_frames

def interpolateHistPMFFrames(all_hist_pmfs, target_num_frames):
    """Interpolate history PMF frames to a target number of frames.
    
    Each PMF file may have different number of frames. This function interpolates
    the frame values to achieve a uniform frame count across all files.
    
    Args:
        all_hist_pmfs (list[list[np.array]]): list of history PMFs, each is a list of frames
        target_num_frames (int): target number of frames to interpolate to
    
    Returns:
        list[list[np.array]]: interpolated history PMFs with uniform frame count
    """
    interpolated = []
    
    for hist_pmf in all_hist_pmfs:
        num_frames = len(hist_pmf)
        
        if num_frames == target_num_frames:
            interpolated.append(hist_pmf)
            continue
        
        if num_frames == 0:
            interpolated.append([])
            continue
        
        # Create interpolated frames
        new_frames = []
        for target_idx in range(target_num_frames):
            # Calculate the source frame index (floating point)
            source_idx = target_idx * (num_frames - 1) / (target_num_frames - 1) if target_num_frames > 1 else 0
            
            # Get the two neighboring frames for interpolation
            lower_idx = int(source_idx)
            upper_idx = min(lower_idx + 1, num_frames - 1)
            
            # Interpolation weight
            weight = source_idx - lower_idx
            
            lower_frame = hist_pmf[lower_idx]
            upper_frame = hist_pmf[upper_idx]
            
            if weight == 0 or lower_idx == upper_idx:
                # No interpolation needed
                new_frames.append(lower_frame.copy())
            else:
                # Interpolate y values (x values should be the same across frames)
                # But to be safe, we interpolate upper_frame to lower_frame's x-coordinates
                upper_interp = interpolate.interp1d(
                    upper_frame[:, 0], upper_frame[:, 1], fill_value="extrapolate"
                )
                
                new_frame = lower_frame.copy()
                new_frame[:, 1] = (1 - weight) * lower_frame[:, 1] + weight * upper_interp(lower_frame[:, 0])
                new_frames.append(new_frame)
        
        interpolated.append(new_frames)
    
    return interpolated

def mergeHistPMF(all_hist_pmfs):
    """Merge multiple history PMF files frame-by-frame.
    
    Each input is a list of frames (from different PMF windows).
    The function first interpolates all to have the same number of frames,
    then merges each frame using the standard mergePMF function.
    
    Args:
        all_hist_pmfs (list[list[np.array]]): list of history PMFs, 
            each is a list of frames from one window
    
    Returns:
        list[np.array]: merged history PMF (list of merged frames)
    """
    if not all_hist_pmfs:
        return []
    
    # Find the maximum number of frames
    max_frames = max(len(hist_pmf) for hist_pmf in all_hist_pmfs)
    
    if max_frames == 0:
        return []
    
    # Interpolate all to have the same number of frames
    interpolated = interpolateHistPMFFrames(all_hist_pmfs, max_frames)
    
    # Merge frame-by-frame
    merged_frames = []
    for frame_idx in range(max_frames):
        # Collect the same frame from all windows
        frames_to_merge = [hist_pmf[frame_idx] for hist_pmf in interpolated if hist_pmf]
        
        if frames_to_merge:
            merged = mergePMF(frames_to_merge)
            merged_frames.append(merged)
    
    return merged_frames

def writeHistPMF(histPmfFile, frames):
    """Write a History PMF file.
    
    Args:
        histPmfFile (str): path to the output history PMF file
        frames (list[np.array]): list of PMF frames to write
    """
    with open(histPmfFile, 'w') as f:
        for frame_idx, frame in enumerate(frames):
            # Write frame header (similar to original format)
            f.write(f'# 1\n')
            if len(frame) > 0:
                x_min = frame[0, 0]
                x_max = frame[-1, 0]
                dx = (x_max - x_min) / (len(frame) - 1) if len(frame) > 1 else 0.1
                f.write(f'#  {x_min:.14e}  {dx:.14e}         {len(frame)}  0\n')
            f.write('\n')
            
            # Write data
            for row in frame:
                f.write(f'  {row[0]:.14e}   {row[1]:.14e}\n')
            
            # Add blank line between frames
            f.write('\n')

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