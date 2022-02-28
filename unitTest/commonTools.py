import os

def checkExistence(files):
    """check the existence of a series of files,
       print unexisted files

    Args:
        files (str): path of files
    """
    exist = True
    for name in files:
        if not os.path.exists(name):
            print(f'Error! {name} does not exist!')
            exist = False
    print('File existence check finished!')
    return exist
    
def compareTwoFiles(refFile, targetFile):
    """compare two files, print different lines

    Args:
        refFile (str): path of reference file
        targetFile (str): path of target file
        
    Returns:
        bool: whether the two files are the same
    """
    
    ref = open(refFile, 'r')
    target = open(targetFile, 'r')
    
    same = True
    
    lineCount = 0
    while True:
        refLine = ref.readline()
        targetLine = target.readline()
        
        if not refLine:
            break
        
        if refLine.strip() != targetLine.strip():
            #print(f'Warning! Lines {lineCount} of {refFile} and {targetFile} are not the same!')
            #print(f'The reference is {refLine}')
            #print(f'The target is {targetLine}')
            same = False
        
        lineCount += 1
        
    #print(f'Comparison between {refFile} and {targetFile} finished!')
    
    ref.close()
    target.close()
    
    return same
    
def checkAllFiles(refPath, targetPath, files):
    """check the existence of all the generated files, compare them with the references one by one,
       print all the differences

    Args:
        refPath (str): path of reference folder'
        targetPath (str): path of target folder'
        files (str): files to compare, relative to '{targetPath}/BFEE' and '{refPath}'
    """
    
    allFilesExist = checkExistence([f'{targetPath}/BFEE/{name}' for name in files])
    if allFilesExist:
        for name in files:
            #print(f'Comparing {targetPath}/BFEE/{name} with {refPath}/{name}')
            if not compareTwoFiles(f'{refPath}/{name}', f'{targetPath}/BFEE/{name}'):
                print(f'Warning, {targetPath}/BFEE/{name} and {refPath}/{name} are not the same!')
    else:
        print(f'Error! some files do not exist!')
    print(f'Checking finished!')
