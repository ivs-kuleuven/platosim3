def fileMatch(fileList, stringList):
    """
    fileMatch(fileList, stringList)
    
    returns the list of files from fileList in which all strings from stringList can be found, in any order
    Case sensitive
    """
    import fnmatch
    while (len(stringList)>1):
        fileList = fileMatch(fileList,[stringList.pop()])
    return [file for file in fileList if fnmatch.fnmatch(file,'*'+stringList[0]+"*")]


def fileSelect(stringList, location="./", listOrder=0):
    """
    fileSelect(stringList, location="./", listOrder=0)

    Returns a list of all files in 'location' with name matching every string in the list
    If listOrder is True, the ordering is forced to be identical to the one in stringList
    """
    import os, fnmatch
    allfiles = os.listdir(location)
    if listOrder:
        pattern = "*"
        i = 0
        while i < len(stringList): 
            pattern += stringList[i] + "*"
            i += 1
        return [file for file in allfiles if fnmatch.fnmatch(file,pattern)]
    else:
        return fileMatch(allfiles, stringList)

