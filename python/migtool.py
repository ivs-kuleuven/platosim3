"""
YAML Migration Tool

usage: migtool <old-inputfile.yaml> <new-inputfile.yaml>

Use Python 3.X to process.
"""
import argparse
import yaml

from collections import OrderedDict


# If we also want to preserve the comments from the YAML file, we will have to dig deep into PyYAML,
# or use a fork of that package called ruamel.yaml [https://pypi.python.org/pypi/ruamel.yaml] which was
# especially build for that.




# Print more info on changes and fixes to the commandline.

verbose = False

# The number of spaces used to indent subkeys when printing or writing to a file

tabSize = 4



class OrderedDictWithInsert(OrderedDict):
    def insert(self, after, newKey, newValue):
        """
        Insert the new key, value pair after the given key (after).
        
        Parameters:
            after    is the existing key after which the new key shall be inserted, if None prepend the new key
            newKey   is the new key to be inserted
            newValue is the new value to be inserted
        """
        items = self.copy().items()
        self.clear()
        if after == None:
            self.update({newKey: newValue})
        for key, value in items:
            self.update({key: value})
            if key == after:
                self.update({newKey: newValue})





def ordered_load(stream, Loader=yaml.Loader, object_pairs_hook=OrderedDictWithInsert):
    """
    Load the YAML file into an Ordered dictionary, i.e. keep the original order from the YAML inputfile.

    Parameters:
        stream  is the source where the YAML input is loaded from
        Loader  is the class to be used for reading the YAML file
        object_pairs_hook is used to map the key, value pairs into an OrderedDict
    """
    class OrderedLoader(Loader):
        pass

    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))

    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, construct_mapping)

    return yaml.load(stream, OrderedLoader)


def ordered_dump(data, stream=None, Dumper=yaml.Dumper, **kwds):
    """
    Dump the (ordered) dictionary into a YAML file.
    
    Parameters:
        data    is an Ordered Dictionary
        stream  is the destination where to stream the YAML output, if None the output is returned
        Dumper  is the class to be used to dump the data into a YAML stream
    """
    class OrderedDumper(Dumper):
        pass

    def _dict_representer(dumper, data):
        return dumper.represent_mapping(
            yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
            data.items())

    OrderedDumper.add_representer(OrderedDict, _dict_representer)

    return yaml.dump(data, stream, OrderedDumper, **kwds)







class YesNo(object):
    """
    Custom format class to make sure boolean values are represented as Yes / No instead of True / False.
    """
    def __init__(self, value):
        self.value = value
        
    def __format__(self, format):
        if (format == 'yesno'):
            if isinstance(self.value, bool):
                return "Yes" if self.value else "No"
        return "{}".format(self.value)




# The following classes define an action to be done on the key, value pairs
# The actions currently implemented are to print to screen or write to a file.

class ActionBase(object):
    def __init__(self, fileObject=None):
        self.fileObject = fileObject
        pass
    def action(self):
        pass

class PrintToScreen(ActionBase):
    def action(self, indent, key, parentKey, value, formatString=None):
        # when the parent key is the top key, add a blank line
        if len(parentKey.split('.')) == 1:
            print ()
        if formatString:
            print (formatString.format(indent, key, parentKey, value))
        else:
            if value != None:
                print ("{}{}: {:yesno}".format(indent, key, YesNo(value)))
            else:
                print ("{}{}:".format(indent, key))

# Key, value pairs are written to the fileObject in a way that closely matches the YAML inputfile
# Values are aligned at a certain column (valuesColumn) and boolean values are written as Yes/No.

# The column at which the values should be written in the YAML file in order to align them properly.
valuesColumn = 45

class WriteToFile(ActionBase):
    def action(self, indent, key, parentKey, value, formatString=None):
        if len(parentKey.split('.')) == 1:
            self.fileObject.write("\n")

        keyString = "{}{}: ".format(indent, key)
        if value != None:
            valueString = "{:yesno}".format(YesNo(value))
            n = valuesColumn - len(keyString)
            self.fileObject.write(keyString)
            self.fileObject.write(" "*n)
            self.fileObject.write(valueString)
        else:
            self.fileObject.write(keyString)
        self.fileObject.write("\n")








# Operations on YAML files

def load_yaml(filename):
    """
    Load a YAML file into an ordered dictionary.
    """
    fileStream = open(filename)
    data = ordered_load(fileStream, yaml.SafeLoader)
    return data

def save_yaml_OLD(filename, data):
    lines = ordered_dump(data, Dumper=yaml.SafeDumper, default_flow_style=False, indent=tabSize)
    lines = lines.split("\n")
    with open(filename, "w") as f:
        for line in lines:
            f.write(line + "\n")
    
def save_yaml(filename, data):
    """
    Save an ordered dictionary as a YAML file.
    """
    with open(filename, "w") as f:
        
        f.write("# PlatoSim3 configuration file\n")
        f.write("---\n")
        
        traverseDict(data, WriteToFile(fileObject=f))









def constructNewParentKey(parentKey, key):
    """
    Convenience function for constructing a dot-separated fully qualified key.
    """
    if parentKey:
        newParentKey = "{}.{}".format(parentKey, key)
    else:
        newParentKey = key
    return newParentKey
    





# Operations on dictionaries

def traverseDict(myDict, action=PrintToScreen(), parentKey=None, indent=""):
    """
    Traverse through the given dictionary and perform the given action at each step
    """
    for key, value in myDict.items():
        
        newParentKey = constructNewParentKey(parentKey, key)
        if isinstance(value, OrderedDict):
            action.action(indent, key, newParentKey, None)
            traverseDict(value, action, newParentKey, indent+" "*tabSize)
        else:
            action.action(indent, key, newParentKey, value)

def removeKeyFromDict(myDict, key):
    """
    Remove a key and it's value(s) permanantly from the dictionary.

    Parameters:        
      myDict is a dictionary (OrderedDict)
      key    is a fully qualified key, i.e. dot-separated key up to the root key
    """
    keys = key.split('.')
    newKey = keys.pop()
    for k in keys:
        myDict = myDict[k]
    myDict.pop(newKey)
    
def addKeyToDict(myDict, key, value):
    """
    Add a key to the dictionary. The new key will always be added to the end of the dictionary.

    Parameters:        
      myDict is a dictionary (OrderedDict)
      key    is a fully qualified key, i.e. dot-separated key up to the root key
      value  is any value for the given key
    """
    keys = key.split('.')
    newKey = keys.pop()
    for k in keys:
        myDict = myDict[k]
    myDict[newKey] = value

def insertKeyToDict(myDict, after, key, value):
    """
    Insert a key to the dictionary after the given key.

    Parameters:        
      myDict is a dictionary (OrderedDict)
      afteris the key after which the new key should be inserted
      key    is a fully qualified key, i.e. dot-separated key up to the root key
      value  is any value for the given key
    """
    keys = key.split('.')
    newKey = keys.pop()
    for k in keys:
        myDict = myDict[k]
    myDict.insert(after, newKey, value)

def replaceKeyInDict(myDict, key, value):
    myDict[key] = value

def findKeys(myDict, searchedKey, parentKey=None, foundKeys=None):
    """
    Return a list of fully qualified keys where the last key matches the searchedKey.
    """
    if not foundKeys:
        foundKeys = []
    for key, value in myDict.items():
        if key == searchedKey:
            foundKeys.append("{}.{}".format(parentKey, key))
            continue
        if isinstance(value, OrderedDict):
            foundKeys = findKeys(value, searchedKey, constructNewParentKey(parentKey, key), foundKeys)

    return foundKeys



def whatIsOldInYaml(oldDict, newDict, parentKey=None, root=None):

    for key, value in oldDict.copy().items():

        newParentKey = constructNewParentKey(parentKey, key)

        if key not in newDict:
            # Check if the key has moved
            foundKeys = findKeys(root, key)

            if not foundKeys:
                removeKeyFromDict(oldDict, key)
                if verbose:
                    print ("DONE - Obsolete key {} has been removed.".format(newParentKey))
            else:
                print ("CHECK - Key {} might have been moved to:".format(newParentKey))
                for foundKey in foundKeys:
                    print ("    {}".format(foundKey))
                # Action might be to move the key, but we don't know if the foundKey already 
                # exists in newDict
            continue
        
        if isinstance(value, OrderedDict):
            if isinstance(newDict[key], OrderedDict):
                whatIsOldInYaml(value, newDict[key], newParentKey, root)
            else:
                print ("CHECK - Key {} no longer contains sub-keys.".format(newParentKey))
        else:
            if isinstance(newDict[key], OrderedDict):
                # This option is handled in the whatIsNewInYaml part.
                # print "DONE - Key {} now contains sub-keys {}.".format(newParentKey, newDict[key].keys())
                pass


def whatIsNewInYaml(oldDict, newDict, parentKey=None, root=None):
    previousKey = None
    for key, value in newDict.copy().items():

        newParentKey = constructNewParentKey(parentKey, key)

        if key not in oldDict:
            if verbose:
                print ("DONE - New key {} has been added.".format(newParentKey))
                if isinstance(value, OrderedDict):
                    traverseDict(value, parentKey=newParentKey, indent=" "*tabSize)
            insertKeyToDict(oldDict, previousKey, key, value)
            previousKey = key
            continue

        if isinstance(value, OrderedDict):
            if isinstance(oldDict[key], OrderedDict):
                whatIsNewInYaml(oldDict[key], value, newParentKey, root)
            else:
                if verbose:
                    print ("DONE - New subkeys for {} have been added.".format(newParentKey))
                    traverseDict(value, parentKey=newParentKey, indent=" "*tabSize)
                replaceKeyInDict(oldDict, key, value)

        else:
            if value != oldDict[key]:
                print ("CHECK - Value changed for {} from {} to {}".format(newParentKey, oldDict[key], value))

        previousKey = key

            


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("oldFile", type=str, help="The old YAML file")
    parser.add_argument("newFile", type=str, help="The new YAML file")
    args = parser.parse_args()
    return args


 
def main():
    """
    newFile = "/Users/rik/Desktop/inputfiles-alphaRelease.yaml"
    oldFile = "/Users/rik/Desktop/inputfiles-master.yaml"
    """

    args = parse_arguments()

    newDict= load_yaml(args.newFile)
    oldDict = load_yaml(args.oldFile)

    whatIsOldInYaml(oldDict, newDict, root=newDict)

    whatIsNewInYaml(oldDict, newDict, root=newDict)

    save_yaml(args.oldFile.strip(".yaml")+"-changed.yaml", oldDict)




if __name__ == '__main__':
    main()
    
