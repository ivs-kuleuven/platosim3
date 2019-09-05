
import ast
import yaml


class YAMLconfigurable(object):
    
    def __init__(self, YAMLFileName):

        self.yamlFileName = YAMLFileName
        
        with open(YAMLFileName, 'r') as stream:
            try:
                self.yamlDocument = yaml.load(stream)
            except yaml.YAMLError as exc:    
                print(exc)




    def getYamlConfiguration(self):

        """
        PURPOSE: Return the YAML configuration as a dictionary.

        OUTPUT:
            - YAML configuration as a dictionary
        """

        return self.yamlDocument





    def __contains__(self, key):

        """
        PURPOSE: Return True if the input parameter (key) is known/exists; False otherwise.

        INPUT:
            - key: string containing the parameter name or "Group/ParameterName" combination

        OUTPUT:
            - Boolean indicating whether or not the given key exists in the configuration file
        """

        if key.find('/') == -1:
            nodeNames = [key]
        else:
            nodeNames = key.split("/")

        node = self.yamlDocument

        for nodeName in nodeNames:
            print("> {}, {}".format(nodeName, type(node)))
            try:
                node = node[nodeName]
            except:
                return False

        return True





    def __getitem__(self, key):

        """
        PURPOSE: Return the value of the input parameter (key).

        INPUT: 
            - key: a string containing the parameter name or "Group/ParameterName" combination

        OUTPUT:
            - Value of the parameter
        """
        
        # Split the path into node names
        # E.g. "PSF/MappedGaussian/Sigma" into ["PSF", "MappedGaussian", "Sigma"]

        if key.find('/') == -1:

            print ("usage: the given parameter name (key) should include the group name of the group that contains the parameter.")
            print ("       E.g in 'Camera/PlateScale', Camera is the group, PlateScale is the parameter.")

            return None
        else:
            nodeNames = key.split("/")

        # Navigate to the deepest node, starting from the document root

        node = self.yamlDocument 
        for nodeName in nodeNames:

            if nodeName in node:
                node = node[nodeName]
            else:
                print("ERROR: The group '{}' was not found in the yaml inputfile '{}'.".format(key, self.configurationFilename))
                return None

        # node is a string, so cast it to its proper value

        try:
            value = ast.literal_eval(node)
        except ValueError:
            value = node
        
        # Return the value of the deepest node

        return value





    def __setitem__(self, key, item):
        
        """
        PURPOSE: Update a specific node.

        INPUT:
            - key: string with parent node name and node name seperated by a slash
            - item: string with the new node value, if not a string the value is converted using str()

        OUTPUT: 
            - True if node could be updated; False otherwise
        """
        
        # Ensure that the given item is a string

        item = str(item)

        # Split the path into node names
        # E.g. "PSF/MappedGaussian/Sigma" into ["PSF", "MappedGaussian", "Sigma"]

        if key.find('/') == -1:

            print ("usage: the given parameter name (key) should include the group name of the group that contains the parameter.")
            print ("       E.g in 'Camera/PlateScale', Camera is the group, PlatScale is the parameter.")
            return None
        
        else:
            nodeNames = key.split("/")

        # Check whether the parent node is in the document. If not, complain

        if nodeNames[0] not in self.yamlDocument:

            print("Error: no node with the name {0} found in input yaml file".format(nodeNames[0]))
            return False

        # If there is only 1 node in the path, we're finished after setting its value

        if len(nodeNames) == 1:

            self.yamlDocument[nodeNames[0]] = item
            return True

        # If we arrive here, there are at least 2 node in the path, check if 2nd parent node exists

        if nodeNames[1] not in self.yamlDocument[nodeNames[0]]:

            print("Error: no node with the name {0} found in input yaml file".format(nodeNames[0]+"/"+nodeNames[1]))
            return False

        # If there are only 2 nodes in the path, we're finished after setting its value

        if len(nodeNames) == 2:

            self.yamlDocument[nodeNames[0]][nodeNames[1]] = item
            return True

        # If we arrive here, there are at least 3 nodes in the path, check if 3rd parent node exists

        if nodeNames[2] not in self.yamlDocument[nodeNames[0]][nodeNames[1]]:

            print("Error: no node with the name {0} found in input yaml file".format(nodeNames[0]+"/"+nodeNames[1]+"/"+nodeNames[2]))
            return False

        # If there are only 3 nodes in the path, we're finished after setting its value

        if len(nodeNames) == 3:

            self.yamlDocument[nodeNames[0]][nodeNames[1]][nodeNames[2]] = item
            return True

        print("Error: detected more than 3 nodes in the path {0}".format(key))

        return False








    def writeYamlConfigurationFile(self, filename):

        """
        PURPOSE: Write the modified configuration to output file location. 
        
        INPUT:
            - filename: Filename
        """

        with open(filename, 'w') as outfile:
            outfile.write( pyaml.dump(self.yamlDocument, indent=4, width=120) )





    def __str__(self):

        """
        PURPOSE: Return a listing of all configuration settings.
                 If a parameter value has been updated for this run, [updated] will be printed after the value.
        """

        root = self.yamlDocument
        msg = "YAML Configuration:\n"
        msg += pyaml.dump(root, indent=4)

        return msg
