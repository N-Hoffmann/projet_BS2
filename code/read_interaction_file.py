import pandas
import numpy
import os

def read_interaction_file_dict(input_file):
    """Reads an input interaction network
    Returns a dictionnary containing nodes as key and interacting nodes as values

    Parameters
    ----------
    input_file : .txt
        Path to .txt file containing two interacting nodes on each line

    Returns
    -------
    dict_node : dict
        Dictionnary containing nodes as key and interaction nodes as values
    """
    df = pandas.read_csv(input_file, sep = None, engine = 'python', skiprows=1, header=None)
    dict_node = {}
    for key in sorted(numpy.unique(numpy.concatenate(df.values))):
        dict_node[key] = []
    for i in range(len(df.index)):
        if df.loc[i,1] not in dict_node[df.loc[i,0]]:
            dict_node[df.loc[i,0]].append(df.loc[i,1])
        if df.loc[i,0] not in dict_node[df.loc[i,1]]:
            dict_node[df.loc[i,1]].append(df.loc[i,0])
    return dict_node

def read_interaction_file_list(input_file):
    """Reads an input interaction network
    Returns a list containing interacting nodes as tuples

    Parameters
    ----------
    input_file : .txt
        Path to .txt file containing two interacting nodes on each line
    
    Returns
    -------
    list_node : list
        List containing interacting nodes as tuples
    """
    df = pandas.read_csv(input_file, sep = None, engine = 'python', skiprows=1, header=None)
    list_node = []
    for i in range(len(df.index)):
        list_node.append((df.loc[i,0], df.loc[i,1]))
    sorted_nodes = sorted(list(set(tuple(sorted(l)) for l in list_node)))
    print(len(sorted_nodes))
    return sorted_nodes

def read_interaction_file_mat(input_file):
    """Reads an interaction file and creates an adjency matrix

    Parameters
    ----------
    input_file : .txt
        Path to .txt file containing two interacting nodes on each line 

    Returns
    -------
    adj_matrix : array
        Adjency matrix
    list_node : list
        Sorted list of nodes used to create the adjency matrix
    """
    int_dict = read_interaction_file_dict(input_file)
    list_node = sorted(int_dict.keys())
    adj_matrix = numpy.zeros((len(list_node), len(list_node)), dtype=int)
    for i in range(len(list_node)):
        for node in int_dict[list_node[i]]:
            adj_matrix[i][list_node.index(node)] = 1
    return adj_matrix, list_node

def read_interaction_file(input_file):
    """Reads an interaction network file and returns a dictionnary with the nodes as keys,
    a list of all the interactions, the adjency matrix and a sorted list of nodes
    This main function calls and returns all previous function

    Parameters
    ----------
    input_file : .txt
        Path to .txt file containing two interacting nodes on each line

    Returns
    -------
    d_int : dict
        Dictionnary containing nodes as key and interacting nodes as values
    l_int : list
        List containing interacting nodes as tuples
    m_int : array
        Adjency matrix
    l_som : list
        Sorted list of nodes
    """
    d_int = read_interaction_file_dict(input_file)
    l_int = read_interaction_file_list(input_file)
    m_int, l_som = read_interaction_file_mat(input_file)
    return d_int, l_int, m_int, l_som

def is_interaction_file(input_file):
    """Checks if input file is in a correct interaction file format
    Checks if input file is not empty
    Checks if first line is an integer and is the correct number of connectors
    Checks if there are only two columns in input file

    Parameters
    ----------
    input_file : .txt file
        Path to .txt file containing two interacting nodes on each line

    Returns
    -------
    True : boolean
        True if input file is a correct interaction file
    """
    if os.stat(input_file).st_size == 0:
        return False
    with open(input_file) as input:
        firstline = input.readlines()[0].rstrip()
    if firstline.isdigit() == False:
        return False
    if int(firstline) != len(read_interaction_file_list(input_file)):
        return False
    df = pandas.read_csv(input_file, sep = None, engine = 'python',skiprows=1, header=None)
    if len(df.columns) != 2:
        return False
    return True