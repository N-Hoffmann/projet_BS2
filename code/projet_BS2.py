import pandas
import numpy
from copy import deepcopy

def read_interaction_file_dict(input_file):
    """Reads an input interaction network
    Returns a dictionnary containing nodes as key and interacting nodes as values

    Parameters
    ----------
    input_file : .csv
        Path to .csv file containing two interacting nodes on each line
    """
    df = pandas.read_csv(input_file,skiprows=1, header=None)
    dict_node = {}
    for key in set(numpy.concatenate(df.values)):
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
    input_file : .csv
        Path to .csv file containing two interacting nodes on each line
    """
    df = pandas.read_csv(input_file,skiprows=1, header=None)
    list_node = []
    for i in range(len(df.index)):
        if (df.loc[i,0], df.loc[i,1]) or (df.loc[i,1],df.loc[i,0]) not in list_node:
            list_node.append((df.loc[i,0], df.loc[i,1]))
    return list_node

def read_interaction_file_mat(input_file):
    """Reads an interaction file and createss an adjency matrix

    Parameters
    ----------
    input_file : .csv
        Path to .csv file containing two interacting nodes on each line 

    Returns
    -------
    adj_matrix
        Adjency matrix
    ord_node
        Sorted list of nodes used to create the adjency matrix
    """
    df = pandas.read_csv(input_file,skiprows=1, header=None)
    list_node = set(numpy.concatenate(df.values))
    ord_node = sorted(list(deepcopy(list_node)))
    int_dict = read_interaction_file_dict(input_file)
    adj_matrix = numpy.zeros((len(ord_node), len(ord_node)), dtype=int)
    for i in range(len(ord_node)):
        for key in int_dict[ord_node[i]]:
            adj_matrix[i][ord_node.index(key)] = 1
    return adj_matrix, ord_node