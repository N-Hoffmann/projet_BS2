import pandas
import numpy

# df.loc[row, column]
def read_interaction_file_dict(input_file):
    """Returns a dictionnary with keys as nodes and values as interacting nodes from an input interaction network

    Args:
        input (csv file): Path to .csv file containing on each line two interacting nodes
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
    return(dict_node)
#print(read_interaction_file_dict('example.csv'))

# def read_interaction_file_dict(input_file):
#     list_interactions= []
#     dict_node = {}
#     with open(input_file,"r") as file:
#         lines = file.readlines()
#     for line in lines:
#         line = line.rsplit('\n')[0].split(',')
#         list_interactions.append(line)
#     for key in [item for line in list_interactions for item in line]:
#         dict_node[key] = []
#     for i in range(1,12):
#         if list_interactions[i][1] not in dict_node[list_interactions[i][0]]:
#             dict_node[list_interactions[i][0]] += (list_interactions[i][1])
#         if list_interactions[i][0] not in dict_node[list_interactions[i][1]]:
#             dict_node[list_interactions[i][1]] += (list_interactions[i][0])
#     print(dict_node)

def read_interaction_file_list(input_file):
    df = pandas.read_csv(input_file,skiprows=1, header=None)
    list_node = []
    for i in range(len(df.index)):
        if (df.loc[i,0], df.loc[i,1]) or (df.loc[i,1],df.loc[i,0]) not in list_node:
            list_node.append((df.loc[i,0], df.loc[i,1]))
    return(list_node)
#print(read_interaction_file_list('example.csv'))

def read_interaction_file_mat(input_file):
    