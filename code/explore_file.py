import numpy
import pandas
import csv
import read_interaction_file as rif

def count_vertices(file):
    """Returns the number of vertices in a graph interaction file

    Parameters
    ----------
    file : .txt
        Path to .txt file containing two interacting nodes on each line

    Returns
    -------
    int
        Number of vertices in file
    """
    df = pandas.read_csv(file, sep = None, engine = 'python', skiprows=1, header=None)
    return len(df.index)

def count_edges(file):
    """Returns the number of edges in a graph interaction file

    Parameters
    ----------
    file : .txt
        Path to .txt file containing two interacting nodes on each line

    Returns
    -------
    int
        Number of edges in file
    """
    df = pandas.read_csv(file, sep = None, engine = 'python', skiprows=1, header=None)
    return len(set(numpy.concatenate(df.values)))

def clean_interactome(filein,fileout):
    """Reads and interaction file and writes a cleaned new interaction file
    Removes redundant interactions and homodimers

    Parameters
    ----------
    filein : .txt
        Path to .txt file containing two interacting nodes on each line
    fileout : .txt
        Path to new cleaned .txt file
    """
    list_out = rif.read_interaction_file_list(filein)
    for intrctn in list_out:
        if intrctn[0] == intrctn[1]:
            list_out.remove(intrctn)
    with open(filein, 'r') as file:
        delim = (csv.Sniffer().sniff(file.readlines()[1])).delimiter
    with open(fileout,'w', newline ="") as output:
        output.write(str(len(list_out)) + "\n")
        writer = csv.writer(output, delimiter=delim)
        for intrctn in list_out:
            writer.writerow(intrctn)

clean_interactome("../example_files/toy_example_dirty.txt","fileout.txt")
def get_degree(file,prot):
    """Returns the number of the degree of interactions of a node(prot)

    Parameters
    ----------
    file : .txt
        Path to .txt file containing two interacting nodes on each line
    prot : str
        Name of a protein in the network

    Returns
    -------
    int
        Degree of interactions of node "prot"
    """
    dict_nodes = rif.read_interaction_file_dict(file)
    return len(dict_nodes[prot])

def get_max_degree(file):
    """Returns the nodes with the highest degree and its corresponding degree

    Parameters
    ----------
    file : file : .txt
        Path to .txt file containing two interacting nodes on each line

    Returns
    -------
    max_nodes : list
        List of nodes with the highest degree
    max_count : int
        Highest degree
    """
    dict_nodes = rif.read_interaction_file_dict(file)
    max_count = max(len(v) for v in dict_nodes.values())
    max_nodes = [k for k, v in dict_nodes.items() if len(v) == max_count]
    return max_nodes, max_count

def get_ave_degree(file):
    """Returns the average degree of nodes in the interaction network

    Parameters
    ----------
    file : .txt
        Path to .txt file containing two interacting nodes on each line

    Returns
    -------
    int
        Average degree of nodes in the network
    """
    dict_nodes = rif.read_interaction_file_dict(file)
    list_degree = []
    for node in dict_nodes:
        list_degree.append(len(dict_nodes[node]))
    return round(sum(list_degree) / len(list_degree), 1)

def count_degree(file,deg):
    """Returns the number of nodes in network with a degree equal to deg

    Parameters
    ----------
    file : .txt
        Path to .txt file containing two interacting nodes on each line
    deg : int
        Degree

    Returns
    -------
    int
        Number of nodes with a degree equal to deg
    """
    dict_nodes = rif.read_interaction_file_dict(file)
    exact_count = 0
    for node in dict_nodes:
        if len(dict_nodes[node]) == deg:
            exact_count += 1
    return exact_count

def histogram_degree(file,dmin,dmax):
    """Prints an histogram with number of nodes having each degree
    Degrees are in the range of dmin,dmax

    Parameters
    ----------
    file : .csv
        Path to .csv file containing two interacting nodes on each line
    dmin : int
        Lower bound of degree range
    dmax : int
        Upper bound of degree range
    """
    dict_degrees = {}
    for deg in range(dmin,dmax+1):
        cnt_degree = count_degree(file,deg)
        dict_degrees[deg] = cnt_degree
        print(deg,"","*"*cnt_degree)
    return dict_degrees

# def histogram_degree(file,dmin,dmax):
#     """Prints an histogram with number of nodes having each degree
#     Degrees are in the range of dmin,dmax (dmax included)

#     Parameters
#     ----------
#     file : .txt
#         Path to .txt file containing two interacting nodes on each line
#     dmin : int
#         Lower bound of degree range
#     dmax : int
#         Upper bound of degree range
#     """
#     for deg in range(dmin,dmax+1):
#         print(deg,"","*"*count_degree(file,deg))