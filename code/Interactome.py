import pandas
import numpy
import csv
import os
import random

class Interactome:
    """This class builds an Interactome object that represents a protein-protein 
    interaction network (PPI)

    Attributes
    ----------
    self.input : Array
        Array built from the input PPI file, used to build other attributes
    
    self.dict : dict
        Dictionnary containing proteins as key and interacting proteins as values
    
    self.list : list
        List containing interacting proteins as tuples
    
    self.protein : list
        List of all proteins in the network
    """
    def __init__(self,input_file):
        if self.is_interaction_file(input_file) == False:
            raise ValueError("The input file is not in a valid PPI network file format")
        self.input = pandas.read_csv(input_file, sep = None, engine = 'python', skiprows=1, header=None)
        self.dict = self.build_dict()
        self.list = self.build_list()
        self.protein = list(self.dict.keys())

    def build_dict(self):
        """
        Returns a dictionnary containing nodes as key and interacting nodes as values

        Parameters
        ----------
        input_file : str
        Path to .txt file containing two interacting nodes on each line

        Returns
        -------
        dict_node : dict
        Dictionnary containing nodes as key and interacting nodes as values
        """
        dict_node = {}
        for key in sorted(numpy.unique(numpy.concatenate(self.input.values))):
            dict_node[key] = []
        for i in range(len(self.input.index)):
            if self.input.loc[i,1] not in dict_node[self.input.loc[i,0]]:
                dict_node[self.input.loc[i,0]].append(self.input.loc[i,1])
            if self.input.loc[i,0] not in dict_node[self.input.loc[i,1]]:
                dict_node[self.input.loc[i,1]].append(self.input.loc[i,0])
        return dict_node

    def build_list(self):
        """Reads an input interaction network
        Returns a list containing interacting nodes as tuples

        Parameters
        ----------
        input_file : str
        Path to .txt file containing two interacting nodes on each line

        Returns
        -------
        list_node : list
        List containing interacting nodes as tuples
        """
        list_node = []
        for i in range(len(self.input.index)):
            list_node.append((self.input.loc[i,0], self.input.loc[i,1]))
        sorted_nodes = sorted(list(set(tuple(sorted(l)) for l in list_node)))
        return sorted_nodes
    
    def build_matrix(self):
        """Reads an interaction file and creates an adjency matrix

        Parameters
        ----------
        input_file : str
        Path to .txt file containing two interacting nodes on each line 

        Returns
        -------
        adj_matrix : array
        Adjency matrix
        list_node : list
        Sorted list of nodes used to create the adjency matrix
        """
        list_node = sorted(self.dict.keys())
        adj_matrix = numpy.zeros((len(list_node), len(list_node)), dtype=int)
        for i in range(len(list_node)):
            for node in self.dict[list_node[i]]:
                adj_matrix[i][list_node.index(node)] = 1
        return adj_matrix, list_node

    def count_edges(self):
        """Returns the number of edges in a graph interaction file

        Parameters
        ----------
        file : str
        Path to .txt file containing two interacting nodes on each line

        Returns
        -------
        int
        Number of edges in file
        """
        return len(self.input)
    
    def count_vertices(self):
        """Returns the number of vertices in a graph interaction file

        Parameters
        ----------
        file : str
        Path to .txt file containing two interacting nodes on each line

        Returns
        -------
        int
        Number of vertices in file
        """
        return len(self.dict.keys())

    def clean_interactome(self, filein, fileout):
        """Reads an interaction file and writes a cleaned new interaction file
        Removes redundant interactions and homodimers

        Parameters
        ----------
        filein : str
        Path to .txt file containing two interacting nodes on each line
        fileout : str
        Path to new cleaned .txt file
        """
        list_out = self.list
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

    def is_interaction_file(self,input_file):
        """Checks if input file is in a correct interaction file format
        Checks if input file is not empty
        Checks if first line is an integer and is the correct number of connectors
        Checks if there are only two columns in input file

        Parameters
        ----------
        input_file : str
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
        df = pandas.read_csv(input_file, sep = None, engine = 'python',skiprows=1, header=None)
        if int(firstline) != len(df):
            return False
        if len(df.columns) != 2:
            return False
        return True

    def get_degree(self, prot):
        """Returns the number of the degree of interactions of a node(prot)

        Parameters
        ----------
        file : str
        Path to .txt file containing two interacting nodes on each line
        prot : str
        Name of a protein in the network

        Returns
        -------
        int
            Degree of interactions of node "prot"
        """
        return len(self.dict[prot])
    
    def get_max_degree(self):
        """Returns the nodes with the highest degree and its corresponding degree

        Parameters
        ----------
        file : file : str
        Path to .txt file containing two interacting nodes on each line

        Returns
        -------
        max_nodes : list
            List of nodes with the highest degree
        max_count : int
            Highest degree
        """
        max_count = max(len(v) for v in self.dict.values())
        max_nodes = [k for k, v in self.dict.items() if len(v) == max_count]
        return max_nodes, max_count

    def get_ave_degree(self):
        """Returns the average degree of nodes in the interaction network

        Parameters
        ----------
        file : str
        Path to .txt file containing two interacting nodes on each line

        Returns
        -------
        int
            Average degree of nodes in the network
        """
        list_degree = []
        for node in self.dict:
            list_degree.append(len(self.dict[node]))
        return round(sum(list_degree) / len(list_degree), 1)

    def count_degree(self,deg):
        """Returns the number of nodes in network with a degree equal to deg

        Parameters
        ----------
        file : str
        Path to .txt file containing two interacting nodes on each line
        deg : int
        Degree

        Returns
        -------
        int
            Number of nodes with a degree equal to deg
        """
        exact_count = 0
        for node in self.dict:
            if len(self.dict[node]) == deg:
                exact_count += 1
        return exact_count

    def histogram_degree(self,dmin,dmax):
        """Prints an histogram with number of nodes having each degree
        Degrees are in the range of dmin,dmax

        Parameters
        ----------
        file : str
        Path to .csv file containing two interacting nodes on each line
        dmin : int
            Lower bound of degree range
        dmax : int
            Upper bound of degree range
        """
        dict_degrees = {}
        for deg in range(dmin,dmax+1):
            cnt_degree = self.count_degree(deg)
            dict_degrees[deg] = cnt_degree
            print(deg,"","*"*cnt_degree)
        return dict_degrees

    def density(self):
        """Gives the density of the graph

        Returns
        -------
        int
            The density of the graph
        """
        return (2*self.count_edges()) / (self.count_vertices()*(self.count_vertices()-1))

    def clustering(self,prot):
        """Gives the local clustering score of a protein (vertix)

        Parameters
        ----------
        prot : str
            The protein (vertix) for which we want to know its local clustering score

        Returns
        -------
        int
            Local clustering score for the protein
        """
        cnt_nghbr = len(self.dict[prot])
        if cnt_nghbr == 1:
            return 0
        else:
            list_nghbr = self.dict[prot]
            len_interactions = len([tup for tup in self.list if tup[0] in list_nghbr and tup[1] in list_nghbr])
            return (2*len_interactions) / (cnt_nghbr * (cnt_nghbr-1))

    def grapher(self,p):
        """Randomly generates vertexes between vertices using the Erdős-Rényi model

        Parameters
        ----------
        p : float
            Probability of a vertex existing between two vertices

        Returns
        -------
        ergraph_dict : dict
            Dictionary

        Raises
        ------
        ValueError
            p must be between 0 and 1
        """
        if p > 1 or p < 0:
            raise ValueError("p must be between 0 and 1")
        ergraph_dict = {}
        for key in self.protein:
            ergraph_dict[key] = []
        for protein_1 in self.protein:
            for protein_2 in self.protein:
                if protein_1 == protein_2:
                    continue
                if random.random() < p:
                    if protein_2 not in ergraph_dict[protein_1]:
                        ergraph_dict[protein_1].append(protein_2)
                    if protein_1 not in ergraph_dict[protein_2]:
                        ergraph_dict[protein_2].append(protein_1)
        return ergraph_dict

    def graphba(self):
        """Randomly generates edges between vertices using the Barabási–Albert model

        Returns
        -------
        bagraph_dict
            Dictionary
        """
        bagraph_dict = {}
        for key in self.protein:
            bagraph_dict[key] = []

        random_vertex = random.sample(self.protein, k=2)
        bagraph_dict[random_vertex[0]].append(random_vertex[1])
        bagraph_dict[random_vertex[1]].append(random_vertex[0])

        total_degree = 0
        for protein in bagraph_dict:
            total_degree += len(bagraph_dict[protein])

        for protein_1 in random.sample(self.protein, len(self.protein)):
            for protein_2 in self.protein:
                if protein_1 == protein_2:
                    continue
                if random.random() < (len(bagraph_dict[protein_1]) / total_degree):
                    if protein_2 not in bagraph_dict[protein_1]:
                        bagraph_dict[protein_1].append(protein_2)
                    if protein_1 not in bagraph_dict[protein_2]:
                        bagraph_dict[protein_2].append(protein_1)

            total_degree = 0
            for protein in bagraph_dict:
                total_degree += len(bagraph_dict[protein])

        return bagraph_dict

    def find_cc(self):
        """Finds connected components in the graph

        Returns
        -------
        list
            List of connected components in the graph
        """
        visited = {}
        for protein in self.protein:
            visited[protein] = False
        cc_list=[]

        def parse_nghbr(to_visit,edge,visited):
            visited[edge] = True
            to_visit.append(edge)
            for adj_edge in self.dict[edge]:
                if visited[adj_edge] == False:
                    to_visit = parse_nghbr(to_visit,adj_edge,visited)
            return to_visit

        for protein in self.protein:
            if visited[protein] == False:
                to_visit = []
                cc_list.append(parse_nghbr(to_visit,protein,visited))
        return cc_list

    def count_cc(self):
        """Finds connected components in the graph and returns their size and the total
        number of connected components

        Returns
        -------
        int
            Total number of CCs and their size
        """
        cc_list = self.find_cc()
        len_list = []
        for cc in cc_list:
            len_list.append(len(cc))
        return len(cc_list), len_list

    def write_cc(self,fileout):
        """Writes the different connected components of the graph in an output file
        Each line represents a connected component, with its size as the firs element and its 
        vertices as the second

        Parameters
        ----------
        fileout : str
            Path to outfile
        """
        cc_list = self.find_cc()
        with open(fileout,'w') as output:
            for i in range(len(cc_list)):
                line = len(cc_list[i]), cc_list[i]
                line = " ".join(map(str,line))
                output.write(line +"\n")

    def extract_cc(self, prot):
        """Returns all the vertices of the connected component where prot is

        Parameters
        ----------
        prot : str
            Protein for which we search the vertices of its connected components

        Returns
        -------
        list
            List of all the vertices of the connected component
        """
        cc_list = self.find_cc()
        for cc in cc_list:
            if prot in cc:
                return cc

    def compute_cc(self):
        """Computes the connected component group for each protein in the graph and returns the 
        group index in a list

        Returns
        -------
        lcc : list
            List containing the connect component group index for each protein in the graph
        """
        cc_list = self.find_cc()
        lcc = [0]*len(self.protein)
        for prot in self.protein:
            for cc in cc_list:
                if prot in cc:
                    lcc[self.protein.index(prot)] = cc_list.index(cc)
        return lcc

if __name__ == "__main__":
    ppi = Interactome("../example_files/toy2.txt")
    import networkx as nx
    import matplotlib.pyplot as plt
    # ppi_rand = ppi.grapher(0.3)
    # g=nx.Graph(ppi_rand)
    # ppi_ba = ppi.graphba()
    # g=nx.Graph(ppi_ba)
    # nx.draw(g, with_labels = True)
    # plt.draw()
    # plt.show()
    g=nx.Graph(ppi.dict)
    nx.draw(g, with_labels = True)
    plt.draw()
    plt.show()
    print(ppi.find_cc())
    print(ppi.count_cc())
    print(ppi.extract_cc("A"))
    print(ppi.compute_cc())
    print(ppi.protein)
    #ppi.write_cc("fileout.txt")