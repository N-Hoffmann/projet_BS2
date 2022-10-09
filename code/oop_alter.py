import pandas
import numpy
import csv

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
        input = pandas.read_csv(input_file, sep = None, engine = 'python', skiprows=1, header=None)
        dict_node = {}
        list_node = []

        for key in sorted(numpy.unique(numpy.concatenate(input.values))):
            dict_node[key] = []
        for i in range(len(input.index)):
            list_node.append((input.loc[i,0], input.loc[i,1]))
            if input.loc[i,1] not in dict_node[input.loc[i,0]]:
                dict_node[input.loc[i,0]].append(input.loc[i,1])
            if input.loc[i,0] not in dict_node[input.loc[i,1]]:
                dict_node[input.loc[i,1]].append(input.loc[i,0])
        
        self.dict = dict_node
        self.list = sorted(list(set(tuple(sorted(l)) for l in list_node)))
        self.protein = list(self.dict.keys())
    
    def get_matrix(self):
        list_node = sorted(self.dict.keys())
        adj_matrix = numpy.zeros((len(list_node), len(list_node)), dtype=int)
        for i in range(len(list_node)):
            for node in self.dict[list_node[i]]:
                adj_matrix[i][list_node.index(node)] = 1
        return adj_matrix, list_node

    def count_edges(self):
        return len(self.list)
    
    def count_vertices(self):
        return len(self.dict.keys())

    def clean_interactome(self, filein, fileout):
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

    def get_degree(self, prot):
        return len(self.dict[prot])
    
    def get_max_degree(self):
        max_count = max(len(v) for v in self.dict.values())
        max_nodes = [k for k, v in self.dict.items() if len(v) == max_count]
        return max_nodes, max_count

    def get_ave_degree(self):
        list_degree = []
        for node in self.dict:
            list_degree.append(len(self.dict[node]))
        return round(sum(list_degree) / len(list_degree), 1)

    def count_degree(self,deg):
        exact_count = 0
        for node in self.dict:
            if len(self.dict[node]) == deg:
                exact_count += 1
        return exact_count

    def histogram_degree(self,dmin,dmax):
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
        list_nghbr = self.dict[prot]
        len_interactions = len([tup for tup in self.list if tup[0] in list_nghbr and tup[1] in list_nghbr])
        return len_interactions / ((cnt_nghbr * (cnt_nghbr-1)) / 2)