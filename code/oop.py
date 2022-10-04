import pandas
import numpy
import csv

class Interactome:
    def __init__(self,input_file):
        self.input = pandas.read_csv(input_file, sep = None, engine = 'python', skiprows=1, header=None)
        self.dict = self.build_dict()
        self.list = self.build_list()
        self.protein = list(self.dict.keys())

    def build_dict(self):
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
        list_node = []
        for i in range(len(self.input.index)):
            list_node.append((self.input.loc[i,0], self.input.loc[i,1]))
        sorted_nodes = sorted(list(set(tuple(sorted(l)) for l in list_node)))
        return sorted_nodes
    
    def build_matrix(self):
        list_node = sorted(self.dict.keys())
        adj_matrix = numpy.zeros((len(list_node), len(list_node)), dtype=int)
        for i in range(len(list_node)):
            for node in self.dict[list_node[i]]:
                adj_matrix[i][list_node.index(node)] = 1
        return adj_matrix, list_node

    def count_edges(self):
        return len(self.input)
    
    def count_vertices(self):
        return len(self.dict.keys())

    def clean(self, filein, fileout):
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