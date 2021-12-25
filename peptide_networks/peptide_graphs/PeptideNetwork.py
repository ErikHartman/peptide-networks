import pandas as pd
import numpy as np
import networkx as nx
from networkx.algorithms import community
import argparse
import matplotlib.pyplot as plt
import re


class PeptideNetwork():
       
    def __init__(self,edge_list):
        """
        Creates network from given edge list

        Args:
            edge_list (Pandas DataFrame): edge list with the columns: from, to, distance, area and accession.
        """
        self.edge_list = edge_list
        self.edge_list['weight'] = 1/(1+edge_list['distance'])
        
    def create_network(self):
        """ Create network using the given edge list. 
        Source = 'from'
        Target = 'to'
        Edge attribute = 'weight'
        """
        self.G = nx.from_pandas_edgelist(
            self.edge_list,
            source='from',
            target='to',
            edge_attr="weight")
                
    def get_communitiess(self):
        """ Returns a communities object using the Girvan-Newman method.
        """
        return community.girvan_newman(self.G)
    
    def subset_edge_list_with_threshold(self,threshold):
        """ 
        Subsets the edge list based on distance threshold
        """
        subset_edge_list = self.edge_list[ (self.edge_list['distance'] <= threshold)]
        self.edge_list = subset_edge_list
    
    def subset_edge_list_with_protein(self, proteins):
        """
        Subsets the edge list based on wanted proteins
        """
        subset_edge_list = self.edge_list[self.edge_list['accession'].isin(proteins)]
        self.edge_list = subset_edge_list
        
    def print_distance_histogram(self,save_path):
        """ Saves a distance histogram to the given save path """
        x = self.edge_list['distance']
        plt.hist(x, density=True, bins=100)
        plt.title(f'Distance histogram {save_path}')
        plt.xlabel('distance')
        plt.savefig(f'{save_path}')
        
        
    def plot_network(self, save_path):
        """ Plots the network """
        pos = nx.spring_layout(self.G, weight="weight")
        degree_centrality = nx.degree_centrality(self.G)
        degree_centrality_color = np.fromiter(degree_centrality.values(), float)
        nx.draw(self.G, pos=pos, node_size=20, with_labels=False, width=1, alpha=0.6, node_color=degree_centrality_color)
        plt.title(f'Network')
        plt.savefig(f'{save_path}')
         
    def degree_analysis(self, save_path):
        """ Saves simple degree analysis to the given save path """
        degree_sequence = sorted([d for n, d in self.G.degree()], reverse=True)
        fig, ax= plt.subplots(2,1, figsize=(8, 8))
        ax[0].plot(degree_sequence, "b-", marker="o")
        ax[0].set_title("Degree Rank Plot")
        ax[0].set_ylabel("Degree")
        ax[0].set_xlabel("Rank")
        ax[1].bar(*np.unique(degree_sequence, return_counts=True))
        ax[1].set_title("Degree histogram")
        ax[1].set_xlabel("Degree")
        ax[1].set_ylabel("# of Nodes")
        fig.tight_layout()
        save_path = save_path.split('/')[-2] + '_' + re.sub('\D', '', save_path)
        plt.savefig(f'{save_path}')

# ------------------------------------------------------------------------------ #


def main(args):
    filepath = args.filepath
    edge_list = pd.read_csv(filepath)
    proteins = ['HBB_HUMAN', 'FIBA_HUMAN']
    threshold = 6
    G = PeptideNetwork(edge_list)
    G.subset_edge_list_with_protein(proteins)
    G.subset_edge_list_with_threshold(threshold)
    G.create_network()
    G.plot_network('findings/peptide_graphs/network_blosum_34.jpg')
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create edge list from file based on distance matrix')
    parser.add_argument("filepath", type=str, help="Path to file")
    args = parser.parse_args()
    main(args)