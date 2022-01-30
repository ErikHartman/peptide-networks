import pandas as pd
import numpy as np
import networkx as nx
from networkx.algorithms import community
import argparse
import matplotlib.pyplot as plt
import re
import itertools
import math


class PeptideNetwork():
       
    def __init__(self,edge_list):
        """
        Creates network from given edge list

        Args:
            edge_list (Pandas DataFrame): edge list with the columns: from, to, distance, area and accession.
        """
        
        edge_list['weight'] = 1/(1+edge_list['distance'])
        self.edge_list = edge_list
        
        
    def create_network(self):
        """ Create network using the given edge list. 
        Source = 'from'
        Target = 'to'
        Edge attribute = 'weight'
        """
        self.edge_list = self.edge_list[ (self.edge_list['distance']>0) ]
        self.edge_list.dropna(inplace=True)
        self.G = nx.from_pandas_edgelist(
            self.edge_list,
            source='from',
            target='to',
            edge_attr=["weight"])
                
    def get_communities(self):
        """ 
        Returns a communities object using the Girvan-Newman method.
        """
        return community.girvan_newman(self.G)
    
    def normalize_area(self):
        """ normalizes the area on total sample intensity """
        total_area = sum(self.edge_list['area_from'].values)
        order = int(math.log(total_area, 2))
        nf = 2**order
        total_area = total_area  / nf
        self.edge_list['area_from'] = self.edge_list['area_from']/total_area
        self.edge_list['area_to'] = self.edge_list['area_to']/total_area
        
        
    def get_connected_components(self):
        """
        Returns a sorted list of connected components.
        """
        return sorted(nx.connected_components(self.G), key=len)
    
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
        subset_edge_list = self.edge_list[self.edge_list['accession_from'].isin(proteins)]
        subset_edge_list = subset_edge_list[subset_edge_list['accession_to'].isin(proteins)]
        self.edge_list = subset_edge_list
        
 
    def get_edge_list(self):
        return self.edge_list
    
    def concatenate_edgelist(self, new_edgelist):
        self.edge_list = pd.concat([self.edge_list, new_edgelist])
        self.edge_list = self.edge_list.groupby(by=['from','to']).mean()
        self.edge_list.reset_index(inplace=True)
        self.edge_list['weight'] = 1/(1+self.edge_list['distance'])
    
    def get_network(self):
        return self.G
    
    def get_degree_sequence(self):
        return sorted([d for n, d in self.G.degree()], reverse=True)
        
    def plot_network(self, save_path):
        """ Plots the network. Color proportional to node color. """
        plt.clf()
        pos = nx.spring_layout(self.G, weight="weight")
        degree_centrality = nx.degree_centrality(self.G)
        degree_centrality_color = np.fromiter(degree_centrality.values(), float)
        nx.draw(self.G, pos=pos, node_size=20, with_labels=False, width=1, alpha=0.6, node_color=degree_centrality_color)
        plt.title(f'Network')
        plt.savefig(f'{save_path}')
        
    def plot_largest_community(self, save_path):
        """ Plots the largest community in the network. Color proportional to node color. """
        plt.clf()
        G_comm =  self.get_connected_components()
        G_comm = self.G.subgraph(G_comm)
        print(G_comm)
        pos = nx.spring_layout(G_comm, weight="weight")
        degree_centrality = nx.degree_centrality(G_comm)
        degree_centrality_color = np.fromiter(degree_centrality.values(), float)
        nx.draw(G_comm, pos=pos, node_size=20, with_labels=False, width=1, alpha=0.6, node_color=degree_centrality_color)
        plt.title(f'Network')
        plt.savefig(f'{save_path}')
         
    
    def global_centrality_analysis(self):
        """ Computes degree centrality, eigenvector centrality and closeness centrality for network """
        deg_cen = nx.degree_centrality(self.G)
        eig_cen = nx.eigenvector_centrality(self.G)
        close_cen = nx.closeness_centrality(self.G)
        return deg_cen, eig_cen, close_cen

    def get_community_positions(self, protein):
        """ Returns starts and end positions of peptides partitioned in communities """
        """ Need to implement this for all proteins in proteins list """
        def get_start(peptide_sequence, protein_sequence):
            start_pos = protein_sequence.find(peptide_sequence)
            return start_pos
        def get_end(peptide_sequence, protein_sequence):
            return get_start(peptide_sequence, protein_sequence) +len(peptide_sequence)
    
        uniprot_df = pd.read_csv('data/human_proteome.gz')
        protein_sequence = uniprot_df[ (uniprot_df['trivname'] == protein) ]
        protein_sequence = protein_sequence['seq'].values[0]
        comp = self.get_communities()
        starts = []
        ends = []
        areas = []
        for communities in itertools.islice(comp, 1):
            sequences = list(sorted(c) for c in communities)
            for seq in sequences:
                starts_temp = []
                ends_temp = []
                areas_temp = []
                for s in seq:
                    starts_temp.append(get_start(s, protein_sequence))
                    ends_temp.append(get_end(s, protein_sequence))
                    try:
                        areas_from = self.edge_list[(self.edge_list['from'] == s)]['area_from'].values[0]
                        areas_temp.append(areas_from)
                    except IndexError:
                        areas_to = self.edge_list[(self.edge_list['to'] == s)]['area_to'].values[0]
                        areas_temp.append(areas_to)
                starts.append(starts_temp)
                ends.append(ends_temp)
                areas.append(areas_temp)
        return starts, ends, areas
    
    def toString(self):
        print(f"Network with {len(self.G.nodes())} nodes and {self.G.number_of_edges()} edges")
        
# ------------------------------------------------------------------------------ #

def plot_community_positions(starts, ends, threshold,save_path):
    plt.clf()
    data = pd.DataFrame(data=zip(starts,ends))
    data.columns = ['starts','ends']
    data['n'] = data['starts'].apply(lambda x: len(x))
    data = data[( data['n'] >= threshold)]
    data['communities'] = data.index
    new_data = []
    for i in data.itertuples():
        comm = i[0]
        lst1 = i[1]
        lst2 = i[2]
        size = i[3]
        for col1,col2 in zip(lst1,lst2):
            new_data.append([col1, col2, comm, np.sqrt(size)])
    df_output = pd.DataFrame(data =new_data, columns=['start','end','communities','size'])
    df_output = df_output[ (df_output['start'] > 0 ) ] 
    df_output = df_output[ (df_output['end'] > 0 ) ] 
    df = df_output
    df = df.groupby('communities').agg(['mean', 'std'])
    df['width'] = df['end']['mean'] - df['start']['mean']
    df['x'] = df['start']['mean'] + (df['end']['mean']-df['start']['mean'])/2
    df['x_sd'] = (df['start']['std'] + df['end']['std'])/2
    df.columns = ['_'.join(col).rstrip('_') for col in df.columns.values]
    print(df)
    plt.bar(x=df['x'].values,height = df['size_mean'].values,width=df['width'].values , xerr = df['x_sd'].values, alpha=0.5, color='green') 
    plt.xlabel('Sequence')
    plt.ylabel('sqrt(peptides in community)')
    plt.savefig(save_path)
    
def plot_distance_histogram(edge_list, save_path):
    plt.clf()
    edge_list = edge_list[(edge_list['distance'] <= 20)] # remove infinities
    """ Saves a distance histogram to the given save path """
    x = edge_list['distance']
    plt.hist(x, density=True, bins=100)
    plt.title(f'Distance histogram {save_path}')
    plt.xlabel('distance')
    plt.savefig(f'{save_path}')
    
def plot_degree_analysis(degree_sequence, save_path):
    """ Saves simple degree analysis to the given save path """
    plt.clf()
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
    plt.savefig(f'{save_path}')
    
def plot_connected_component_sizes(connected_components, save_path):
    sizes = [len(cc) for cc in connected_components]
    print(sizes)
    x=range(len(sizes))
    print(x)
    plt.clf()
    plt.bar(x=x, height=sizes)
    plt.savefig(f'{save_path}')

    
def main(args):
    filepath = args.filepath
    edge_list = pd.read_csv(filepath)

    
    proteins = ['HBB_HUMAN']
    threshold = 3
    
    G = PeptideNetwork(edge_list)
    G.subset_edge_list_with_protein(proteins)
    G.subset_edge_list_with_threshold(threshold)
    G.create_network()
    cc = G.get_connected_components()
    plot_connected_component_sizes(cc, 'test2.jpg')

    #plot_distance_histogram(edge_list, 'findings/peptide_graphs/biophysical_34.jpg')
    #plot_degree_analysis(G.get_degree_sequence(), 'findings/peptide_graphs/degree_biophysical_34.jpg')
    
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create edge list from file based on distance matrix')
    parser.add_argument("filepath", type=str, help="Path to file")
    args = parser.parse_args()
    main(args)