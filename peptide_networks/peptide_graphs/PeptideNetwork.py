import pandas as pd
import numpy as np
import networkx as nx
from networkx.algorithms import community
import argparse
import matplotlib.pyplot as plt
import re
import itertools


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
        
 
    def get_edge_list(self):
        return self.edge_list
    
    def get_network(self):
        return self.G
    
    def get_degree_sequence(self):
        return sorted([d for n, d in self.G.degree()], reverse=True)
        
    def plot_network(self, save_path):
        plt.clf()
        """ Plots the network. Color proportional to node color. """
        pos = nx.spring_layout(self.G, weight="weight")
        degree_centrality = nx.degree_centrality(self.G)
        degree_centrality_color = np.fromiter(degree_centrality.values(), float)
        nx.draw(self.G, pos=pos, node_size=20, with_labels=False, width=1, alpha=0.6, node_color=degree_centrality_color)
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
        for communities in itertools.islice(comp, 1):
            sequences = list(sorted(c) for c in communities)
            for seq in sequences:
                starts.append([get_start(s, protein_sequence) for s in seq])
                ends.append([get_end(s, protein_sequence) for s in seq])
        return starts, ends
        
# ------------------------------------------------------------------------------ #

def plot_position_variance(starts, ends, threshold):
    plt.clf()
    data = pd.DataFrame(data=zip(starts,ends))
    data.columns = ['starts','ends']
    data['n'] = data['starts'].apply(lambda x: len(x))
    data = data[( data['n'] >= threshold)]
    data['communities'] = data.index
    new_data = []
    for i in data.itertuples():
        lst1 = i[1]
        lst2 = i[2]
        comm = i[0]
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
    plt.savefig('test.jpg')
    
def plot_distance_histogram(edge_list, save_path):
    plt.clf()
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

    
def main(args):
    filepath = args.filepath
    edge_list = pd.read_csv(filepath)

    plot_distance_histogram(edge_list, 'findings/peptide_graphs/levenshtein_31.jpg')
    # proteins = ['HBB_HUMAN']
    # threshold = 2
    # G = PeptideNetwork(edge_list)
    # G.subset_edge_list_with_protein(proteins)
    # G.subset_edge_list_with_threshold(threshold)
    # G.create_network()
    # G.plot_network('test2.jpg')
    




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create edge list from file based on distance matrix')
    parser.add_argument("filepath", type=str, help="Path to file")
    args = parser.parse_args()
    main(args)