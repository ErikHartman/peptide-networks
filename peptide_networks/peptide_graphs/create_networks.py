import pandas as pd
import numpy as np
import networkx as nx
import argparse
import matplotlib.pyplot as plt
import re


def create_network_from_edge_list(edge_list):
    edge_list['weight'] = 1/(edge_list['distance'] + 1)
    G = nx.from_pandas_edgelist(
    edge_list,
    source='from',
    target='to',
    edge_attr="weight")
    return G


def plot_network(G, filepath):
    """ want to create color mapper based on accession value"""
    save_path = filepath.split('/')[-2] + '_' + re.sub('\D', '', filepath)
    pos = nx.spring_layout(G, weight="weight")
    degree_centrality = nx.degree_centrality(G)
    degree_centrality_color = np.fromiter(degree_centrality.values(), float)
    nx.draw(G, pos=pos, node_size=20, with_labels=False, width=1, alpha=0.6, node_color=degree_centrality_color)
    plt.title(f'Network {save_path}')
    plt.savefig(f'findings/peptide_graphs/network_{save_path}.jpg')
    return 

def degree_analysis(G, filepath):
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
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
    save_path = filepath.split('/')[-2] + '_' + re.sub('\D', '', filepath)
    plt.savefig(f'findings/peptide_graphs/degree_{save_path}.jpg')

def print_edge_histogram(edge_list, filepath):
    save_path = filepath.split('/')[-2] + '_' + re.sub('\D', '', filepath)
    x = edge_list['distance']
    plt.hist(x, density=True, bins=100)
    plt.title(f'Distance histogram {save_path}')
    plt.xlabel('distance')
    plt.savefig(f'findings/peptide_graphs/density_{save_path}.jpg')
    return

def subset_edge_list_with_threshold(edge_list, threshold):
    subset_edge_list = edge_list[ (edge_list['distance'] <= threshold)]
    return subset_edge_list
    
def subset_edge_list_with_protein(edge_list, proteins):
    subset_edge_list = edge_list[edge_list['accession'].isin(proteins)]
    return subset_edge_list
    

def main(args):
    filepath = args.filepath
    edge_list = pd.read_csv(filepath)
    
    proteins = ['HBB_HUMAN', 'FIBA_HUMAN']
    edge_list = subset_edge_list_with_protein(edge_list, proteins)
    edge_list = subset_edge_list_with_threshold(edge_list, 6)
    G = create_network_from_edge_list(edge_list)
    plot_network(G, filepath)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create edge list from file based on distance matrix')
    parser.add_argument("filepath", type=str, help="Path to file")
    args = parser.parse_args()
    main(args)