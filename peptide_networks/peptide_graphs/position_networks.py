from unittest import result
import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
import matplotlib.pyplot as plt
import os
from PeptideNetwork import PeptideNetwork
import networkx as nx


def main(args):
    """
    Create master networks and get pos
    create networks for each edge list by subsetting the master network
    
    plot networks
    
    """
    dirpath = args.dirpath
    edge_list = pd.read_csv(dirpath + os.listdir(dirpath)[0])
    protein = ['FIBA_HUMAN']
    threshold = 3
    
    master_network = PeptideNetwork(edge_list)
    master_network.subset_edge_list_with_protein(protein)
    master_network.subset_edge_list_with_threshold(threshold)
    networks = {'file':[], 'G':[]}
    networks['file'].append(os.listdir(dirpath)[0])
    networks['G'].append(master_network)
    for file in tqdm(os.listdir(dirpath)[1:]):
        edge_list = pd.read_csv(dirpath + file)
        G = PeptideNetwork(edge_list)
        G.subset_edge_list_with_protein(protein)
        G.subset_edge_list_with_threshold(threshold)
        edge_list = G.get_edge_list()
        master_network.concatenate_edgelist(edge_list)
        networks['file'].append(file)
        networks['G'].append(G)
    
    master_network.create_network()
    print(master_network.get_edge_list())
    master_pos = nx.spring_layout(master_network.get_network(), weight='weight') 
    
    """
    Create an empty copy of the master network.
    Create a network for a file
    
    If the file network contains an edge which is not in empty copy
        insert edge
    """
    fig, axs  = plt.subplots(2,3, figsize=(20,20))
    for i, ax in zip(range(len(networks['file'])), axs.ravel()):
        G_file = networks['G'][i]
        print(networks['file'][i])
        G_file.create_network()
        G_file.toString()
        G_file = G_file.get_network()
        G_ec = nx.create_empty_copy(master_network.get_network(), with_data=True)
        added_edges = 0
        for edge in G_file.edges(data=True):
                e = edge[0:2]
                weight = edge[2]
                if not G_ec.has_edge(*e):
                    added_edges += 1
                    G_ec.add_edge(*e, weight=weight['weight'])
        print("added edges: " , added_edges)
        node_colors = []
        node_sizes = []
        nr_of_nodes = 0
        for dc in nx.degree_centrality(G_file).values():
                node_colors.append(dc)
        for node in G_ec.nodes():
                if len(list(G_ec.edges(node))) > 0:
                        
                        node_sizes.append(0.5)
                        nr_of_nodes += 1
 
        print("added nodes: " , nr_of_nodes)
        options = {
                "edge_cmap": plt.cm.winter,
                "with_labels": False,
                "cmap":plt.cm.hot,
                "node_color": node_colors,
                "alpha":0.8,
                "linewidths":0.5,
                "width":0.7,
                "node_size":15,

                }
        nx.draw(G_file, pos=master_pos, ax=ax, **options)
        ax.set_title(networks['file'][i])
        ax.annotate(f'Added edges: {added_edges}', xy=(0.7, 0.1), xycoords='axes fraction')
        ax.annotate(f'Number of nodes: {nr_of_nodes}', xy=(0.7, 0.15), xycoords='axes fraction')
    plt.savefig('test.jpg')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create edge list from file based on distance matrix')
    parser.add_argument("dirpath", type=str, help="Path to dir")
    args = parser.parse_args()
    main(args)