import pandas as pd
import numpy as np
import networkx as nx
import argparse
import matplotlib.pyplot as plt


def create_network_from_edge_list(edge_list):
    G = nx.from_pandas_edgelist(
    edge_list,
    edge_attr=["distance"],
    create_using=nx.MultiGraph())
    return G

def print_edge_histogram(edge_list):
    x = edge_list['distance']
    plt.hist(x, density=True, bins=30)
    plt.show()
    return

def subset_edge_list(edge_list):
    return edge_list
    
    

def main(args):
    filepath = args.filepath
    edge_list = pd.read_excel(filepath)
    
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create edge list from file based on distance matrix')
    parser.add_argument("filepath", type=str, help="Path to file")
    args = parser.parse_args()
    main()