import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
import matplotlib.pyplot as plt
import os
from PeptideNetwork import PeptideNetwork



def main(args):
    """
    Here we want to 
    
    for each edge list
        for threshold in range
            create networks
            get measure
        
    plot results
    """
    dirpath = args.dirpath
    results_dict = {'file':[], 'nr':[]}
    for file in tqdm(os.listdir(dirpath)):
        edge_list = pd.read_csv(dirpath+file)
        G = PeptideNetwork(edge_list)
        nr = []
        for threshold in reversed(range(1, 12)):
            G.subset_edge_list_with_threshold(threshold)
            G.create_network()
            measure = len(G.get_connected_components())
            print(threshold, measure)
            nr.append(measure)
        results_dict['file'].append(file)
        results_dict['nr'].append(nr)

    fig, axs = plt.subplots(3,2, figsize=(10,10))
    dict_size = len(results_dict['file'])
    for i, ax in zip(range(dict_size), axs.ravel()):
        file = results_dict['file'][i]
        y = results_dict['nr'][i]
        ax.plot(y)
        ax.set_ylabel('nr subgraphs')
        ax.set_xlabel('threshold')
        ax.set_title(file)
    plt.tight_layout()
    plt.savefig('subgraphs_levenshtein.jpg')
    
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create edge list from file based on distance matrix')
    parser.add_argument("dirpath", type=str, help="Path to dir")
    args = parser.parse_args()
    main(args)