import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
import matplotlib.pyplot as plt
import os
from PeptideNetwork import PeptideNetwork

def plot_community_positions(results_dict, threshold, save_path):
    plt.clf()
    fig, axs = plt.subplots(3,2, figsize=(10,10))
    dict_size = len(results_dict['file'])
    for i, ax in zip(range(0, dict_size), axs.ravel()):
        starts = results_dict['starts'][i]
        ends = results_dict['ends'][i]
        file = results_dict['file'][i]
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
        ax.bar(x=df['x'].values,height = df['size_mean'].values,width=df['width'].values , xerr = df['x_sd'].values, alpha=0.5, color='green') 
        ax.set_xlabel('Sequence')
        ax.set_ylabel('sqrt(peptides in community)')
        ax.set_title(file)
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
    
def plot_degree_analysis(results_dict, save_path):
    """ Saves simple degree analysis to the given save path """
    plt.clf()
    dict_size = len(results_dict['file'])
    fig, axs= plt.subplots(dict_size,2, figsize=(12, 12))
    
    for i in range(0,dict_size):
        degree_sequence = results_dict['degree_sequences'][i]
        file = results_dict['file'][i]
        axs[i,0].plot(degree_sequence, "b-", marker="o")
        axs[i,0].set_title("Degree Rank Plot" + file)
        axs[i,0].set_ylabel("Degree")
        axs[i,0].set_xlabel("Rank")
        axs[i,1].bar(*np.unique(degree_sequence, return_counts=True))
        axs[i,1].set_title("Degree histogram" + file)
        axs[i,1].set_xlabel("Degree")
        axs[i,1].set_ylabel("# of Nodes")
    plt.tight_layout()
    plt.savefig(f'{save_path}')
    
def plot_connected_component_sizes(results_dict, save_path):
    plt.clf()
    fig, axs = plt.subplots(3,2, figsize=(10,10))
    dict_size = len(results_dict['file'])
    for i, ax in zip(range(0, dict_size), axs.ravel()):
        connected_components = results_dict['connected_component_sizes'][i]
        sizes = [len(cc) for cc in connected_components]
        x=range(len(sizes))
        ax.bar(x=x, height=sizes, label = max(sizes))
        ax.legend()
        ax.set_title(results_dict['file'][i])
    plt.suptitle('connected component size distribution')
    plt.savefig(f'{save_path}')

    
def main(args):
    """
    Here we want to:
    results_dict = {}
    for each edge list:
        Create network
        subset
        generate results
        save in results dict
    plot results
    """
    dirpath = args.dirpath
    proteins = ['HBB_HUMAN']
    threshold = 3
    results_dict = {'file':[], 'connected_component_sizes':[], 'degree_sequences' :[], 'starts':[], 'ends':[]}
    for file in tqdm(os.listdir(dirpath)):
        edge_list = pd.read_csv(dirpath+file)
        G = PeptideNetwork(edge_list)
        G.subset_edge_list_with_protein(proteins)
        G.subset_edge_list_with_threshold(threshold)
        G.create_network()
        cc = G.get_connected_components()
        degree_sequence = G.get_degree_sequence()
        starts, ends = G.get_community_positions(proteins[0])
        results_dict['file'].append(file)
        results_dict['connected_component_sizes'].append(cc)
        results_dict['degree_sequences'].append(degree_sequence)
        results_dict['starts'].append(starts)
        results_dict['ends'].append(ends)
    #plot_degree_analysis(results_dict, 'test')
    plot_community_positions(results_dict,2, 'test.jpg')
    #plot_connected_component_sizes(results_dict, 'FIBA_levenshtein.jpg')
    
    #plot_distance_histogram(edge_list, 'findings/peptide_graphs/biophysical_34.jpg')
    #plot_degree_analysis(G.get_degree_sequence(), 'findings/peptide_graphs/degree_biophysical_34.jpg')
    
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create edge list from file based on distance matrix')
    parser.add_argument("dirpath", type=str, help="Path to dir")
    args = parser.parse_args()
    main(args)