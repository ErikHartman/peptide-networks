# -*- coding: utf-8 -*-
"""
Created on Fri Dec 24 13:15:07 2021

@author: matti
"""

from hash_util import PeptideSequence, import_clean_data
import pandas as pd
import time
import glob

class ParentString(object):
    def __init__(self, loc):
        self.loc = loc
        self.children = []
    
    def addChild(self,child):
        self.children.append(child)
    
    def getChildren(self):
        return self.children 
    
    def getLoc(self):
        return self.loc

class GroupSubPeptides(object):
    def __init__(self, path):
        self.parents = []
        self.df = import_clean_data(path)
        
        sizes = self.df.index.get_level_values(0).unique()
        self.hash_map = {key: {} for key in sizes}
        self.leftedge_chain_map = {key: {} for key in range(7,max(sizes) + 1)}
        self.rightedge_chain_map = {key: {} for key in range(7,max(sizes) + 1)}
        
        for index, row in self.df.iloc[::-1].iterrows():
            if not (self.hash_map.get(index[0]).get(index[1]) is None):
                self.hash_map[index[0]][index[1]].addChild(index)
            else:
                parent = ParentString(index)
                
                if self.edgechains(row['Peptide'], parent):
                    self.parents.append(parent)
                    
                for size in [x for x in sizes if x < index[0]]:
                    self.hash_all_of_size(row['Peptide'], size, parent)
                
    
    def hash_all_of_size(self, sequence, size, parent):
        hashstring = PeptideSequence('')
        index = 0
        sequence_len = len(sequence)
        while index < size:
            hashstring.insert(sequence[index])
            index += 1
        while index < sequence_len:
            self.add_hash(hashstring.getHash(), size, parent)
            index += 1
    
    def add_hash(self, hash_value, size, parent):
        self.hash_map[size][hash_value] = parent
        
    def edgechains(self, sequence, parent):
        hashstring = PeptideSequence(sequence[:7])
        size = 7
        return_value = True
        
        while size < len(sequence):
            hashstring.insert(sequence[size])
            size += 1
            return_value = return_value and self.add_edge(hashstring.getHash(), size, parent, left=True)
        while size > 7:
            hashstring.pop_front()
            size -= 1
            return_value = return_value and self.add_edge(hashstring.getHash(), size, parent, left=False)
            
        return return_value
    
    def add_edge(self, hash_value, size, parent, left = True):
        check_map = self.rightedge_chain_map
        update_map = self.leftedge_chain_map
        if not left:
            check_map, update_map = update_map, check_map
            
        if not (check_map.get(size).get(hash_value) is None):
            check_map[size][hash_value].addChild(parent)
            return False
        elif update_map.get(size).get(hash_value) is None:
            update_map[size][hash_value] = parent
            return True
            
    
if __name__ == '__main__':
    n_parents = []
    n_daughters = []
    n_peptides = []
    nms = []
    t1 = time.time()
    
    for path in glob.glob('../../data/clean_data/*/*'):
        grouped = GroupSubPeptides(path)
        print(grouped.df.shape)
        n_parents.append(len(grouped.parents))
        n_daughters.append(grouped.df.shape[0] - len(grouped.parents))
        n_peptides.append(grouped.df.shape[0])
        nms.append(path.split('/')[0])
        print(f'Done in {time.time() - t1}')
        t1 = time.time()
    
    datafrm = pd.DataFrame({'file': nms, 
                            'n_peptides':n_peptides, 
                        'n_daughters' : n_daughters, 
                        'n_parents': n_parents})
    print(datafrm.head())
    datafrm.to_csv('parent_specs.csv')