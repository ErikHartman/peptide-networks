# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 09:20:59 2021

@author: matti
"""

import urllib
from bs4 import BeautifulSoup
import pandas as pd
import tarfile
with tarfile.open('../data/pig_prot/sus-bar_annotation.tar.gz', "r:*") as tar:
    csv_path1 = tar.getnames()[1]
    csv_path2 = tar.getnames()[3]
    df1 = pd.read_csv(tar.extractfile(csv_path1), 
                      names=['name','hash','cluster','ontology',
                               'term id','description'], sep="\t")
    df2 = pd.read_csv(tar.extractfile(csv_path2), names=['name','hash','cluster','ontology',
                               'term id','description'], 
                         sep="\t")
    df = pd.concat([df1,df2])

print(df.shape)
df.drop_duplicates(subset='name', inplace=True)
print(df.shape)


import requests as r
import sys
from Bio import SeqIO
from io import StringIO

e = 0

def findseq(cID):
    try:
        if cID[:7] != 'ENSSSCP':
            baseUrl="http://www.uniprot.org/uniprot/"
            currentUrl=baseUrl+cID+".fasta?version=1"
            response = r.post(currentUrl)
            cData=''.join(response.text)
            
            Seq=StringIO(cData)
            return list(SeqIO.parse(Seq,'fasta'))[0].seq
        else:
            server = "https://rest.ensembl.org"
            ext = "/sequence/id/" + cID + "?"
            req = r.get(server+ext, headers={ "Content-Type" : "text/plain"})
             
            if not req.ok:
              print(cID)
              return ""
             
            return req.text
    except:
        print(cID)
        return ""
    
df['sequence'] = df['name'].apply(findseq)
df.drop(df['sequence'] == "", inplace=True)
print(df.shape)
df.to_csv('../data/pig_prot/all_sequences.csv')
print(e)