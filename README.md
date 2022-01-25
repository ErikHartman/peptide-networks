# peptide-networks

## Installing python package

To run pytest from other folder one needs to run:

pip install -e .

in source folder.

---

## Folder structure

- /data
- /peptide_networks
  - /experiment 1
  - /experiment 2
- /test
  - /experiment 1
  - /experiment 2
- /findings
  - /experiment 1
  - /experiment 2

---

# Status

## Mattias

1. Proteasprediktion
   Lite frygt 30 % av alla peptider är direkta derivat av andra peptider i populationen.
   Ca 70 % av alla peptider kan härledas till andra peptider.

## Erik

---

Notes:

- There are **three** main ways of computing the distance in a peptide graph: string distance (e.g., levenshtein distance), evolutionary substitutional matrices (e.g., BLOSUM62) or biophysical properties.
- Need to find a way to compute the optimal thresholding and optimal way of creating edge-lists.
- I think BLOSUM and substitional matrices are bad. Let's go with biophysical and levenshtein.

Discuss: biophysical properties, what to look for in networks, saving findings.

Att göra:

- Börja skapa nätverk från subsets.
- Kvantifiera skillnader mellan inf och ninf (centrality measures, connectivity, assortativity etc)

---

<code>python create_edge_lists.py -filepath -matrix [matrix choices: blosum, levenshtein, biophysical]</code>
<br>
<code>python create_networks.py -filepath (for edge list)</code>

Graph for positional partitioning works fairly well.

![example graph](/findings/peptide_graphs/lev_HBB_31_network.jpg)
![example graph](/findings/peptide_graphs/lev_HBB_31_position.jpg)
