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

Notes: There are **three** main ways of computing the distance in a peptide graph: string distance (e.g., levenshtein distance), evolutionary substitutional matrices (e.g., BLOSUM62) or biophysical properties.

Discuss: biophysical properties, creating edge lists (make it quicker), what to look for in networks, layout for graph (make class?)

---

<code>python create_edge_lists.py -filepath -matrix [matrix choices: blosum, levenshtein, biophysical]</code>
<br>
<code>python create_networks.py -filepath (for edge list)</code>

Rn the graphs look pretty messy...

![example graph](/findings/peptide_graphs/network_blosum_34.jpg)
