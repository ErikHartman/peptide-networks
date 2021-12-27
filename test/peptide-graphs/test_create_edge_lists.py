import pytest
from peptide_networks.peptide_graphs import create_edge_list



def test_levenshtein_distance():
    seq1 = "AAA"
    seq2 = "AAA"
    seq3 = "AAB"
    seq4 = "ABC"
    assert create_edge_list.levenshtein_distance(seq1, seq2) == 0
    assert create_edge_list.levenshtein_distance(seq3,seq4) == 2

def test_blosum_distance():
    seq1 = "AAA"
    seq2 = "AAG"
    assert create_edge_list.blosum_distance(seq1, seq2) == 8.0