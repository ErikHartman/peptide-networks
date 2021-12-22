from attr import has
import pytest
from peptide_networks.preprocessing import hash_util

@pytest.mark.parametrize("string_in, expected_result", [
    ("", 0),
    ("A", 0),
    ("R", 1),
    ("AA", 0),
    ("RA", 20),
    ("RAAAAAAAAAAAAAAA", 14321255926290448385),
])
def test_HashString_init(string_in, expected_result):
    hash_string = hash_util.HashString(string_in)
    assert hash_string.getHash() == expected_result

def test_hashString_init2():
    assert hash_util.HashString().getHash() == 0

@pytest.mark.parametrize("inp, expected_result", [
    (0, 1),
    (1, 2),
    (3, 8),
    (64, 1),
])
def test_HashString_po(inp, expected_result):
    hs = hash_util.HashString()
    assert hs._po(2, inp) == expected_result

@pytest.mark.parametrize("inp, result", [
    ('RY', True),
    ('YR', False),
    ('RYR', False),
    ('YRY', False),
    ('', False)
])
def test_HashString_eq(inp, result):
    assert (hash_util.HashString('RY') == hash_util.HashString(inp)) is result

def test_HashString_insert():
    hs = hash_util.HashString('W')
    hs.insert('R')
    assert hs.getString() == 'WR' and hs.getSize() == 2 and hs == hash_util.HashString('WR')

def test_HashString_pop():
    hs1 = hash_util.HashString('QR')
    hs2 = hash_util.HashString('WQR')
    hs3 = hash_util.HashString('QRT')
    hs2.pop_front()
    hs3.pop_back()
    assert (hs2.getSize() == 2) and (hs3.getSize() == 2) and (hs1 == hs2) and (hs1 == hs3)