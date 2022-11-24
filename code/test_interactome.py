import pytest
from pathlib import Path
import interactome
import matplotlib.pyplot as plt

test_file = Path(__file__).parent / '../example_files/toy_example.txt'
test_ppi = interactome.Interactome(test_file)

def test_class():
    assert isinstance(test_ppi, interactome.Interactome)

def test_attributes():
    assert test_ppi.dict == {'A': ['B', 'C'], 'B': ['A', 'C', 'D'], 'C': ['A', 'B'], \
        'D': ['B', 'E', 'F'], 'E': ['D'], 'F': ['D']}
    assert test_ppi.list == [('A', 'B'), ('A','C'), ('B', 'C'), ('B', 'D'), ('D', 'E'), ('D', 'F')]
    assert test_ppi.protein == ['A', 'B', 'C', 'D', 'E', 'F']

def test_methods():
    assert test_ppi.count_edges() == 6
    assert test_ppi.count_vertices() == 6
    assert test_ppi.get_degree("B") == 3
    assert test_ppi.get_degree("E") == 1
    assert test_ppi.get_max_degree() == (["B","D"], 3)
    assert test_ppi.get_ave_degree() == 2.0
    assert test_ppi.count_degree(2) == 2
    assert test_ppi.density() == 0.4
    assert test_ppi.clustering("B") == 0.3333333333333333
    assert test_ppi.clustering("F") == 0
    assert test_ppi.count_cc() == (1, [6])

def test_cc():
    test_file = Path(__file__).parent / '../example_files/toy2.txt'
    test_ppi = interactome.Interactome(test_file)
    assert test_ppi.count_cc() == (3, [6, 2, 3])
    assert test_ppi.find_cc() == [['A', 'G', 'C', 'D', 'E', 'F'], ['B', 'H'], ['I', 'J', 'K']]
    assert test_ppi.extract_cc("B") == ['B', 'H']