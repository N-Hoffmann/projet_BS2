import pytest
import os
from pathlib import Path
import read_interaction_file as rif
import explore_file as ef

explore_file = pytest.importorskip("explore_file")
test_file = Path(__file__).parent / '../example_files/toy_example.txt'

def test_count_vertices():
    assert ef.count_vertices(test_file) == 6

def test_count_edges():
    assert ef.count_edges(test_file) == 6

def test_clean_interactome():
    # Find a better way to create temp fileout.txt ?
    # Find a better way to assert ? tried filecmp
    ef.clean_interactome(Path(__file__).parent / "../example_files/toy_example_dirty.txt","test_out.txt")
    assert rif.read_interaction_file_list("../example_files/toy_example.txt") \
        == rif.read_interaction_file_list("test_out.txt") 
    os.remove("test_out.txt")

def test_get_degree():
    assert ef.get_degree(test_file,"D") == 3
    assert ef.get_degree(test_file,"A") == 2
    assert ef.get_degree(test_file,"E") == 1

def test_get_max_degree():
    assert ef.get_max_degree(test_file) == (['B', 'D'], 3) or (['D', 'B'], 3)

def test_get_ave_degree():
    assert ef.get_ave_degree(test_file) == 2.0

def test_count_degree():
    assert ef.count_degree(test_file,1) == 2
    assert ef.count_degree(test_file,2) == 2
    assert ef.count_degree(test_file,3) == 2

def test_histogram_degree():
    assert ef.histogram_degree(test_file,1,4) == {1: 2, 2: 2, 3: 2, 4: 0}