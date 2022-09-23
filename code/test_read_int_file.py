import pytest
import numpy
import read_interaction_file as rif

read_interaction_file = pytest.importorskip("read_interaction_file")

test_file = '../example_files/toy_example.txt'

def test_dict():
    result = rif.read_interaction_file_dict(test_file)
    assert type(result) == dict
    assert len(result) == 6
    assert result == {'D': ['B', 'E', 'F'], 'B': ['A', 'C', 'D'], 'F': ['D'], \
        'A': ['B', 'C'], 'C': ['A', 'B'], 'E': ['D']}

def test_list():
    result = rif.read_interaction_file_list(test_file)
    assert type(result) == list
    assert len(result) == 6
    assert ('B', 'C') or ('C', 'B') in result

def test_matrix():
    matrix, node_list = rif.read_interaction_file_mat(test_file)
    expected_array = [[0, 1, 1, 0, 0, 0],
       [1, 0, 1, 1, 0, 0],
       [1, 1, 0, 0, 0, 0],
       [0, 1, 0, 0, 1, 1],
       [0, 0, 0, 1, 0, 0],
       [0, 0, 0, 1, 0, 0]]
    assert numpy.array_equal(matrix, expected_array) == True
    assert node_list == ['A','B','C','D','E','F']

def test_read_interaction_file():
    d_int, l_int, m_int, l_som = rif.read_interaction_file(test_file)
    assert d_int == {'D': ['B', 'E', 'F'], 'B': ['A', 'C', 'D'], 'F': ['D'], \
        'A': ['B', 'C'], 'C': ['A', 'B'], 'E': ['D']}
    assert l_int == [('A', 'B'), ('A', 'C'), ('B', 'C'), ('B', 'D'), ('D', 'E'), ('D', 'F')]
    expected_array = [[0, 1, 1, 0, 0, 0],
       [1, 0, 1, 1, 0, 0],
       [1, 1, 0, 0, 0, 0],
       [0, 1, 0, 0, 1, 1],
       [0, 0, 0, 1, 0, 0],
       [0, 0, 0, 1, 0, 0]]
    assert numpy.array_equal(m_int, expected_array) == True
    assert l_som == ['A', 'B', 'C', 'D', 'E', 'F']

def test_is_interaction_file():
    assert rif.is_interaction_file(test_file) == True
    assert rif.is_interaction_file("../example_files/empty_file.csv") == False
    assert rif.is_interaction_file("../example_files/wrong_first_line.csv") == False
    assert rif.is_interaction_file("../example_files/wrong_columns.csv") == False