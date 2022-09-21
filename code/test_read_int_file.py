import pytest
import read_interaction_file as rif
import numpy

read_interaction_file = pytest.importorskip("read_interaction_file")

test_file = '../data/example.csv'

def test_dict():
    result = rif.read_interaction_file_dict(test_file)
    assert type(result) == dict
    assert len(result) == 9
    assert result["A"] == ['D', 'B', 'C', 'F']
    assert result["E"] == ['F', 'D']

def test_list():
    result = rif.read_interaction_file_list(test_file)
    assert type(result) == list
    assert len(result) == 11
    assert ('A', 'D') or ('D', 'A') in result

def test_matrix():
    matrix, node_list = rif.read_interaction_file_mat(test_file)
    expected_array = [[0, 1, 1, 1, 0, 1, 0, 0, 0],
       [1, 0, 1, 0, 0, 0, 0, 1, 0],
       [1, 1, 0, 0, 0, 0, 1, 0, 1],
       [1, 0, 0, 0, 1, 0, 0, 1, 0],
       [0, 0, 0, 1, 0, 1, 0, 0, 0],
       [1, 0, 0, 0, 1, 0, 0, 0, 0],
       [0, 0, 1, 0, 0, 0, 0, 0, 0],
       [0, 1, 0, 1, 0, 0, 0, 0, 0],
       [0, 0, 1, 0, 0, 0, 0, 0, 0]]
    assert numpy.array_equal(matrix, expected_array) == True
    assert node_list == ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'X','Y']

def test_read_interaction_file():
    d_int, l_int, m_int, l_som = rif.read_interaction_file(test_file)

def test_is_interaction_file():
    assert rif.is_interaction_file(test_file) == True
    assert rif.is_interaction_file("../data/empty_file.csv") == False
    assert rif.is_interaction_file("../data/wrong_first_line.csv") == False
    assert rif.is_interaction_file("../data/wrong_columns.csv") == False