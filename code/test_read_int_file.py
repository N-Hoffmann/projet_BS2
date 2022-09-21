import pytest
import read_interaction_file as rif

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
    assert ('A', 'D') in result

def test_matrix():
    pass