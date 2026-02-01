"""Unit tests for deli.cube.count module"""

import pytest

from deli.decode.count import (
    _build_graph,
    _build_graph_semi_directional,
    _count_connected_components,
    _get_dna_hamming_neighbors,
    corrected_count,
    dedup_count,
    raw_count,
)


def test_get_dna_hamming_neighbors():
    """Test _get_dna_hamming_neighbors function."""
    seq = "ACGT"
    neighbors = _get_dna_hamming_neighbors(seq)
    expected_neighbors = [
        "CCGT",
        "GCGT",
        "TCGT",  # Change at position 0
        "AAGT",
        "AGGT",
        "ATGT",  # Change at position 1
        "ACAT",
        "ACCT",
        "ACTT",  # Change at position 2
        "ACGA",
        "ACGC",
        "ACGG",  # Change at position 3
    ]
    assert set(neighbors) == set(expected_neighbors)


def test_build_graph():
    """Test _build_graph function."""
    umis = {"ACGT": 10, "CCGT": 5, "GCGT": 3}
    graph = _build_graph(umis)
    expected_graph = {"ACGT": ["CCGT", "GCGT"], "CCGT": ["ACGT", "GCGT"], "GCGT": ["ACGT", "CCGT"]}

    for key in expected_graph:
        assert set(graph[key]) == set(expected_graph[key])

    # a harder one
    umis = {"ATAGATC": 10, "AGATAGC": 5, "AGATAGT": 3, "AGATAAT": 2, "ATCGATC": 6, "CCCCCCC": 1}
    graph = _build_graph(umis)
    expected_graph = {
        "ATAGATC": ["ATCGATC"],
        "AGATAGC": ["AGATAGT"],
        "AGATAGT": ["AGATAGC", "AGATAAT"],
        "AGATAAT": ["AGATAGT"],
        "ATCGATC": ["ATAGATC"],
        "CCCCCCC": [],
    }

    for key in expected_graph:
        assert set(graph[key]) == set(expected_graph[key])


def test_build_graph_semi_directional():
    """Test _build_graph_semi_directional function."""
    umis = {"ACGT": 10, "CCGT": 5, "GCGT": 3}
    graph = _build_graph_semi_directional(umis, degree=2)
    expected_graph = {"ACGT": ["GCGT"], "CCGT": [], "GCGT": ["ACGT"]}

    for key in expected_graph:
        assert set(graph[key]) == set(expected_graph[key])

    # a harder one
    umis = {"ATAGATC": 10, "AGATAGC": 5, "AGATAGT": 3, "AGATAAT": 1, "ATCGATC": 2, "CCCCCCC": 1}
    graph = _build_graph_semi_directional(umis)
    expected_graph = {
        "ATAGATC": ["ATCGATC"],
        "AGATAGC": [],
        "AGATAGT": ["AGATAAT"],
        "AGATAAT": ["AGATAGT"],
        "ATCGATC": ["ATAGATC"],
        "CCCCCCC": [],
    }

    for key in expected_graph:
        assert set(graph[key]) == set(expected_graph[key])


def test_count_connected_components():
    """Test _count_connected_components function."""
    graph = {"ACGT": ["CCGT"], "CCGT": ["ACGT"], "GCGT": []}
    count = _count_connected_components(graph)
    assert count == 2


def test_raw_count():
    """Test raw_count function."""
    umis = {"ACGT": 10, "CCGT": 5, "GCGT": 3}
    assert raw_count(umis) == 18


def test_dedup_count():
    """Test dedup_count function."""
    umis = {"ACGT": 10, "CCGT": 5, "GCGT": 3}
    assert dedup_count(umis) == 3


def test_corrected_count_cluster():
    """Test corrected_count function with 'cluster' method."""
    umis = {"ACGT": 10, "CCGT": 5, "GCGT": 3}
    count = corrected_count(umis, method="cluster")
    assert count == 1


def test_corrected_count_directional():
    """Test corrected_count function with 'directional' method."""
    umis = {"ACGT": 10, "CCGT": 5, "GCGT": 3}
    count = corrected_count(umis, method="directional")
    assert count == 2


def test_corrected_count_invalid_method():
    """Test corrected_count function with invalid method."""
    umis = {"ACGT": 10, "CCGT": 5, "GCGT": 3}
    with pytest.raises(ValueError):
        corrected_count(umis, method="invalid")
