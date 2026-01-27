"""Functions for handling barcode counts based on UMI sequences"""

from typing import Literal

from numba import njit


@njit
def _get_dna_hamming_neighbors(seq) -> set[str]:
    """
    Given a DNA sequence, return all sequences with Hamming distance of 1

    Parameters
    ----------
    seq: str

    Returns
    -------
    list[str]
    """
    seqs = set()
    for i, let in enumerate(seq):
        for new in ["A", "G", "C", "T"]:
            if new == let:
                continue
            seqs.add(seq[:i] + new + seq[i + 1 :])
    return seqs


def _build_graph(umis: dict[str, int]) -> dict[str, list[str]]:
    """
    Build a graph of UMIs

    There will be one Node for each unique UMI that stores the number of times it
    was observed. Edges are drawn if the hamming distance between two UMIs is 1

    Notes
    -----
    Expected input is the dict of all unique UMIs observed for a given compound ID

    Parameters
    ----------
    umis: dict[str, int]
        Map of UMI sequence to its count

    Returns
    -------
    dict[str, list[str]]
        Adjacency representation of the UMI graph, keys are the node, values are the list of neighboring nodes
    """
    graph: dict[str, list[str]] = {umi: [] for umi in umis}
    known_umis = tuple(umis.keys())
    for i, umi in enumerate(umis):
        neighbors = _get_dna_hamming_neighbors(umi)
        for known_umi in known_umis[:i]:
            if known_umi in neighbors:
                graph[umi].append(known_umi)
                graph[known_umi].append(umi)
    return graph


def _build_graph_semi_directional(umis: dict[str, int], degree: int = 2) -> dict[str, list[str]]:
    """
    Build a "semi-directional" graph of UMIs

    There will be one Node for each unique UMI that stores the number of times it
    was observed. Edges are drawn only drawn if the hamming distance between two UMIs is 1
    AND the count of the source UMI is at least <degree> times that of the target UMI.

    "Semi-directional" means that edges are undirectional, but the edge condition itself if directional.
    If either possible direction meets the condition, the edge is drawn for both directions.
    The result is an undirected graph.

    Notes
    -----
    Expected input is the dict of all unique UMIs observed for a given compound ID

    Parameters
    ----------
    umis: dict[str, int]
        Map of UMI sequence to its count
    degree: int
        The degree of reduction in counts required to draw a directed edge.
        Edge is only drawn if count(source) >= (degree * count(target)) + 1

    Returns
    -------
    dict[str, list[str]]
        Adjacency representation of the UMI graph, keys are the node, values are the list of neighboring nodes
    """
    graph: dict[str, list[str]] = {umi: [] for umi in umis}
    known_umis = tuple(umis.keys())
    for i, umi in enumerate(umis):
        neighbors = _get_dna_hamming_neighbors(umi)
        for known_umi in known_umis[:i]:
            if known_umi in neighbors:
                if (umis[umi] >= (degree * umis[known_umi]) + 1) or (umis[known_umi] >= (degree * umis[umi]) + 1):
                    graph[umi].append(known_umi)
                    graph[known_umi].append(umi)
    return graph


def _count_connected_components(graph: dict[str, list[str]]) -> int:
    """
    Count the number of weakly connected components in a directed graph using DFS.

    Parameters
    ----------
    graph : dict[str, list[str]]
        Adjacency list representation of the graph.

    Returns
    -------
    int
        The number of connected components.
    """
    visited = set()
    count = 0
    for node in graph:
        if node not in visited:
            count += 1
            stack = [node]
            visited.add(node)
            while stack:
                current = stack.pop()
                for neighbor in graph[current]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)
    return count


def raw_count(umis: dict[str, int]) -> int:
    """
    Get the raw count for a set of UMIs

    'Raw count' is simply the total number of UMIs observed, including duplicates

    Parameters
    ----------
    umis: dict[str, int]
        Map of UMI sequence to its count

    Returns
    -------
    int
        The raw count of UMIs
    """
    return sum(umis.values())


def dedup_count(umis: dict[str, int]) -> int:
    """
    Get the dedup count for a set of UMIs

    'Dedup count' is simply the number of unique UMIs observed

    Parameters
    ----------
    umis: dict[str, int]
        Map of UMI sequence to its count

    Returns
    -------
    int
        The degenerate count of unique UMIs
    """
    return len(umis)


def corrected_count(umis: dict[str, int], method: Literal["cluster", "directional"] = "directional") -> int:
    """
    Get the error-corrected count for a set of UMIs

    'Error-corrected count' is the number of unique UMIs observed after
    collapsing UMIs that are likely sequencing / PCR errors

    Notes
    -----
    This works by implementing a set-up similar to UMI Tools "directional graph" approach [1]_

    References
    ----------
    .. [1] Smith, T., Heger, A., & Sudbery, I. (2017). UMI-tools: modeling sequencing errors in Unique Molecular
       Identifiers to improve quantification accuracy. Genome Research, 27(3), 491–499.
       https://doi.org/10.1101/gr.209601.116

    Parameters
    ----------
    umis: dict[str, int]
        Map of UMI sequence to its count
    method: Literal["cluster", "directional"], default "directional"
        The method to use for building the UMI graph.
        - `cluster`: undirected edges between UMIs with Hamming distance of 1
        - `directional`: directed edges between UMIs with Hamming distance of 1
                         only if count(source) >= 2 * count(target) - 1

    Returns
    -------
    int
        The error-corrected count of UMIs
    """
    if method == "directional":
        graph = _build_graph_semi_directional(umis)  # get the graph of UMIs
    elif method == "cluster":
        graph = _build_graph(umis)  # get the graph of UMIs
    else:
        raise ValueError(f"Unrecognized method '{method}' for corrected_count. Options are ['cluster', 'directional']")

    return _count_connected_components(graph)
