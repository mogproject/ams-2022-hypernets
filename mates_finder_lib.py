import numpy as np
import hypernetx as hnx
import networkx as nx
from networkx.algorithms.bipartite.matrix import from_biadjacency_matrix
import scipy
from collections import defaultdict
from enum import Enum
from typing import Mapping, List, Tuple
import matplotlib.pyplot as plt


class MateType(Enum):
    CLIQUE_MATE = 0
    LINE_MATE = 1
    GRAM_MATE = 2


def weighted_clique_expansion(A: np.matrix) -> np.matrix:
    return A @ A.T


def weighted_line_graph(A: np.matrix) -> np.matrix:
    return A.T @ A


def get_bit(x: int, i: int) -> int:
    return (x >> i) & 1


def get_adj(n: int, m: int, index: int) -> np.matrix:
    # (n,m,index) maps to a unique n x m matrix
    assert n * m < 64 and 0 <= index < (1 << n * m)

    return np.matrix([[get_bit(index, j * n + i) for j in range(m)] for i in range(n)])


def next_perm(n: int, m: int, index: int) -> int:
    if index + 1 == 1 << (n * m):
        return -1  # no next permutation

    # columns
    cols = [(index >> (n * j)) % (1 << n) for j in range(m)]
    j = m - 1
    while j >= 1:
        if cols[j - 1] != cols[j]:
            break
        index -= cols[j] << (n * j)
        j -= 1
    return index + (1 << (n * j))


def m2s(A: np.matrix) -> str:
    """Converts a matrix to a string."""
    # Assumption: each element is one-digit integer
    return ' '.join(''.join(str(A[i, j]) for j in range(A.shape[1])) for i in range(A.shape[0]))


def is_isomorphic(H1: np.matrix, H2: np.matrix) -> bool:
    """
    Checks if two hypergraphs are isomorphic after removing isolates.
    """
    G1 = from_biadjacency_matrix(scipy.sparse.csr_matrix(H1))
    G2 = from_biadjacency_matrix(scipy.sparse.csr_matrix(H2))
    G1.remove_nodes_from(list(nx.isolates(G1)))
    G2.remove_nodes_from(list(nx.isolates(G2)))
    return nx.faster_could_be_isomorphic(G1, G2) and nx.is_isomorphic(G1, G2)


def get_x_mates(n: int, m: int, mate_type: MateType) -> Mapping[str, List[Tuple[int, int, int]]]:
    ret = defaultdict(list)
    index = 0

    while index >= 0:
        A = get_adj(n, m, index)

        if mate_type == MateType.CLIQUE_MATE:
            key = m2s(weighted_clique_expansion(A))
        elif mate_type == MateType.LINE_MATE:
            key = m2s(weighted_line_graph(A))
        elif mate_type == MateType.GRAM_MATE:
            key = m2s(weighted_clique_expansion(A)) + \
                ',' + m2s(weighted_line_graph(A))
        else:
            raise ValueError('unknown mate type')

        # isomorphic to existing hypergraphs?
        if not any(is_isomorphic(A, get_adj(nn, mm, idx)) for nn, mm, idx in ret[key]):
            ret[key] += [(n, m, index)]

        index = next_perm(n, m, index)
    return ret


def get_cliquemates(n: int, m_max: int) -> Mapping[str, List[Tuple[int, int, int]]]:
    return get_x_mates(n, m_max, MateType.CLIQUE_MATE)


def get_linemates(n_max: int, m: int) -> Mapping[str, List[Tuple[int, int, int]]]:
    return get_x_mates(n_max, m, MateType.LINE_MATE)


def get_grammates(n: int, m: int) -> Mapping[str, List[Tuple[int, int, int]]]:
    return get_x_mates(n, m, MateType.GRAM_MATE)


def get_largest_cliquemates(n: int, m_max: int) -> List[Tuple[int, int, int]]:
    d = get_cliquemates(n, m_max)
    return max(d.values(), key=len)


def get_largest_linemates(n_max: int, m: int) -> List[Tuple[int, int, int]]:
    d = get_linemates(n_max, m)
    return max(d.values(), key=len)


def get_largest_grammates(n: int, m: int) -> List[Tuple[int, int, int]]:
    d = get_grammates(n, m)
    return max(d.values(), key=len)


def print_nontrivial_cliquemates(n: int, m_max: int, threshold: int = 2) -> None:
    d = get_cliquemates(n, m_max)
    for xs in d.values():
        if len(xs) >= threshold:
            print(xs)


def print_nontrivial_linemates(n_max: int, m: int, threshold: int = 2) -> None:
    d = get_cliquemates(n_max, m)
    for xs in d.values():
        if len(xs) >= threshold:
            print(xs)


def print_nontrivial_grammates(n: int, m: int, threshold: int = 2) -> None:
    d = get_grammates(n, m)
    for xs in d.values():
        if len(xs) >= threshold:
            print(xs)


def visualize_adj(A: np.matrix):
    a = {}
    for i in range(A.shape[1]):
        s = np.nonzero(A[:, i])[0].tolist()
        if s:
            a[f'e{i}'] = s
    H = hnx.Hypergraph(a)
    hnx.drawing.draw(H)
