"""Optimal one-to-one assignment (Hungarian / Kuhn-Munkres, O(n^3)).

Implemented here rather than pulled from SciPy so the pipeline has no heavy
numeric dependency. `max_weight_matching` takes a rectangular profit matrix and
returns the row->column pairing that maximises the total score, with every row
and every column used at most once.
"""
from __future__ import annotations

from typing import List, Sequence, Tuple

INF = float("inf")


def _hungarian(cost: List[List[float]]) -> List[int]:
    """Minimise total cost for a square matrix. Returns col index per row."""
    n = len(cost)
    if n == 0:
        return []
    # 1-indexed potentials, standard JV formulation.
    u = [0.0] * (n + 1)
    v = [0.0] * (n + 1)
    p = [0] * (n + 1)      # p[col] = row assigned to col
    way = [0] * (n + 1)

    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = [INF] * (n + 1)
        used = [False] * (n + 1)
        while True:
            used[j0] = True
            i0, delta, j1 = p[j0], INF, 0
            for j in range(1, n + 1):
                if used[j]:
                    continue
                cur = cost[i0 - 1][j - 1] - u[i0] - v[j]
                if cur < minv[j]:
                    minv[j], way[j] = cur, j0
                if minv[j] < delta:
                    delta, j1 = minv[j], j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while j0:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1

    assignment = [-1] * n
    for j in range(1, n + 1):
        if p[j]:
            assignment[p[j] - 1] = j - 1
    return assignment


def max_weight_matching(profit: Sequence[Sequence[float]],
                        min_profit: float = 0.0) -> List[Tuple[int, int]]:
    """Pair rows with columns to maximise total profit.

    `profit` may be rectangular. Pairs scoring at or below `min_profit` are
    dropped, so a row with nothing worth matching stays unmatched rather than
    being forced onto a poor column.
    """
    n_rows = len(profit)
    n_cols = len(profit[0]) if n_rows else 0
    if n_rows == 0 or n_cols == 0:
        return []

    size = max(n_rows, n_cols)
    # Pad to square with zero profit, and convert to a cost matrix.
    cost = [[0.0] * size for _ in range(size)]
    for i in range(n_rows):
        for j in range(n_cols):
            cost[i][j] = -float(profit[i][j])

    assignment = _hungarian(cost)
    pairs = []
    for i, j in enumerate(assignment):
        if i < n_rows and 0 <= j < n_cols and profit[i][j] > min_profit:
            pairs.append((i, j))
    return pairs
