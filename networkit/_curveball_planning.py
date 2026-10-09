"""Conditional ordinary-Curveball planning, without sampling or decisions.

Copyright 2026 DaniilKi contributors. SPDX-License-Identifier: Apache-2.0
Adapted from degree_null.planning / planning_core, Curveball Error Budget,
published v0.3.0 commit a4858a43354b9e26064ac11d3705b88dec8d4f42.
License and attribution:
_licenses/curveball-planner/LICENSE and NOTICE. No sampler code is copied.
"""
from dataclasses import dataclass
from fractions import Fraction
from math import comb
import re


@dataclass(frozen=True)
class CurveballTradePlan:
    """Immutable attempted-trade count, conditional premises and count-cap result.

    This is a planning record, not an Algorithm or a sampling certificate.
    ``withinTradeLimit`` compares a count cap, not time, memory or scientific
    adequacy. It is None when no cap was supplied. ``degrees`` retains labels
    by position, including isolates; disconnected realizations are included.
    """

    degrees: tuple[int, ...]
    epsilon: Fraction
    attemptedTrades: int
    unorderedVertexPairs: int
    stateCountLog2Ceiling: int
    binaryBlocks: int
    forcedUnique: bool
    maxTrades: int | None

    @property
    def withinTradeLimit(self):
        return None if self.maxTrades is None else self.attemptedTrades <= self.maxTrades

    @property
    def assumptions(self):
        """Premises needed to interpret a nontrivial plan as a TV budget."""
        return (
            "ordinary uniform unordered vertex-pair heat-bath kernel; "
            "nonnegative spectrum and spectral gap at least 1/binom(n,2)",
            "independent ideal uniform pairs and fixed-size exclusive-neighbor subsets",
            "all attempted trades, including unchanged outcomes, count",
            "uniform labeled simple undirected realizations of the exact degrees; "
            "no connectedness restriction",
            "implementation/kernel correspondence and PRNG quality are not certified",
        )


def _validateDegrees(degrees):
    if not isinstance(degrees, (list, tuple)):
        raise TypeError("degrees must be a finite list or tuple")
    n = len(degrees)
    if n > 1000:
        raise ValueError("this bounded planner supports at most 1000 vertices")
    if any(type(d) is not int or not 0 <= d < n for d in degrees) or sum(degrees) % 2:
        raise ValueError("invalid simple-undirected degree vector")
    # Havel-Hakimi graphicality check; no graph or samples are constructed.
    remaining = list(degrees)
    while remaining:
        remaining.sort(reverse=True)
        d = remaining.pop(0)
        if d == 0:
            break
        if d > len(remaining) or remaining[d - 1] == 0:
            raise ValueError("degree vector is not graphical")
        for i in range(d):
            remaining[i] -= 1
    return tuple(degrees)


def _forcedUnique(degrees):
    # Sufficient isolated/universal peeling test; a false result is inconclusive.
    remaining = list(degrees)
    while remaining:
        if 0 in remaining:
            remaining.remove(0)
        elif len(remaining) - 1 in remaining:
            remaining.remove(len(remaining) - 1)
            remaining = [d - 1 for d in remaining]
        else:
            return False
    return True


def curveballTradePlan(degrees, epsilon="1/10", *, maxTrades=None):
    """Plan sufficient *attempted* ordinary trades under explicit ideal premises.

    Parameters
    ----------
    degrees : list[int] or tuple[int, ...]
        Graphical labeled simple-undirected degrees, including isolates, at
        most 1000 entries. Each position is a label; its order is preserved.
    epsilon : fractions.Fraction or str
        Exact per-output total-variation allowance, strictly between 0 and
        1/2. Strings are bounded decimals or integer fractions. Binary floats
        and scientific-notation strings are refused. Fraction numerator and
        denominator are limited to 4096 bits each.
    maxTrades : int or None, optional
        Nonnegative attempted-trade cap. An insufficient cap is reported by
        ``withinTradeLimit=False`` without shortening the requested plan.

    Returns
    -------
    CurveballTradePlan
        Exact count and premises. No backend, graph, PRNG or inference is run.

    Notes
    -----
    This conditional bound assumes the specified ordinary heat-bath operator
    has nonnegative spectrum and gap at least 1/binom(n,2); the helper does
    not prove these premises or certify NetworKit's implementation/PRNG.
    For C=binom(n,2), m=sum(degrees)/2, B=binom(C,m), choose upward-rounded
    u=ceil(log2(B)), v=max(0,ceil(log2(1/(2*epsilon)))), k=ceil(u/2)+v, t=C*k.
    Exact integers verify B <= 2**u and B <= 4*epsilon**2*2**(2*k).
    Unique realizations certified by isolated/universal peeling need no
    trades. This says nothing about usefulness of a significance test.
    This budget is not for GlobalCurveball rounds, successful switches,
    directed/weighted graphs, matrix kernels or connected-only ensembles.
    Integer cost can be large even within the vertex bound; no runtime/RSS
    or benefit over a heuristic schedule is promised.
    """
    labeledDegrees = _validateDegrees(degrees)
    if maxTrades is not None and (type(maxTrades) is not int or maxTrades < 0):
        raise ValueError("maxTrades must be a nonnegative integer or None")
    if type(epsilon) is str:
        pattern = r"(?:[0-9]{1,200}(?:/[0-9]{1,200})?|[0-9]{0,200}\.[0-9]{1,200})"
        if len(epsilon) > 401 or re.fullmatch(pattern, epsilon) is None:
            raise ValueError("epsilon must be a bounded exact decimal or fraction")
    elif type(epsilon) is not Fraction:
        raise TypeError("epsilon must be a Fraction or exact decimal/fraction string")
    epsilon = Fraction(epsilon)
    if not 0 < epsilon < Fraction(1, 2):
        raise ValueError("epsilon must lie strictly between 0 and 1/2")
    if max(epsilon.numerator.bit_length(), epsilon.denominator.bit_length()) > 4096:
        raise ValueError("epsilon numerator/denominator exceed the 4096-bit planning limit")
    n = len(labeledDegrees)
    pairs = n * (n - 1) // 2
    unique = _forcedUnique(labeledDegrees)
    u = blocks = trades = 0
    if not unique:
        m = sum(labeledDegrees) // 2
        bound = comb(pairs, m)
        u = (bound - 1).bit_length()
        v = max(0, epsilon.denominator.bit_length() - (2 * epsilon.numerator).bit_length())
        if 2 * epsilon.numerator * (1 << v) < epsilon.denominator:
            v += 1
        blocks = (u + 1) // 2 + v
        if bound * epsilon.denominator**2 > 4 * epsilon.numerator**2 * (1 << (2 * blocks)):
            raise ArithmeticError("exact upward-rounding check failed")
        trades = pairs * blocks
    return CurveballTradePlan(labeledDegrees, epsilon, trades, pairs, u, blocks, unique, maxTrades)
