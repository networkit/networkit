.. Copyright 2026 DaniilKi contributors. SPDX-License-Identifier: Apache-2.0

Conditional planning for ordinary Curveball
------------------------------------------

``curveballTradePlan`` computes a conservative attempted-trade budget without
running a sampler or making an inferential decision. It is for the ordinary
``Curveball`` kernel paired with ``CurveballUniformTradeGenerator``. It does
not give a budget for ``GlobalCurveball`` rounds, successful edge switches,
directed/weighted graphs, binary matrices or connected-only ensembles.

For example, plan for a labeled six-cycle using only synthetic data::

    from fractions import Fraction
    import networkit as nk

    plan = nk.randomization.curveballTradePlan(
        [2] * 6, Fraction(1, 19900), maxTrades=314)
    assert plan.attemptedTrades == 315
    assert plan.withinTradeLimit is False  # the count is not truncated

``degrees`` includes isolates and preserves their positions as labels. The
bounded helper checks graphicality and supports at most 1000 vertices.
Large integer planning costs can still be substantial; the returned count
is not a runtime, memory or performance prediction. Only exact ``Fraction``
or bounded decimal/fraction text tolerances are accepted. Fraction numerator
and denominator are limited to 4096 bits each.

The mathematical interpretation is conditional on a specific ordinary
uniform-pair heat-bath operator: nonnegative spectrum and spectral gap at
least ``1/binom(n,2)``, independent ideal uniform pair/subset draws, and all
attempted outcomes counted, including unchanged graphs. The helper neither
proves that operator premise nor certifies NetworKit's backend or seeded PRNG.
The state space includes every simple undirected labeled realization of the
degree vector, including disconnected graphs.

With ``C=binom(n,2)``, ``m=sum(degrees)/2`` and ``B=binom(C,m)``, choose
``u=ceil(log2(B))``, ``v=max(0,ceil(log2(1/(2*epsilon))))``,
``k=ceil(u/2)+v`` and ``t=C*k``. Exact integer checks establish the upward
rounding in ``B <= 4*epsilon**2*2**(2*k)``. Under the stated premises the
binary bound is ``TV <= 0.5*sqrt(B)*2**(-k)``. Isolated/universal peeling is
a sufficient unique-realization check; such a state needs zero trades.
Uniqueness does not make a useful significance test.

Planning a batch of independent outputs requires a separately chosen
per-output allowance. For example, a union bound divides a chosen batch
allowance by its fixed number of outputs. The helper does not select a
test statistic, tail, significance threshold, batch or failure policy.
No accuracy or runtime improvement over heuristic schedules is asserted.

If an application separately decides to execute a plan, existing ordinary
Curveball can consume successive bounded trade lists. The budget counts
*attempted trades*, not affected edges; all draws must be retained. Never
allocate one list of every planned trade blindly, truncate an insufficient
count cap or change to global rounds while retaining this interpretation.
Generating pairs requires at least two vertices; a zero-trade unique plan
does not call the trade generator.

This planner adapts the Apache-2.0 Curveball Error Budget project's existing
conditional arithmetic; it is not a new sampler or new mixing theorem:
https://github.com/DaniilKi/curveball-error-budget/tree/a4858a43354b9e26064ac11d3705b88dec8d4f42
This is the published v0.3.0 planner source; the earlier kernel arithmetic has
historical lineage at ``9319b8289617a6a08aeb528000bddff4da7ca121``.
The conditional operator input is credited there to OpenAI/math family 131,
commit ``fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb``. Ordinary Curveball itself
is the existing algorithm of Carstens, Berger and Strona, and Carstens et al.
(ESA 2018), distinct from the planner's conditional budgeting corollary.
Component license and attribution accompany the installed helper under
``networkit/_licenses/curveball-planner/``; no original license is replaced.
Package metadata identifies both MIT and Apache-2.0 components. These are
different component licenses, not a choice of license for every file.
