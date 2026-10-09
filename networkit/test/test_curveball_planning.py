"""Small deterministic controls for the planner, not sampler-law tests.

Copyright 2026 DaniilKi contributors. SPDX-License-Identifier: Apache-2.0
See networkit/_licenses/curveball-planner/LICENSE and NOTICE.
"""
import itertools
import math
import unittest
from dataclasses import FrozenInstanceError
from fractions import Fraction

from networkit._curveball_planning import curveballTradePlan


class TestCurveballPlanning(unittest.TestCase):
    def testSyntheticSixCycle(self):
        plan = curveballTradePlan([2] * 6, Fraction(1, 19900))
        self.assertEqual(plan.attemptedTrades, 315)
        self.assertEqual(plan.binaryBlocks, 21)
        self.assertEqual(plan.stateCountLog2Ceiling, 13)
        self.assertEqual(plan.unorderedVertexPairs, 15)
        self.assertFalse(plan.forcedUnique)
        self.assertEqual(plan.degrees, (2,) * 6)
        self.assertIn("not certified", plan.assumptions[-1])

    def testIndependentIntegerInequalities(self):
        for degrees in ([1] * 4, [2] * 5, [2] * 6, [3] * 6):
            for epsilon in (Fraction(1, 4), Fraction(1, 10), Fraction(1, 16), Fraction(1, 19900)):
                plan = curveballTradePlan(degrees, epsilon)
                bound = math.comb(math.comb(len(degrees), 2), sum(degrees) // 2)
                self.assertLessEqual(bound, 1 << plan.stateCountLog2Ceiling)
                self.assertLessEqual(bound * epsilon.denominator**2,
                                     4 * epsilon.numerator**2 * (1 << (2 * plan.binaryBlocks)))
                self.assertEqual(plan.attemptedTrades, plan.unorderedVertexPairs * plan.binaryBlocks)

    def testExhaustiveSmallGraphicalityAndUniqueClaims(self):
        # Enumerate every labeled n<=5 graph independently of Havel-Hakimi.
        for n in range(6):
            edges = list(itertools.combinations(range(n), 2))
            realizations = {}
            for mask in range(1 << len(edges)):
                degrees = tuple(sum(bool(mask & (1 << k)) for k, e in enumerate(edges) if v in e)
                                for v in range(n))
                realizations[degrees] = realizations.get(degrees, 0) + 1
            for degrees in itertools.product(range(n), repeat=n):
                if degrees not in realizations:
                    with self.assertRaises(ValueError):
                        curveballTradePlan(degrees)
                else:
                    plan = curveballTradePlan(degrees)
                    if plan.forcedUnique:
                        self.assertEqual(realizations[degrees], 1)
                        self.assertEqual(plan.attemptedTrades, 0)

    def testIsolatesAndLabelsRetained(self):
        degrees = [0, 2, 2, 2, 2, 0]
        plan = curveballTradePlan(degrees)
        self.assertEqual(plan.degrees, tuple(degrees))
        self.assertEqual(degrees, [0, 2, 2, 2, 2, 0])
        self.assertFalse(plan.forcedUnique)

    def testUniqueRealizations(self):
        for degrees in ([], [0], [0] * 6, [5] * 6, [5, 1, 1, 1, 1, 1], [0, 1, 1]):
            plan = curveballTradePlan(degrees)
            self.assertTrue(plan.forcedUnique)
            self.assertEqual(plan.attemptedTrades, 0)

    def testCountCapDoesNotTruncate(self):
        plan = curveballTradePlan([2] * 6, "1/19900", maxTrades=314)
        self.assertFalse(plan.withinTradeLimit)
        self.assertEqual(plan.attemptedTrades, 315)
        self.assertTrue(curveballTradePlan([2] * 6, "1/19900", maxTrades=315).withinTradeLimit)
        self.assertIsNone(curveballTradePlan([2] * 6).withinTradeLimit)

    def testImmutableRecordAndInputCopy(self):
        degrees = [2] * 6
        plan = curveballTradePlan(degrees)
        degrees[0] = 0
        self.assertEqual(plan.degrees, (2,) * 6)
        with self.assertRaises(FrozenInstanceError):
            plan.attemptedTrades = 0

    def testExactToleranceFormats(self):
        self.assertEqual(curveballTradePlan([2] * 6, ".1"), curveballTradePlan([2] * 6, "1/10"))
        for epsilon in (0.1, True, 1, None):
            with self.assertRaises(TypeError):
                curveballTradePlan([2] * 6, epsilon)
        for epsilon in ("0", "1/2", "1", "-1/10", "1e-9", "1/0", "9" * 402):
            with self.assertRaises((ValueError, ZeroDivisionError)):
                curveballTradePlan([2] * 6, epsilon)

    def testInputErrors(self):
        for degrees in ([True, 1], [1.0, 1], [-1, 1], [2, 0], [1, 0], [3, 3, 1, 1], [0] * 1001):
            with self.assertRaises(ValueError):
                curveballTradePlan(degrees)
        with self.assertRaises(TypeError):
            curveballTradePlan(iter([2] * 6))
        for cap in (-1, True, 1.5):
            with self.assertRaises(ValueError):
                curveballTradePlan([2] * 6, maxTrades=cap)

    def testTighterToleranceCannotReduceCount(self):
        counts = [curveballTradePlan([2] * 6, e).attemptedTrades
                  for e in ("1/4", "1/10", "1/100", "1/19900")]
        self.assertEqual(counts, sorted(counts))

    def testBinaryRoundingBoundaries(self):
        for k in range(1, 20):
            for denominator in ((1 << k) - 1, 1 << k, (1 << k) + 1):
                if denominator <= 2:
                    continue
                epsilon = Fraction(1, denominator)
                expected = 0
                while 2 * (1 << expected) < denominator:
                    expected += 1
                plan = curveballTradePlan([2] * 6, epsilon)
                self.assertEqual(plan.binaryBlocks, 7 + expected)
        with self.assertRaises(ValueError):
            curveballTradePlan([2] * 6, Fraction(1, 1 << 4096))


class TestCurveballPlanningPublicExport(unittest.TestCase):
    def testPublicRandomizationExport(self):
        # Requires the rebuilt randomization extension, unlike pure helper tests.
        import networkit.randomization as randomization
        self.assertIs(randomization.curveballTradePlan, curveballTradePlan)
        self.assertEqual(randomization.curveballTradePlan([2] * 6, "1/19900").attemptedTrades, 315)
