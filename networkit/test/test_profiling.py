#!/usr/bin/env python3
import unittest

import networkit as nk
from networkit.profiling.profiling import Config, Profile


class TestProfilingConfig(unittest.TestCase):
	def testCreateConfigUnknownPresetRaisesValueError(self):
		with self.assertRaises(ValueError) as ctx:
			Config.createConfig("typo")
		msg = str(ctx.exception)
		self.assertIn("unknown preset", msg)
		self.assertIn("typo", msg)

	def testCreateConfigKnownPresets(self):
		for preset in ("complete", "minimal", "default"):
			cfg = Config.createConfig(preset)
			self.assertIsInstance(cfg, Config)


class TestProfileCreate(unittest.TestCase):
	def testCreateOnOneAndTwoNodeGraphs(self):
		# Issue #1498: KPathCentrality scores are NaN on these graphs, and
		# Profile.create used to divide by zero inside ranked().
		for n in (1, 2):
			with self.subTest(n=n):
				self.assertIsInstance(Profile.create(nk.Graph(n)), Profile)

	def testCreateWithPartitionMeasures(self):
		# Timings for assortativity and centralization only exist for node centralities;
		# partition-only configs used to raise UnboundLocalError, and a partition after
		# a centrality reported that centrality's timings.
		G = nk.Graph(3)
		G.addEdge(0, 1)
		G.addEdge(1, 2)
		G.addEdge(0, 2)
		for measures in (["Partition.Communities"], ["Centrality.Degree", "Partition.Communities"]):
			with self.subTest(measures=measures):
				cfg = Config()
				for measure in measures:
					cfg.setMeasure(measure)
				profile = Profile.create(G, config=cfg)
				self.assertEqual(profile.getElapsedTime("Partition.Communities")[3:], (0.0, 0.0))


if __name__ == "__main__":
	unittest.main()
