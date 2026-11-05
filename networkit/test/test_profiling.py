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


if __name__ == "__main__":
	unittest.main()
