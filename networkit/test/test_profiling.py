#!/usr/bin/env python3
import unittest

from networkit.profiling.profiling import Config


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


if __name__ == "__main__":
	unittest.main()
