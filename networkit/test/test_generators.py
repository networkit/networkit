#!/usr/bin/env python3
import unittest
import networkit as nk
import math

class TestGenerators(unittest.TestCase):
	def setUp(self):
		pass

	def test_GeometricInhomogenous(self):
		# Test that lazy initialization works on all methods
		assert(len(nk.generators.GeometricInhomogenousGenerator(100, 10).getPositions()) == 100)
		assert(len(nk.generators.GeometricInhomogenousGenerator(100, 10).getWeights()) == 100)
		assert(nk.generators.GeometricInhomogenousGenerator(100, 10).generate().numberOfNodes() == 100)
		assert(nk.generators.GeometricInhomogenousGenerator(100, 10).generateKeepingInput().numberOfNodes() == 100)

		# Test that generator produces graphs as requested
		errorFound = False
		for n in [100, 200, 500]:
			for avgDeg in [n / 10, n / 5]:
				for ple in [2.5, 3.0, 4.0]:
					for alpha in [1.5, 10, math.inf]:
						for dim in [1, 2, 5]:
							try:
								gen = nk.generators.GeometricInhomogenousGenerator(n, avgDeg, ple, alpha, dim)

								positions = gen.getPositions()
								weights = gen.getWeights()

								assert(len(positions) == n)
								assert(len(weights) == n)

								for pt in positions:
									assert(len(pt) == dim)

								graphs = []
								graphs.append( (gen.generateKeepingInput(), "Keep") )

								assert(gen.getPositions() == positions)
								assert(gen.getWeights() == weights)

								graphs.append( (gen.generate(), "NotKeep") )

								assert(len(gen.getPositions()) == 0)
								assert(len(gen.getPositions()) == 0)

								gen2 = nk.generators.GeometricInhomogenousGenerator.fromPoints(positions, weights, alpha)
								graphs.append( (gen2.generate(), "FromPSNoScale") )

								# now scale weights and let generate scale them back
								scaled_weights = [2*x for x in weights]
								gen3 = nk.generators.GeometricInhomogenousGenerator.fromPoints(positions, scaled_weights, alpha, avgDeg)
								graphs.append( (gen3.generate(), "FromPSScale") )

								for G, label in graphs:
									assert G.numberOfNodes() == n, label
									ad = 2 * G.numberOfEdges() / n
									assert ad > avgDeg / 2, label + ": AverageDegree generated %f" % ad
									assert ad < 2 * avgDeg, label + ": AverageDegree generated %f" % ad

							except AssertionError as e:
								print("n=", n, " avgDeg=", avgDeg, " ple=", ple, " alpha=", alpha, " dim=", dim)
								print(" -> Error: ", e)
								errorFound = True

		if errorFound:
			raise "An assertion failed"

if __name__ == "__main__":
	unittest.main()
