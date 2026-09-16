#!/usr/bin/env python3

import time
import unittest

import networkit as nk

class TestSubgraphIsomorphism(unittest.TestCase):

	def setUp(self):
		# Pattern: Edge (0,1)
		self.arc = nk.Graph(2)
		self.arc.addEdge(0, 1)

		# Target: Square {0, 1, 2, 3}
		self.square = nk.Graph(4)
		self.square.addEdge(0, 1)
		self.square.addEdge(1, 2)
		self.square.addEdge(2, 3)
		self.square.addEdge(3, 0)

		# Target: Square {0, 1, 2, 3} with diagonal (0,2)
		self.diagonal = nk.Graph(5)
		self.diagonal.addEdge(0, 1)
		self.diagonal.addEdge(1, 2)
		self.diagonal.addEdge(2, 3)
		self.diagonal.addEdge(3, 0)
		self.diagonal.addEdge(0, 2)

		self.expected = {
			(0, 1),
			(1, 2),
			(2, 3),
			(3, 0),
			(1, 0),
			(2, 1),
			(3, 2),
			(0, 3),
		}

	def testInstantiationOfAbstractBaseClass(self):
		with self.assertRaises(RuntimeError):
			nk.isomorphism.SubgraphIsomorphism(self.arc, self.square)

	def testNoMatch(self):
		target = nk.Graph(2)

		vf2 = nk.isomorphism.VF2(self.arc, target)
		vf2.run()

		self.assertFalse(vf2.hasMatch())
		self.assertEqual(vf2.numberOfMatches(), 0)

	def testMatchReporting(self):
		vf2 = nk.isomorphism.VF2(self.arc, self.square)
		vf2.setStoreMatches(False)
		vf2.run()

		self.assertEqual(vf2.numberOfWorkers(), 1)
		self.assertTrue(vf2.hasMatch())
		self.assertEqual(vf2.numberOfMatches(), 8)

		with self.assertRaises(RuntimeError):
			vf2.getMatches()

		vf2.setStoreMatches(True)
		vf2.run()

		self.assertEqual({tuple(match) for match in vf2.getMatches()}, self.expected)

		# Check that running again gives the same matches
		vf2.run()

		self.assertEqual({tuple(match) for match in vf2.getMatches()}, self.expected)

	def testIsoVsMono(self):
		# Due to the diagonal edge in the target there should be no induced matching
		vf2_iso = nk.isomorphism.VF2(self.square, self.diagonal, nk.isomorphism.Semantics.INDUCED)
		vf2_iso.run()
		self.assertFalse(vf2_iso.hasMatch())

		vf2_mono = nk.isomorphism.VF2(self.square, self.diagonal, nk.isomorphism.Semantics.MONOMORPHISM)
		vf2_mono.run()
		self.assertEqual(vf2_mono.numberOfMatches(), 8)

	def testNodeAndEdgeLabels(self):
		vf2 = nk.isomorphism.VF2(self.arc, self.square)

		vf2.setNodeLabels([0, 1], [0, 1, 2, 1])
		vf2.run()
		self.assertEqual(vf2.numberOfMatches(), 2)

		self.arc.indexEdges()
		self.square.indexEdges()
		vf2_2 = nk.isomorphism.VF2(self.arc, self.square)
		vf2_2.setNodeLabels([0, 1], [0, 1, 2, 1])
		vf2_2.setEdgeLabels([0], [0, 0, 0, 1])
		vf2_2.run()
		self.assertEqual(vf2_2.numberOfMatches(), 1)

		# Reset all node labels to wildcards so now only edge labels matter
		vf2_2.setNodeLabels([nk.none, nk.none], [nk.none, nk.none, nk.none, nk.none])
		vf2_2.run()
		self.assertEqual(vf2_2.numberOfMatches(), 6)

	def testMatchCap(self):
		vf2_noCap = nk.isomorphism.VF2(self.arc, self.square)
		vf2_cap0 = nk.isomorphism.VF2(self.arc, self.square, maxMatches=0)
		vf2_cap3 = nk.isomorphism.VF2(self.arc, self.square, maxMatches=3)

		vf2_noCap.run()
		vf2_cap0.run()
		vf2_cap3.run()

		self.assertEqual(vf2_noCap.numberOfMatches(), 8)
		self.assertEqual(vf2_cap0.numberOfMatches(), 8)
		self.assertEqual(vf2_cap3.numberOfMatches(), 3)

	def testDifferentAlgos(self):
		vf2 = nk.isomorphism.VF2(self.arc, self.square)
		vf2.run()

		ri = nk.isomorphism.RI(self.arc, self.square)
		ri.run()

		riPar = nk.isomorphism.ParallelRI(self.arc, self.square)
		riPar.run()

		self.assertTrue(vf2.hasMatch())
		self.assertEqual(vf2.numberOfMatches(), 8)

		self.assertTrue(ri.hasMatch())
		self.assertEqual(ri.numberOfMatches(), 8)

		self.assertTrue(riPar.hasMatch())
		self.assertEqual(riPar.numberOfMatches(), 8)

		self.assertEqual({tuple(match) for match in vf2.getMatches()}, self.expected)

		self.assertEqual({tuple(match) for match in ri.getMatches()}, self.expected)

		self.assertEqual({tuple(match) for match in riPar.getMatches()}, self.expected)

	def testRIVariants(self):
		ri = nk.isomorphism.RI(self.arc, self.square, variant=nk.isomorphism.Variant.RI)
		ri.run()

		ri_ds = nk.isomorphism.RI(self.arc, self.square, variant=nk.isomorphism.Variant.RI_DS)
		ri_ds.run()

		self.assertTrue(ri.hasMatch())
		self.assertEqual(ri.numberOfMatches(), 8)

		self.assertTrue(ri_ds.hasMatch())
		self.assertEqual(ri_ds.numberOfMatches(), 8)

		self.assertEqual({tuple(match) for match in ri.getMatches()}, self.expected)

		self.assertEqual({tuple(match) for match in ri_ds.getMatches()}, self.expected)

	def testVF2Callback(self):
		matches = []

		def callback(match):
			matches.append(match)

		vf2 = nk.isomorphism.VF2(self.arc, self.square)
		vf2.setCallback(callback)
		vf2.setStoreMatches(True)
		vf2.run()

		self.assertEqual(len(matches), 8)
		self.assertEqual({tuple(match) for match in matches}, self.expected)
		# SubgraphIsomorphism used with a callback does not store the matches
		with self.assertRaises(RuntimeError):
			vf2.getMatches()

	def testRICallback(self):
		matches = []

		def callback(match):
			matches.append(match)

		ri = nk.isomorphism.RI(self.arc, self.square)
		ri.setStoreMatches(True)
		ri.setCallback(callback)
		ri.run()

		self.assertEqual(len(matches), 8)
		self.assertEqual({tuple(match) for match in matches}, self.expected)
		# SubgraphIsomorphism used with a callback does not store the matches
		with self.assertRaises(RuntimeError):
			ri.getMatches()

	def testParallelRICallback(self):
		matches = []

		def callback(workerId, match):
			matches.append((workerId, match))

		riPar = nk.isomorphism.ParallelRI(self.arc, self.square)
		riPar.setParallelCallback(callback)
		riPar.run()

		self.assertEqual(len(matches), 8)
		for workerId, match in matches:
			self.assertGreaterEqual(workerId, 0)
			self.assertLess(workerId, riPar.numberOfWorkers())
		self.assertEqual({tuple(match) for workerId, match in matches}, self.expected)
		# SubgraphIsomorphism used with a callback does not store the matches
		with self.assertRaises(RuntimeError):
			riPar.getMatches()

	def testSerialCallbackWithParallelRI(self):
		matches = []
		inFlight = 0
		maxInFlight = 0

		def callback(match):
			nonlocal inFlight, maxInFlight
			inFlight += 1
			maxInFlight = max(maxInFlight, inFlight)
			# Sleeping releases the GIL, which would let an overlapping call start
			time.sleep(0.001)
			matches.append(tuple(match))
			inFlight -= 1

		riPar = nk.isomorphism.ParallelRI(self.arc, self.square)
		riPar.setCallback(callback)
		riPar.run()

		# The workers of ParallelRI take turns at a callback set with setCallback
		self.assertEqual(maxInFlight, 1)
		self.assertEqual(len(matches), 8)
		self.assertEqual(set(matches), self.expected)

	def testParallelCallbackWithSequentialAlgorithms(self):
		matches_vf2 = []
		matches_ri = []

		def callback_vf2(workerId, match):
			matches_vf2.append((workerId, tuple(match)))

		def callback_ri(workerId, match):
			matches_ri.append((workerId, tuple(match)))

		vf2 = nk.isomorphism.VF2(self.arc, self.square)
		vf2.setParallelCallback(callback_vf2)
		vf2.run()

		ri = nk.isomorphism.RI(self.arc, self.square)
		ri.setParallelCallback(callback_ri)
		ri.run()

		# The sequential algorithms report every match as worker 0
		for matches in (matches_vf2, matches_ri):
			self.assertEqual(len(matches), 8)
			self.assertEqual({workerId for workerId, match in matches}, {0})
			self.assertEqual({match for workerId, match in matches}, self.expected)

	def testCallbackTypeMismatches(self):
		matches_seq = []
		matches_par = []

		def callback_seq(match):
			matches_seq.append(match)

		def callback_par(workerId, match):
			matches_par.append((workerId, match))

		vf2 = nk.isomorphism.VF2(self.arc, self.square)
		with self.assertRaises(TypeError):
			vf2.setParallelCallback(callback_seq)

		riPar = nk.isomorphism.ParallelRI(self.arc, self.square)
		with self.assertRaises(TypeError):
			riPar.setCallback(callback_par)

		with self.assertRaises(TypeError):
			vf2.setCallback(callback_par)

		with self.assertRaises(TypeError):
			riPar.setParallelCallback(callback_seq)

	def testCallbackCalledExactlyOncePerMatch(self):
		callback_count = 0

		def callback(match):
			nonlocal callback_count
			callback_count += 1

		vf2 = nk.isomorphism.VF2(self.arc, self.square)
		vf2.setCallback(callback)
		vf2.run()

		self.assertEqual(callback_count, 8)

	def testSetCallbackAfterRun(self):
		matches = []

		def callback(match):
			matches.append(match)

		vf2 = nk.isomorphism.VF2(self.arc, self.square)
		vf2.run()

		# No callback was set before run so there should not be anything in matches
		self.assertEqual(len(matches), 0)

		vf2.setCallback(callback)
		vf2.run()

		# Callback has been set before run so matches should be filled
		self.assertEqual(len(matches), 8)
		self.assertEqual({tuple(match) for match in matches}, self.expected)

	def testUseLastSetCallback(self):
		matches1 = []
		matches2 = []

		def callback1(match):
			matches1.append(match)

		def callback2(match):
			matches2.append(match)

		vf2 = nk.isomorphism.VF2(self.arc, self.square)
		vf2.setCallback(callback1)
		vf2.setCallback(callback2)
		vf2.run()

		# callback2 was most recently set so it should be filled, callback1 should be empty
		self.assertEqual(len(matches1), 0)
		self.assertEqual(len(matches2), 8)
		self.assertEqual({tuple(match) for match in matches2}, self.expected)

	def testUseLastSetCallbackOfEitherForm(self):
		matches_serial = []
		matches_parallel = []

		def callback_serial(match):
			matches_serial.append(tuple(match))

		def callback_parallel(workerId, match):
			matches_parallel.append(tuple(match))

		vf2 = nk.isomorphism.VF2(self.arc, self.square)
		vf2.setParallelCallback(callback_parallel)
		vf2.setCallback(callback_serial)

		# Nothing references the replaced callback anymore, so it must never be called
		del callback_parallel
		vf2.run()

		self.assertEqual(len(matches_parallel), 0)
		self.assertEqual(len(matches_serial), 8)
		self.assertEqual(set(matches_serial), self.expected)

	def testCallbackCanBeCallableObject(self):
		class Callback:
			def __init__(self):
				self.matches = []

			def __call__(self, match):
				self.matches.append(tuple(match))

		callback = Callback()

		vf2 = nk.isomorphism.VF2(self.arc, self.square)
		vf2.setCallback(callback)
		vf2.run()

		self.assertEqual(len(callback.matches), 8)
		self.assertEqual(set(callback.matches), self.expected)

	def testDifferentCallbacksDoNotShareState(self):
		matches_1 = []
		matches_2 = []

		def callback1(match):
			matches_1.append(tuple(match))

		def callback2(match):
			matches_2.append(tuple(match))

		vf2_1 = nk.isomorphism.VF2(self.arc, self.square)
		vf2_2 = nk.isomorphism.VF2(self.arc, self.square)

		vf2_1.setCallback(callback1)
		vf2_2.setCallback(callback2)

		vf2_1.run()
		vf2_2.run()

		self.assertEqual(len(matches_1), 8)
		self.assertEqual(len(matches_2), 8)
		self.assertEqual(matches_1, matches_2)

	def testCallbackKeepsObjectAlive(self):
		class Callback:
			def __init__(self):
				self.count = 0

			def __call__(self, match):
				self.count += 1

		callback = Callback()

		vf2 = nk.isomorphism.VF2(self.arc, self.square)
		vf2.setCallback(callback)

		del callback
		vf2.run()

		self.assertEqual(vf2.numberOfMatches(), 8)

	def testCallbackRaisesOnFirstMatchSeq(self):
		matches = []

		def callback(match):
			matches.append(tuple(match))
			raise RuntimeError("sequential callback failed")

		ri = nk.isomorphism.RI(self.arc, self.square)
		ri.setCallback(callback)

		with self.assertRaises(RuntimeError):
			ri.run()

		self.assertEqual(len(matches), 1)

	def testCallbackRaisesOnFirstMatchPar(self):
		matches = []

		def callback(workerId, match):
			matches.append((workerId, tuple(match)))
			raise RuntimeError("parallel callback failed")

		riPar = nk.isomorphism.ParallelRI(self.arc, self.square)
		riPar.setParallelCallback(callback)

		with self.assertRaises(RuntimeError):
			riPar.run()

		self.assertGreaterEqual(len(matches), 1)

	def testCallbackRaisesBaseException(self):
		matches = []

		def callback_seq(match):
			matches.append(tuple(match))
			raise SystemExit("callback exited")

		def callback_par(workerId, match):
			callback_seq(match)

		# SystemExit does not derive from Exception, but it has to stop the search as well
		ri = nk.isomorphism.RI(self.arc, self.square)
		ri.setCallback(callback_seq)

		with self.assertRaisesRegex(RuntimeError, "callback exited"):
			ri.run()

		self.assertEqual(len(matches), 1)

		riPar = nk.isomorphism.ParallelRI(self.arc, self.square)
		riPar.setParallelCallback(callback_par)

		with self.assertRaisesRegex(RuntimeError, "callback exited"):
			riPar.run()

	def testParallelCallbackResultsMatchSerialResults(self):
		matches_ri = []
		matches_riPar = []

		def callback_ri(match):
			matches_ri.append(tuple(match))

		def callback_riPar(workerId, match):
			matches_riPar.append(tuple(match))

		ri = nk.isomorphism.RI(self.arc, self.square)
		ri.setCallback(callback_ri)
		ri.run()

		riPar = nk.isomorphism.ParallelRI(self.arc, self.square)
		riPar.setParallelCallback(callback_riPar)
		riPar.run()

		self.assertEqual(set(matches_ri), set(matches_riPar))
		self.assertEqual(len(matches_ri), len(matches_riPar))

if __name__ == "__main__":
	unittest.main()