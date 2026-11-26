#!/usr/bin/env python3
import unittest
import os

import networkit as nk

class TestStructures(unittest.TestCase):

	def setUp(self):
		self.G = nk.readGraph("input/jazz.graph", nk.Format.METIS)
		self.G.indexEdges()

	def testCoverBounds(self):
		co = nk.Cover(4)
		co.extend() #cover of size 5
		co.setUpperBound(4)
		self.assertEqual(co.upperBound(), 4)
		self.assertEqual(co.lowerBound(), 0)

	def testCoverSubsetsSizes(self):
		co = nk.Cover(3)
		co.setUpperBound(3)
		co.addToSubset(0,0)
		co.addToSubset(1,1)
		co.addToSubset(2,2)
		self.assertEqual(co.numberOfElements(), 3)
		self.assertEqual(co.numberOfSubsets(), 3)
		self.assertSetEqual(co.getSubsetIds(), {0,1,2})
		self.assertListEqual(co.subsetSizes(), [1,1,1])
		self.assertDictEqual(co.subsetSizeMap(), {0: 1, 1: 1, 2: 1})
		
	def testCoverSubsetsRemove(self):
		co = nk.Cover(2)
		co.setUpperBound(2)
		co.addToSubset(0,0)	
		self.assertSetEqual(co.getMembers(0), {0})
		co.removeFromSubset(0,0)
		self.assertFalse(co.contains(0))

	def testCoverSubsetsMerge(self):	
		co = nk.Cover(3)
		co.setUpperBound(3)
		co.addToSubset(0,0)
		co.addToSubset(1,1)
		co.mergeSubsets(0,1) #0+1 -> 2
		self.assertTrue(co.inSameSubset(0,1))
		
	def testCoverSingletons(self):	
		co = nk.Cover(3)
		co.allToSingletons()
		co.moveToSubset(0,1) 
		co.toSingleton(0)
		self.assertEqual(co.numberOfElements(),co.numberOfSubsets())

	def testPartitionBounds(self):
		pa = nk.Partition(4)
		pa.extend() #partition of size 5
		pa.setUpperBound(5)
		self.assertEqual(pa.upperBound(),5)
		self.assertEqual(pa.lowerBound(),0)
		
	def testPartitionSubsetsSizes(self):
		pa = nk.Partition(3)
		pa.setUpperBound(3)	
		pa.addToSubset(0,0)
		pa.addToSubset(1,1)
		pa.addToSubset(2,2)
		self.assertEqual(pa.numberOfSubsets(), 3)
		self.assertSetEqual(pa.getSubsetIds(), {0,1,2})		

	def testPartitionSubsetsMerge(self):
		pa = nk.Partition(2)
		pa.setUpperBound(2)	
		pa.addToSubset(0,0)
		pa.addToSubset(1,1)
		pa.mergeSubsets(0,1) #0+1 -> 2
		self.assertTrue(pa.inSameSubset(0,1))
		self.assertTrue(pa.contains(0))
		
	def testPartitionSingletons(self):
		pa = nk.Partition(3)	
		pa.allToSingletons()
		pa.moveToSubset(1,0)
		pa.toSingleton(1)
		self.assertEqual(pa.numberOfElements(), pa.numberOfSubsets())

	def testPartitionEquality(self):
		p1 = nk.Partition(3)
		p1.setUpperBound(2)	
		p1.addToSubset(0,0)
		p1.addToSubset(1,1)
		p2 = nk.Partition(3)
		p2.setUpperBound(2)
		p2.addToSubset(1,0)
		p2.addToSubset(0,1)
		self.assertTrue(p1==p2)
		p1.addToSubset(0,2)
		p2.addToSubset(0,2)
		self.assertFalse(p1==p2)
		
	def testMoveToSubset(self):
		p1 = nk.Partition(3)
		p1.setUpperBound(2)	
		p1.addToSubset(0,0)
		p1.addToSubset(0,1)
		p1.moveToSubset(1,1)
		self.assertEqual(p1.subsetOf(0),0)
		self.assertEqual(p1.subsetOf(1),1)

	def testFailedAssertionsPartition(self):
		p = nk.Partition(3)
		p.setUpperBound(2)
		p.addToSubset(0,0)
		self.assertEqual(p.getVector(), [0, nk.none, nk.none])
		self.assertEqual(p.upperBound(), 2)

		# cannot get or set subset of element that is out of bounds
		with self.assertRaises(IndexError):
			p.subsetOf(3)
		with self.assertRaises(IndexError):
			p[3]
		with self.assertRaises(IndexError):
			p[3] = 0

		# cannot set subset id that is out of bounds
		with self.assertRaises(IndexError):
			p[1] = 99

		# cannot add or move element that is out of bounds
		with self.assertRaises(IndexError):
			p.addToSubset(0,3)
		with self.assertRaises(IndexError):
			p.moveToSubset(0,3)

		# cannot add or move element to subset with out of bounds id
		with self.assertRaises(IndexError):
			p.addToSubset(99,1)
		with self.assertRaises(IndexError):
			p.moveToSubset(99,1)

		# cannot add element that has already been assigned
		with self.assertRaises(ValueError):
			p.addToSubset(1,0)
		# cannot move element that has not yet been assigned
		with self.assertRaises(ValueError):
			p.moveToSubset(0,1)

		# cannot make a singleton of element that is out of bounds
		with self.assertRaises(IndexError):
			p.toSingleton(3)

		# cannot merge subsets if one of their ids is out of bounds
		with self.assertRaises(IndexError):
			p.mergeSubsets(0,99)

		# cannot check if two elements belong to the same subset if one of them is out of bounds or unassigned
		with self.assertRaises(IndexError):
			p.inSameSubset(0,3)
		with self.assertRaises(ValueError):
			p.inSameSubset(0,1)

		# cannot get the members of a subset if its id is out of bounds
		with self.assertRaises(IndexError):
			p.getMembers(99)

		# can reduce upperBound as long as it stays above the largest assigned subset id
		p.setUpperBound(1)
		self.assertEqual(p.upperBound(), 1)

		# cannot reduce upperBound below the largest assigned subset id
		p.setUpperBound(2)
		p.addToSubset(1,1)
		with self.assertRaises(ValueError):
			p.setUpperBound(1)

	def testFailedAssertionsCover(self):
		c = nk.Cover(3)
		c.setUpperBound(2)
		c.addToSubset(0,0)
		self.assertEqual(c.subsetsOf(0), {0})
		self.assertEqual(c.subsetsOf(1), set())
		self.assertEqual(c.subsetsOf(2), set())
		self.assertEqual(c.upperBound(), 2)

		# cannot get subset of element that is out of bounds
		with self.assertRaises(IndexError):
			c.subsetsOf(3)

		# cannot add, move or remove element that is out of bounds
		with self.assertRaises(IndexError):
			c.addToSubset(0,3)
		with self.assertRaises(IndexError):
			c.moveToSubset(0,3)
		with self.assertRaises(IndexError):
			c.removeFromSubset(0,3)

		# cannot add, move or remove element to subset with out of bounds id
		with self.assertRaises(IndexError):
			c.addToSubset(99,1)
		with self.assertRaises(IndexError):
			c.moveToSubset(99,1)
		with self.assertRaises(IndexError):
			c.removeFromSubset(99, 1)

		# cannot make a singleton of element that is out of bounds
		with self.assertRaises(IndexError):
			c.toSingleton(3)

		# cannot merge subsets if one of their ids is out of bounds
		with self.assertRaises(IndexError):
			c.mergeSubsets(0,99)

		# cannot check if two elements belong to the same subset if one of them is out of bounds or unassigned
		with self.assertRaises(IndexError):
			c.inSameSubset(0,3)
		with self.assertRaises(ValueError):
			c.inSameSubset(0,1)

		# cannot get the members of a subset if its id is out of bounds
		with self.assertRaises(IndexError):
			c.getMembers(99)

		# can reduce upperBound as long as it stays above the largest assigned subset id
		c.setUpperBound(1)
		self.assertEqual(c.upperBound(), 1)

		# cannot reduce upperBound below the largest assigned subset id
		c.setUpperBound(2)
		c.addToSubset(1,1)
		with self.assertRaises(ValueError):
			c.setUpperBound(1)


	def testPartitionEqualityDifferentStructure(self):
		"""Equal-length partitions with different grouping are not equal."""
		singles = nk.Partition(3)
		singles.allToSingletons()

		one = nk.Partition(3)
		one.setUpperBound(1)
		one.addToSubset(0, 0)
		one.addToSubset(0, 1)
		one.addToSubset(0, 2)
		self.assertFalse(singles == one)

		# Same singleton structure, different subset ids → equal (id-independent)
		a = nk.Partition(3)
		a.setUpperBound(10)
		a.addToSubset(9, 0)
		a.addToSubset(7, 1)
		a.addToSubset(4, 2)
		b = nk.Partition(3)
		b.allToSingletons()
		self.assertTrue(a == b)

		# Same multiset of sizes [2,1] but different element assignment → unequal
		pair_left = nk.Partition(3)
		pair_left.setUpperBound(2)
		pair_left.addToSubset(0, 0)
		pair_left.addToSubset(0, 1)
		pair_left.addToSubset(1, 2)
		pair_right = nk.Partition(3)
		pair_right.setUpperBound(2)
		pair_right.addToSubset(0, 0)
		pair_right.addToSubset(1, 1)
		pair_right.addToSubset(1, 2)
		self.assertFalse(pair_left == pair_right)

		# Unassigned (none) vs assigned at the same index is not equal.
		left = nk.Partition(2)
		left.addToSubset(0, 0)
		right = nk.Partition(2)
		right.addToSubset(0, 1)
		self.assertFalse(left == right)

		# Same assignment, both unassigned at the same index, different ids.
		same_a = nk.Partition(2)
		same_a.setUpperBound(5)
		same_a.addToSubset(4, 0)
		same_b = nk.Partition(2)
		same_b.addToSubset(0, 0)
		self.assertTrue(same_a == same_b)
		self.assertTrue(nk.Partition(3) == nk.Partition(3))

	def testSetItemReassignsAndUnassigns(self):
		p = nk.Partition(3)
		p.allToSingletons()
		self.assertEqual(p.getVector(), [0, 1, 2])
		p[0] = 1
		self.assertEqual(p.getVector(), [1, 1, 2])
		p[0] = nk.none
		self.assertEqual(p.getVector(), [nk.none, 1, 2])

if __name__ == "__main__":
	unittest.main()
