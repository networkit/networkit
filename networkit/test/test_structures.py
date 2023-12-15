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
		pa.addToSubset(1,0)
		pa.toSingleton(1)
		self.assertEqual(pa.numberOfElements(), pa.numberOfSubsets())

	def testPartitionEquality(self):
		p1 = nk.Partition(3)	
		p1.addToSubset(0,0)
		p1.addToSubset(1,1)
		p2 = nk.Partition(3)
		p2.addToSubset(0,0)
		p2.addToSubset(1,1)
		self.assertTrue(p1==p2)
		p1.extend()
		self.assertFalse(p1==p2)

if __name__ == "__main__":
	unittest.main()
