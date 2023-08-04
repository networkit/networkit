#!/usr/bin/env python3
import numpy as np
import os
import random
import unittest

import networkit as nk

class TestEngineering(unittest.TestCase):

	def setUp(self):
		self.L = nk.readGraph("input/looptest1.gml", nk.Format.GML) #without self-loops

	def testScaling(self):
		
		def createAlgo(input=self.L):
			return nk.centrality.Betweenness(input)
		
		self.L.indexEdges()
		Data1 = nk.engineering.strongScaling(createAlgo, inargs={}, threadSequence=[1,2,4], outPath="./scaling.data")
		Data2 = nk.engineering.weakScaling(createAlgo, inargs={}, inputTitles=["Test1", "Test2"], threadSequence=[1,2], inputSequence=[self.L, self.L], outPath="./scaling.data")
		self.assertEqual(len(Data1), 3)
		self.assertEqual(len(Data2), 2)
		
	def testThreads(self):
		nk.engineering.setNumberOfThreads(8)
		self.assertEqual(nk.engineering.getMaxNumberOfThreads(), 8)
		self.assertLessEqual(nk.engineering.getCurrentNumberOfThreads(), 8)	

	def testLogLevel(self):
		nk.engineering.setLogLevel("ERROR")
		self.assertEqual("ERROR", nk.engineering.getLogLevel())	

if __name__ == "__main__":
	unittest.main()
