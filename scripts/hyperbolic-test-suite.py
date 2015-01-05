#! /usr/bin/env python3

from networkit import *
import time
import json
import math
import psutil

cost = []

processors = psutil.NUM_CPUS

stamp = time.time()

cpucount = 2 #TODO: properly import number of cpus here
balancelist = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995]
capacitylist = range(7, 13)
threadlist = range(0, int(math.log2(cpucount)+1))
scalelist = range(10, 26)# upper limit depends on platform. (Maybe) lower on phi, definitely higher on vsmp
iterations = 10
# Forget about openmp chunk size and schedule. They are not adjustable on the python level and the optimal value probably doesn't depend much on the other
# parameters.

# All out parameter study:
# For each node size, make one pass over all parameters to find an optimum while keeping everything else at default, then make a second pass with 
# updated parameters to see whether the optimum persists. Use all threads. Make seperate runs with parallel/sequential quadtree construction manually.

setNumberOfThreads(processors)
scaleindex = 0
for scale in scalelist:
	n = 2**scale
	gen = generators.HyperbolicGenerator(n)#TODO: adapt m somehow
	gen.setTheoreticalSplit(True)
	minBindex = 0
	minCindex = 0
	bindex = 0
	for balance in balancelist:
		gen.setBalance(balance)
		cindex = 0
		for capacityscale in capacitylist:
			capacity = 2**capacityscale
			gen.setLeafCapacity(capacity)
			print("Starting ", iterations, " runs with", n, " nodes, capacity ", capacity, " and balance ", balance)
			start = time.time()
			for i in range(0, iterations):
				G = gen.generate()
			end = time.time()
			cost.append((end-start)/5)
			if (end-start)/5 < cost[scaleindex * len(balancelist) * len(capacitylist) + minBindex * len(capacitylist) + minCindex]:
				minBindex = bindex
				minCindex = cindex
			with open(str(stamp)+'-cost.json', 'w') as j:
				json.dump(cost, j)
			print("Finished iterations in ", end - start, " seconds")
			cindex = cindex + 1
			sys.stdout.flush()
		bindex = bindex + 1
	print("Optimal performance at balance ", balancelist[minBindex], " and capacity ", 2**capacitylist[minCindex], ".")
	scaleindex = scaleindex + 1
#TODO: find optimal value for each n, second pass

threadcost = []
# Weak scaling
minlength = min(len(scalelist, threadlist))
for i in range(0,minlength):
	scale = scalelist[i]
	threads = 2**threadlist[i]
	setNumberOfThreads(threads)
	n = 2**scale
	gen = generators.HyperbolicGenerator(n)
	print("Starting ", iterations, " runs with", n, " nodes, and ", threads, " threads.")
	start = time.time()
	for i in range(0, iterations):
		G = gen.generate()
	end = time.time()
	threadcost.append((end-start)/5)
	print("Finished iterations in ", end - start, " seconds")
	
# Strong scaling
for scale in scalelist:
	n = 2**scale
	gen = generators.HyperbolicGenerator(n)
	for th in threadlist:
	threads = 2**threadlist[i]
	setNumberOfThreads(threads)
	print("Starting ", iterations, " runs with", n, " nodes, and ", threads, " threads.")
	start = time.time()
	for i in range(0, iterations):
		G = gen.generate()
	end = time.time()
	threadcost.append((end-start)/5)
	print("Finished iterations in ", end - start, " seconds")

