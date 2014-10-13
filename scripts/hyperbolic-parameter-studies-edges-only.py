#! /usr/bin/env python3.4

from networkit import *
import time
import json

cost = []
numberOfEdges = []

alphalist = [0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6] # creating lists to dump in json
alphalist.reverse()
stretchlist = [0.6,0.8,1.0,1.2,1.4,1.6]
stretchlist.reverse()
factorlist = list(range(11,21))
maglist = list(range(6,9))
stamp = time.time()

# dumping parameters

with open(str(stamp)+'-alpha.json', 'w') as j:
	json.dump(alphalist, j)
with open(str(stamp)+'-stretch.json', 'w') as j:
	json.dump(stretchlist, j)
with open(str(stamp)+'-mag.json', 'w') as j:
	json.dump(maglist, j)

factor = 10

for magnitude in maglist:
	n = 10**magnitude
	for alpha in alphalist:
		for stretch in stretchlist:
			print("Starting run with", n, " nodes, alpha ", alpha, "factor ", factor/10, " and stretch ", stretch/10)
			start = time.time()
			G = generators.HyperbolicGenerator(n,factor/10,alpha,stretch/10).generate()
			end = time.time()
			m = G.numberOfEdges()
				
			cost.append(end - start)
			numberOfEdges.append(m)

			with open(str(stamp)+'-cost-'+str(n)+'.json', 'w') as j:
				json.dump(cost, j)
			with open(str(stamp)+'-edges-'+str(n)+'.json', 'w') as j:
				json.dump(numberOfEdges, j)

			print("Finished run in ", end - start, " seconds, resulting in ", G.numberOfEdges(), " edges. AvgDegree: ", 2*m / n)
