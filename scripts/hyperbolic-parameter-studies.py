#! /usr/bin/env python3.4

from networkit import *
import time
import json

cost = []
numberOfEdges = []
powerLaw = []
diameters = []
clustCoeff = []
degeneracy = []
degreeAss = []
powerLawExponents = []
ncoms = []
mods = []

alphalist = [0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0] # creating lists to dump in json
alphalist.reverse()
stretchlist = list(range(5,21))
factorlist = list(range(1,11))
maglist = list(range(3,8))
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
			fit = properties.powerLawFit(G)
			m = G.numberOfEdges()
			diameter = properties.Diameter.estimatedDiameterRange(G)[0]
			clusteringCoefficient = properties.clustering(G)
			degen = properties.degeneracy(G)
			degAss = properties.degreeAssortativity(G)
			gamma = properties.powerLawExponent(G)
			plm = community.PLM()
			coms = plm.run(G)
			communities = coms.numberOfSubsets()
			modularity = community.Modularity().getQuality(coms, G)
				
			cost.append(end - start)
			numberOfEdges.append(m)
			powerLaw.append(fit[1])
			diameters.append(diameter)
			clustCoeff.append(clusteringCoefficient)
			degeneracy.append(degen)
			degreeAss.append(degAss)
			powerLawExponents.append(gamma)
			ncoms.append(communities)
			mods.append(modularity)

			#with open(str(stamp)+'-result-'+str(n)+'.csv', 'a') as f:
			#	f.write(str(n) + "," + str(alpha) + "," + str(factor/10) + "," + str(end-start) + "," + str(m) + "," + str(fit[1]) + "\n")

			with open(str(stamp)+'-cost-'+str(n)+'.json', 'w') as j:
				json.dump(cost, j)
			with open(str(stamp)+'-edges-'+str(n)+'.json', 'w') as j:
				json.dump(numberOfEdges, j)
			with open(str(stamp)+'-powerLaw-'+str(n)+'.json', 'w') as j:
				json.dump(powerLaw, j)
			with open(str(stamp)+'-diameter-'+str(n)+'.json', 'w') as j:
				json.dump(diameters, j)
			with open(str(stamp)+'-clusterCoefficients-'+str(n)+'.json', 'w') as j:
				json.dump(clustCoeff, j)
			with open(str(stamp)+'-degeneracy-'+str(n)+'.json', 'w') as j:
				json.dump(degeneracy, j)
			with open(str(stamp)+'-degAssortativity-'+str(n)+'.json', 'w') as j:
				json.dump(degreeAss, j)
			with open(str(stamp)+'-powerLawExponent-'+str(n)+'.json', 'w') as j:
				json.dump(powerLawExponents, j)
			with open(str(stamp)+'-ncoms-'+str(n)+'.json', 'w') as j:
				json.dump(ncoms, j)
			with open(str(stamp)+'-modularities-'+str(n)+'.json', 'w') as j:
				json.dump(mods, j)

			print("Finished run in ", end - start, " seconds, resulting in ", G.numberOfEdges(), " edges. AvgDegree: ", 2*m / n, ", powerLaw:" , fit[1], ", mod:", modularity)

