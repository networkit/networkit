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

maglist = list(range(3,7))

stamp = time.time()

# dumping parameters

with open(str(stamp)+'-mag.json', 'w') as j:
	json.dump(maglist, j)

stretch = 10

for magnitude in maglist:
	n = 10**magnitude
	print("Starting run with", n, " nodes")
	start = time.time()
	G = generators.DorogovtsevMendesGenerator(n).generate()
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

	with open(str(stamp)+'-cost.json', 'w') as j:
		json.dump(cost, j)
	with open(str(stamp)+'-edges.json', 'w') as j:
		json.dump(numberOfEdges, j)
	with open(str(stamp)+'-powerLaw.json', 'w') as j:
		json.dump(powerLaw, j)
	with open(str(stamp)+'-diameter.json', 'w') as j:
		json.dump(diameters, j)
	with open(str(stamp)+'-clusterCoefficients.json', 'w') as j:
		json.dump(clustCoeff, j)
	with open(str(stamp)+'-degeneracy.json', 'w') as j:
		json.dump(degeneracy, j)
	with open(str(stamp)+'-degAssortativity.json', 'w') as j:
		json.dump(degreeAss, j)
	with open(str(stamp)+'-powerLawExponent.json', 'w') as j:
		json.dump(powerLawExponents, j)
	with open(str(stamp)+'-ncoms.json', 'w') as j:
		json.dump(ncoms, j)
	with open(str(stamp)+'-modularities.json', 'w') as j:
		json.dump(mods, j)

	print("Finished run in ", end - start, " seconds, resulting in ", G.numberOfEdges(), " edges. AvgDegree: ", 2*m / n, ", powerLaw:" , fit[1], ", mod:", modularity)

