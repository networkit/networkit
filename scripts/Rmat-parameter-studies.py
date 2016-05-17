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

edgeFactorlist = list(range(1,200))
alist = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
scalelist = list(range(10,21))
stamp = time.time()

# dumping parameters

with open(str(stamp)+'-scale.json', 'w') as j:
	json.dump(scalelist, j)
with open(str(stamp)+'-eF.json', 'w') as j:
	json.dump(edgeFactorlist, j)
with open(str(stamp)+'-alist.json', 'w') as j:
	json.dump(alist, j)

factor = 10
sum = 0

for scale in scalelist:
	n = 2**scale
	for ef in edgeFactorlist:
		for a in alist:
			sum = a
			for bscaled in range(0,int(min(a*10,(1-sum)*10))+1):
				b = bscaled/10
				sum = a+b
				for cscaled in range(0,int(min(a*10,(1-sum)*10))+1):
					c = cscaled/10
					sum = a+b+c
					if sum > 1 or 1-sum > a:
						continue
					d = 1-sum
					print("Starting run with", n, " nodes and", ef*n, "edges, a=", a, ",b=", b, ", c=",c," and d=",d)
					start = time.time()
					G = generators.RmatGenerator(scale,ef, a,b,c,d).generate()
					end = time.time()
					fit = properties.powerLawFit(G)
					m = G.numberOfEdges()
					diameter = properties.Diameter.estimatedDiameterRange(G)[0]
					print("Diameter:", diameter) 
					clusteringCoefficient = properties.clustering(G)
					print("cc:", clusteringCoefficient) 
					#degen = properties.degeneracy(G)
					#print("K-core:", degen) 
					degAss = properties.degreeAssortativity(G)
					print("Degree Assortativity:", degAss) 
					gamma = properties.powerLawExponent(G)
					print("PowerLaw Exponent:", gamma) 
					plm = community.PLM()
					coms = plm.run(G)
					communities = coms.numberOfSubsets()
					print("Communities:", communities) 
					modularity = community.Modularity().getQuality(coms, G)
					print("Modularity:", modularity) 
				
					cost.append(end - start)
					numberOfEdges.append(m)
					powerLaw.append(fit[1])
					diameters.append(diameter)
					clustCoeff.append(clusteringCoefficient)
					#degeneracy.append(degen)
					degreeAss.append(degAss)
					powerLawExponents.append(gamma)
					ncoms.append(communities)
					mods.append(modularity)

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
					#with open(str(stamp)+'-degeneracy-'+str(n)+'.json', 'w') as j:
					#	json.dump(degeneracy, j)
					with open(str(stamp)+'-degAssortativity-'+str(n)+'.json', 'w') as j:
						json.dump(degreeAss, j)
					with open(str(stamp)+'-powerLawExponent-'+str(n)+'.json', 'w') as j:
						json.dump(powerLawExponents, j)
					with open(str(stamp)+'-ncoms-'+str(n)+'.json', 'w') as j:
						json.dump(ncoms, j)
					with open(str(stamp)+'-modularities-'+str(n)+'.json', 'w') as j:
						json.dump(mods, j)

					print("Finished run in ", end - start, " seconds, resulting in ", G.numberOfEdges(), " edges. AvgDegree: ", 2*m / n, ", powerLaw:" , fit[1])
