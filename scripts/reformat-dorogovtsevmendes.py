#! /usr/bin/env python3.4

import json
import sys
import os.path

stamp = sys.argv[1]

# load parameters
with open(str(stamp)+'-mag.json', 'r') as j:
	maglist = json.load(j)

with open(str(stamp)+'-cost.json', 'r') as j:
	cost = json.load(j)
with open(str(stamp)+'-edges.json', 'r') as j:
	edges = json.load(j)
with open(str(stamp)+'-powerLaw.json', 'r') as j:
	powerLaw = json.load(j)
with open(str(stamp)+'-clusterCoefficients.json', 'r') as j:
	clusterCoeff = json.load(j)
with open(str(stamp)+'-degAssortativity.json', 'r') as j:
	degreeAss = json.load(j)
with open(str(stamp)+'-powerLawExponent.json', 'r') as j:
	powerLawExponents = json.load(j)
with open(str(stamp)+'-diameter.json', 'r') as j:
	diameters = json.load(j)
with open(str(stamp)+'-degeneracy.json', 'r') as j:
	degeneracy = json.load(j)
with open(str(stamp)+'-ncoms.json', 'r') as j:
	ncoms = json.load(j)
with open(str(stamp)+'-modularities.json', 'r') as j:
	mods = json.load(j)


if not os.path.isfile(str(stamp)+'-cost.json'):
	print ("File ", str(stamp),'-cost.json', 'does not exist. aborting.')
print("Reformating", str(stamp),".")

i = 0	
with open(str(stamp)+'-table-edges.dat', 'w') as e:
	with open(str(stamp)+'-table-cost.dat', 'w') as c:
		with open(str(stamp)+'-table-powerlaw.dat', 'w') as p:
			with open(str(stamp)+'-table-clusterCoeff.dat', 'w') as cc:
				with open(str(stamp)+'-table-assortativity.dat', 'w') as a:
					with open(str(stamp)+'-table-powerLawExp.dat', 'w') as exp:			
						with open(str(stamp)+'-table-diameter.dat', 'w') as diam:
							with open(str(stamp)+'-table-degeneracy.dat', 'w') as degen:
								with open(str(stamp)+'-table-ncoms.dat', 'w') as nc:
									with open(str(stamp)+'-table-modularities.dat', 'w') as mo:
										with open(str(stamp)+'-table-comsize.dat', 'w') as cs:
											for mag in maglist:
												n = 10**mag
													if (i < len(edges)):
														e.write(str(n) + '\t' + str(edges[i]) + '\n')
														c.write(str(n) + '\t' + str(cost[i]) + '\n')
														p.write(str(n) + '\t' + str(powerLaw[i]) + '\n')
														cc.write(str(n) + '\t' + str(clusterCoeff[i]) + '\n')
														a.write(str(n) + '\t' + str(degreeAss[i]) + '\n')
														exp.write(str(n) + '\t' + str(powerLawExponents[i]) + '\n')
														diam.write(str(n) + '\t' + str(diameters[i]) + '\n')
														degen.write(str(n) + '\t' + str(degeneracy[i]) + '\n')
														nc.write(str(n) + '\t' + str(ncoms[i]) + '\n')
														mo.write(str(n) + '\t' + str(mods[i]) + '\n')
														cs.write(str(n) + '\t' + str(int(n/ncoms[i])) + '\n')
														i = i+1
													else:
														break
												else:
													e.write('\n')
													c.write('\n')
													p.write('\n')
													cc.write('\n')
													a.write('\n')
													exp.write('\n')
													diam.write('\n')
													degen.write('\n')
													nc.write('\n')
													mo.write('\n')
													cs.write('\n')
			
				
