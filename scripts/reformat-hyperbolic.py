#! /usr/bin/env python3.4

import json
import sys
import os.path

stamp = sys.argv[1]

# load parameters
with open(str(stamp)+'-alpha.json', 'r') as j:
	alphalist = json.load(j)
with open(str(stamp)+'-stretch.json', 'r') as j:
	stretchlist = json.load(j)
with open(str(stamp)+'-mag.json', 'r') as j:
	maglist = json.load(j)

i = 0
for mag in maglist:
	# load data
	n = 10**mag
	if not os.path.isfile(str(stamp)+'-cost-'+str(n)+'.json'):
		print ("File ", str(stamp),'-cost-',str(n),'.json', 'does not exist. aborting.')
		break
	print("Reformating", str(stamp),",",str(n),".")
	with open(str(stamp)+'-cost-'+str(n)+'.json', 'r') as j:
		cost = json.load(j)
	with open(str(stamp)+'-edges-'+str(n)+'.json', 'r') as j:
		edges = json.load(j)
	with open(str(stamp)+'-powerLaw-'+str(n)+'.json', 'r') as j:
		powerLaw = json.load(j)
	with open(str(stamp)+'-clusterCoefficients-'+str(n)+'.json', 'r') as j:
		clusterCoeff = json.load(j)
	with open(str(stamp)+'-degAssortativity-'+str(n)+'.json', 'r') as j:
		degreeAss = json.load(j)
	with open(str(stamp)+'-powerLawExponent-'+str(n)+'.json', 'r') as j:
		powerLawExponents = json.load(j)
	with open(str(stamp)+'-diameter-'+str(n)+'.json', 'r') as j:
		diameters = json.load(j)
	with open(str(stamp)+'-degeneracy-'+str(n)+'.json', 'r') as j:
		degeneracy = json.load(j)
	with open(str(stamp)+'-ncoms-'+str(n)+'.json', 'r') as j:
		ncoms = json.load(j)
	with open(str(stamp)+'-modularities-'+str(n)+'.json', 'r') as j:
		mods = json.load(j)

	with open(str(stamp)+'-table-edges-'+str(n)+'.dat', 'w') as e:
		with open(str(stamp)+'-table-cost-'+str(n)+'.dat', 'w') as c:
			with open(str(stamp)+'-table-powerlaw-'+str(n)+'.dat', 'w') as p:
				with open(str(stamp)+'-table-clusterCoeff-'+str(n)+'.dat', 'w') as cc:
					with open(str(stamp)+'-table-assortativity-'+str(n)+'.dat', 'w') as a:
						with open(str(stamp)+'-table-powerLawExp-'+str(n)+'.dat', 'w') as exp:			
							with open(str(stamp)+'-table-diameter-'+str(n)+'.dat', 'w') as diam:
								with open(str(stamp)+'-table-degeneracy-'+str(n)+'.dat', 'w') as degen:
									with open(str(stamp)+'-table-ncoms-'+str(n)+'.dat', 'w') as nc:
										with open(str(stamp)+'-table-modularities-'+str(n)+'.dat', 'w') as mo:
											with open(str(stamp)+'-table-comsize-'+str(n)+'.dat', 'w') as cs:
												for alpha in alphalist:
													for stretch in stretchlist:
														if (i < len(edges)):
															e.write(str(alpha) + '\t' + str(stretch/10) + '\t' + str(edges[i]) + '\n')
															c.write(str(alpha) + '\t' + str(stretch/10) + '\t' + str(cost[i]) + '\n')
															p.write(str(alpha) + '\t' + str(stretch/10) + '\t' + str(powerLaw[i]) + '\n')
															cc.write(str(alpha) + '\t' + str(stretch/10) + '\t' + str(clusterCoeff[i]) + '\n')
															a.write(str(alpha) + '\t' + str(stretch/10) + '\t' + str(degreeAss[i]) + '\n')
															exp.write(str(alpha) + '\t' + str(stretch/10) + '\t' + str(powerLawExponents[i]) + '\n')
															diam.write(str(alpha) + '\t' + str(stretch/10) + '\t' + str(diameters[i]) + '\n')
															degen.write(str(alpha) + '\t' + str(stretch/10) + '\t' + str(degeneracy[i]) + '\n')
															nc.write(str(alpha) + '\t' + str(stretch/10) + '\t' + str(ncoms[i]) + '\n')
															mo.write(str(alpha) + '\t' + str(stretch/10) + '\t' + str(mods[i]) + '\n')
															cs.write(str(alpha) + '\t' + str(stretch/10) + '\t' + str(int(n/ncoms[i])) + '\n')
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
			
				
