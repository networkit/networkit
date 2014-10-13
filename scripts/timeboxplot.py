#! /usr/bin/env python3

import json
import sys
import os.path

stamp = sys.argv[1]

# load parameters
with open(str(stamp)+'-mag.json', 'r') as j:
	maglist = json.load(j)

classlist = ['a','b','c','d','e']
magindex = 0
offset = 0

with open(str(stamp)+'-table-eps-avg.dat', 'w') as eps:
	with open(str(stamp)+'-table-scatter.dat', 'w') as scatter:
		scatter.write('x\ty\tlabel\n')
		eps.write('x\ty\tlabel\n')
		for mag in maglist:
			# load data
			n = 10**mag
			if not os.path.isfile(str(stamp)+'-cost-'+str(n)+'.json'):
				print ("File ", str(stamp),'-cost-',str(n),'.json', 'does not exist. aborting.')
				break
			with open(str(stamp)+'-cost-'+str(n)+'.json', 'r') as j:
				cost = json.load(j)
			with open(str(stamp)+'-edges-'+str(n)+'.json', 'r') as j:
				edges = json.load(j)
			print("Reformating", str(stamp),",",str(n),".")
			with open(str(stamp)+'-cost-'+str(n)+'.json', 'r') as j:
				cost = json.load(j)
			with open(str(stamp)+'-edges-'+str(n)+'.json', 'r') as j:
				edges = json.load(j)

			for i in range(offset,len(edges)):
				scatter.write(str(edges[i]) + '\t' + str(cost[i]) + '\t' + classlist[magindex] + '\n')
				eps.write(str(edges[i]/n) + '\t' + str(edges[i]/cost[i]) + '\t' + classlist[magindex] + '\n')
			offset = len(edges) #workaround for a bug in the parameter studies code: lists are not emptied when going to higher magnitude
			scatter.write('\n')
			eps.write('\n')
			magindex = magindex + 1

