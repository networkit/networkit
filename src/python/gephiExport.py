import os
from structures import Partition

def writeCSV(values, fpath, column_name):
	f = open(fpath, "w+")
	f.write("id,{0}\n".format(column_name))
	for i in range(0,len(values)):
		f.write("{0},{1}\n".format(i,values[i]))
	f.close()
