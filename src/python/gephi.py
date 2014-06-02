""" TODO: module docstring """

import os

def exportNodeValues(values, fpath, column_name):
	""" TODO: docstring """
	f = open(fpath, "w+")
	f.write("id,{0}\n".format(column_name))
	for i in range(0,len(values)):
		f.write("{0},{1}\n".format(i,values[i]))
	f.close()
