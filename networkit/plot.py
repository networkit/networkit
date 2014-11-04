from networkit import *
import matplotlib.pyplot as plt
import numpy as np

try:
	import seaborn
	seaborn.set_style("whitegrid")
except ImportError as importError:
	print("WARNING: module 'seaborn' is not installed, but recommended")
	print(importError)


def degreeDistribution(G):
	"""Plots the degree distribution of the given network."""
	dd = properties.degreeDistribution(G)
	plt.yscale("symlog")
	plt.xscale("log")
	plt.xlabel("nodes")
	plt.ylabel("degree")
	plt.plot(dd)

# TODO: hop plot
