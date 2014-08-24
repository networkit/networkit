from _NetworKit import SimmelianBackboneParametric, SimmelianBackboneNonParametric, MultiscaleBackbone, LocalSimilarityBackbone, SimmelianMultiscaleBackbone, ChibaNishizekiTriangleCounter, SimmelianJaccardAttributizer, GlobalThresholdFilter, LocalSimilarityAttributizer, MultiscaleAttributizer, SimmelianOverlapAttributizer, RandomAttributizer, RandomBackbone

def writeCsv(file, G, B):
	f = open(file, "w+")
	f.write("Source,Target,Type,{0}\n".format("Backbone"))
	for edge in G.edges():
		value = 1 if edge in B.edges() else 0
		f.write("{0},{1},{2},{3}\n".format(edge[0],edge[1],"Undirected",value))
	f.close()
