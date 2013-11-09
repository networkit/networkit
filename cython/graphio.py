from _NetworKit import (METISGraphReader, FastMETISGraphReader, DotGraphWriter, EdgeListIO, \
						 LineFileReader, SNAPGraphWriter, ClusteringReader, ClusteringWriter)

class GraphConverter:
    
    def __init__(self, reader, writer):
        self.reader = reader
        self.writer = writer
        
    def convert(self, inPath, outPath):
        G = self.reader.read(inPath)
        self.writer.write(G, outPath)
        
    def __str__(self):
        return "GraphConverter: {0} => {0}".format(self.reader, self.writer)

def getConverter(fromFormat, toFormat):
    
    readers = {"metis": METISGraphReader, "edgelist" : EdgeListIO}    
    writers = {"edgelist": EdgeListIO}
    
    reader = readers[fromFormat]()
    writer = writers[toFormat]()
    
    return GraphConverter(reader, writer)