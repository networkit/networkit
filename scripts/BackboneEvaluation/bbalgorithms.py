from types import *
from NetworKit import *

class bb_SimmelianBackboneNonParametric:
    def getName(self):
        return "SimmelianBackboneNonParametric"

    def getShortName(self):
        return "Simmelian NonParametric"

    def getAlgorithmExpr(self, parameter):
        return "backbones.SimmelianBackboneNonParametric(" + str(parameter) + ")"

    def getAttribute(self, graph):
        chiba = backbones.ChibaNishizekiTriangleCounter()
        triangles = chiba.getAttribute(graph)
        sj = backbones.SimmelianJaccardAttributizer()
        a_sj = sj.getAttribute(graph, triangles)
        return a_sj

    def getBackboneFromAttribute(self, graph, attribute, value):
        gf = backbones.GlobalThresholdFilter(value, True)
        return gf.calculate(graph, attribute)

    def parameterizationType(self):
        return "Float"

    def increasing(self):
        return False

    def requiresWeight(self):
        return False

# -----------------------------------------------------------

class bb_Original:
    def getName(self):
        return "Original Graph"

    def getShortName(self):
        return "Original"

    def getAlgorithmExpr(self, parameter):
        return "None"

    def requiresWeight(self):
        return False

    def getAttribute(self, graph):
        return None

    def getBackboneFromAttribute(self, graph, attribute, value):
        return graph

    def parameterizationType(self):
        return "None"   #Algorithm does not take a parameter as input.

# -----------------------------------------------------------

class bb_SimmelianMultiscaleBackbone:
    def getName(self):
        return "SimmelianMultiscaleBackbone"

    def getShortName(self):
        return "Simmelian Multiscale"

    def getAlgorithmExpr(self, parameter):
        return "backbones.SimmelianMultiscaleBackbone(" + str(parameter) + ")"

    def getAttribute(self, graph):
        chiba = backbones.ChibaNishizekiTriangleCounter()
        triangles = chiba.getAttribute(graph)
        ms = backbones.MultiscaleAttributizer()
        a_ms = ms.getAttribute(graph, triangles)
        return a_ms

    def getBackboneFromAttribute(self, graph, attribute, value):
        gf = backbones.GlobalThresholdFilter(value, False)
        return gf.calculate(graph, attribute)

    def parameterizationType(self):
        return "Float"

    def increasing(self):
        return True

    def requiresWeight(self):
        return False
# -----------------------------------------------------------

class bb_SimmelianBackboneParametric:
    def getName(self):
        return "SimmelianBackboneParametric"

    def getShortName(self):
        return "Simmelian Parametric"

    def getAlgorithmExpr(self, parameter):
        return "backbones.SimmelianBackboneParametric(10, " + str(parameter) + ")"

    def getAttribute(self, graph):
        chiba = backbones.ChibaNishizekiTriangleCounter()
        triangles = chiba.getAttribute(graph)
        so = backbones.SimmelianOverlapAttributizer(10)
        a_so = so.getAttribute(graph, triangles)
        return a_so

    def getBackboneFromAttribute(self, graph, attribute, value):
        gf = backbones.GlobalThresholdFilter(value, True)
        return gf.calculate(graph, attribute)

    def parameterizationType(self):
        return "Int"

    def requiresWeight(self):
        return False

# -----------------------------------------------------------

class bb_LocalSimilarityBackbone:
    def getName(self):
        return "LocalSimilarityBackbone"

    def getShortName(self):
        return "Local Similarity"

    def getAlgorithmExpr(self, parameter):
        return "backbones.LocalSimilarityBackbone(" + str(parameter) + ")"

    def getAttribute(self, graph):
        attributizer = backbones.LocalSimilarityAttributizer()
        a_ls = attributizer.getAttribute(graph, [])
        return a_ls

    def getBackboneFromAttribute(self, graph, attribute, value):
        gf = backbones.GlobalThresholdFilter(value, False)
        return gf.calculate(graph, attribute)

    def parameterizationType(self):
        return "Float"

    def increasing(self):
        return False

    def requiresWeight(self):
        return False

# -----------------------------------------------------------

class bb_MultiscaleBackbone:
    def getName(self):
        return "MultiscaleBackbone"

    def getShortName(self):
        return "Multiscale"

    def getAlgorithmExpr(self, parameter):
        return "backbones.MultiscaleBackbone(" + str(parameter) + ")"

    def getAttribute(self, graph):
        #TODO we might use a precalculated edge attribute for speedup, but that
        # requires writable edge attributes in python.
        return None

    def getBackboneFromAttribute(self, graph, attribute, value):
        msb = backbones.MultiscaleBackbone(value)
        return msb.calculate(graph)

    def parameterizationType(self):
        return "Float"

    def increasing(self):
        return True

    def requiresWeight(self):
        return True

# -----------------------------------------------------------

class bb_RandomBackbone:
    def __init__(self, tag):
        self._tag = tag

    def getName(self):
        return "RandomBackbone " + self._tag

    def getShortName(self):
        return "Random " + self._tag

    def getAlgorithmExpr(self, parameter):
        return "backbones.RandomBackbone(" + str(parameter) + ")"

    def getAttribute(self, graph):
        return None

    def getBackboneFromAttribute(self, graph, attribute, value):
        rb = backbones.RandomBackbone(value)
        return rb.calculate(graph)

    def parameterizationType(self):
        return "Trivial"  #Trivial: No parameterizitation needed.

    def increasing(self):
        return True

    def requiresWeight(self):
        return False
