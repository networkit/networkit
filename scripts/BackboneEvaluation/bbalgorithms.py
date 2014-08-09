from types import *
from NetworKit import *

class bb_SimmelianBackboneNonParametric:
    def getName(self):
        return "SimmelianBackboneNonParametric"

    def getShortName(self, parameter):
        return "SBNP " + str(parameter)

    def getAlgorithmExpr(self, parameter):
        return "backbones.SimmelianBackboneNonParametric(" + str(parameter) + ")"

    def getPrecalcAttribute(self, graph):
        chiba = backbones.ChibaNishizekiTriangleCounter()
        empty = backbones.EdgeAttribute(0)
        triangles = chiba.getAttribute(graph, empty)
        sj = backbones.SimmelianJaccardAttributizer()
        a_sj = sj.getAttribute(graph, triangles)
        return a_sj

    def getPrecalcBackbone(self, graph, attribute, value):
        gf = backbones.GlobalThresholdFilter(value, True)
        return gf.calculate(graph, attribute)

    def parameterType(self):
        return "FloatType"

    def increasing(self):
        return False

# -----------------------------------------------------------

class bb_Original:
    def getName(self):
        return "Original Graph"

    def getShortName(self):
        return "Original"

    def getAlgorithmExpr(self):
        return "OriginalGraph()"

# -----------------------------------------------------------

class bb_SimmelianMultiscaleBackbone:
    def getName(self):
        return "SimmelianMultiscaleBackbone"

    def getShortName(self, parameter):
        return "SMB " + str(parameter)

    def getAlgorithmExpr(self, parameter):
        return "backbones.SimmelianMultiscaleBackbone(" + str(parameter) + ")"

    def getPrecalcAttribute(self, graph):
        chiba = backbones.ChibaNishizekiTriangleCounter()
        empty = backbones.EdgeAttribute(0)
        triangles = chiba.getAttribute(graph, empty)
        ms = backbones.MultiscaleAttributizer()
        a_ms = ms.getAttribute(graph, triangles)
        return a_ms

    def getPrecalcBackbone(self, graph, attribute, value):
        gf = backbones.GlobalThresholdFilter(value, False)
        return gf.calculate(graph, attribute)

    def parameterType(self):
        return "FloatType"

    def increasing(self):
        return True
# -----------------------------------------------------------

class bb_SimmelianBackboneParametric:
    def getName(self):
        return "SimmelianBackboneParametric"

    def getShortName(self, parameter):
        return "SBP " + str(parameter)

    def getAlgorithmExpr(self, parameter):
        return "backbones.SimmelianBackboneParametric(10, " + str(parameter) + ")"

    def getPrecalcAttribute(self, graph):
        chiba = backbones.ChibaNishizekiTriangleCounter()
        empty = backbones.EdgeAttribute(0)
        triangles = chiba.getAttribute(graph, empty)
        so = backbones.SimmelianOverlapAttributizer(10)
        a_so = so.getAttribute(graph, triangles)
        return a_so

    def getPrecalcBackbone(self, graph, attribute, value):
        gf = backbones.GlobalThresholdFilter(value, True)
        return gf.calculate(graph, attribute)

    def parameterType(self):
        return "IntType"

# -----------------------------------------------------------

class bb_LocalSimilarityBackbone:
    def getName(self):
        return "LocalSimilarityBackbone"

    def getShortName(self, parameter):
        return "LSB " + str(parameter)

    def getAlgorithmExpr(self, parameter):
        return "backbones.LocalSimilarityBackbone(" + str(parameter) + ")"

    def getPrecalcAttribute(self, graph):
        empty = backbones.EdgeAttribute(0)
        attributizer = backbones.LocalSimilarityAttributizer()
        a_ls = attributizer.getAttribute(graph, empty)
        return a_ls

    def getPrecalcBackbone(self, graph, attribute, value):
        gf = backbones.GlobalThresholdFilter(value, True)
        return gf.calculate(graph, attribute)

    def parameterType(self):
        return "FloatType"

    def increasing(self):
        return False
