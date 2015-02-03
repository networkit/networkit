from types import *
from networkit import *

class S_SimmelianBackboneNonParametric:

    def getShortName(self):
        return "Simmelian NonParametric"

    def getAlgorithm(self):
        return sparsification.SimmelianBackboneNonParametric()

    def requiresWeight(self):
        return False

    def parameterizationType(self):
        return "Default"

# -----------------------------------------------------------

class OriginalAlgorithm(sparsification.Sparsifier):

    def getAttribute(self, G):
        return None

    def _getSparsifiedGraph(self, G, parameter, attribute):
        return G

    def _getParameterizationAlgorithm(self):
        return sparsification.SimpleParameterization()

class S_Original:

    def getShortName(self):
        return "Original"

    def requiresWeight(self):
        return False

    def getAlgorithm(self):
        return OriginalAlgorithm()

    def parameterizationType(self):
        return "None"

# -----------------------------------------------------------

class S_SimmelianMultiscale:

    def getShortName(self):
        return "Simmelian Multiscale"

    def getAlgorithm(self):
        return sparsification.SimmelianMultiscaleBackbone()

    def requiresWeight(self):
        return False

    def parameterizationType(self):
        return "Default"
# -----------------------------------------------------------

class S_SimmelianBackboneParametric:

    def __init__(self, maxRank):
        self.maxRank = maxRank

    def getShortName(self):
        return "Simmelian Parametric"

    def getAlgorithm(self):
        return sparsification.SimmelianBackboneParametric(self.maxRank)

    def requiresWeight(self):
        return False

    def parameterizationType(self):
        return "Default"

# -----------------------------------------------------------

class S_LocalSimilarity:

    def getShortName(self):
        return "Local Similarity"

    def getAlgorithm(self):
        return sparsification.LocalSimilarityBackbone()

    def requiresWeight(self):
        return False

    def parameterizationType(self):
        return "Default"

# -----------------------------------------------------------

class S_Multiscale:

    def getShortName(self):
        return "Multiscale"

    def getAlgorithm(self):
        return sparsification.MultiscaleBackbone()

    def requiresWeight(self):
        return True

    def parameterizationType(self):
        return "Default"

# -----------------------------------------------------------

class S_Random:
    def __init__(self, tag):
        self._tag = tag

    def getShortName(self):
        return ("Random " + self._tag).strip()

    def getAlgorithm(self):
        return sparsification.RandomBackbone()

    def requiresWeight(self):
        return False

    def parameterizationType(self):
        return "Default"

# -----------------------------------------------------------

class S_ForestFire:

    def __init__(self, tag, pf, tber):
        self.tag = tag
        self.pf = pf
        self.tber = tber

    def getShortName(self):
        return ("ForestFire " + self.tag).strip()

    def getAlgorithm(self):
        return sparsification.ForestFireBackbone(self.pf, self.tber)

    def requiresWeight(self):
        return False

    def parameterizationType(self):
        return "Default"

# -----------------------------------------------------------

class S_LocalDegree:

    def getShortName(self):
        return "Local Degree"

    def getAlgorithm(self):
        return sparsification.LocalDegreeBackbone()

    def requiresWeight(self):
        return False

    def parameterizationType(self):
        return "Default"
