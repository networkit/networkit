from types import *
from networkit import *

class bb_SimmelianBackboneNonParametric:
    
    def getShortName(self):
        return "Simmelian NonParametric"

    def getAlgorithm(self):
        return backbones.SimmelianBackboneNonParametric()

    def requiresWeight(self):
        return False
    
    def parameterizationType(self):
        return "Default"

# -----------------------------------------------------------

class OriginalAlgorithm(backbones.BackboneAlgorithm):
    
    def getAttribute(self, G):
        return None

    def _getBackbone(self, G, parameter, attribute):
        return G
 
    def _getParameterizationAlgorithm(self):
        return backbones.SimpleParameterization()

class bb_Original:
    
    def getShortName(self):
        return "Original"

    def requiresWeight(self):
        return False
    
    def getAlgorithm(self):
        return OriginalAlgorithm()
    
    def parameterizationType(self):
        return "None"

# -----------------------------------------------------------

class bb_SimmelianMultiscale:
    
    def getShortName(self):
        return "Simmelian Multiscale"

    def getAlgorithm(self):
        return None #TODO

    def requiresWeight(self):
        return False
    
    def parameterizationType(self):
        return "Default"
# -----------------------------------------------------------

class bb_SimmelianBackboneParametric:

    def getShortName(self):
        return "Simmelian Parametric"
    
    def getAlgorithm(self):
        return backbones.SimmelianBackboneParametric()

    def requiresWeight(self):
        return False
    
    def parameterizationType(self):
        return "Default"

# -----------------------------------------------------------

class bb_LocalSimilarity:

    def getShortName(self):
        return "Local Similarity"
    
    def getAlgorithm(self):
        return backbones.LocalSimilarityBackbone()

    def requiresWeight(self):
        return False
    
    def parameterizationType(self):
        return "Default"

# -----------------------------------------------------------

class bb_Multiscale:
    
    def getShortName(self):
        return "Multiscale"

    def getAlgorithm(self):
        return backbones.MultiscaleBackbone()
    
    def requiresWeight(self):
        return True
    
    def parameterizationType(self):
        return "Default"

# -----------------------------------------------------------

class bb_Random:
    def __init__(self, tag):
        self._tag = tag

    def getShortName(self):
        return ("Random " + self._tag).strip()

    def getAlgorithm(self):
        return backbones.RandomBackbone()

    def requiresWeight(self):
        return False
    
    def parameterizationType(self):
        return "Default"

# -----------------------------------------------------------

class bb_ForestFire:

    def __init__(self, tag, pf, tber):
        self.tag = tag
        self.pf = pf
        self.tber = tber

    def getShortName(self):
        return ("ForestFire " + self.tag).strip()

    def getAlgorithm(self):
        return backbones.ForestFireBackbone(pf, tber)
    
    def requiresWeight(self):
        return False
    
    def parameterizationType(self):
        return "Default"

# -----------------------------------------------------------

class bb_LocalDegree:

    def getShortName(self):
        return "Local Degree"

    def getAlgorithm(self):
        return backbones.LocalDegreeBackbone()
    
    def requiresWeight(self):
        return False
    
    def parameterizationType(self):
        return "Default"