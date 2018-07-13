"""
This module provides graph generators that produce synthetic networks according to various models.
"""

__author__ = "Christian Staudt"

# extension imports
from _NetworKit import BarabasiAlbertGenerator, PubWebGenerator, ErdosRenyiGenerator, ClusteredRandomGraphGenerator, DorogovtsevMendesGenerator, DynamicPubWebGenerator, DynamicPathGenerator, ChungLuGenerator, HyperbolicGenerator, MocnikGenerator, MocnikGeneratorBasic, DynamicHyperbolicGenerator, HavelHakimiGenerator, DynamicDorogovtsevMendesGenerator, RmatGenerator, DynamicForestFireGenerator, RegularRingLatticeGenerator, WattsStrogatzGenerator, PowerlawDegreeSequence, EdgeSwitchingMarkovChainGenerator, EdgeSwitchingMarkovChainGenerator as ConfigurationModelGenerator, LFRGenerator

from . import graphio

import subprocess
import os
import tempfile
import scipy

class BTERReplicator:
	"""
	Wrapper class that calls the BTER graph generator implementation in
	FEASTPACK from http://www.sandia.gov/~tgkolda/feastpack/ using GNU
	Octave.

	Note that BTER needs the rng method which is unavailable in Octave, but
	the call in bter.m can be easily replaced.
	"""
	matlabname = 'octave'
	matlabScript = """
	addpath('{0}');
	filename = 'bter_input.mat';
	load(filename);
	addpath('{1}');
	tStart = tic;
	[ccd,gcc] = ccperdeg(G);
	nd = accumarray(nonzeros(sum(G,2)),1);
	nd = nd * {2};
	tFit = toc(tStart);
	tStart = tic;
	[E1,E2] = bter(nd,ccd,'verbose',false,'blowup',10);
	tGenerate = toc(tStart);
	G_bter = bter_edges2graph(E1,E2);
	save('-v7', '{3}', 'G_bter', 'tFit', 'tGenerate');
	exit;
	"""
	feastpackPath = "."


	@classmethod
	def setPaths(cls, feastpackPath):
		cls.feastpackPath = feastpackPath

	def __init__(self, G, scale=1):
		self.G = G
		self.scale = scale

	def generate(self):
		with tempfile.TemporaryDirectory() as tmpdir:
			scriptPath = os.path.join(tmpdir, "bter_wrapper.m")
			tempFileOut = os.path.join(tmpdir, 'bter_output.mat')
			tempFileIn = os.path.join(tmpdir, 'bter_input.mat')
			# write MATLAB script
			with open(scriptPath, 'w') as matlabScriptFile:
				matlabScriptFile.write(self.matlabScript.format(tmpdir, self.feastpackPath, self.scale, tempFileOut))
			graphio.writeMat(self.G, tempFileIn)
			subprocess.call([self.matlabname, '-qf', scriptPath])
			G_bter = graphio.readMat(tempFileOut, key='G_bter')
			matlabObject = scipy.io.loadmat(tempFileOut)
			self.t_fit = matlabObject["tFit"][0][0]
			self.t_generate = matlabObject["tGenerate"][0][0]
			return G_bter

	@classmethod
	def fit(cls, G, scale=1):
		return cls(G, scale)
