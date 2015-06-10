/*
 * LAMGBenchmark.cpp
 *  Adaptation of LaplacianSolverBenchmark by Daniel Hoske
 *  Created on: May 27, 2015
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "LAMGBenchmark.h"
#include "BenchGraphs.h"

#include <iostream>
#include <fstream>
#include <future>
#include <atomic>
#include <algorithm>

#include "../../properties/ConnectedComponents.h"
#include "../../auxiliary/Timer.h"
#include "../../auxiliary/StringTools.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"
#include "../../io/MatrixMarketReader.h"
#include "../GaussSeidelRelaxation.h"

using namespace std;

namespace NetworKit {

enum SmootherType {
	GAUSS_SEIDEL
};

/** Maximum solve time */
constexpr int MAX_CONVERGENCE_TIME = 10 * 60 * 1000;

/** Floating point epsilon to use in comparisons. */
constexpr double EPSILON = 1e-9;

struct Benchmark {
	string name;

	// LAMG settings
	count setupTries;
	count solveTriesPerSetup;
	double residual;
	SmootherType smootherType;

	vector<string> graphs;
};

const vector<Benchmark> BENCHS = { // available benchmarks for different graph types
	{
		"Grids",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		GRIDS
	},
	{
		"Barabasi",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		BARABASI
	},
	{
		"Laplace",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		LAPLACE
	},
	{
		"Dimacs Numerics",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		DIMACS_NUMERICS
	},
	{
		"Dimacs Sparse",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		DIMACS_SPARSE
	},
	{
		"Dimacs Clustering",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		DIMACS_CLUSTERING
	},
	{
		"Dimacs Street",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		DIMACS_STREET_NETWORKS
	},
	{
		"Walshaw",
		10,
		10,
		1e-6,
		GAUSS_SEIDEL,
		WALSHAW
	}
};

Vector randZeroSum(const Graph& G, size_t seed) {
	mt19937 rand(seed);
	auto rand_value = uniform_real_distribution<double>(-1.0, 1.0);
	ConnectedComponents con(G);
	count n = G.numberOfNodes();
	con.run();
	Partition comps = con.getPartition();

	/* Fill each component randomly such that its sum is 0 */
	Vector b(n, 0.0);

	for (int id : comps.getSubsetIds()) {
		auto indexes = comps.getMembers(id);
		assert(!indexes.empty());
		double sum = 0.0;
		for (auto entry : indexes) {
			b[entry] = rand_value(rand);
			sum += b[entry];
		}
		b[*indexes.begin()] -= sum;
	}

	return b;
}

Vector randVector(count dimension) {
	Vector randVector(dimension);
	for (index i = 0; i < dimension; ++i) {
		randVector[i] = 2.0 * Aux::Random::probability() - 1.0;
	}

	// introduce bias
	for (index i = 0; i < dimension; ++i) {
		randVector[i] = randVector[i] * randVector[i];
	}

	return randVector;
}

bool isSymmetric(const CSRMatrix& A) {
	bool output = true;
	A.forNonZeroElementsInRowOrder([&] (index i, index j, edgeweight w) {
		if (A(j, i) != w) {
			output = false;
		}
	});
	return output;
}

bool isLaplacian(const CSRMatrix& A) {
	if (!isSymmetric(A)) {
		return false;
	}

	/* Criterion: \forall_i \sum_j A_ij = 0  */
	vector<double> row_sum(A.numberOfRows());
	atomic<bool> right_sign(true);
	A.parallelForNonZeroElementsInRowOrder([&] (node i, node j, double value) {
		if (i != j && value > EPSILON) {
			right_sign = false;
		}
		row_sum[i] += value;
	});

	return right_sign && all_of(row_sum.begin(), row_sum.end(), [] (double val) {return abs(val) < EPSILON;});
}

Graph laplacianToGraph(const CSRMatrix &matrix) {
	Graph G(max(matrix.numberOfRows(), matrix.numberOfColumns()), true, true);
	matrix.forNonZeroElementsInRowOrder([&](node u, node v, double val) {
		if (v != u) {
			G.addEdge(u, v, -val);
		}
	});

	return G;
}

Graph matrixToGraph(const CSRMatrix &matrix) {
	Graph G(max(matrix.numberOfRows(), matrix.numberOfColumns()), true, true);
	matrix.forNonZeroElementsInRowOrder([&](node u, node v, double val) {
		G.addEdge(u, v, val);
	});

	return G;
}

Graph readGraph(const string& path) {
	using namespace Aux::StringTools;
	if (ends_with(path, ".graph")) {
		METISGraphReader reader;
		Graph G = reader.read(path);
		CSRMatrix M = CSRMatrix::adjacencyMatrix(G);
		if (isLaplacian(M)) {
			INFO("already laplacian!!");
		}
		return reader.read(path);
	} else if (ends_with(path, ".mtx")) {
		INFO("reading ", path);
		MatrixMarketReader reader;
		CSRMatrix M = reader.read(path);
		if (isLaplacian(M)) {
			return laplacianToGraph(M);
		} else {
			return matrixToGraph(M);
		}
	} else {
		throw std::runtime_error("unknown graph format " + path);
	}
}

void writeVector(const Vector &vector, const string &filename) {
	ofstream outStream;
	outStream.precision(16);
	outStream.open(filename);
	for (index i = 0; i < vector.getDimension() - 1; ++i) {
		outStream << vector[i] << ", ";
	}
	outStream << vector[vector.getDimension()-1];
	outStream.close();
}

void writeProblemVectors(const Vector &b, const Vector &x, const string graphFilepath) {
	string filepath = graphFilepath.substr(0, graphFilepath.find_last_of("."));
	stringstream ss;
	ss << filepath << "_x";
	writeVector(x, ss.str());

	ss.str("");
	ss << filepath << "_b";
	writeVector(b, ss.str());
}

void writeBenchmarkResults(string texContent, string filename) {
	ofstream output;
	output.open(filename);
	output << texContent;
	output.close();
}

template <typename T>
string numprint(T val) {
	stringstream ss;
	ss << setprecision(16);
	ss << "\\numprint{" << val << "}";
	return ss.str();
}

string printLatexDocumentHeader() {
	stringstream header;
	header << "\\documentclass[11pt,a4paper]{article} \n" <<
					"\\usepackage{tabularx} \n" <<
					"\\usepackage{numprint} \n" <<
					"\\usepackage[margin=1in]{geometry} \n\n" <<
					"\\begin{document} \n\n" <<
					"\\\nprounddigits{2}";

	time_t time_cur = std::time(nullptr);
	tm* time_local = localtime(&time_cur); // not reentrant
	char time_str[80];
	strftime(&time_str[0], sizeof(time_str), "%Y-%m-%d %H:%M:%S", time_local);

	header << "\\section*{Benchmark " << Aux::toStringF("%s", time_str) << "}";

	return header.str();
}

string printLatexDocumentFooter() {
	string footer = "\\end{document}";
	return footer;
}


string printBenchHeader(const Benchmark &bench) {
	stringstream ss;
	ss << "\\begin{table}[h] \n";
	ss << "\\begin{tabularx}{\\textwidth}{l r r r r r r} \n";
	ss << "Graph & size & nonzeros & setup & solve & iterations & residual \\\\ \\hline \n";

	return ss.str();
}

string printBenchFooter(const Benchmark &bench) {
	stringstream ss;
	ss << "\\end{tabularx} \n";
	ss << "\\caption{" << bench.graphs[0].substr(0, bench.graphs[0].find_last_of("/")) << "} \n";
	ss << "\\end{table} \n";

	return ss.str();
}

string printTableRow(const vector<LAMGSolverStatus> &solverStati, double avgSetupTime, double avgSolveTime) {
	stringstream ss;
	ss.precision(16);

	count numConverged = 0;
	double numIters = 0;
	double finalResidual = 0.0;
	for (LAMGSolverStatus status : solverStati) {
		numConverged += status.converged;
		numIters += status.numIters;
		finalResidual += status.residual;
	}

	ss << numprint(avgSetupTime) << " & " << numprint(avgSolveTime) << " & " << numprint(numIters / (double) solverStati.size()) << " & " << numprint(finalResidual / (double) solverStati.size()) << " \\\\ \n";

	return ss.str();
}

void outputPlotData(const LAMGSolverStatus &status, const string &filename) {
	ofstream output;
	output.open(filename);
	output << setprecision(16);
	std::vector<double> residualHistory = status.residualHistory;
	for (int i = 0; i < residualHistory.size(); ++i) {
		output << i << "\t" << residualHistory[i] << endl;
	}

	output.close();
}


string benchmark(const CSRMatrix &matrix, const Vector &initialX, const Vector &b, const string &graphPath, count numSetups, count numSolvesPerSetup, const Smoother &smoother) {
	Aux::Timer tSetup;
	Aux::Timer tSolve;

	count setupTime = 0;
	count solveTime = 0;
	vector<LAMGSolverStatus> solverStati(numSetups * numSolvesPerSetup);

	for (index i = 0; i < numSetups; ++i) {
		MultiLevelSetup setup(smoother);
		LevelHierarchy hierarchy;

		tSetup.start();
			setup.setup(matrix, hierarchy);
			SolverLamg solver(hierarchy, smoother);
		tSetup.stop();
		setupTime += tSetup.elapsedMilliseconds();

		for (index j = 0; j < numSolvesPerSetup; ++j) {
			Vector x = initialX;
			solverStati[i * numSolvesPerSetup + j].maxConvergenceTime = MAX_CONVERGENCE_TIME;

			tSolve.start();
				solver.solve(x, b, solverStati[i * numSolvesPerSetup + j]);
			tSolve.stop();
			solveTime += tSolve.elapsedMilliseconds();
		}
	}

	double avgSetupTime = (double) setupTime / (double) numSetups;
	double avgSolveTime = (double) solveTime / (double) (numSetups*numSolvesPerSetup);

	stringstream ss;
	ss << graphPath.substr(0, graphPath.find_last_of(".")) << ".plot";
	outputPlotData(solverStati[0], ss.str());

	return printTableRow(solverStati, avgSetupTime, avgSolveTime);
}

string benchmark(const Benchmark &bench) {
	/* Read all given graphs in parallel */
	vector<Graph> graphs;
	vector<future<Graph>> graph_reads;
	for (const auto& filename : bench.graphs) {
		graph_reads.emplace_back(async(launch::async, [&] {
			Graph G = readGraph(filename);
			G.setName(filename);
			return G;
		}));
	}

	for (auto& read : graph_reads) {
		graphs.emplace_back(read.get());
	}
	sort(graphs.begin(), graphs.end(), [&] (const Graph& G1, const Graph& G2) {
		return G1.numberOfNodes() < G2.numberOfNodes();
	});

	/* Convert formats */
	vector<CSRMatrix> L;

	vector<Vector> b;
	vector<Vector> x;
	for (index i = 0; i < bench.graphs.size(); ++i) {
		auto& G = graphs[i];
		CSRMatrix LCur = CSRMatrix::graphLaplacian(G);
		L.emplace_back(LCur);

		Vector bCur = randZeroSum(G, 12345);
		b.emplace_back(bCur);

		Vector xCur = randVector(G.numberOfNodes());
		x.emplace_back(xCur);

		writeProblemVectors(bCur, xCur, graphs[i].getName());
	}

	string output = "";
	output += printBenchHeader(bench);

	Smoother *smoother;
	switch(bench.smootherType) {
		case GAUSS_SEIDEL:
			smoother = new GaussSeidelRelaxation();
			break;
		default:
			smoother = new GaussSeidelRelaxation();
			break;
	}

	stringstream ss;
	for (index i = 0; i < graphs.size(); ++i) {
		string graphFile = graphs[i].getName().substr(bench.graphs[i].find_last_of("/") + 1, bench.graphs[i].length());
		output += graphFile.substr(0, graphFile.find_last_of("."));
		ss << " & " << L[i].numberOfRows() << " x " << L[i].numberOfColumns() << " & ";
		ss << L[i].nnz() << " & ";
		output += ss.str();

		INFO(graphFile);

		ConnectedComponents con(graphs[i]);
		con.run();
		if (con.numberOfComponents() == 1) { // LAMG solver currently only supports connected graphs
			output += benchmark(L[i], x[i], b[i], graphs[i].getName(), bench.setupTries, bench.solveTriesPerSetup, *smoother);
		} else {
			output += " -  & - & - & - \\\\ \n";
		}
		ss.str("");
	}

	output += printBenchFooter(bench);
	delete smoother;

	return output;
}

TEST_F(LAMGBenchmark, bench) {
	string texContent = printLatexDocumentHeader();
	stringstream ss;

	Benchmark bench = BENCHS[0]; // grids
	texContent += benchmark(bench);
	ss << bench.name;

	texContent += printLatexDocumentFooter();

	ss << ".tex";
	string filename = "benchmark_" + ss.str();

	writeBenchmarkResults(texContent, filename);
}

} /* namespace NetworKit */
