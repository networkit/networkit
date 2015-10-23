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
#include <iomanip>

#include "../../components/ConnectedComponents.h"
#include "../../auxiliary/Timer.h"
#include "../../auxiliary/StringTools.h"
#include "../../graph/Graph.h"
#include "../../io/METISGraphReader.h"
#include "../../io/MatrixMarketReader.h"
#include "../GaussSeidelRelaxation.h"
#include "../LAMG/Lamg.h"

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

	vector<Instance> instances;
};

const vector<Benchmark> BENCHS = { // available benchmarks for different graph types
	{ // 0
		"Grids",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		GRIDS
	},
	{ // 1
		"Barabasi",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		BARABASI
	},
	{ // 2
		"Laplace",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		LAPLACE
	},
	{ // 3
		"Dimacs Numerics",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		DIMACS_NUMERICS
	},
	{ // 4
		"Dimacs Sparse",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		DIMACS_SPARSE
	},
	{ // 5
		"Dimacs Clustering",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		DIMACS_CLUSTERING
	},
	{ // 6
		"Dimacs Street",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		DIMACS_STREET_NETWORKS
	},
	{ // 7
		"Facebook100",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		FACEBOOK100
	},
	{ // 8
		"Walshaw",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		WALSHAW
	},
	{ // 9
		"Snap",
		10,
		5,
		1e-6,
		GAUSS_SEIDEL,
		SNAP
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

Graph sddToLaplacian(const CSRMatrix& A) {
  assert(CSRMatrix::isSDD(A));
  count n = A.numberOfColumns();

  /* Compute row sum and excess */
  vector<double> abs_row_sum(n);
  A.parallelForNonZeroElementsInRowOrder([&] (node i, node j, double value) {
    if (i != j) {
      abs_row_sum[i] += abs(value);
    }
  });

  vector<double> excess(n);
#pragma omp parallel for
  for (index i = 0; i < n; ++i) {
     excess[i] = A(i, i) - abs_row_sum[i];
  }

  /* Set lower and upper diagonal elements */
  Graph G(2*n, true);
  for (index i = 0; i < n; ++i) {
    G.addEdge(i, i + n, excess[i]/2);
  }

  /* Distribute original entries of A into the 2 times larger graph G */
  A.forNonZeroElementsInRowOrder([&] (node i, node j, double value) {
    if (i < j) {
      /* Off-diagonals */
      if (value < 0) {
        G.addEdge(i,     j,     -value);
        G.addEdge(i + n, j + n, -value);
      } else {
        G.addEdge(i + n,     j, value);
        G.addEdge(    i, j + n, value);
      }
    }
  });

  return G;
}

Graph readGraph(Instance& instance) {
	using namespace Aux::StringTools;
	if (ends_with(instance.path, ".graph") || ends_with(instance.path, ".metis")) {
		METISGraphReader reader;
		Graph G = reader.read(instance.path);
		if (instance.type == MATRIX) {
			CSRMatrix A = CSRMatrix::adjacencyMatrix(G);
			if (CSRMatrix::isLaplacian(A)) {
				instance.type = LAPLACIAN;
			} else if (CSRMatrix::isSDD(A)) {
				instance.type = SDD;
			} else {
				INFO(instance.path, " is unknown");
				instance.type = UNKNOWN;
			}
		}

		return G;
	} else if (ends_with(instance.path, ".mtx")) {
		MatrixMarketReader reader;
		INFO("reading ", instance.path);
		CSRMatrix M = reader.read(instance.path);
		INFO("done");
		if (CSRMatrix::isLaplacian(M)) {
			instance.type = LAPLACIAN;
		} else if (CSRMatrix::isSDD(M)) {
			instance.type = SDD;
		} else {
			instance.type = UNKNOWN;
		}

		return CSRMatrix::matrixToGraph(M);
	} else {
		throw std::runtime_error("unknown graph format " + instance.path);
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
string numprint(T val, bool scientific = false) {
	stringstream ss;
	ss << setprecision(16);
	if (scientific)	ss.setf(std::ios_base::scientific, std::ios_base::floatfield);
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
					"\\nprounddigits{2}\n";

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
	ss << "\\caption{" << bench.instances[0].path.substr(0, bench.instances[0].path.find_last_of("/")) << "} \n";
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

	ss << numprint(avgSetupTime) << " & " << numprint(avgSolveTime) << " & " << numprint(numIters / (double) solverStati.size()) << " & " << numprint(finalResidual / (double) solverStati.size(), true) << " \\\\ \n";

	return ss.str();
}

string printTableRow(const SolverStatus &status, double avgSetupTime, double avgSolveTime) {
	stringstream ss;
	ss.precision(16);

	double numIters = status.numIters;
	double finalResidual = status.residual;

	ss << numprint(avgSetupTime) << " & " << numprint(avgSolveTime) << " & " << numprint(numIters) << " & " << numprint(finalResidual, true) << " \\\\ \n";

	return ss.str();
}

void outputPlotData(const LAMGSolverStatus &status, const string &filename) {
	ofstream output;
	output.open(filename);
	output << setprecision(16);
	std::vector<double> residualHistory = status.residualHistory;
	for (index i = 0; i < residualHistory.size(); ++i) {
		output << i << "\t" << residualHistory[i] << endl;
	}

	output.close();
}


string benchmark(const CSRMatrix &matrix, const Vector &initialX, const Vector &b, const string &graphPath, count numSetups, count numSolvesPerSetup, const Smoother &smoother, double desiredResidual) {
	Aux::Timer tSetup;
	Aux::Timer tSolve;

	count setupTime = 0;
	count solveTime = 0;

	Lamg lamg(desiredResidual);
	SolverStatus status;
	for (index i = 0; i < numSetups; ++i) {
		tSetup.start();
			lamg.setup(matrix);
		tSetup.stop();
		setupTime += tSetup.elapsedMilliseconds();

		for (index j = 0; j < numSolvesPerSetup; ++j) {
			Vector x = initialX;
			tSolve.start();
				status = lamg.solve(b, x, MAX_CONVERGENCE_TIME);
			tSolve.stop();
			solveTime += tSolve.elapsedMilliseconds();
		}
	}

	double avgSetupTime = (double) setupTime / (double) numSetups;
	double avgSolveTime = (double) solveTime / (double) (numSetups*numSolvesPerSetup);

	stringstream ss;
	ss << graphPath.substr(0, graphPath.find_last_of(".")) << ".plot";
	//outputPlotData(solverStati[0], ss.str());

	return printTableRow(status, avgSetupTime, avgSolveTime);
}

string benchmark(Benchmark &bench) {
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
	for (index i = 0; i < bench.instances.size(); ++i) {
		Graph G = readGraph(bench.instances[i]);
		CSRMatrix L;
		switch(bench.instances[i].type) {
			case LAPLACIAN:
				L = CSRMatrix::adjacencyMatrix(G);
				break;
			case SDD:
				L = CSRMatrix::adjacencyMatrix(sddToLaplacian(CSRMatrix::adjacencyMatrix(G)));
				break;
			case NETWORK:
				L = CSRMatrix::graphLaplacian(G);
				break;
			default:
				INFO(bench.instances[i].path, " is unknown matrix. Will not solve this!");
				continue;
		}

		string path = bench.instances[i].path;
		G.setName(path);
		Vector b = randZeroSum(CSRMatrix::matrixToGraph(L), 12345);
		Vector x = randVector(L.numberOfColumns());

		INFO("L = ", L.numberOfRows(), "x", L.numberOfColumns(), "; x = ", x.getDimension(), "; b = ", b.getDimension());

		writeProblemVectors(b, x, G.getName());

		string graphFile = G.getName().substr(path.find_last_of("/") + 1, path.length());
		output += graphFile.substr(0, graphFile.find_last_of("."));
		ss << " & " << L.numberOfRows() << " x " << L.numberOfColumns() << " & ";
		ss << L.nnz() << " & ";
		output += ss.str();

		INFO(graphFile);

//		Graph Gcheck = matrixToGraph(L);
//		ConnectedComponents con(Gcheck);
//		con.run();
//		if (con.numberOfComponents() == 1) { // LAMG solver currently only supports connected graphs
			output += benchmark(L, x, b, G.getName(), bench.setupTries, bench.solveTriesPerSetup, *smoother, bench.residual);
//			INFO(output);
//		} else {
//			INFO("not connected");
//			Partition p = con.getPartition();
//			output += benchmarkDisconnected(G, p, x, b, G.getName(), bench.setupTries, bench.solveTriesPerSetup, *smoother, bench.residual);
//
//			//output += " -  & - & - & - \\\\ \n";
//		}
		ss.str("");
	}

	output += printBenchFooter(bench);
	delete smoother;

	return output;
}

TEST_F(LAMGBenchmark, bench) {
	string texContent = printLatexDocumentHeader();
	stringstream ss;

	Benchmark bench = BENCHS[8]; // snap
	texContent += benchmark(bench);
	ss << bench.name;

	texContent += printLatexDocumentFooter();

	ss << ".tex";
	string filename = "benchmark_" + ss.str();

	writeBenchmarkResults(texContent, filename);
}

} /* namespace NetworKit */
