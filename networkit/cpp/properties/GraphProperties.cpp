/*
 */

#include "GraphProperties.h"
#include "../auxiliary/Log.h"


namespace NetworKit {

std::vector<count> GraphProperties::degreeDistribution(const Graph& G) {
	std::vector<count> distribution;
	if (G.isDirected()) {
		auto maxDegree = minMaxDegreeDirected(G).second;
		distribution.assign(maxDegree.first+maxDegree.second+1, 0);
		G.parallelForNodes([&](node v){
			count i = G.degreeOut(v) + G.degreeIn(v);
			#pragma omp atomic
			distribution[i]++;
		});
		// workaround as now minmaxdegree for combined degree is implemented yet.
		count i = distribution.size()-1;
		while (i > 0 && distribution[i] == 0) --i;
		distribution.resize(i+1);
	} else {
		count maxDegree = minMaxDegree(G).second;
		distribution.assign(maxDegree+1, 0);
		G.parallelForNodes([&](node v){
			count i = G.degree(v);
			#pragma omp atomic
			distribution[i]++;
		});
	}
	return distribution;
}


std::vector<double> GraphProperties::localClusteringCoefficients(const Graph& G) {
	count n = G.upperNodeIdBound();
	std::vector<double> numerator(n); //
	std::vector<double> denominator(n); // $\deg(u) \cdot ( \deg(u) - 1 )$
	std::vector<double> coefficient(n); // $c(u) := \frac{2 \cdot |E(N(u))| }{\deg(u) \cdot ( \deg(u) - 1)}$

	G.balancedParallelForNodes([&](node u){
		count edgeCount = 0;
		G.forEdgesOf(u, [&](node u, node v) {
			G.forEdgesOf(v, [&](node v, node w){
				if (G.hasEdge(u, w)) {
					edgeCount += 1;
				}
			});
		});

		numerator[u] = edgeCount; // factor 2 is omitted because each edge has been counted twice
	});

	G.balancedParallelForNodes([&](node u){
		denominator[u] = G.degree(u) * (G.degree(u) - 1);
	});

	G.balancedParallelForNodes([&](node u){
		if (denominator[u] == 0.0) {
			coefficient[u] = 0.0;
		} else {
			coefficient[u] = numerator[u] / denominator[u];
		}
	});

	return coefficient;
}

std::vector<double> GraphProperties::localClusteringCoefficientPerDegree(const Graph& G) {

	std::vector<count> degDist = degreeDistribution(G);
	std::vector<double> coefficient;
	std::vector<double> perDegree(degDist.size(), 0.0);

	if (G.numberOfNodes() > 1 ) {
		coefficient = localClusteringCoefficients(G);

		G.forNodes([&](node u){
			perDegree[G.degree(u)] += coefficient[u];
		});

		// get the average local clustering coefficient for nodes of each degreee
		for (index i = 2; i < degDist.size(); ++i) {
			if (degDist[i] == 0) {
				perDegree[i] = 0.0; // TODO: should this be -1
			} else {
				perDegree[i] = perDegree[i] / (double) degDist[i];
			}
		}
	}



	// allows to avoid the situation, when local clustering coefficient is calculated for 0-1 degree nodes.
	// These nodes are warranted not to be triangle centers, thus we avoid calculating the coefficients for the,
	degDist[0] = 0;
	if (G.numberOfNodes() > 0 ) degDist[1] = 0;

	return perDegree;
}

double GraphProperties::averageLocalClusteringCoefficient(const Graph& G) {
	std::vector<double> coefficients = GraphProperties::localClusteringCoefficients(G);
	double sum = 0.0;
	for (double c : coefficients) {
		sum += c;
	}
	double avg = sum / G.numberOfNodes();
	return avg;
}

std::pair<count, count> GraphProperties::minMaxDegree(const Graph& G) {
	assert(!G.isDirected());
	count min = G.numberOfNodes();
	count max = 0;

	G.forNodes([&](node v){
		count d = G.degree(v);
		if (d < min) {
			min = d;
		}
		if (d > max) {
			max = d;
		}
	});

	return std::pair<count, count>(min, max);
}

std::pair<std::pair<count,count>, std::pair<count,count>> GraphProperties::minMaxDegreeDirected(const Graph& G) {
	assert(G.isDirected());
	count minIn = G.numberOfNodes();
	count minOut = G.numberOfNodes();
	count maxIn = 0;
	count maxOut = 0;

	G.forNodes([&](node v){
		count d = G.degreeIn(v);
		if (d < minIn) {
			minIn = d;
		}
		if (d > maxIn) {
			maxIn = d;
		}
		d = G.degreeOut(v);
		if (d < minOut) {
			minOut = d;
		}
		if (d > maxIn) {
			maxOut = d;
		}
	});

	return {{minIn,minOut},{maxIn,maxOut}};
}

std::vector< count > GraphProperties::degreeSequence(const NetworKit::Graph &G) {
	std::vector<count> sequence(G.upperNodeIdBound());

	G.parallelForNodes([&](node v) {
		sequence[v] = G.degree(v);
	});

	return sequence;
}

double GraphProperties::averageDegree(const Graph& G) {

	count n = G.numberOfNodes();

	count degSum = G.parallelSumForNodes([&](node v){
		return G.degree(v);
	});

	double avgDeg = degSum / (double) n;
	return avgDeg;
}

double GraphProperties::degreeAssortativitySlower(const Graph& G, bool useWeights) {
	// note: a parallel implementation would rather follow Newman's book, p. 267
	// however, parallelism introduces inaccuracies due to numerical issues in reduction

	double r = 0.0; // result
	double A = 0.0; // accumulates degree products
	double B = 0.0; // accumulates degree sums
	double C = 0.0; // accumulates sum of degree squares

	double degu = 0.0; // temp storage for degree(u)
	double degv = 0.0; // temp storage for degree(v)

	double halfVolume = 0.0; // if needed, halfVolume accumulates the total edge weight of the graph (such a routine exists, but is not called for performance reasons)

	// iterate over edges and accumulate
	if (G.isWeighted() && useWeights) {
		G.forEdges([&](node u, node v, edgeweight ew) {
			degu = G.weightedDegree(u);
			degv = G.weightedDegree(v);
			A += degu * degv;
			B += degu + degv;
			C += degu*degu + degv*degv;
			halfVolume += ew;
		});
	}
	else {
		G.forEdges([&](node u, node v) {
			degu = G.degree(u);
			degv = G.degree(v);
			A += degu * degv;
			B += degu + degv;
			C += degu*degu + degv*degv;
		});

		halfVolume = G.numberOfEdges();
	}

	double volume = 2.0 * halfVolume;
	A = A / halfVolume;
	B = B / volume;
	B = B*B;
	C = C / volume;

	TRACE("A: ", A, ", B: ", B, ", C: ", C);

	assert(C != B);
	r = (A - B) / (C - B);
	return r;
}



double GraphProperties::degreeAssortativity(const Graph& G, bool useWeighted) {
	double r = 0.0; // result

	double S1 = 0.0; // accumulates degrees
	double S2 = 0.0; // accumulates squared degrees
	double S3 = 0.0; // accumulates cubed degress
	double Se = 0.0; // accumulates degree products

	double deg = 0.0; // temp storage for degree
	double sqr = 0.0;  // temp storage square of degree

	// iterate over edges and accumulate
	if (G.isWeighted() && useWeighted) {
		return degreeAssortativitySlower(G, useWeighted);
	}
	else {
		TRACE("start sum");
		G.forEdges([&](node u, node v) { // parallel iteration leads to crash
			Se += G.degree(u) * G.degree(v);
		});
		TRACE("end sum");

		G.forNodes([&](node u) {
			deg = G.degree(u);
			S1 += deg;
			sqr = deg * deg;
			S2 += sqr;
			S3 += sqr * deg;
		});
	}
	Se = 2.0 * Se;

	assert(S1 * S3 != S2 * S2);
	r = (S1 * Se - S2 * S2) / (S1 * S3 - S2 * S2);
	return r;
}

double GraphProperties::degreeAssortativityDirected(const Graph& G, bool alpha, bool beta) {
	assert(G.isDirected());
	count alphaAvg = 0;
	count betaAvg = 0;
	G.forEdges([&](node u, node v){
		alphaAvg += (alpha) ? G.degreeOut(u) : G.degreeIn(u);
		betaAvg += (beta) ? G.degreeOut(v) : G.degreeIn(v);
	});
	alphaAvg /= G.numberOfEdges();
	betaAvg /= G.numberOfEdges();
	double normalize = 1.f/G.numberOfEdges();
	count sum = 0;
	count jAggr = 0;
	count kAggr = 0;
	G.forEdges([&](node u, node v){
		count j = ((alpha) ? G.degreeOut(u) : G.degreeIn(u)) - alphaAvg;
		count k = ((beta) ? G.degreeOut(v) : G.degreeIn(v)) - betaAvg;
		sum += (j * k);
		jAggr += (j*j);
		kAggr += (k*k);
	});
	// multiplication with normalize doesn't change anything; could be avoided?!
	return (normalize*sum)/(std::sqrt(normalize*jAggr)*std::sqrt(normalize*kAggr));

}



} /* namespace NetworKit */
