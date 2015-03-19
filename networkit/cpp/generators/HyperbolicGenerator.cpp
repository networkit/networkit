/*
 * HyperbolicGenerator.cpp
 *
 *  Created on: 20.05.2014
 *      Author: Moritz v. Looz (moritz.looz-corswarem@kit.edu)
 * This generator is a simplified version of the model presented in
 * "Hyperbolic geometry of complex networks. Physical Review E, 82:036106, Sep 2010." by
 *   Dmitri Krioukov, Fragkiskos Papadopoulos, Maksim Kitsak, Amin Vahdat, and Marian Boguna
 *
 */

#include <cstdlib>
#include <random>
#include <math.h>
#include <assert.h>
#include <omp.h>

#include "../graph/GraphBuilder.h"
#include "HyperbolicGenerator.h"
#include "Quadtree/Quadtree.h"
#include "../auxiliary/Random.h"
#include "../auxiliary/ProgressMeter.h"


namespace NetworKit {

HyperbolicGenerator::HyperbolicGenerator() {
	stretch = 1;
	alpha = 1;
	factor = 1;
	nodeCount = 10000;
	initialize();
}

HyperbolicGenerator::HyperbolicGenerator(count n) {
	nodeCount = n;
	alpha = 1;
	factor = 1;
	stretch=1;
	temperature=0;
	initialize();
}

/**
 * Construct a generator for n nodes and m edges
 */
HyperbolicGenerator::HyperbolicGenerator(count n, double avgDegree, double plexp, double T) {
	nodeCount = n;
	if (plexp < 2) throw std::runtime_error("Exponent of power-law degree distribution must be >= 2");
	if (T < 0 || T == 1) throw std::runtime_error("Temperature must be non-negative and not 1.");//Really necessary? Graphs with T=1 can be generated, only their degree is not controllable
	if (avgDegree > n) throw std::runtime_error("Average Degree must be at most n");
	alpha = (plexp-1)/2;
	double R = HyperbolicSpace::hyperbolicAreaToRadius(n);
	double targetR = HyperbolicSpace::getTargetRadius(n, n*avgDegree/2, alpha, T);// 2*log(8*n / (M_PI*(m/n)*2));
	stretch = targetR / R;
	factor = 1;
	temperature=T;
	initialize();
}

void HyperbolicGenerator::initialize() {
	if (temperature == 0) {
		capacity = 1000;
	} else {
		capacity = 10;
	}
	theoreticalSplit = false;
	threadtimers.resize(omp_get_max_threads());
	balance = 0.5;
}

Graph HyperbolicGenerator::generate() {
	return generate(nodeCount, factor, alpha, stretch, temperature);
}

Graph HyperbolicGenerator::generate(count n, double distanceFactor, double alpha, double stretchradius, double T) {
	double R = stretchradius*HyperbolicSpace::hyperbolicAreaToRadius(n);
	vector<double> angles(n);
	vector<double> radii(n);
	double r = HyperbolicSpace::hyperbolicRadiusToEuclidean(R);
	//sample points randomly

	HyperbolicSpace::fillPoints(angles, radii, stretchradius, alpha);
	vector<index> permutation(n);

	index p = 0;
	std::generate(permutation.begin(), permutation.end(), [&p](){return p++;});

	//can probably be parallelized easily, but doesn't bring much benefit
	std::sort(permutation.begin(), permutation.end(), [&angles,&radii](index i, index j){return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);});

	vector<double> anglecopy(n);
	vector<double> radiicopy(n);

	#pragma omp parallel for
	for (index j = 0; j < n; j++) {
		anglecopy[j] = angles[permutation[j]];
		radiicopy[j] = radii[permutation[j]];
	}

	INFO("Generated Points");
	return generate(anglecopy, radiicopy, r, R*distanceFactor, T);
}

double HyperbolicGenerator::expectedNumberOfEdges(count n, double stretch) {
	double R = stretch*HyperbolicSpace::hyperbolicAreaToRadius(n);
	return (8 / M_PI) * n * exp(-R/2)*(n/2);
}

Graph HyperbolicGenerator::generate(const vector<double> &angles, const vector<double> &radii, double r, double thresholdDistance, double T) {
	Aux::Timer timer;
	timer.start();
	index n = angles.size();
	assert(radii.size() == n);
	Quadtree<index> quad(r, theoreticalSplit, alpha, capacity, balance);

	//initialize a graph builder for n nodes and an undirected, unweighted graph with direct swap
	for (index i = 0; i < n; i++) {
		assert(radii[i] < r);
		quad.addContent(i, angles[i], radii[i]);
	}

	quad.trim();
	timer.stop();
	INFO("Filled Quadtree, took ", timer.elapsedMilliseconds(), " milliseconds.");

	return generate(angles, radii, quad, thresholdDistance, T);
}

Graph HyperbolicGenerator::generateCold(const vector<double> &angles, const vector<double> &radii, Quadtree<index> &quad, double thresholdDistance) {
	index n = angles.size();
	assert(radii.size() == n);
	Aux::Timer timer;
	timer.start();
	vector<double> empty;
	GraphBuilder result(n, false, false, directSwap);

	Aux::ProgressMeter progress(n, 10000);
	#pragma omp parallel
	{
		index id = omp_get_thread_num();
		threadtimers[id].start();
		#pragma omp for schedule(guided) nowait
		for (index i = 0; i < n; i++) {
			//get neighbours for node i
			count expectedDegree = (4/M_PI)*n*exp(-HyperbolicSpace::EuclideanRadiusToHyperbolic(radii[i])/2);
			vector<index> near;
			near.reserve(expectedDegree*1.1);
			quad.getElementsInHyperbolicCircle(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), thresholdDistance, near);
			//count realDegree = near.size();
			//std::swap(expectedDegree, realDegree);//dummy statement for debugging
			if (directSwap) {
				std::remove(near.begin(), near.end(), i); //no self loops!
				near.pop_back();//std::remove doesn't remove element but swaps it to the end
				result.swapNeighborhood(i, near, empty, false);
			} else {
				for (index j : near) {
					if (j < i) result.addEdge(i,j);
				}
			}

			if (i % 10000 == 0) {
				#pragma omp critical (progress)
				{
					progress.signal(i);
				}
			}
		}
		threadtimers[id].stop();
	}

	timer.stop();
	INFO("Generating Edges took ", timer.elapsedMilliseconds(), " milliseconds.");
	return result.toGraph(true);
}

Graph HyperbolicGenerator::generate(const vector<double> &angles, const vector<double> &radii, Quadtree<index> &quad, double thresholdDistance, double T) {
	if (T < 0) throw std::runtime_error("Temperature cannot be negative.");
	if (T == 0) return generateCold(angles, radii, quad, thresholdDistance);
	assert(T > 0);
	count n = angles.size();
	assert(radii.size() == n);
	assert(quad.size() == n);

	//now define lambda
	double beta = 1/T;
	auto edgeProb = [beta, thresholdDistance](double distance) -> double {return 1 / (exp(beta*(distance-thresholdDistance)/2)+1);};

	//get Graph
	GraphBuilder result(n, false, false, false);//no direct swap with probabilistic graphs, sorry
	Aux::ProgressMeter progress(n, 10000);
	count totalCandidates = 0;
	#pragma omp parallel for
	for (index i = 0; i < n; i++) {
		vector<index> near;
		totalCandidates += quad.getElementsProbabilistically(HyperbolicSpace::polarToCartesian(angles[i], radii[i]), edgeProb, near);
		for (index j : near) {
			if (j < i) result.addEdge(i, j);
		}

		if (i % 10000 == 0) {
			#pragma omp critical (progress)
			{
				progress.signal(i);
			}
		}

	}
	DEBUG("Candidates tested: ", totalCandidates);
	return result.toGraph(true);

}
}
