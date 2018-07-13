/*
 * MocnikGenerator.cpp
 *
 * Created on: July 7, 2018
 * Author: Franz-Benjamin Mocnik <mail@mocnik-science.net>
 */

#include "MocnikGenerator.h"
#include "../auxiliary/Random.h"
#include <cmath>
#include <unordered_map>

namespace NetworKit {

MocnikGenerator::MocnikGenerator(count dim, count n, double k, bool weighted): dim(dim), weighted(weighted) {
	ns.push_back(n);
	ks.push_back(k);
}

MocnikGenerator::MocnikGenerator(count dim, std::vector<count> ns, double k, bool weighted): dim(dim), ns(ns), weighted(weighted) {
	for (auto &n : ns) {
		ks.push_back(k);
	}
}

MocnikGenerator::MocnikGenerator(count dim, std::vector<count> ns, std::vector<double> ks, bool weighted): dim(dim), ns(ns), ks(ks), weighted(weighted) {
}

MocnikGenerator::MocnikGenerator(count dim, count n, double k, std::vector<double> weighted): dim(dim), weighted(true), relativeWeights(weighted) {
	ns.push_back(n);
	ks.push_back(k);
}

MocnikGenerator::MocnikGenerator(count dim, std::vector<count> ns, double k, std::vector<double> weighted): dim(dim), ns(ns), weighted(true), relativeWeights(weighted) {
	for (auto &n : ns) {
		ks.push_back(k);
	}
}

MocnikGenerator::MocnikGenerator(count dim, std::vector<count> ns, std::vector<double> ks, std::vector<double> weighted): dim(dim), ns(ns), ks(ks), weighted(true), relativeWeights(weighted) {
}

// TOOLS

template<typename KeyType, typename ValueType> std::pair<KeyType,ValueType> get_max(const std::map<KeyType,ValueType> &x) {
	using pairtype=std::pair<KeyType,ValueType>;
	return *std::max_element(x.begin(), x.end(), [](const pairtype &p1, const pairtype &p2) {
		return p1.second < p2.second;
	});
}

// GEOMETRY

// norm of a vector
static inline double norm(std::vector<double> &v, const double &shift) {
	double x = 0;
	for (count j = 0; j < v.size(); j++) {
		x += (v[j] + shift) * (v[j] + shift);
	}
	return sqrt(x);
}

// Euclidean distance between two vectors
static inline double dist(std::vector<double> &v, std::vector<double> &w) {
	double x = 0;
	for (count j = 0; j < v.size(); j++) {
		x += std::pow(v[j] - w[j], 2);
	}
	x = sqrt(x);
	return x;
}

// NODE ARRAY

void MocnikGenerator::initCellArray(MocnikGenerator::State &s, const count &m) {
	s.aMax = m;
	for (count j = 0; j < std::pow(s.aMax, dim); j++) {
		NodePositionMap tmp;
		s.a.push_back(tmp);
	}
}

MocnikGenerator::NodePositionMap MocnikGenerator::getNodes(MocnikGenerator::State &s, const int &i) {
		return s.a[i];
}

const void MocnikGenerator::addNode(MocnikGenerator::State &s, const std::pair<node, std::vector<double>> &p) {
	s.a[toIndex(s, p.second)].push_back(p);
}

const int MocnikGenerator::toIndex(MocnikGenerator::State &s, const std::vector<double> &v) {
	std::vector<int> w;
	for (count j = 0; j < v.size(); j++) {
		w.push_back(fmin(floor(v[j] * s.aMax), s.aMax - 1));
	}
	return toIndex(s, w);
}

const int MocnikGenerator::toIndex(MocnikGenerator::State &s, const std::vector<int> &v) {
	int x = 0;
	for (count j = v.size() - 1; j >= 0 && j < v.size(); j--) {
		x = x * s.aMax + v[j];
	}
	return x;
}

const std::vector<int> MocnikGenerator::fromIndex(MocnikGenerator::State &s, const int &i) {
	std::vector<int> v;
	int i2 = i;
	for (count j = 0; j < dim; j++) {
		int i2New = i2 % s.aMax;
		v.push_back(i2New);
		i2 = (i2 - i2New) / s.aMax;
	}
	return v;
}

const std::vector<int> MocnikGenerator::boxSurface(MocnikGenerator::State &s, const int &i, const int &r) {
	std::vector<std::vector<int>> se;
	if (r == 0) {
		se.push_back(fromIndex(s, i));
	} else {
		std::vector<int> iV = fromIndex(s, i);
		for (count d = 0; d < dim; d++) {
			std::vector<std::vector<int>> v;
			std::vector<int> tmp;
			v.push_back(tmp);
			for (count j = 0; j < d; j++) {
				std::vector<std::vector<int>> w;
				for (int mu = fmax(iV[j] - r + 1, 0); mu <= fmin(iV[j] + r - 1, s.aMax - 1); mu++) {
					for (std::vector<int> vElem : v) {
						std::vector<int> x(vElem);
						x.push_back(mu);
						w.push_back(x);
					}
				}
				v = w;
			}
			std::vector<std::vector<int>> w;
			for (std::vector<int> vElem : v) {
				if (iV[d] - r >= 0) {
					std::vector<int> x(vElem);
					x.push_back(iV[d] - r);
					w.push_back(x);
				}
				if (iV[d] + r < s.aMax) {
					std::vector<int> x(vElem);
					x.push_back(iV[d] + r);
					w.push_back(x);
				}
			}
			v = w;
			for (count j = d + 1; j < dim; j++) {
				std::vector<std::vector<int>> w;
				for (int mu = fmax(iV[j] - r, 0); mu <= fmin(iV[j] + r, s.aMax - 1); mu++) {
					for (std::vector<int> vElem : v) {
						std::vector<int> x(vElem);
						x.push_back(mu);
						w.push_back(x);
					}
				}
				v = w;
			}
			se.insert(se.end(), v.begin(), v.end());
		}
	}
	std::vector<int> seResult;
	for (std::vector<int> v : se) {
		seResult.push_back(toIndex(s, v));
	}
	sort(seResult.begin(), seResult.end());
	seResult.erase(unique(seResult.begin(), seResult.end()), seResult.end());
	return seResult;
}

const std::vector<int> MocnikGenerator::boxVolume(MocnikGenerator::State &s, const int &j, const double &r) {
	int r2 = ceil(r * s.aMax);
	std::vector<std::vector<int>> se;
	std::vector<int> tmp;
	se.push_back(tmp);
	std::vector<int> iV = fromIndex(s, j);
	// find all boxes
	for (count d = 0; d < dim; d++) {
		std::vector<std::vector<int>> w;
		for (int mu = fmax(iV[d] - r2, 0); mu <= fmin(iV[d] + r2, s.aMax - 1); mu++) {
			for (std::vector<int> sElem : se) {
				std::vector<int> x(sElem);
				x.push_back(mu);
				w.push_back(x);
			}
		}
		se = w;
	}
	std::vector<int> seResult;
	for (std::vector<int> v : se) {
		seResult.push_back(toIndex(s, v));
	}
	return seResult;
}

// EDGES TO GRAPH

void MocnikGenerator::addEdgesToGraph(Graph &G, const count &n, const double &k, const double &relativeWeight, const bool &firstRun) {
	// map vector containing the nodes resp. their positions
	MocnikGenerator::State s;
	initCellArray(s, ceil(std::pow(n, 1./dim) / k));

	// add the nodes to the state
	for (int j = 0; j < n; j++) {
		addNode(s, npm[j]);
	}

	// create the edges
	count jMax = std::pow(s.aMax, dim);
#pragma omp parallel for
	for (count j = 0; j < jMax; j++) {
		double dmNew, dm;
		NodePositionMap ns = getNodes(s, j);
		if (ns.empty()) {
			continue;
		}
		// compute the minimal distance from a node to all other nodes
		std::map<node, double> distMins;
		for (auto &it : ns) {
			distMins[it.first] = -1;
		}
		bool nodesFound = false;
		int r = -1;
		while (!nodesFound || r < get_max(distMins).second * s.aMax) {
			r++;
			for (auto &x : boxSurface(s, j, r)) {
				NodePositionMap ns2 = getNodes(s, x);
				for (auto &it : ns) {
					dm = distMins[it.first];
					for (auto &it2 : ns2) {
						if (it.first != it2.first) {
							dmNew = dist(it.second, it2.second);
							if (dmNew < dm || dm == -1) {
								nodesFound = true;
								dm = dmNew;
							}
						}
					}
					distMins[it.first] = dm;
				}
			}
		}
		// add the edges
		double kdMin, d;
		for (auto &it : ns) {
			kdMin = k * distMins[it.first];
			for (auto &x : boxVolume(s, j, kdMin)) {
				for (auto &it2 : getNodes(s, x)) {
					d = dist(it.second, it2.second);
					if (d <= kdMin && it.first != it2.first && (firstRun || !G.hasEdge(it.first, it2.first))) {
						G.addEdge(it.first, it2.first, d * relativeWeight);
					}
				}
			}
		}
	}
}

// GRAPH GENERATION

Graph MocnikGenerator::generate() {
	// create relative weights
	if (relativeWeights.empty()) {
		for (auto &n : ns) {
			relativeWeights.push_back(1);
		}
	}
	
	assert(dim > 0);
	for (auto &n : ns) {
		assert(n > 1);
	}
	for (auto &k : ks) {
		assert(k > 1);
	}
	assert(ns.size() > 0);
	assert(ks.size() == ns.size());
	assert(relativeWeights.size() == ns.size());
	
	// create graph
	Graph G(0, weighted, true);
	
	// create the nodes
	node curr = 0;
	while (curr < *std::max_element(ns.begin(), ns.end())) {
		std::vector<double> v = {};
		for (count j = 0; j < dim; j++) {
			v.push_back(Aux::Random::real());
		}
		// test wheather the new node would be contained in the ball B_{.5}(.5, ..., .5)
		if (norm(v, -.5) < .5) {
			npm.push_back(std::make_pair(G.addNode(), v));
			curr++;
		}
	}
	
	// create the edges
	for (count j = 0; j < ns.size(); j++) {
		addEdgesToGraph(G, ns[j], ks[j], relativeWeights[j], j == 0);
	}

	G.shrinkToFit();
	return G;
}

} /* namespace NetworKit */
