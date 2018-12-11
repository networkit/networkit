/*
 * Volume.cpp
 *
 * Created on: July 7, 2018
 * Author: Franz-Benjamin Mocnik <mail@mocnik-science.net>
 */

#include "Volume.h"
#include "../graph/Graph.h"
#include <unordered_map>

namespace NetworKit {

std::unordered_map<node, double> Volume::nodesWithinDistance(const Graph &G, double r, node n) {
	std::unordered_map<node, double> ms;
	ms.insert(std::make_pair(n, 0));
	std::vector<node> msToCheck;
	msToCheck.push_back(n);
	double r2;
	while (!msToCheck.empty()) {
		std::vector<node> msToCheckNew;
		for (auto &m : msToCheck) {
			for (auto &m2 : G.neighbors(m)) {
				r2 = ms[m] + G.weight(m, m2);
				if (ms.count(m2) == 0) {
					if (r2 <= r) {
						ms[m2] = r2;
						msToCheckNew.push_back(m2);
					}
				} else {
					ms[m2] = fmin(ms[m2], r2);
				}
			}
		}
		msToCheck = msToCheckNew;
	}
	return ms;
}

double Volume::volume(const Graph &G, const double r, const count samples) {
	double x = 0;
	for (count j = 0; j < samples; j++) {
		x += Volume::nodesWithinDistance(G, r, G.randomNode()).size();
	}
	return x / samples;
}

std::vector<double> Volume::volume(const Graph &G, const std::vector<double> rs, const count samples) {
	std::vector<double> xs;
	for (auto &r : rs) {
		xs.push_back(0);
	}
	double rmax = *std::max_element(std::begin(rs), std::end(rs));
	for (count j = 0; j < samples; j++) {
		std::unordered_map<node, double> ms = Volume::nodesWithinDistance(G, rmax, G.randomNode());
		count i = 0;
		for (auto &r : rs) {
			for (auto &it : ms) {
				if (it.second <= r) {
					xs[i] += 1;
				}
			}
			i++;
		}
	}
	std::vector<double> ys;
	for (auto &x : xs) {
		ys.push_back(x / samples);
	}
	return ys;
}

} /* namespace NetworKit */
