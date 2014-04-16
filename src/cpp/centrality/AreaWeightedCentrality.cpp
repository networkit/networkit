/*
 * AreaWeightedCentrality.cpp
 *
 *  Created on: 14.04.2014
 *      Author: Maximilian Vogel
 */

#include "AreaWeightedCentrality.h"
#include "../io/LineFileReader.h"
#include <cmath>
#include <string>

namespace NetworKit {

AreaWeightedCentrality::AreaWeightedCentrality(const Graph& G, bool normalized) : Centrality(G, normalized), cosine(cosine) {
}

void AreaWeightedCentrality::initCoordinates(std::string lat_file, std::string lon_file) {
	LineFileReader reader;
	auto toDouble = [](std::vector<std::string> content) {
		std::vector<double> result;
		for (auto e : content) {
			if (!e.empty()) {
				result.push_back(std::stod(e));
			}
		}
		std::cout << std::endl;
		return result;
	};
	auto toCosine = [](std::vector<double> latitude) {
		std::vector<double> result;
		for (auto val : latitude) {
			result.push_back(std::cos(val * 3.14159265 / 180.0));
		}
		return result;
	};
/*	auto content = reader.read(lat_file);
	INFO("Reading coordinates successful");
	auto latitude = toDouble(content);
	for (index i = 0, end = latitude.size(); i < end; ++i) {
		std::cout << latitude[i] << " ";
	}
	std::cout << std::endl;
	INFO("Converting to double successful");
	cosine = toCosine(latitude);
	INFO("Calculating cosine successful");*/
	cosine = toCosine(toDouble(reader.read(lat_file)));	
}

void AreaWeightedCentrality::run() {
	scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);
	double sum = 0.0f;

	G.parallelForNodes([&](node u) {
		double current = 0.0f;	
		G.forNeighborsOf(u, [&](node v){
			current += cosine[v]; //G.weight(u,v) not necessary at all...
		        //	* std::cos(latitude[v] * 60 / 3.14159265)); //retrieve cos from somewhere else
		});
		scoreData[u] +=	current;
		sum += current;

	});

	if (normalized) {
		G.parallelForNodes([&](node u) {
			scoreData[u] = scoreData[u] / sum;
		});
	}
}
/*
std::vector<double> scores() {
	return scoreData;
}

double score(node u) {
	return scoreData[u];
}*/


} /* namespace NetworKit */

