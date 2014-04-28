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
#include <fstream>

namespace NetworKit {

AreaWeightedCentrality::AreaWeightedCentrality(const Graph& G, bool normalized) : Centrality(G, normalized) {
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

	// read coordinate files.
	this->lat = toDouble(reader.read(lat_file));
	this->lon = toDouble(reader.read(lon_file));

	// set Coordinates to the graph. TODO: new class to read Graph with coordinates.
	/*G.initCoordinates();
	for ( index i = 0, end = lat.size(); i < end; ++i ) {
		G.setCoordinate(i,Point<float>(lat[i],lon[i]));
	}*/

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
	// precompute cosine values.
	cosine = toCosine(this->lat);	
}

void AreaWeightedCentrality::run() {
	scoreData = std::vector<double>(G.upperNodeIdBound(), 0.0);
	double sum = 0.0f;

	G.parallelForNodes([&](node u) {
		double current = 0.0f;	
		G.forNeighborsOf(u, [&](node v){
			current += cosine[v];
		});
		scoreData[u] +=	current;
		sum += current;
	});

	if (normalized) {
		G.parallelForNodes([&](node u) {
			scoreData[u] /= sum;
		});
	}
}

void AreaWeightedCentrality::writeValuesToCSV(std::string file, double scale) {
	auto scaleTo = [&](std::vector<double> val, double factor) {
		std::vector<double> result(val.size());
		for (index i = 0, end = val.size(); i < end; ++i) {
			result[i] = factor * val[i];
		}
		return result;
	};
	auto val = scores();
	if (scale != 0.0f) {
		double factor = scale / ranking()[0].second;
		val = scaleTo(val, factor);
	} 
	std::ofstream output(file, std::ofstream::out);
	for ( index i = 0, end = val.size(); i < end; ++i) {
		output << lon[i] << "," << lat[i] << "," << val[i] << std::endl;
	}
	output.close();
}

void AreaWeightedCentrality::writeClusteringToCSV(Partition values, std::string file) {
	std::ofstream output(file, std::ofstream::out);
	for ( index i = 0, end = values.numberOfElements(); i < end; ++i) {
		output << lon[i] << "," << lat[i] << "," << values[i] << std::endl;
	}
	output.close();
}


/*
std::vector<double> scores() {
	return scoreData;
}

double score(node u) {
	return scoreData[u];
}*/


} /* namespace NetworKit */

