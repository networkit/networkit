/*
 * DynamicDGSParser.cpp
 *
 *  Created on: Jun 17, 2013
 *      Author: forigem
 */

#include "DynamicDGSParser.h"

namespace NetworKit {

std::vector<std::vector<std::string>> localNodeCategories;

DynamicDGSParser::DynamicDGSParser(std::string path) : graphInitialized(false) {
	dgsFile.open(path.c_str(), std::ifstream::in);
}

DynamicDGSParser::~DynamicDGSParser() {

}

void DynamicDGSParser::initializeGraph() {
	if (! dgsFile.is_open()) {
		throw std::runtime_error("DGS input file could not be opened.");
	}
	else {
		DEBUG("Opened DGS file");

		std::string line;
		std::string cookie = "DGS004";

		// handle first line
		std::getline(dgsFile, line);
		if (!line.compare(0, cookie.size(), cookie)) { // compare prefix
			DEBUG("found magic cookie: DGS004");
		} else {
			DEBUG("First line: " << line);
			throw std::runtime_error("This does not seem to be a valid DGS file. Expected magic cookie 'DGS004' in first line");
		}

		// handle second line: optional name of file, number of clock ticks, total number of events
		std::getline(dgsFile, line);

		// Throw away the st0
		std::getline(dgsFile, line);

		graphInitialized = true;
	}
}

void DynamicDGSParser::generate() {
	if (graphInitialized == false) {
		throw std::runtime_error("Can not call generate() before graph was initialized.");
	}
	std::string line;
	bool breakTimeStep = false; // true if breaking from the while loop was due to a time step event

	while (std::getline(dgsFile, line)) {
		std::vector<std::string> split = Aux::StringTools::split(line);
		std::string tag = split[0];

		//std::unordered_map<std::string, node> edgeNames;

		if (tag.compare("st") == 0 && split.size() == 2) { // clock
			Gproxy->timeStep();
			breakTimeStep = true;
			break;

		} else if (tag.compare("an") == 0 && split.size() >= 2) { // add node
			// Get the node name from the input
			std::string nodeName = split[1];
			// Add a node to a graph, mapping it to the node name inside the nodeNames map
			nodeNames[nodeName] = Gproxy->addNode();
			if (split.size() >= 4) { // DGS with ground truth

				std::string categoriesFullString = split[2]; /// Example: category="cond-mat.stat-mech, q-fin.ST"
				std::vector<std::string> categoriesFullStringSplit = Aux::StringTools::split(categoriesFullString, '"');

				std::string categoriesCommaSeparated = categoriesFullStringSplit[1]; // Example: cond-mat.stat-mech, q-fin.ST
				std::vector<std::string> categories = Aux::StringTools::split(categoriesCommaSeparated, ',');

				std::vector<std::string> currentNodeCategories;
				for (std::string category : categories) {
					currentNodeCategories.push_back(category);
				}
				nodeCategories.push_back(currentNodeCategories);
				INFO("NODE CATEGORIES SIZE " << nodeCategories.size());

				std::string dateFullString = split[3]; // Example: date="08-1997"
				INFO("DATE SPLIT " << dateFullString);

				std::vector<std::string> dateFullStringSplit = Aux::StringTools::split(dateFullString, '"');
				INFO("DATE SPLIT1 tt" << dateFullStringSplit.size());

				std::string date = dateFullStringSplit[1];
				INFO("DATE SPLIT2 ");

				nodeDates.push_back(date);
				INFO("DATE SPLIT3 ");

			}

		} else if (tag.compare("ae") == 0 && split.size() >= 4) { // add edge
			std::string edge_from = split[2];
			std::string edge_to = split[3];
			std::string edge_name = split[1];
			Gproxy->addEdge(nodeNames[edge_from], nodeNames[edge_to], 1.0);

		} else if (tag.compare("ce") == 0 && split.size() == 3) { // update edge. Only the "weight" attribute is supported so far

			std::string from_to_edges = split[1];
			std::vector<std::string> edgesSplit = Aux::StringTools::split(from_to_edges, '-');
			std::string edge_from = edgesSplit[0];
			std::string edge_to = edgesSplit[1];

			std::string weight = split[2];
			std::vector<std::string> weightSplit = Aux::StringTools::split(weight, '=');
			double weightValue = atoi(weightSplit[1].c_str() );

			Gproxy->setWeight(nodeNames[edge_from], nodeNames[edge_to], weightValue);

		} else if (tag.compare("dn") == 0 && split.size() == 2) {
			std::string nodeName = split[1];
			node deleteNode = nodeNames[nodeName];
			// Delete the nodes only if there are no edges connected to it
			if (Gproxy->G->degree(deleteNode) == 0) {
				Gproxy->removeNode(deleteNode);
			} else {
				throw std::runtime_error("The node was not deleted, since there are edges attached to it.");
			}

		} else if (tag.compare("de") == 0 && split.size() == 2) {
			std::string from_to_edges = split[1];
			std::vector<std::string> edgesSplit = Aux::StringTools::split(from_to_edges, '-');
			std::string edge_from = edgesSplit[0];
			std::string edge_to = edgesSplit[1];
			node u = nodeNames[edge_from];
			node v = nodeNames[edge_to];

			Gproxy->removeEdge(u, v);
		}

	} // end while
	if (!breakTimeStep) {
		// while loop finished because it hit the end of the file
		throw std::logic_error("reached the end of the .dgs file");
	}



}

void DynamicDGSParser::evaluateClusterings(const Clustering& clustering) {
	INFO("eval clust");


	// Stores the category breakdown for each of the clusterings
	std::vector<std::unordered_map<std::string, int>> clusterMappings;
	clusterMappings.reserve(100000000); //should this be clustering.numberOfClusters()?

	// Stores the mapping between  clustring IDs (seem somewhat random
	// and not sequential) to "normalized" IDSs (0..numberOfClusterings-1)
	std::unordered_map<cluster, int> clusterIDMap;

	cluster normalizedID = 0;
	int clusterIDCounter = 0;

	//Vector of cluster sizes
	std::vector<int> clusterSizes;

	// Iterating through nodes and clusters to which they belong
	clustering.forEntries([&](node v, cluster c){

		// Checking if the specific cluster has been already assigned with a normalized ID
		std::unordered_map<cluster, int>::const_iterator gotClusterID = clusterIDMap.find (c);

		// If not, assigning it one
		if ( gotClusterID == clusterIDMap.end() ) {
			clusterIDMap[c] = clusterIDCounter;
			//clusterSizes.push_back(c.g);
			//c.clusterSize
			clusterIDCounter++;
		}

		// Using the normalized ID to address entries in the vector
		normalizedID = clusterIDMap[c];
		TRACE("normalizedID " << normalizedID);
		TRACE("node categories size " << nodeCategories.size());
		INFO("c " << c);
		TRACE("v " << v);

		// Assuming 1 to 1 node correspondence (VERIFY),
		// getting the categories to which this node (a scientific paper) belongs
		std::vector<std::string> currentNodeCategories = nodeCategories[v];

		for (std::string category : currentNodeCategories) {
			INFO("category " << category);

			INFO("C was fine " << clusterMappings[normalizedID].size());
			clusterMappings[normalizedID].reserve(100000);
			//clusterMappings[normalizedID].insert ( {"dummy_category",1});

			INFO("Dummy ok");
			clusterMappings[normalizedID].find (category);

			std::unordered_map<std::string, int>::const_iterator got = clusterMappings[normalizedID].find (category);
			INFO("got was fine");

			if ( got == clusterMappings[normalizedID].end() ) {
				std::cout << "not found";
				clusterMappings[normalizedID][category] = 1;
			} else {
				std::cout << "Added";
				clusterMappings[normalizedID][category] = clusterMappings[normalizedID][category] + 1;
			}
		}
	});

	std::ofstream fs("clust.csv");
	if(!fs) {
		std::cerr<<"Cannot open the output file" << std::endl;
	} else {
		INFO("FILE OK");

		for (int i=0; i<clustering.numberOfClusters(); i++) {
			fs << "cluster-" << i << '\t';
			for (const auto &pair : clusterMappings[i]) {
				fs << pair.first << '\t' << pair.second << '\t';
			}
			fs << '\n';
		}
	}

	fs.close();
}

} /* namespace NetworKit */
