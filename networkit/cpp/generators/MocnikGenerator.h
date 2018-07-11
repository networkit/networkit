/*
 * MocnikGenerator.h
 *
 * Created on: July 7, 2018
 *	 Author: Franz-Benjamin Mocnik <mail@mocnik-science.net>
 */

#ifndef MOCNIKGENERATOR_H_
#define MOCNIKGENERATOR_H_

#include "StaticGraphGenerator.h"
#include <unordered_map>

namespace NetworKit {

/**
 * @ingroup generators
 * Creates G(d, n, p) graphs.
 */
class MocnikGenerator: public StaticGraphGenerator {
private:
		typedef std::vector<std::pair<node, std::vector<double>>> NodePositionMap;
		typedef std::vector<NodePositionMap> NodeArray;

		struct State {
				NodeArray a;
				int aMax;
		};

		NodePositionMap npm;

		void initNodeArray(State &s, const count &m);
		NodePositionMap getNodes(State &s, const int &i);
		const void addNode(State &s, const std::pair<node, std::vector<double>> &p);
		const int toIndex(State &s, const std::vector<double> &v);
		const int toIndex(State &s, const std::vector<int> &v);
		const std::vector<int> fromIndex(State &s, const int &i);
		const std::vector<int> boxSurface(State &s, const int &i, const int &r);
		const std::vector<int> boxVolume(State &s, const int &j, const double &r);
		void addEdgesToGraph(Graph &G, const count &n, const double &k, const double &relativeWeight, const bool &firstRun);

protected:
		count dim;
		std::vector<count> ns;
		std::vector<double> ks;
		bool weighted;
		std::vector<double> relativeWeights;

public:
		/**
		 * Creates random spatial graphs according to the Mocnik model.
		 *
		 * Please cite the following publications, in which you will find a
		 * description of the model:
		 * 
		 * Franz-Benjamin Mocnik, Andrew Frank: "Modelling Spatial Structures",
		 * Proceedings of the 12th Conference on Spatial Information Theory (COSIT),
		 * 2015, pages 44-64. doi: 10.1007/978-3-319-23374-1_3
		 *
		 * Improved algorithm.
		 *
		 * @param dim	Dimension of the space.
		 * @param n	Number of nodes in the graph; or a list containing the numbers
		 *	   of nodes in each layer in case of a hierarchical model.
		 * @param k	Density parameter, determining the ratio of edges to nodes; in
		 *	   case of a hierarchical model, also a list of density parameters can be
		 *	   provided.
		 * @param weighted	Determines whether weights should be added to the edges;
		 *	   in case of a hierarchical model, also a list of relative weights can be
		 *	   provided.
		 */
		MocnikGenerator(count dim, count n, double k, bool weighted=false);
		MocnikGenerator(count dim, std::vector<count> ns, double k, bool weighted=false);
		MocnikGenerator(count dim, std::vector<count> ns, std::vector<double> ks, bool weighted=false);
		MocnikGenerator(count dim, count n, double k, std::vector<double> weighted);
		MocnikGenerator(count dim, std::vector<count> ns, double k, std::vector<double> weighted);
		MocnikGenerator(count dim, std::vector<count> ns, std::vector<double> ks, std::vector<double> weighted);

		virtual Graph generate();
};

} /* namespace NetworKit */
#endif /* MOCNIKGENERATOR_H_ */
