/*
 * Overlapper.h
 *
 *  Created on: 30.10.2012
 *      Author: cls
 */

#ifndef OVERLAPPER_H_
#define OVERLAPPER_H_

#include <set>
#include <vector>


namespace EnsembleClustering {

// TODO: import from module
class Clustering;
class Graph;

class Overlapper {

public:
	Overlapper();
	virtual ~Overlapper();

	virtual Clustering run(std::set<Clustering> clusterings, Graph G) = 0;
};

} /* namespace EnsembleClustering */
#endif /* OVERLAPPER_H_ */
