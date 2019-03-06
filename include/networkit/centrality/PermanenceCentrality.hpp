#ifndef PERMANENCECENTRALITY_H
#define PERMANENCECENTRALITY_H

#include "Centrality.hpp"
#include "../structures/Partition.hpp"

namespace NetworKit {

// TODO documentation completely missing!
class PermanenceCentrality : public Algorithm {
public:
	PermanenceCentrality(const Graph &G, const Partition &P);
	void run();
	double getPermanence(node u);
	double getIntraClustering(node u);
private:
	const Graph &G;
	const Partition &P;
	std::vector<index> inBegin;
	std::vector<node> inEdges;
	std::vector<bool> marker;
};


}

#endif // PERMANENCECENTRALITY_H
