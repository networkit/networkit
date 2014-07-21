#include "BidirectionalBFS.h"
#include "../auxiliary/Log.h"


namespace NetworKit {

BidirectionalBFS::BidirectionalBFS(const Graph& G) : G(G), forward(G,0), backward(G,0) {

}

void BidirectionalBFS::run(node s, node t) {
	TRACE("bidirectional BFS run method called with ", s, " and ", t);
	forward.init(s);
	backward.init(t);
	std::vector<node> junction;
	bool forwardNext = true;
	double minDist = std::numeric_limits<edgeweight>::max();

	auto mustContinue = [](double minDist, double forwardMin, double BackwardMin) {
		return (minDist > (forwardMin + BackwardMin));
	};

	while (!forward.isFinished() && !backward.isFinished() && mustContinue(minDist,forward.getCurrentMin(),backward.getCurrentMin())) {
		bool path_found = false;
		node u;
		if (forwardNext) {
			u = forward.settleNext();
			TRACE("forward settleNext resulted in node ",u);
			if (backward.wasVisited(u)) path_found = true;
		} else {
			u = backward.settleNext();
			TRACE("backward settleNext resulted in node ",u);
			if (forward.wasVisited(u)) path_found = true;
		}
		forwardNext = !forwardNext;

		if (path_found) {
			double distance = forward.extractDistance(u) + backward.extractDistance(u);
			if (distance < minDist) {
				TRACE("new minimal path found");
				minDist = distance;
				junction = {u};
			} else if (distance == minDist) {
				TRACE("new node with same distance found");
				junction.push_back(u);
			}

		}
	}
}

}
