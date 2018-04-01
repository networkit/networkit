	#ifndef MAHMOODYGROUPBETWEENNESS_H_
	#define MAHMOODYGROUPBETWEENNESS_H_

	#include "../graph/Graph.h"
	#include "../base/Algorithm.h"

namespace NetworKit{


class MahmoodyGroupBetweenness{

public:
	/** Constructs the GroupBetweenness-class for a given undirected graph @a G. 
	* @param setSize Size of the set of nodes.
	* @þaram epsilon Determines the accuracy of the approx. 
	*/
	MahmoodyGroupBetweenness(const Graph& g,count groupSize,double epsilon);

	/**
	* Approximately computes a set of nodes with maximum groupbetweenness. Based on the algorithm of Mahmoody,Tsourakakis and Upfal.
	*/
	void run();

	/**
	*Returns a vector of nodes containing the set of nodes with approximated maximum groupbetweenness.
	*/
	std::vector<node> groupMaxBetweenness();

protected:
	const Graph& G;
	bool hasRun = false;
	std::vector<node> maxGroup;
	count groupSize;
	double epsilon;
		
};

inline std::vector<node> MahmoodyGroupBetweenness::groupMaxBetweenness()
{
	if (!hasRun) {
		throw std::runtime_error("run method has not been called");
	}
	return maxGroup;
}

}/* namespace NetworKit */

#endif /* MAHMOODYGROUPBETWEENNESS_H_ */
