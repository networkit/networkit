#ifndef PREFIXJACCARDCOEFFICIENT_H
#define PREFIXJACCARDCOEFFICIENT_H

#include "../graph/Graph.h"

namespace NetworKit {

template <typename AttributeT>
class PrefixJaccardCoefficient {
private:
	const Graph& G;
	const std::vector<AttributeT>& inAttribute;
	std::vector<double> outAttribute;
	bool hasRun;

public:
	PrefixJaccardCoefficient(const Graph& G, const std::vector<AttributeT>& attribute);
	void run();
	std::vector<double> getAttribute();
};

}

#endif // PREFIXJACCARDCOEFFICIENT_H
