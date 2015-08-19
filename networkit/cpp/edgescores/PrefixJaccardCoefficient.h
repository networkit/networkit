#ifndef PREFIXJACCARDCOEFFICIENT_H
#define PREFIXJACCARDCOEFFICIENT_H

#include "../edgescores/EdgeScore.h"

namespace NetworKit {

template <typename AttributeT>
class PrefixJaccardCoefficient : public EdgeScore<double> {

public:
	PrefixJaccardCoefficient(const Graph& G, const std::vector<AttributeT>& attribute);
	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;
	virtual void run() override;


private:
	const std::vector<AttributeT>& inAttribute;

};

}

#endif // PREFIXJACCARDCOEFFICIENT_H
