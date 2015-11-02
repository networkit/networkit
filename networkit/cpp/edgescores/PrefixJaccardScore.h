#ifndef PREFIXJACCARDSCORE_H
#define PREFIXJACCARDSCORE_H

#include "../edgescores/EdgeScore.h"

namespace NetworKit {

template <typename AttributeT>
class PrefixJaccardScore : public EdgeScore<double> {

public:
	PrefixJaccardScore(const Graph& G, const std::vector<AttributeT>& attribute);
	virtual double score(edgeid eid) override;
	virtual double score(node u, node v) override;
	virtual void run() override;


private:
	const std::vector<AttributeT>& inAttribute;

};

}

#endif // PREFIXJACCARDSCORE_H
