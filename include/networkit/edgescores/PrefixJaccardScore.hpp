#ifndef NETWORKIT_EDGESCORES_PREFIX_JACCARD_SCORE_HPP_
#define NETWORKIT_EDGESCORES_PREFIX_JACCARD_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

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

#endif // NETWORKIT_EDGESCORES_PREFIX_JACCARD_SCORE_HPP_
