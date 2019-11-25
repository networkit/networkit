
#ifndef NETWORKIT_SPARSIFICATION_SCAN_STRUCTURAL_SIMILARITY_SCORE_HPP_
#define NETWORKIT_SPARSIFICATION_SCAN_STRUCTURAL_SIMILARITY_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

    class SCANStructuralSimilarityScore : public EdgeScore<double> {

    public:
        SCANStructuralSimilarityScore(const Graph& G, const std::vector<count>& triangles);
        virtual void run() override;
        virtual double score(edgeid eid) override;
        virtual double score(node u, node v) override;

    private:
        const std::vector<count>& triangles;

    };

} // namespace NetworKit

#endif // NETWORKIT_SPARSIFICATION_SCAN_STRUCTURAL_SIMILARITY_SCORE_HPP_
