/*
 * ChibaNishizekiQuadrangleEdgeScore.h
 *
 *  Created on: 18.11.2014
 *      Author: Michael Hamann, Gerd Lindner
 */

#ifndef NETWORKIT_EDGESCORES_CHIBA_NISHIZEKI_QUADRANGLE_EDGE_SCORE_HPP_
#define NETWORKIT_EDGESCORES_CHIBA_NISHIZEKI_QUADRANGLE_EDGE_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

class ChibaNishizekiQuadrangleEdgeScore : public EdgeScore<count> {

public:
    ChibaNishizekiQuadrangleEdgeScore(const Graph& G);
    virtual count score(edgeid eid) override;
    virtual count score(node u, node v) override;
    virtual void run() override;
};

} // namespace NetworKit

#endif // NETWORKIT_EDGESCORES_CHIBA_NISHIZEKI_QUADRANGLE_EDGE_SCORE_HPP_
