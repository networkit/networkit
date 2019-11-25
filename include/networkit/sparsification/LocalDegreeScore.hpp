/*
 * LocalDegreeScore.h
 *
 *  Created on: 28.08.2014
 *      Author: Gerd Lindner
 */

#ifndef NETWORKIT_SPARSIFICATION_LOCAL_DEGREE_SCORE_HPP_
#define NETWORKIT_SPARSIFICATION_LOCAL_DEGREE_SCORE_HPP_

#include <networkit/edgescores/EdgeScore.hpp>

namespace NetworKit {

/**
 * Local Degree sparsification method.
 * See 'Structure-Preserving Sparsification of Social Networks' by Lindner, Staudt, Hamann.
 */
class LocalDegreeScore : public EdgeScore<double> {

public:

    LocalDegreeScore(const Graph& G);
    virtual void run() override;
    virtual double score(edgeid eid) override;
    virtual double score(node u, node v) override;

};

}
/* namespace NetworKit */
#endif // NETWORKIT_SPARSIFICATION_LOCAL_DEGREE_SCORE_HPP_
