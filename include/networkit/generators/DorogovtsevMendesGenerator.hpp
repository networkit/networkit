/*
* DorogovtsevMendesGenerator.h
*
*  Created on: 27.05.2014
*      Author: Christian Staudt
*/

#ifndef NETWORKIT_GENERATORS_DOROGOVTSEV_MENDES_GENERATOR_HPP_
#define NETWORKIT_GENERATORS_DOROGOVTSEV_MENDES_GENERATOR_HPP_

#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

/**
 * @ingroup generators
 */
class DorogovtsevMendesGenerator: public StaticGraphGenerator {

public:
    /**
    * TODO:
    *
    * @param nNodes 	number of nodes in target graph
    */
    DorogovtsevMendesGenerator(count nNodes);

    virtual Graph generate();

protected:
        count nNodes;

};

} /* namespace NetworKit */
#endif // NETWORKIT_GENERATORS_DOROGOVTSEV_MENDES_GENERATOR_HPP_
