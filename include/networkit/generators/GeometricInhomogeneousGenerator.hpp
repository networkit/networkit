/*
 * GeometricInhomogeneousGenerator.hpp
 *
 *  Created on: 26.10.2023
 *      Author: Christopher Weyand, Manuel Penschuck <networkit@manuel.jetzt>
 */

#ifndef GEOMETRICINHOMOGENEOUS_H_
#define GEOMETRICINHOMOGENEOUS_H_

#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

class GeometricInhomogeneousGenerator: public StaticGraphGenerator {
public:
    /**
     * @param[in] n Number of nodes
     * @param[in] avgDegree The desired average degree with 0 < avgDegree < n-1
     * @param[in] powerlawExp Target exponent of power-law distribution with powerlawExp > 2.0
     * @param[in] temperature Temperature adds noise to the instance with 0 <= T <= 1
     * @param[in] dimension Dimension of the underlying geometry with 1 <= dimension <= 5
     */
    GeometricInhomogeneousGenerator(count n, double avgDegree, double powerlawExp=3, double temperature=0.0, unsigned dim=1);

    // Add virtual destructor
    virtual ~GeometricInhomogeneousGenerator() = default;

    /// @return Graph to be generated according to parameters specified in constructor
    Graph generate() override;

private:
    count n;
    double avgDegree;
    double powerlawExp;
    double alpha;
    unsigned dim;
};

} // NetworKit

#endif // GEOMETRICINHOMOGENEOUS_H_