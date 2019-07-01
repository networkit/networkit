/*
 * GeometricInhomogenousGenerator.h
 *
 *  Created on: 09.05.2019
 *      Author: Manuel Penschuck <networkit@manuel.jetzt>
 */

#ifndef GEOMETRICINHOMOGENOUS_H_
#define GEOMETRICINHOMOGENOUS_H_

#include <limits>
#include <vector>

#include <networkit/generators/StaticGraphGenerator.hpp>

namespace NetworKit {

class GeometricInhomogenousGenerator: public StaticGraphGenerator {
public:
    using Coordinate = std::vector<double>;

    /**
     * @param[in] n Number of nodes
     * @param[in] avgDegree The desired average degree with 0 < avgDegree < n-1
     * @param[in] powerlawExp Target exponent of power-law distribution with powerlawExp > 2.0
     * @param[in] alpha parameter (1/temperature) with alpha > 1
     * @param[in] dimension Dimension of the underlying geometry with 1 <= dimension <= 5
     */
    GeometricInhomogenousGenerator(count n, double avgDegree, double powerlawExp=3, double alpha=std::numeric_limits<double>::infinity(), unsigned dim=1);

    /**
     * Construct generator from weights that are then scale to match avgDegree and T.
     *
     * @param[in] points Coordinates of points
     * @param[in] weights Unscaled weights (assumed to be powerlaw distributed)
     * @param[in] avgDegree The desired average degree with 0 < avgDegree < n-1
     * @param[in] alpha parameter (1/temperature) with alpha > 1
     *
     * @warning points and weights are moved into the container. The we're not using
     * rvalue refs because Cython does not handle them.
     */
    GeometricInhomogenousGenerator(std::vector<Coordinate> points, std::vector<double> weights, double avgDegree, double alpha);

    /**
     * Construct generator from *already scaled* weights.
     *
     * @param[in] points Coordinates of points
     * @param[in] weights *Scaled* weights
     * @param[in] alpha parameter (1/temperature) with alpha > 1
     *
     * @warning points and weights are moved into the container. The we're not using
     * rvalue refs because Cython does not handle them.
     */
    GeometricInhomogenousGenerator(std::vector<Coordinate> points, std::vector<double> weights, double alpha);

    // Add virtual destructor
    virtual ~GeometricInhomogenousGenerator() = default;

    /// @return Graph to be generated according to parameters specified in constructor freeing memory by deleting the input point set.
    Graph generate() override;

    /// @return Graph to be generated according to parameters specified in constructor keeping the input point set.
    Graph generateKeepingInput();

    /**
     * @return  Weights used to generate the graph
     * @warning The data is destroyed if generate is called with keep_input = false (default).
     */
    const std::vector<double>& weights() const noexcept {
        return pointWeights;
    }

    /**
     * @return  Point positions used to generate the graph
     * @warning The data is destroyed if generate is called with keep_input = false (default).
     */
    const std::vector<Coordinate>& positions() const noexcept {
        return pointPositions;
    }

private:
    double alpha;

    std::vector<Coordinate> pointPositions;
    std::vector<double> pointWeights;

    void checkInputParameters() const;
};

} // NetworKit

#endif // GEOMETRICINHOMOGENOUS_H_
