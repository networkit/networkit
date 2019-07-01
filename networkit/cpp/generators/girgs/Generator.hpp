/*
 * Generator.hpp
 *
 *  Created on: 03. May 2019
 *      Author: Christopher Weyand <Christopher.Weyand@hpi.de>, Manuel Penschuck <networkit@manuel.jetzt>
 *
 * Code is adopted from https://github.com/chistopher/girgs
 *
 */

#ifndef GENERATORS_GIRGS_GENERATOR_H_
#define GENERATORS_GIRGS_GENERATOR_H_

#include <vector>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace girgs {

/**
 * @brief
 *  The weights are sampled according to a power law distribution between [1, n)
 *
 * @param n
 *  The size of the graph. Should match with size of positions.
 * @param ple
 *  The power law exponent to sample the new weights. Should be 2.0 to h3.0.
 *
 * @return
 *  The weights according to the desired distribution.
 */
std::vector<double> generateWeights(int n, double ple, bool parallel = true);

/**
 * @brief
 *  Samples d dimensional coordinates for n points on a torus \f$[0,1)^d\f$.
 *
 * @param n
 *  Size of the graph.
 * @param dimension
 *  Dimension of the geometry.
 *
 * @return
 *  The positions on a torus. All inner vectors have the same length.
 */
std::vector<std::vector<double>> generatePositions(int n, int dimension, bool parallel = true);

/**
 * @brief
 *  Scales all weights so that the expected average degree equals desiredAvgDegree.
 *  Implemented as binary search over an estimation function.
 *
 * @bug
 *  For \f$\alpha > 10\f$ we use the estimation for threshold graphs due to numerical difficulties.
 *  This leads to slightly higher degrees than desired.
 *  Also I experienced inaccurate results for \f$9 \leq \alpha < 10\f$.
 *
 * @param weights
 *  The weights to be modified.
 * @param desiredAvgDegree
 *  The desired average degree.
 * @param dimension
 *  Dimension of the underlying geometry. Should equal the dimensionality of the positions.
 * @param alpha
 *  Parameter of the algorithm. Should be the same as for the generation process.
 *
 * @return
 *  The scaling s applied to all weights.
 *  The constant c hidden in the theta of the edge probabilities is \f$s^\alpha\f$ for \f$\alpha < \infty\f$
 *  and \f$s^{1/d}\f$ in the threshold case.
 */
double scaleWeights(std::vector<double>& weights, double desiredAvgDegree, int dimension, double alpha);

/**
 * @brief
 *  Samples edges according to weights and positions.
 *  An edge between node u and v is formed with probability \f$ \left(\frac{w_u w_v / W}{|| x_u - x_v ||^d}\right)^\alpha \f$ or 1.0 if the term exceeds 1.0.
 *
 * @param weights
 *  Power law distributed weights.
 * @param positions
 *  Positions on a torus. All inner vectors should have the same length indicating the dimension of the torus.
 * @param alpha
 *  Edge probability parameter.
 *
 * @return
 *  NetworKit graph object
 */
Graph generateEdges(std::vector<double>& weights, std::vector<std::vector<double>>& positions, double alpha, bool keep_input = true);

} // namespace girgs
} // namespace NetworKit

#endif // GENERATORS_GIRGS_GENERATOR_H_
