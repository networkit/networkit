/*
 * BitManipulation.h
 *
 *  Created on: 03. May 2019
 *      Author: Christopher Weyand <Christopher.Weyand@hpi.de>, Manuel Penschuck <networkit@manuel.jetzt>
 */

#ifndef GENERATORS_GIRGS_WEIGHT_SCALING_H_
#define GENERATORS_GIRGS_WEIGHT_SCALING_H_

#include <vector>

namespace NetworKit {
namespace girgs {

double estimateWeightScaling(const std::vector<double> &weights, double desiredAvgDegree, int dimension, double alpha);

double estimateWeightScalingThreshold(const std::vector<double> &weights, double desiredAvgDegree, int dimension);

} // namespace girgs
} // namespace NetworKit

#endif // GENERATORS_GIRGS_WEIGHT_SCALING_H_
