/*
 * ROCMetric.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "ROCMetric.h"

namespace NetworKit {

std::pair<std::vector<double>, std::vector<double>> ROCMetric::generatePoints() {
  std::pair<std::vector<double>, std::vector<double>> points;
  points.first.reserve(thresholds.size());
  points.second.reserve(thresholds.size());

  INFO("falsePositives.size() = ", falsePositives.size());

  for (index i = 0; i < thresholds.size(); ++i) {
    double falsePositiveRatio = 1.0 * falsePositives.at(i) / numNegatives;
    if (points.first.size() == 0 || points.first.back() < falsePositiveRatio) {
      points.first.push_back(falsePositiveRatio);
      points.second.push_back(1.0 * truePositives.at(i) / numPositives);
    }
  }
  return points;
}

} // namespace NetworKit