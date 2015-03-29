/*
 * ROCMetric.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "ROCMetric.h"

namespace NetworKit {

void ROCMetric::generatePointsImpl() {
  generatedPoints.first.clear();
  generatedPoints.second.clear();
  for (index i = 0; i < thresholds.size(); ++i) {
    double falsePositiveRatio = 1.0 * falsePositives.at(i) / numNegatives;
    if (generatedPoints.first.size() == 0 || generatedPoints.first.back() < falsePositiveRatio) {
      generatedPoints.first.push_back(falsePositiveRatio);
      generatedPoints.second.push_back(1.0 * truePositives.at(i) / numPositives);
    }
  }
}

} // namespace NetworKit