/*
 * PrecisionRecallMetric.cpp
 *
 *  Created on: 21.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "PrecisionRecallMetric.h"

namespace NetworKit {

std::pair<std::vector<double>, std::vector<double>> PrecisionRecallMetric::generatePoints() {
  std::pair<std::vector<double>, std::vector<double>> points;
  points.first.reserve(thresholds.size());
  points.second.reserve(thresholds.size());
  for (index i = 0; i < thresholds.size(); ++i) {
    double recall = 1;
    double precision = 1;
    if (truePositives.at(i) > 0 || falseNegatives.at(i) > 0) {
      recall = 1.0 * truePositives.at(i) / (truePositives.at(i) + falseNegatives.at(i));
      //INFO("Recall = ", truePositives.at(i), " / (", truePositives.at(i), " + ", falseNegatives.at(i), ") = ", recall, ".");
    } else {
      //INFO("Recall = 1");
    }
    if (truePositives.at(i) > 0 || falsePositives.at(i) > 0) {
      precision = 1.0 * truePositives.at(i) / (truePositives.at(i) + falsePositives.at(i));
      //INFO("Precision = ", truePositives.at(i), " / (", truePositives.at(i), " + ", falsePositives.at(i), ") = ", precision, ".");
    } else {
      //INFO("Precision = 1");
    }
    if (points.first.size() == 0 || points.first.back() < recall) {
      points.first.push_back(recall);
      points.second.push_back(precision);
    }
  }
  return points;
}

} // namespace NetworKit