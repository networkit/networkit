/*
 * PrecisionRecallMetric.cpp
 *
 *  Created on: 21.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "PrecisionRecallMetric.h"

namespace NetworKit {

void PrecisionRecallMetric::generatePointsImpl() {
  generatedPoints.first.clear();
  generatedPoints.second.clear();
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
    if (generatedPoints.first.size() == 0 || generatedPoints.first.back() < recall) {
      generatedPoints.first.push_back(recall);
      generatedPoints.second.push_back(precision);
    }
  }
  INFO("thresholds.size() = ", thresholds.size());
  INFO("truePositives.size() = ", truePositives.size());
  INFO("falseNegatives.size() = ", falseNegatives.size());
  INFO("falsePositives.size() = ", falsePositives.size());
}

} // namespace NetworKit