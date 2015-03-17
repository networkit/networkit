/*
 * EvaluationCurve.cpp
 *
 *  Created on: 17.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "EvaluationCurve.h"

namespace NetworKit {

EvaluationCurve::EvaluationCurve(const Graph& testGraph, std::vector<LinkPredictor::node_dyad_score_pair> dyadScorePairs)
    : testGraph(testGraph), dyadScorePairs(dyadScorePairs) {
}

std::pair<std::vector<double>, std::vector<double>> EvaluationCurve::getPoints() const {
  return generatedPoints;
}

double EvaluationCurve::areaUnderCurve() const {
  double sum = 0;
  for (index i = 0; i < generatedPoints.first.size() - 1; ++i) {
    // Trapezoid rule
    sum += 0.5 * (generatedPoints.first.at(i + 1) - generatedPoints.first.at(i))
           * (generatedPoints.second.at(i) + generatedPoints.second.at(i + 1));
  }
  return sum;
}

} // namespace NetworKit