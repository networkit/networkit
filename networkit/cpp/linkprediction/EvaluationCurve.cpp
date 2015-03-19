/*
 * EvaluationCurve.cpp
 *
 *  Created on: 17.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "EvaluationCurve.h"

namespace NetworKit {

EvaluationCurve::EvaluationCurve() : testGraph(nullptr) {
}

EvaluationCurve::EvaluationCurve(const Graph& testGraph, std::vector<LinkPredictor::node_dyad_score_pair> predictions)
    : testGraph(&testGraph), predictions(predictions) {
}

void EvaluationCurve::setTestGraph(const Graph& newTestGraph) {
  testGraph = &newTestGraph;
}

void EvaluationCurve::setPredictions(std::vector<LinkPredictor::node_dyad_score_pair> newPredictions) {
  predictions = newPredictions;
}

void EvaluationCurve::generatePoints() {
  if (testGraph == nullptr) {
    throw std::logic_error("Set testGraph first.");
  } else if (predictions.size() == 0) {
    throw std::logic_error("Set predictions first.");
  }
  generatePointsImpl();
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