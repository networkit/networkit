/*
 * EvaluationMetric.cpp
 *
 *  Created on: 17.03.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include "EvaluationMetric.h"

namespace NetworKit {

EvaluationMetric::EvaluationMetric() : testGraph(nullptr) {
}

EvaluationMetric::EvaluationMetric(const Graph& testGraph, std::vector<LinkPredictor::node_dyad_score_pair> predictions)
    : testGraph(&testGraph), predictions(predictions) {
}

void EvaluationMetric::setTestGraph(const Graph& newTestGraph) {
  testGraph = &newTestGraph;
}

void EvaluationMetric::setPredictions(std::vector<LinkPredictor::node_dyad_score_pair> newPredictions) {
  predictions = newPredictions;
}

void EvaluationMetric::generatePoints(count numThresholds) {
  if (testGraph == nullptr) {
    throw std::logic_error("Set testGraph first.");
  } else if (predictions.size() == 0) {
    throw std::logic_error("Set predictions first.");
  }
  if (predictions.size() < numThresholds || numThresholds == 0) {
    this->numThresholds = predictions.size() + 1;
  } else {
    this->numThresholds = numThresholds;
  }
  calculateStatisticalMeasures();
  generatePointsImpl();
}

std::pair<std::vector<double>, std::vector<double>> EvaluationMetric::getPoints() const {
  return generatedPoints;
}

double EvaluationMetric::areaUnderCurve() const {
  double sum = 0;
  for (index i = 0; i < generatedPoints.first.size() - 1; ++i) {
    // Trapezoid rule
    sum += 0.5 * (generatedPoints.first.at(i + 1) - generatedPoints.first.at(i))
           * (generatedPoints.second.at(i) + generatedPoints.second.at(i + 1));
  }
  return sum;
}

void EvaluationMetric::calculateStatisticalMeasures() {
  truePositives.clear();
  falsePositives.clear();
  trueNegatives.clear();
  falseNegatives.clear();
  setPositivesAndNegatives();
  for (index i = 0; i < numThresholds; ++i) {
    count test = predictions.size() * (1.0 * i / (numThresholds - 1));
    addStatisticsForThresholdIndex(test);
  }
}

void EvaluationMetric::addStatisticsForThresholdIndex(index thresholdIndex) {
  count tmpTruePositives = 0;
  count tmpFalsePositives = 0;
  count tmpTrueNegatives = 0;
  count tmpFalseNegatives = 0;
  #pragma omp parallel for
  for (index i = 0; i < predictions.size(); ++i) {
    bool hasEdge = testGraph->hasEdge(predictions[i].first.first, predictions[i].first.second);
    if (i < thresholdIndex) { // predicted as an edge
      if (hasEdge) {
        #pragma omp atomic
        tmpTruePositives++;
      } else {
        #pragma omp atomic
        tmpFalsePositives++;
      }
    } else {
      if (hasEdge) {
        #pragma omp atomic
        tmpFalseNegatives++;
      } else {
        #pragma omp atomic
        tmpTrueNegatives++;
      }
    }
  }
  truePositives.push_back(tmpTruePositives);
  falsePositives.push_back(tmpFalsePositives);
  trueNegatives.push_back(tmpTrueNegatives);
  falseNegatives.push_back(tmpFalseNegatives);
}

void EvaluationMetric::setPositivesAndNegatives() {
  numPositives = 0;
  #pragma omp parallel for
  for (index i = 0; i < predictions.size(); ++i) {
    if (testGraph->hasEdge(predictions[i].first.first, predictions[i].first.second)) {
      #pragma omp atomic
      numPositives++;
    }
  }
  numNegatives = predictions.size() - numPositives;
}

} // namespace NetworKit