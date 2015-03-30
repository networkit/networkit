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

std::pair<std::vector<double>, std::vector<double>> EvaluationMetric::getPoints() const {
  return generatedPoints;
}

void EvaluationMetric::generatePoints(count numThresholds) {
  if (testGraph == nullptr) {
    throw std::logic_error("Set testGraph first.");
  } else if (predictions.size() == 0) {
    throw std::logic_error("Set predictions first.");
  }
  // Don't overshoot with the number of thresholds.
  if (predictions.size() < numThresholds || numThresholds == 0) {
    numThresholds = predictions.size() + 1;
  }
  thresholds.clear();
  thresholds.reserve(numThresholds);
  for (index i = 0; i < numThresholds; ++i) {
    // Percentile calculation through nearest rank method.
    thresholds.push_back(std::ceil(predictions.size() * (1.0 * i / (numThresholds - 1))));
  }
  calculateStatisticalMeasures();
  generatePointsImpl();
}

double EvaluationMetric::areaUnderCurve() const {
  double sum = 0;
  for (index i = 0; i < generatedPoints.first.size() - 1; ++i) {
    // Trapezoid rule
    sum += 0.5 * (generatedPoints.first[i + 1] - generatedPoints.first[i])
           * (generatedPoints.second[i] + generatedPoints.second[i + 1]);
  }
  return sum;
}

void EvaluationMetric::calculateStatisticalMeasures() {
  #pragma omp parallel sections
  {
    #pragma omp section
    setTrueAndFalsePositives();
    #pragma omp section
    setTrueAndFalseNegatives();
    #pragma omp section
    setPositivesAndNegatives();
  }
}

void EvaluationMetric::setTrueAndFalsePositives() {
  count tmpTruePositives = 0;
  count tmpFalsePositives = 0;
  truePositives.clear();
  falsePositives.clear();
  std::vector<index>::iterator thresholdIt = thresholds.begin();
  for (index i = 0; i < predictions.size(); ++i) {
    // Check in the beginning to catch threshold = 0 as well.
    if (thresholdIt != thresholds.end() && i == *thresholdIt) {
      truePositives.push_back(tmpTruePositives);
      falsePositives.push_back(tmpFalsePositives);
      ++thresholdIt;
    }
    bool hasEdge = testGraph->hasEdge(predictions[i].first.first, predictions[i].first.second);
    if (hasEdge) {
      tmpTruePositives++;
    } else {
      tmpFalsePositives++;
    }
  }
  if (thresholdIt != thresholds.end()) {
    truePositives.push_back(tmpTruePositives);
    falsePositives.push_back(tmpFalsePositives);
  }
}

void EvaluationMetric::setTrueAndFalseNegatives() {
  count tmpTrueNegatives = 0;
  count tmpFalseNegatives = 0;
  trueNegatives.clear();
  falseNegatives.clear();
  std::vector<index>::reverse_iterator thresholdIt = thresholds.rbegin();
  INFO("thresholds.size() = ", thresholds.size());
  INFO("thresholds.front() = ", thresholds.front());
  INFO("thresholds.back() = ", thresholds.back());
  for (index i = predictions.size(); i > 0; --i) {
    // Check in the beginning to catch threshold = predictions.size().
    if (thresholdIt != thresholds.rend() && i == *thresholdIt) {
      trueNegatives.push_back(tmpTrueNegatives);
      falseNegatives.push_back(tmpFalseNegatives);
      if (*thresholdIt != 0)
        ++thresholdIt;
    }
    bool hasEdge = testGraph->hasEdge(predictions[i - 1].first.first, predictions[i - 1].first.second);
    if (hasEdge) {
      tmpFalseNegatives++;
    } else {
      tmpTrueNegatives++;
    }
  }
  if (thresholdIt != thresholds.rend()) {
    trueNegatives.push_back(tmpTrueNegatives);
    falseNegatives.push_back(tmpFalseNegatives);
  }
  ++thresholdIt;
  std::reverse(trueNegatives.begin(), trueNegatives.end());
  std::reverse(falseNegatives.begin(), falseNegatives.end());
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