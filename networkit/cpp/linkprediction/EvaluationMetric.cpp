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
  // Determine the correct number of thresholds to use.
  if (predictions.size() < numThresholds || numThresholds == 0) {
    this->numThresholds = predictions.size() + 1;
  } else {
    this->numThresholds = numThresholds;
  }
  std::vector<index> thresholdIndices;
  thresholdIndices.reserve(this->numThresholds);
  for (index i = 0; i < this->numThresholds; ++i) {
    INFO("Added threshold ", std::ceil(predictions.size() * (1.0 * i / (this->numThresholds - 1))));
    thresholdIndices.push_back(std::ceil(predictions.size() * (1.0 * i / (this->numThresholds - 1))));
  }

  #pragma omp parallel sections
  {
    #pragma omp section
    setTrueAndFalsePositives(thresholdIndices);
    #pragma omp section
    setTrueAndFalseNegatives(thresholdIndices);
    #pragma omp section
    setPositivesAndNegatives();
  }
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
    count thresholdIndex = std::ceil(predictions.size() * (1.0 * i / (numThresholds - 1)));
    addStatisticsForThresholdIndex(thresholdIndex);
  }
}

void EvaluationMetric::setTrueAndFalsePositives(std::vector<index> thresholdIndices) {
  count tmpTruePositives = 0;
  count tmpFalsePositives = 0;
  std::vector<index>::iterator thresholdIt = thresholdIndices.begin();
  for (index i = 0; i < predictions.size(); ++i) {
    // Check in the beginning to catch threshold = 0 as well.
    if (thresholdIt != thresholdIndices.end() && i == *thresholdIt) {
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
  if (thresholdIt != thresholdIndices.end()) {
    truePositives.push_back(tmpTruePositives);
    falsePositives.push_back(tmpFalsePositives);
  }
}

void EvaluationMetric::setTrueAndFalseNegatives(std::vector<index> thresholdIndices) {
  count tmpTrueNegatives = 0;
  count tmpFalseNegatives = 0;
  std::vector<index>::reverse_iterator thresholdIt = thresholdIndices.rbegin();
  for (index i = predictions.size(); i > 0; --i) {
    index correctedIndex = i - 1;
    // Check in the beginning to catch threshold = predictions.size().
    if (thresholdIt != thresholdIndices.rend() && correctedIndex <= *thresholdIt) {
      trueNegatives.push_back(tmpTrueNegatives);
      falseNegatives.push_back(tmpFalseNegatives);
      if (*thresholdIt != 0)
        ++thresholdIt;
    }
    bool hasEdge = testGraph->hasEdge(predictions[correctedIndex].first.first, predictions[correctedIndex].first.second);
    if (hasEdge) {
      tmpFalseNegatives++;
    } else {
      tmpTrueNegatives++;
    }
  }
  if (thresholdIt != thresholdIndices.rend()) {
    trueNegatives.push_back(tmpTrueNegatives);
    falseNegatives.push_back(tmpFalseNegatives);
  }
  std::reverse(trueNegatives.begin(), trueNegatives.end());
  std::reverse(falseNegatives.begin(), falseNegatives.end());
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