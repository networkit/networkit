/*
 * MaxentStressGTest.cpp
 *
 *  Created on: Apr 19, 2016
 *      Author: Michael
 */

#include <gtest/gtest.h>
#include <vector>
#include <string>

#include "../../graph/Graph.h"
#include "../Point.h"

#include "../../io/METISGraphReader.h"
#include "../../io/METISGraphWriter.h"
#include "../../graph/Graph.h"
#include "../../components/ConnectedComponents.h"
#include "../MaxentStress.h"

#include "../../numerics/LAMG/Lamg.h"
#include "../../numerics/ConjugateGradient.h"
#include "../../numerics/Preconditioner/IdentityPreconditioner.h"
#include "../../numerics/Preconditioner/DiagonalPreconditioner.h"

#include "../../community/PLM.h"
#include "../../clique/MaxClique.h"

#include "../../sparsification/LocalDegreeScore.h"
#include "../../sparsification/RandomEdgeScore.h"
#include "../../sparsification/GlobalThresholdFilter.h"

#include "../../auxiliary/Timer.h"
#include "../../auxiliary/Random.h"

#include <iostream>
#include <unordered_map>
#include <random>

#include "../PivotMDS.h"

#include <stdio.h>

namespace NetworKit {

class MaxentStressGTest : public testing::Test {};

TEST_F(MaxentStressGTest, benchMaxentStressCoordinatesLAMG) {

    std::vector<std::string> graphFiles = {"input/airfoil1.graph"};
	METISGraphReader reader;

	for (const std::string& graphFile : graphFiles) {
		Graph graph = reader.read(graphFile);

		double runtime = 0;
		double fullStress = 0;
		double maxentStress = 0;
		Aux::Random::setSeed(Aux::Random::integer(), false);

		Aux::Timer t;
		t.start();
		PivotMDS pivotMds(graph, 2, 30);
		pivotMds.run();
        MaxentStress maxentStressAlgo(graph, 2, pivotMds.getCoordinates(), 1, 0.001, MaxentStress::LinearSolverType::LAMG);
		maxentStressAlgo.run();
		t.stop();

		runtime = t.elapsedMicroseconds();
		if (graph.numberOfNodes() < 1e5) {
			maxentStressAlgo.scaleLayout();
			fullStress = maxentStressAlgo.fullStressMeasure();
			maxentStress = maxentStressAlgo.maxentMeasure();
		}


		runtime /= 1000;

		INFO(graphFile, "\t", maxentStress, "\t", fullStress, "\t", runtime);
	}
}

TEST_F(MaxentStressGTest, benchMaxentStressConjGradIdPrecAlgebraicDistance) {
    std::vector<std::string> graphFiles = {"input/airfoil1.graph"};
    METISGraphReader reader;

    for (const std::string& graphFile : graphFiles) {
        Graph graph = reader.read(graphFile);

        double runtime = 0;
        double fullStress = 0;
        double maxentStress = 0;
        Aux::Random::setSeed(Aux::Random::integer(), false);

        Aux::Timer t;
        t.start();
        MaxentStress maxentStressAlgo(graph, 2, 1, 0.001, MaxentStress::LinearSolverType::CONJUGATE_GRADIENT_IDENTITY_PRECONDITIONER, true, MaxentStress::GraphDistance::ALGEBRAIC_DISTANCE);
        maxentStressAlgo.run();
        t.stop();

        runtime = t.elapsedMicroseconds();
        if (graph.numberOfNodes() < 1e5) {
            maxentStressAlgo.scaleLayout();
            fullStress = maxentStressAlgo.fullStressMeasure();
            maxentStress = maxentStressAlgo.maxentMeasure();
        }


        runtime /= 1000;

    INFO(graphFile, "\t", maxentStress, "\t", fullStress, "\t", runtime);
    }
}

TEST_F(MaxentStressGTest, benchMaxentStressConjGradDiagPrecond) {
    std::vector<std::string> graphFiles = {"input/airfoil1.graph"};
    METISGraphReader reader;

    for (const std::string& graphFile : graphFiles) {
        Graph graph = reader.read(graphFile);

        double runtime = 0;
        double fullStress = 0;
        double maxentStress = 0;
        Aux::Random::setSeed(Aux::Random::integer(), false);

        Aux::Timer t;
        t.start();
        MaxentStress maxentStressAlgo(graph, 2, 1, 0.001, MaxentStress::LinearSolverType::CONJUGATE_GRADIENT_DIAGONAL_PRECONDITIONER);
        maxentStressAlgo.run();
        t.stop();

        runtime = t.elapsedMicroseconds();
        if (graph.numberOfNodes() < 1e5) {
            maxentStressAlgo.scaleLayout();
            fullStress = maxentStressAlgo.fullStressMeasure();
            maxentStress = maxentStressAlgo.maxentMeasure();
        }


        runtime /= 1000;

        INFO(graphFile, "\t", maxentStress, "\t", fullStress, "\t", runtime);
    }
}

TEST_F(MaxentStressGTest, benchMaxentStressCoordConjGradIdPrecond) {
    std::vector<std::string> graphFiles = {"input/airfoil1.graph"};
    METISGraphReader reader;

    for (const std::string& graphFile : graphFiles) {
        Graph graph = reader.read(graphFile);

        double runtime = 0;
        double fullStress = 0;
        double maxentStress = 0;
        Aux::Random::setSeed(Aux::Random::integer(), false);

        Aux::Timer t;
        t.start();
        PivotMDS pivotMds(graph, 2, 30);
        pivotMds.run();
        MaxentStress maxentStressAlgo(graph, 2, pivotMds.getCoordinates(), 1, 0.001, MaxentStress::LinearSolverType::CONJUGATE_GRADIENT_IDENTITY_PRECONDITIONER, false, MaxentStress::GraphDistance::ALGEBRAIC_DISTANCE);
        maxentStressAlgo.run();
        t.stop();

        runtime = t.elapsedMicroseconds();
        if (graph.numberOfNodes() < 1e5) {
            maxentStressAlgo.scaleLayout();
            fullStress = maxentStressAlgo.fullStressMeasure();
            maxentStress = maxentStressAlgo.maxentMeasure();
        }


        runtime /= 1000;

        INFO(graphFile, "\t", maxentStress, "\t", fullStress, "\t", runtime);
    }
}

} /* namespace NetworKit */
