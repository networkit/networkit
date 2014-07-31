/*
 * BackboneBenchmark.cpp
 *
 *  Created on: 31.07.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#include "BackboneBenchmark.h"
#include "../Backbones.h"

namespace NetworKit {

BackboneBenchmark::BackboneBenchmark() {
	this->n = 250;
	INFO("n = " , this->n);
}

BackboneBenchmark::~BackboneBenchmark() {
}

TEST_F(BackboneBenchmark, completeGraphSimmelianBackboneParametric) {
	int64_t n = this->n;
	Aux::Timer runtime;

	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	runtime.start();
	SimmelianBackboneParametric algSBP(10, 5);
	Graph B = algSBP.calculate(G, EdgeAttribute());

	runtime.stop();
	INFO("[DONE] completeGraphSimmelianBackboneParametric (" , runtime.elapsed().count() , " ms)");
}

TEST_F(BackboneBenchmark, completeGraphSimmelianBackboneNonParametric) {
	int64_t n = this->n;
	Aux::Timer runtime;

	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	runtime.start();
	SimmelianBackboneNonParametric algSBNP(0.5);
	Graph B = algSBNP.calculate(G, EdgeAttribute());

	runtime.stop();
	INFO("[DONE] SimmelianBackboneNonParametric (" , runtime.elapsed().count() , " ms)");
}

TEST_F(BackboneBenchmark, completeGraphMultiscaleBackbone) {
	int64_t n = this->n;
	Aux::Timer runtime;

	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	runtime.start();
	MultiscaleBackbone algMB(0.5);
	Graph B = algMB.calculate(G, EdgeAttribute());

	runtime.stop();
	INFO("[DONE] MultiscaleBackbone (" , runtime.elapsed().count() , " ms)");
}

TEST_F(BackboneBenchmark, completeGraphLocalSimilarityBackbone) {
	int64_t n = this->n;
	Aux::Timer runtime;

	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	runtime.start();
	LocalSimilarityBackbone algLSB(0.5);
	Graph B = algLSB.calculate(G, EdgeAttribute());

	runtime.stop();
	INFO("[DONE] LocalSimilarityBackbone (" , runtime.elapsed().count() , " ms)");
}

TEST_F(BackboneBenchmark, completeGraphSimmelianMultiscaleBackbone) {
	int64_t n = this->n;
	Aux::Timer runtime;

	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);

	runtime.start();
	SimmelianMultiscaleBackbone algSMB(0.5);
	Graph B = algSMB.calculate(G, EdgeAttribute());

	runtime.stop();
	INFO("[DONE] SimmelianMultiscaleBackbone (" , runtime.elapsed().count() , " ms)");
}

} /* namespace NetworKit */

#endif /*NOGTEST */
