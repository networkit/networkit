/*
 * BackboneBenchmark.cpp
 *
 *  Created on: 31.07.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#include "BackboneBenchmark.h"
#include "../../edgeattributes/ChibaNishizekiTriangleCounter.h"
#include "../../edgeattributes/TriangleCounter.h"
#include "../SimmelianJaccardAttributizer.h"
#include "../SimmelianOverlapAttributizer.h"
#include "../MultiscaleAttributizer.h"
#include "../LocalSimilarityAttributizer.h"
#include "../RandomEdgeAttributizer.h"
#include "../GlobalThresholdFilter.h"
#include "../../io/METISGraphReader.h"
#include "../../auxiliary/Log.h"

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
	G.indexEdges();

	runtime.start();

	ChibaNishizekiTriangleCounter counter(G);
	std::vector<count> counts = counter.getAttribute();

	SimmelianOverlapAttributizer attributizer(G, counts, 10);
	auto attribute = attributizer.getAttribute();

	runtime.stop();
	INFO("[DONE] completeGraphSimmelianBackboneParametric (" , runtime.elapsed().count() , " ms)");
}

TEST_F(BackboneBenchmark, completeGraphSimmelianBackboneNonParametric) {
	int64_t n = this->n;
	Aux::Timer runtime;

	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);
	G.indexEdges();

	runtime.start();

	ChibaNishizekiTriangleCounter counter(G);
	std::vector<count> counts = counter.getAttribute();

	SimmelianJaccardAttributizer attributizer(G, counts);
	auto attribute = attributizer.getAttribute();

	runtime.stop();
	INFO("[DONE] SimmelianBackboneNonParametric (" , runtime.elapsed().count() , " ms)");
}

TEST_F(BackboneBenchmark, completeGraphMultiscaleBackbone) {
	int64_t n = this->n;
	Aux::Timer runtime;

	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);
	G.indexEdges();

	runtime.start();

	std::vector<double> weight(G.upperEdgeIdBound());
	G.forEdges([&](node u, node v, edgeid eid) {
		weight[eid] = G.weight(u, v);
	});

	MultiscaleAttributizer attributizer(G, weight);
	auto attribute = attributizer.getAttribute();

	runtime.stop();
	INFO("[DONE] MultiscaleBackbone (" , runtime.elapsed().count() , " ms)");
}

TEST_F(BackboneBenchmark, completeGraphLocalSimilarityBackbone) {
	int64_t n = this->n;
	Aux::Timer runtime;

	GraphGenerator graphGen;
	Graph G = graphGen.makeCompleteGraph(n);
	G.indexEdges();

	runtime.start();

	ChibaNishizekiTriangleCounter counter(G);
	std::vector<count> triangles = counter.getAttribute();

	LocalSimilarityAttributizer attributizer(G, triangles);
	auto attribute = attributizer.getAttribute();

	runtime.stop();
	INFO("[DONE] LocalSimilarityBackbone (" , runtime.elapsed().count() , " ms)");
}

TEST_F(BackboneBenchmark, backboneBenchmarkGraphFile) {
	std::string path = "";

	std::cout << "[INPUT] .graph file path >" << std::endl;
	std::getline(std::cin, path);
	Aux::Timer runtime;

	// --------- IO
	std::cout << "[BEGIN] reading graph: " << path << std::endl;
	runtime.start();
	METISGraphReader reader;
	Graph g = reader.read(path);
	runtime.stop();
	std::cout << "[DONE] reading graph " << runtime.elapsedTag() << std::endl;

	// --------- Edge indexing
	std::cout << "[BEGIN] edge indexing: " << std::endl;
	runtime.start();
	g.indexEdges();
	runtime.stop();
	std::cout << "[DONE] edge indexing " << runtime.elapsedTag() << std::endl;

	// --------- Triangle counting
	std::cout << "[BEGIN] triangle counting: " << std::endl;
	runtime.start();
	ChibaNishizekiTriangleCounter oldTriangleAttributizer(g);
	std::vector<count> oldTriangles = oldTriangleAttributizer.getAttribute();
	runtime.stop();
	std::cout << "[DONE] Chiba Nishizeki triangle counting " << runtime.elapsedTag() << std::endl;

	// --------- Triangle counting
	std::cout << "[BEGIN] triangle counting: " << std::endl;
	runtime.start();
	TriangleCounter triangleAttributizer(g);
	std::vector<count> triangles = triangleAttributizer.getAttribute();
	runtime.stop();
	std::cout << "[DONE] triangle counting " << runtime.elapsedTag() << std::endl;

	// --------- Multiscale
	std::cout << "[BEGIN] multiscale attribute: " << std::endl;
	runtime.start();
	MultiscaleAttributizer multiscaleAttributizer(g, std::vector<double>(triangles.begin(), triangles.end()));
	std::vector<double> multiscale = multiscaleAttributizer.getAttribute();
	runtime.stop();
	std::cout << "[DONE] multiscale attribute " << runtime.elapsedTag() << std::endl;

	std::cout << "[BEGIN] global filter (multiscale attribute): " << std::endl;
	runtime.start();
	GlobalThresholdFilter filter(g, multiscale, 0.5, false);
	Graph b = filter.calculate();
	runtime.stop();
	std::cout << "[DONE] global filter (multiscale attribute) " << runtime.elapsedTag() << std::endl;

	// --------- Simmelian Backbone (Jaccard)
	std::cout << "[BEGIN] Simmelian Jaccard attribute: " << std::endl;
	runtime.start();
	SimmelianJaccardAttributizer jaccardAttributizer(g, triangles);
	std::vector<double> jaccard = jaccardAttributizer.getAttribute();
	runtime.stop();
	std::cout << "[DONE] Simmelian Jaccard attribute " << runtime.elapsedTag() << std::endl;

	std::cout << "[BEGIN] global filter (simmelian jaccard attribute): " << std::endl;
	runtime.start();
	GlobalThresholdFilter filter2(g, jaccard, 0.5, true);
	b = filter2.calculate();
	runtime.stop();
	std::cout << "[DONE] global filter (simmelian jaccard attribute) " << runtime.elapsedTag() << std::endl;

	// --------- Simmelian Backbone (Overlap)
	std::cout << "[BEGIN] Simmelian Overlap attribute: " << std::endl;
	runtime.start();
	SimmelianOverlapAttributizer overlapAttributizer(g, triangles, 10);
	std::vector<double> overlap = overlapAttributizer.getAttribute();
	runtime.stop();
	std::cout << "[DONE] Simmelian Overlap attribute " << runtime.elapsedTag() << std::endl;

	std::cout << "[BEGIN] global filter (simmelian overlap attribute): " << std::endl;
	runtime.start();
	GlobalThresholdFilter filter3(g, overlap, 5, true);
	b = filter3.calculate();
	runtime.stop();
	std::cout << "[DONE] global filter (simmelian overlap attribute) " << runtime.elapsedTag() << std::endl;

	// --------- Local similarity Backbone
	std::cout << "[BEGIN] Local Similarity attribute: " << std::endl;
	runtime.start();
	LocalSimilarityAttributizer localSimAttributizer(g, triangles);
	std::vector<double> minExponent = localSimAttributizer.getAttribute();
	runtime.stop();
	std::cout << "[DONE] Local Similarity attribute " << runtime.elapsedTag() << std::endl;

	std::cout << "[BEGIN] global filter (local similarity attribute): " << std::endl;
	runtime.start();
	GlobalThresholdFilter filter4 (g, minExponent, 0.37, true);
	b = filter4.calculate();
	runtime.stop();
	std::cout << "[DONE] global filter (local similarity attribute) " << runtime.elapsedTag() << std::endl;
}

} /* namespace NetworKit */

#endif /*NOGTEST */
