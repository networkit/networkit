/*
 * BackboneBenchmark.cpp
 *
 *  Created on: 31.07.2014
 *      Author: Gerd Lindner
 */

#ifndef NOGTEST

#include "BackboneBenchmark.h"
#include "../ChibaNishizekiTriangleCounter.h"
#include "../SimmelianJaccardAttributizer.h"
#include "../SimmelianOverlapAttributizer.h"
#include "../MultiscaleAttributizer.h"
#include "../LocalSimilarityAttributizer.h"
#include "../RandomAttributizer.h"
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

	SimmelianOverlapAttributizer attributizer(10);
	auto attribute = attributizer.getAttribute(G, std::vector<int>());

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

	SimmelianJaccardAttributizer attributizer;
	auto attribute = attributizer.getAttribute(G, std::vector<int>());

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

	MultiscaleAttributizer attributizer;
	std::vector<double> weight(G.upperEdgeIdBound());
	G.forEdges([&](node u, node v, edgeid eid) {
		weight[eid] = G.weight(u, v);
	});
	auto attribute = attributizer.getAttribute(G, weight);

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

	LocalSimilarityAttributizer attributizer;
	auto attribute = attributizer.getAttribute(G, std::vector<int>());

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
	ChibaNishizekiTriangleCounter triangleAttributizer;
	std::vector<int> triangles = triangleAttributizer.getAttribute(g, std::vector<int>(0));
	runtime.stop();
	std::cout << "[DONE] triangle counting " << runtime.elapsedTag() << std::endl;

	// --------- Multiscale
	std::cout << "[BEGIN] multiscale attribute: " << std::endl;
	runtime.start();
	MultiscaleAttributizer multiscaleAttributizer;
	std::vector<double> multiscale = multiscaleAttributizer.getAttribute(g, std::vector<double>(triangles.begin(), triangles.end()));
	runtime.stop();
	std::cout << "[DONE] multiscale attribute " << runtime.elapsedTag() << std::endl;

	std::cout << "[BEGIN] global filter (multiscale attribute): " << std::endl;
	runtime.start();
	GlobalThresholdFilter filter(0.5, false);
	Graph b = filter.calculate(g, multiscale);
	runtime.stop();
	std::cout << "[DONE] global filter (multiscale attribute) " << runtime.elapsedTag() << std::endl;

	// --------- Simmelian Backbone (Jaccard)
	std::cout << "[BEGIN] Simmelian Jaccard attribute: " << std::endl;
	runtime.start();
	SimmelianJaccardAttributizer jaccardAttributizer;
	std::vector<double> jaccard = jaccardAttributizer.getAttribute(g, triangles);
	runtime.stop();
	std::cout << "[DONE] Simmelian Jaccard attribute " << runtime.elapsedTag() << std::endl;

	std::cout << "[BEGIN] global filter (simmelian jaccard attribute): " << std::endl;
	runtime.start();
	filter = GlobalThresholdFilter(0.5, true);
	b = filter.calculate(g, jaccard);
	runtime.stop();
	std::cout << "[DONE] global filter (simmelian jaccard attribute) " << runtime.elapsedTag() << std::endl;

	// --------- Simmelian Backbone (Overlap)
	std::cout << "[BEGIN] Simmelian Overlap attribute: " << std::endl;
	runtime.start();
	SimmelianOverlapAttributizer overlapAttributizer(10);
	std::vector<double> overlap = overlapAttributizer.getAttribute(g, triangles);
	runtime.stop();
	std::cout << "[DONE] Simmelian Overlap attribute " << runtime.elapsedTag() << std::endl;

	std::cout << "[BEGIN] global filter (simmelian overlap attribute): " << std::endl;
	runtime.start();
	filter = GlobalThresholdFilter(5, true);
	b = filter.calculate(g, overlap);
	runtime.stop();
	std::cout << "[DONE] global filter (simmelian overlap attribute) " << runtime.elapsedTag() << std::endl;

	// --------- Local similarity Backbone
	std::cout << "[BEGIN] Local Similarity attribute: " << std::endl;
	runtime.start();
	LocalSimilarityAttributizer localSimAttributizer;
	std::vector<double> minExponent = localSimAttributizer.getAttribute(g, std::vector<int>());
	runtime.stop();
	std::cout << "[DONE] Local Similarity attribute " << runtime.elapsedTag() << std::endl;

	std::cout << "[BEGIN] global filter (local similarity attribute): " << std::endl;
	runtime.start();
	filter = GlobalThresholdFilter(0.37, true);
	b = filter.calculate(g, minExponent);
	runtime.stop();
	std::cout << "[DONE] global filter (local similarity attribute) " << runtime.elapsedTag() << std::endl;
}

} /* namespace NetworKit */

#endif /*NOGTEST */
