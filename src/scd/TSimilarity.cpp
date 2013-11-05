/*
 * TSimilarity.cpp
 *
 *  Created on: 14.09.2013
 *      Author: cls, Yassine Marrakchi
 */

#include "TSimilarity.h"

namespace NetworKit {

TSimilarity::TSimilarity(const Graph& G): G(&G) {
}

TSimilarity::~TSimilarity() {
}


TDummySimilarity::TDummySimilarity(const Graph& G) : TSimilarity(G) {
}

TDummySimilarity::~TDummySimilarity() {
}

double TDummySimilarity::getValue(node u, node v) {
	return 1;
}

TNodesSimilarity::TNodesSimilarity(const Graph& G) : TSimilarity(G) {
}

TNodesSimilarity::~TNodesSimilarity() {
}

double TNodesSimilarity::getValue(node u, node v) {
	double similarity = 1 / log(G->degree(u));
	this->G->forNeighborsOf(u, [&](node x) {
		if (G->hasEdge(v, x)) {
			similarity = similarity + (1 / log(G->degree(x)));
		}
	});
	if (!(G->hasEdge(v, v))) {
		similarity = similarity + 1 / log(G->degree(v));
	}
	return similarity;
}
} /* namespace NetworKit */
