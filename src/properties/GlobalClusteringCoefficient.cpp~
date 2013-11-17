/*
 * ClusteringCoefficient.cpp
 *
 *  Created on: 08.04.2013
 *      Author: cls
 */

#include "GlobalClusteringCoefficient.h"

namespace NetworKit {

GlobalClusteringCoefficient::GlobalClusteringCoefficient() {
}

GlobalClusteringCoefficient::~GlobalClusteringCoefficient() {
}

float GlobalClusteringCoefficient::run(Graph& G) {

	int numerator=0; 
	int denominator=0;
	int numerator_vorher=0;
	
	G.forNodes([&](node u)
	{
		if(G.degree(u)>=2) {

			numerator_vorher=numerator;
		
			denominator+=(G.degree(u))*(G.degree(u)-1)/2;
			G.forEdgesOf(u,[&](node u, node w)
				{	
					
					G.forEdgesOf(w,[&](node w, node v)
					{
						if(v!=u)
						{
							if(G.hasEdge(v,u)==1)
								numerator++;// Jedes Dreieck wird von drei unterschiedlichen Knoten gezählt, 										    // deshalb Faktor 3 weggelassen.
						}
					});
				});
		numerator=(numerator-numerator_vorher)/2+numerator_vorher;// Für jeden Knoten wird jedes Dreieck doppelt gezählt.
									  // Das muss rausgerechnet werden.
		}
	});
	float coefficient = (float) numerator / denominator;
	return coefficient;
}
}	





























