/*
 * MaxentStress.cpp
 *
 *  Created on: 22.01.2014
 *      Author: Henning
 */

#include "MaxentStress.h"

namespace NetworKit {

MaxentStress::MaxentStress(Point<float> bottomLeft, Point<float> topRight):
				Layouter(bottomLeft, topRight)
{

}

MaxentStress::~MaxentStress() {

}

void MaxentStress::draw(Graph& G) {


	// Paper: q = 0 for many graphs, q = 0.8 for graphs with many degree-1 nodes
	// alpha: initially 1, then in each iteration alpha := 0.3 * alpha

	float q = 0;
	float alpha = 1.0; // FIXME: make alpha dependent on iteration

	G.forNodes([&](node u) {
		Point<float> uPoint(G.getCoordinate(u, 0), G.getCoordinate(u, 1));
		Point<float> attractiveForce(0.0, 0.0);
		Point<float> repulsiveForce(0.0, 0.0);

		G.forNodes([&](node v) {
			if (u < v) { // only unordered pairs
				Point<float> vPoint(G.getCoordinate(v, 0), G.getCoordinate(v, 1));
				float diffX = uPoint[0] - vPoint[0];
				float diffY = uPoint[1] - vPoint[1];
				Point<float> diffVec(diffX, diffY);
				float len = diffVec.length();

				if (G.hasEdge(u, v)) {
					// sum over all node pairs in S
					// $\frac{1}{\rho_i} \sum_{i,j \in S} w_{ij} * (x_j + d_{ij} \frac{x_i - x_j}{\Vert x_i - x_j \Vert})$

					float dist = 1.0; // FIXME: dist = ... (algebraic distance between u and v)
					diffVec.scale(dist / len);
					diffVec += vPoint;
					attractiveForce += diffVec.scale(G.weight(u, v));
				}
				else {
					// traverse remaining vertices not in neighborhood for repulsive forces
					// sum over all node pairs not in S
					// $\frac{\alpha}{\rho_i} \sum_{i,j \notin S} w_{ij} \frac{x_i - x_j}{\Vert x_i - x_j \Vert^{q+2}}$
					// TODO: approximation (e.g. Barnes-Hut)

					float denom = 1.0 / pow(len, q+2);
					diffVec.scale(denom);
					repulsiveForce += diffVec.scale(G.weight(u, v));
				}
			}
		});

		// apply forces to node u
		float rhoInv = 1.0 / (float) G.weightedDegree(u);
		attractiveForce.scale(rhoInv);
		repulsiveForce.scale(alpha * rhoInv);
	});
}

} /* namespace NetworKit */
