/*
 * PubWebGenerator.cpp
 *
 *  Created on: Apr 10, 2013
 *      Author: Henning
 */

#include "PubWebGenerator.h"

namespace NetworKit {

// FIXME: change to constants
#define THRSH_CENTER_ITER 10
#define ONE_DIM_SCALE 0.0f
#define MAX_RAND_VAL 1000
#define BENEFIT_INC 0.25f
#define PI 3.141f
#define MAX_DENSE_AREA_RADIUS 0.25f
#define MIN_MAX_DENSE_AREA_FACTOR 5 // ###
#define PART_CHANGE_FACTOR 0.0001f
#define INIT_LOAD 100.0f

struct circle {
	float x;
	float y;
	float rad;
};


void PubWebGenerator::moveNodeIntoUnitSquare(float& x, float& y) {
	auto move([&](float& z) {
		if (z > 1.0f) {
			z -= 1.0f;
		} else if (z < 0.0f) {
			z += 1.0f;
		}
	});

	move(x);
	move(y);
}

float PubWebGenerator::squaredDistanceInUnitTorus(float x1, float y1, float x2, float y2) {
	auto adjustForUnitTorus([&](float& z) {
		if (z > 0.5) {
			z = 1.0 - z;
		}
		else if (z < -0.5) {
			z = z + 1.0;
		}
	});

	float distx = x1 - x2;
	float disty = y1 - y2;
	adjustForUnitTorus(distx);
	adjustForUnitTorus(disty);

	return (distx * distx + disty * disty);
}

PubWebGenerator::PubWebGenerator(count numNodes, count numberOfDenseAreas,
		float neighborhoodRadius, count maxNumberOfNeighbors) :
		n(numNodes), numDenseAreas(numberOfDenseAreas), neighRad(
				neighborhoodRadius), maxNeigh(maxNumberOfNeighbors) {
	// TODO Auto-generated constructor stub

}

PubWebGenerator::~PubWebGenerator() {
	// TODO Auto-generated destructor stub
}

void PubWebGenerator::determineNeighbors(Graph& g) {
	auto isValidEdge([&](node u, node v, float squaredDistance) {
		return ((squaredDistance <= neighRad * neighRad)
				&& (g.degree(u) <= maxNeigh)
				&& (g.degree(v) <= maxNeigh));
	});

	g.forNodePairs([&](node u, node v) { // FIXME: quadratic loop!
		float sqrDist = squaredDistanceInUnitTorus(g.getCoordinate(u, 0), g.getCoordinate(u, 1), g.getCoordinate(v, 0), g.getCoordinate(v, 1));
		if (isValidEdge(u, v, sqrDist)) {
			g.addEdge(u, v);
		}
	});
}


Graph PubWebGenerator::generate() {
	// init
	Graph g(0);
	count dims = 2;
	std::vector<float> coordinates;
	Aux::RandomProbability randGen;

	// choose area sizes
	std::vector<circle> denseAreaXYR(numDenseAreas);
	for (index area = 0; area < numDenseAreas; ++area) {
		// anti-quadratic probability distribution
		float f = randGen.randomFloat() * MIN_MAX_DENSE_AREA_FACTOR + 1.0f;
		denseAreaXYR[area].rad = (MAX_DENSE_AREA_RADIUS * f * f)
				/ (MIN_MAX_DENSE_AREA_FACTOR * MIN_MAX_DENSE_AREA_FACTOR);
	}

	// compute number of nodes per cluster, each cluster has approx. same density
	float f = 0.0;
	for (index i = 0; i < numDenseAreas; ++i) {
		f += pow(denseAreaXYR[i].rad, 1.5);
	}
	f = ((float) n * ((float) numDenseAreas / ((float) numDenseAreas + 2.0f)))
			/ f; // TODO: better formula?

	std::vector<count> numPerArea(numDenseAreas);
	for (index i = 0; i < numDenseAreas; ++i) {
		numPerArea[i] = roundf(
				f * sqrt(denseAreaXYR[i].rad) * denseAreaXYR[i].rad);
	}

	// fill dense areas
	index t = 0;
	for (index area = 0; area < numDenseAreas; ++area) {
		// choose center randomly, ensure complete cluster is within (0,1) without modifications
		denseAreaXYR[area].x = randGen.randomFloat();
		denseAreaXYR[area].y = randGen.randomFloat();

		for (index j = 0; j < numPerArea[area]; ++j, ++t) {
			// compute random angle between [0, 2pi) and distance between [0, width/2]
			float angle = randGen.randomFloat() * 2.0 * PI;
			float dist = randGen.randomFloat() * denseAreaXYR[area].rad;

			// compute coordinates and adjust them
			float x = denseAreaXYR[area].x + cosf(angle) * dist;
			float y = denseAreaXYR[area].y + sinf(angle) * dist;
			moveNodeIntoUnitSquare(x, y);

			// create vertex with this coordinate
			g.addNode();
			coordinates.push_back(x);
			coordinates.push_back(y);
		}
	}

	// randomly spread remaining vertices over whole area
	for (; t < n; ++t) {
		g.addNode();
		coordinates.push_back(randGen.randomFloat());
		coordinates.push_back(randGen.randomFloat());
	}

	// insert coordinates into graph
	g.initCoordinates(dims);
	for (index v = 0; v < n; ++v) {
		g.setCoordinate(v, 0, coordinates[v * dims]);
		g.setCoordinate(v, 1, coordinates[v * dims + 1]);
	}

	// determine neighbors before adjusting the coordinates
	determineNeighbors(g);

	// adjust coordinates for postscript output
	g.forNodes([&](node u) {
		g.setCoordinate(u, 0, 1000.0 * g.getCoordinate(u, 0));
		g.setCoordinate(u, 1, 1000.0 * g.getCoordinate(u, 1));
	});

	return g;
}

} /* namespace NetworKit */

#if 0

// pLeave: probability that a node leaves the network
// fGrowth: network growth factor
// no update of numNeigh, nb, ew, degree necessary because this will be recalculated anyway by next DetermineNeighbors call
void UpdatePubwebGraph(Graph *const g, int *const part, int P,
		float** load, float** pload, float delta,
		float pLeave, float fGrowth, int vMax)
{
	int i, j, n, s, t, p, neigh;
	float f, angle, dist, x, y, internalDegree;

	// delete each vertex with probability pLeave
	for (i = 0; i < g->V; ++i) {
		if (FloatRand() <= pLeave) { // vertex will be deleted
			internalDegree = InternalDegree(g, i, part);

			// send remaining load to (internal) neighbors
			if (internalDegree > 0.0f) {
				for (p = 0; p < P; ++p) {
					for (j = 0; j < g->numNeigh[i]; ++j) {
						neigh = g->nb[i][j];
						if (part[i] == part[neigh]) {
							load[p][neigh] += load[p][i] * g->ew[i][j] / internalDegree;
							pload[p][neigh] += pload[p][i] * g->ew[i][j] / internalDegree;
						}
					}
				}
			}
			else { // no internal neighbors
				for (p = 0; p < P; ++p) {
					for (j = 0; j < g->numNeigh[i]; ++j) {
						load[p][g->nb[i][j]] += load[p][i] * g->ew[i][j] / g->degree[i];
						pload[p][g->nb[i][j]] += pload[p][i] * g->ew[i][j] / g->degree[i];
					}
				}
			}

			// make data structures contiguous again
			if (i < g->V - 1) {
				MemoryMove(&g->xy[2*i], &g->xy[2*(i+1)], float, 2 * (g->V - i - 1));
				MemoryMove(&g->bw[i], &g->bw[i+1], float, g->V - i - 1);
				for (p = 0; p < P; ++p) {
					MemoryMove(&load[p][i], &load[p][i+1], float, g->V - i - 1);
					MemoryMove(&pload[p][i], &pload[p][i+1], float, g->V - i - 1);
				}
			}
			--(g->V);
		}
	}

	// add nodes according to growth factor, randomize value in range [95 % .. 105 %)
	n = g->V * fGrowth * (0.95 + (FloatRand() / 10.0));
	if (g->V + n > vMax)
	n = vMax - g->V;

	// compute number of nodes per cluster, each cluster has approx. same density
	for (i = 0, f = 0.0; i < g->numDenseAreas; ++i)
	f += (sqrt(g->denseAreaXYR[3*i+2]) * g->denseAreaXYR[3*i+2]);
	f = ((float) n * ((float) g->numDenseAreas / ((float) g->numDenseAreas + 2.0f))) / f;// ###

	// add nodes to dense areas
	for (i = 0, t = g->V; i < g->numDenseAreas; ++i) {
		s = roundf(f * sqrt(g->denseAreaXYR[3*i+2]) * g->denseAreaXYR[3*i+2]);
		for (j = 0; j < s; ++j, ++t) {
			// berechne nun zufaellig den Winkel zwischen [0,2pi) und den Abstand zwischen [0,width/2]
			angle = FloatRand() * 2.0 * PI;
			dist = FloatRand() * g->denseAreaXYR[3*i+2];
			// fuege den Knoten mit der Koordinate in g ein
			x = g->denseAreaXYR[3*i] + cosf(angle) * dist;
			y = g->denseAreaXYR[3*i+1] + sinf(angle) * dist;

			AdjustWrapAroundNode(x, y);
			g->xy[t*2] = x;
			g->xy[t*2 + 1] = y;

			// initial load
			for (p = 0; p < P; ++p) {
				load[p][t] = 0.0f;
				pload[p][t] = 0.0f; // delta;
			}

			// initial partition
			part[t] = rand() % P;
		}
	}

	// randomly spread remaining vertices over whole area
	for (; t < g->V + n; ++t) {
		g->xy[t*2] = FloatRand();
		g->xy[t*2 + 1] = FloatRand();
	}

	g->V += n;
}

// factor: factor of total nodes to move
void MovePubweb(Graph *const g, int *const part, int P,
		float *const *const load, float *const *const pload,
		float factor)
{
	int i, j, n, s, t, p, e, neigh;
	float f, angle, dist, x, y;

	// choose nodes according to factor, randomize value in range [95 % .. 105 %)
	n = g->V * factor * (0.95 + (FloatRand() / 10.0));

	// compute number of nodes per cluster, each cluster has approx. same density
	for (i = 0, f = 0.0; i < g->numDenseAreas; ++i)
	f += (sqrt(g->denseAreaXYR[3*i+2]) * g->denseAreaXYR[3*i+2]);
	f = ((float) n * ((float) g->numDenseAreas / ((float) g->numDenseAreas + 2.0f))) / f;// ###

	// move nodes to dense areas
	for (i = 0; i < g->numDenseAreas; ++i) {
		// determine how many nodes to move to this dense area
		s = roundf(f * sqrt(g->denseAreaXYR[3*i+2]) * g->denseAreaXYR[3*i+2]);
		for (j = 0; j < s; ++j, --n) {
			// choose node to move randomly
			t = IntRand(g->V);

			// berechne nun zufaellig den Winkel zwischen [0,2pi) und den Abstand zwischen [0,width/2]
			angle = FloatRand() * 2.0 * PI;
			dist = FloatRand() * g->denseAreaXYR[3*i+2];
			// berechne Koordinaten
			x = g->denseAreaXYR[3*i] + cosf(angle) * dist;
			y = g->denseAreaXYR[3*i+1] + sinf(angle) * dist;
			AdjustWrapAroundNode(x, y);

			// send remaining load to (internal) neighbors
			float internalDegree = InternalDegree(g, t, part);
			if (internalDegree > 0.0f) {
				for (p = 0; p < P; ++p) {
					for (e = 0; e < g->numNeigh[t]; ++e) {
						neigh = g->nb[t][e];
						if (part[t] == part[neigh]) {
							load[p][neigh] += load[p][t] * g->ew[t][e] / internalDegree;
							pload[p][neigh] += pload[p][t] * g->ew[t][e] / internalDegree;
						}
					}
				}
			}
			else { // no internal neighbors
				for (p = 0; p < P; ++p) {
					for (e = 0; e < g->numNeigh[t]; ++e) {
						neigh = g->nb[t][e];
						load[p][neigh] += load[p][t] * g->ew[t][e] / g->degree[t];
						pload[p][neigh] += pload[p][t] * g->ew[t][e] / g->degree[t];
					}
				}
			}

			// Knoten eintragen
			g->xy[t*2] = x;
			g->xy[t*2 + 1] = y;
			part[t] = rand() % P;
			for (int p = 0; p < P; ++p) {
				load[p][t] = 0.0f;
				pload[p][t] = 0.0f;
			}
		}
	}

	// randomly spread remaining vertices over whole area
	for (; n > 0; --n) {
		// choose node to move randomly
		t = IntRand(g->V);

		// send remaining load to (internal) neighbors
		float internalDegree = InternalDegree(g, t, part);
		if (internalDegree > 0.0f) {
			for (p = 0; p < P; ++p) {
				for (e = 0; e < g->numNeigh[t]; ++e) {
					neigh = g->nb[t][e];
					if (part[t] == part[neigh]) {
						load[p][neigh] += load[p][t] * g->ew[t][e] / internalDegree;
						pload[p][neigh] += pload[p][t] * g->ew[t][e] / internalDegree;
					}
				}
			}
		}
		else { // no internal neighbors
			for (p = 0; p < P; ++p) {
				for (e = 0; e < g->numNeigh[t]; ++e) {
					neigh = g->nb[t][e];
					load[p][neigh] += load[p][t] * g->ew[t][e] / g->degree[t];
					pload[p][neigh] += pload[p][t] * g->ew[t][e] / g->degree[t];
				}
			}
		}

		// Knoten eintragen
		g->xy[t*2] = FloatRand();
		g->xy[t*2 + 1] = FloatRand();
		part[t] = rand() % P;

		for (int p = 0; p < P; ++p) {
			load[p][t] = 0.0f;
			pload[p][t] = 0.0f;
		}

	}
}

#endif

