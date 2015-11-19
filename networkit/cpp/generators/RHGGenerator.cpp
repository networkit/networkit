#include <cstdlib>
#include <random>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include <algorithm>
#include "../graph/GraphBuilder.h"
#include "../graph/Graph.h"
#include "RHGGenerator.h"
#include "../auxiliary/Random.h"
#include <set>

using std::vector;

namespace NetworKit {

	RHGGenerator::RHGGenerator() {
		stretch = 1;
		alpha = 1;
		factor = 1;
		nodeCount = 10000;
	}

	RHGGenerator::RHGGenerator(count n) {
		nodeCount = n;
		alpha = 1;
		factor = 1;
		stretch = 1;
	}
	/**
	* Construct a generator for n nodes and m edges
	*/

	RHGGenerator::RHGGenerator(count n, double avgDegree, double plexp) {
		nodeCount = n;
		if (plexp < 2) throw std::runtime_error("Exponent of power-law degree distribution must be >= 2");
		alpha = (plexp-1)/2;
		double R = HyperbolicSpace::hyperbolicAreaToRadius(n);
		double targetR = HyperbolicSpace::getTargetRadius(n, n*avgDegree/2, alpha, 0);// 2*log(8*n / (M_PI*(m/n)*2));
		stretch = targetR / R;
		factor = 1;
		threadtimers.resize(omp_get_max_threads());
	}

	Graph RHGGenerator::generate() {
		return generate(nodeCount, factor, alpha, stretch);
	}

	Graph RHGGenerator::generate(count n, double distanceFactor, double alpha, double stretchradius) {
		double R = stretchradius*HyperbolicSpace::hyperbolicAreaToRadius(n);
		assert(R > 0);
		vector<double> angles(n);
		vector<double> radii(n);

		//sample points randomly
		fillPoints(angles, radii, stretchradius, alpha);
		vector<index> permutation(n);

		index p = 0;
		std::generate(permutation.begin(), permutation.end(), [&p](){return p++;});

		//can probably be parallelized easily, but doesn't bring much benefit
		std::sort(permutation.begin(), permutation.end(), [&angles,&radii](index i, index j){return angles[i] < angles[j] || (angles[i] == angles[j] && radii[i] < radii[j]);});

		vector<double> anglecopy(n);
		vector<double> radiicopy(n);

		#pragma omp parallel for
		for (index j = 0; j < n; j++) {
			anglecopy[j] = angles[permutation[j]];
			radiicopy[j] = radii[permutation[j]];
		}
		INFO("Generated Points");
		return generate(anglecopy, radiicopy, R*distanceFactor);
	}

	Graph RHGGenerator::generate(const vector<double> &angles, const vector<double> &radii, double thresholdDistance) {
		Aux::Timer timer;
		timer.start();
		index n = angles.size();
		assert(radii.size() == n);
		//1.Generate bandRadius'
		vector<double> bandRadius;
		getBandRadius(n, bandRadius, thresholdDistance);
		//2. Initialize empty bands
		vector<vector<Point2D<double>>> bands(bandRadius.size() - 1);
		//3. Put points to bands
		for(index i = 0; i < n; i++){
			for(index j = 0; j < bands.size(); j++){
				if(radii[i] >= bandRadius[j] && radii[i] <= bandRadius[j+1]){
					bands[j].push_back(Point2D<double>(angles[i], radii[i], i));
					break;
				}
			}
		}
		//the bands are already sorted since we sorted angle&radii before
		timer.stop();
		INFO("Filled bands, took ", timer.elapsedMilliseconds(), " milliseconds.");
		return generate(angles, radii, bands, bandRadius, thresholdDistance);
	}


	Graph RHGGenerator::generate(const vector<double> &angles, const vector<double> &radii, const vector<vector<Point2D<double>>> &bands, const vector<double> &bandRadius,
		double thresholdDistance) {

			index n = angles.size();
			assert(radii.size() == n);
			Aux::Timer timer;
			vector<double> empty;
			GraphBuilder result(n, false, false);
			timer.start();

			//1.Extract band angles to use it later without increasing complexity, can create a band class to handle this more elegantly
			vector<vector<double>> bandAngles(bands.size());
			for(index i=0; i < bands.size(); i++){
				for(index j=0; j < bands[i].size(); j++)
				bandAngles[i].push_back(bands[i][j].getX());
			}

			//2.Insert edges
			#pragma omp parallel
			{
				index id = omp_get_thread_num();
				threadtimers[id].start();
				#pragma omp for schedule(guided) nowait
				for (index i = 0; i < n; i++) {
					count expectedDegree = (4/M_PI)*n*exp(-(radii[i])/2);//TODO: adapt for alpha!=1
					vector<index> near;
					near.reserve(expectedDegree*1.1);
					Point2D<double> pointV(angles[i], radii[i], i);
					for(index j = 0; j < bands.size(); j++){
						if(bandRadius[j+1] > radii[i]){
							double minTheta, maxTheta;
							std::tie (minTheta, maxTheta) = getMinMaxTheta(angles[i], radii[i], bandRadius[j], thresholdDistance);
							vector<Point2D<double>> slab;
							getPointsWithinAngles(minTheta, maxTheta, bands[j], bandAngles[j], slab);
							for(index w = 0; w < slab.size(); w++){
								if(getHyperbolicDistance(pointV, slab[w]) <= thresholdDistance){
									if(slab[w].getIndice() != i){
										near.push_back(slab[w].getIndice());
									}
								}
							}
						}
					}
					if (directSwap) {
						auto newend = std::remove(near.begin(), near.end(), i); //no self loops!
						if (newend != near.end()) {
							assert(newend+1 == near.end());
							assert(*(newend)==i);
							near.pop_back();//std::remove doesn't remove element but swaps it to the end
						}
						result.swapNeighborhood(i, near, empty, false);
					} else {
						for (index j : near) {
							if (j >= n) ERROR("Node ", j, " prospective neighbour of ", i, " does not actually exist. Oops.");
							if(radii[j] > radii[i] || (radii[j] == radii[i] && angles[j] < angles[i]))
								result.addHalfEdge(i,j);
						}
					}
				}
				threadtimers[id].stop();
			}
			timer.stop();
			INFO("Generating Edges took ", timer.elapsedMilliseconds(), " milliseconds.");
			return result.toGraph(!directSwap, true);
		}
	}
