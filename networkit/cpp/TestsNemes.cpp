#include "TestsNemes.h"
#include "properties/EffectiveDiameter.h"
#include "centrality/Closeness.h"
#include "centrality/Betweenness.h"
#include "centrality/KPathCentrality.h"
#include "properties/ConnectedComponents.h"
#include "auxiliary/Random.h"
#include "io/METISGraphReader.h"

#include <math.h>
#include <iterator>
#include <stdlib.h>
#include <omp.h>
#include <map>

#include <numeric>
#include "graph/BFS.h"
#include "graph/Dijkstra.h"
#include <ctime>

namespace NetworKit {

double TestsNemes::effectiveDiameterExact() {
	using namespace std;

	cout << "Testing EffectiveDiameterExact:" << endl;

	double begin;
	double end;
	double result;

	//vector<string> testInstances = {"power", "PGPgiantcompo", "memplus", "as-22july06", "cs4", "citationCiteseer"};
	vector<string> testInstances = {"karate", "lesmis", "power"};
	for (auto testInstance : testInstances) {
		METISGraphReader reader;
		Graph g = reader.read("input/" + testInstance + ".graph");
		begin = omp_get_wtime();
		result = EffectiveDiameter::effectiveDiameterExact(g);
		end = omp_get_wtime();
		cout << testInstance << ", effective diameter: " << result << ", time elapsed: " << end - begin << endl;
	}
	return 0;
}

/*
 * do not call this method on the same graph multiple times in a single run
 * since the RNG and hence the result will be the same every time.
 */
double TestsNemes::effectiveDiameter(std::string graph, count k) {
	using namespace std;

	cout << "Testing EffectiveDiameter:" << endl;

	double begin;
	double end;
	double result;

	METISGraphReader reader;
	Graph g = reader.read("input/" + graph + ".graph");
	begin = omp_get_wtime();
	result = EffectiveDiameter::effectiveDiameter(g,0,k,7);
	end = omp_get_wtime();
	cout << graph << ", k=" << k << ", effective diameter: " << result << ", time elapsed: " << end - begin << endl;
	return 0;
}

double TestsNemes::hopPlot() {
	using namespace std;

	cout << "Testing Hop-Plot:" << endl;

	vector<pair<string, count>> testInstances= {
			pair<string, count>("power", 2048),
			pair<string, count>("citationCiteseer", 512),
			pair<string, count>("as-22july06", 1024),
			pair<string, count>("bcsstk30", 1024)
	};

	for (auto testInstance : testInstances) {
		METISGraphReader reader;
		Graph g = reader.read("input/" + testInstance.first + ".graph");
		map<count, double> hopplot = EffectiveDiameter::hopPlot(g,0,testInstance.second,7);
		for (count i = 0; i < hopplot.size(); i++) {
			cout << i << ": " << hopplot[i] << endl;
		}
	}
	return 0;
}

double TestsNemes::closeness() {
	using namespace std;

	double begin;
	double end;

	cout << "Testing Closeness:" << endl;

	vector<string> testInstances = {"power", "PGPgiantcompo", "astro-ph", "memplus", "cs4", "as-22july06", "bcsstk30", "cond-mat-2003", "rgg_n_2_15_s0", "fe_pwt", "cond-mat-2005", "delaunay_n16", "rgg_n_2_16_s0",	"kron_g500-simple-logn16", "luxembourg.osm", "delaunay_n17", "citationCiteseer", "coPapersCiteseer"};

	for (auto testInstance : testInstances) {
		METISGraphReader reader;
		Graph g = reader.read("input/" + testInstance + ".graph");
		begin = omp_get_wtime();
		Closeness close(g);
		close.run();
		end = omp_get_wtime();
		cout << testInstance << ", time elapsed: " << end - begin << endl;
	}
	return 0;
}

double TestsNemes::betweenness() {
	using namespace std;

	double begin;
	double end;

	cout << "Testing Betweenness:" << endl;

	vector<string> testInstances = {"power", "PGPgiantcompo", "memplus", "cs4", "as-22july06", "bcsstk30", "rgg_n_2_15_s0", "fe_pwt", "cond-mat-2005", "delaunay_n16", "rgg_n_2_16_s0", "luxembourg.osm", "delaunay_n17"};

	for (auto testInstance : testInstances) {
		METISGraphReader reader;
		Graph g = reader.read("input/" + testInstance + ".graph");
		begin = omp_get_wtime();
		Betweenness bet1(g);
		bet1.run();
		end = omp_get_wtime();
		cout << testInstance << " (only nodes, no indexing), time elapsed: " << end - begin << endl;

		begin = omp_get_wtime();
		g.indexEdges();
		Betweenness bet2(g,false,false);
		bet2.run();
		end = omp_get_wtime();
		cout << testInstance << " (only nodes, with indexing), time elapsed: " << end - begin << endl;

		begin = omp_get_wtime();
		g.indexEdges();
		Betweenness bet3(g,false,true);
		bet3.run();
		end = omp_get_wtime();
		cout << testInstance << " (nodes and edges, with indexing), time elapsed: " << end - begin << endl;
	}
	return 0;
}


double TestsNemes::kpath() {
	using namespace std;

	double begin;
	double end;
	METISGraphReader reader;

	cout << "Testing K-Path:" << endl;

	vector<double> kappas;
	kappas.push_back(5);
	kappas.push_back(10);
	kappas.push_back(15);
	kappas.push_back(20);

	vector<double> alphas;
	alphas.push_back(-0.1);
	alphas.push_back(0);
	alphas.push_back(0.2);
	alphas.push_back(0.5);

	vector<string> testInstances = {"PGPgiantcompo", "as-22july06"};

	for (auto testInstance : testInstances) {
		Graph g = reader.read("input/" + testInstance + ".graph");
		for (auto kappa : kappas) {
			for (auto alpha : alphas) {
				begin = omp_get_wtime();
				KPathCentrality kpath(g, alpha, kappa);
				kpath.run();
				end = omp_get_wtime();
				cout << testInstance << ", alpha:" << alpha << ", kappa:" << kappa << ", time elapsed: " << end - begin << endl;
			}
		}
		cout << endl;
	}

	vector<string> testInstances2 = {"PGPgiantcompo", "memplus", "cs4", "as-22july06", "cond-mat-2003", "rgg_n_2_15_s0", "fe_pwt", "cond-mat-2005", "delaunay_n16", "rgg_n_2_16_s0", "luxembourg.osm", "delaunay_n17", "citationCiteseer", "coPapersCiteseer"};

	vector<double> kappas2;
	kappas2.push_back(5);
	kappas2.push_back(10);
	kappas2.push_back(15);

	vector<double> alphas2;
	alphas2.push_back(0);
	alphas2.push_back(0.2);
	alphas2.push_back(0.5);

	for (auto testInstance : testInstances2) {
		Graph g = reader.read("input/" + testInstance + ".graph");

		//fixed alpha=0
		for (auto kappa : kappas2) {
			begin = omp_get_wtime();
			KPathCentrality kpath(g, 0, kappa);
			kpath.run();
			end = omp_get_wtime();
			cout << testInstance << ", alpha: 0, kappa:" << kappa << ", time elapsed: " << end - begin << endl;
		}
		cout << endl;

		//fixed alpha=0.2
		for (auto kappa : kappas2) {
			begin = omp_get_wtime();
			KPathCentrality kpath(g, 0.2, kappa);
			kpath.run();
			end = omp_get_wtime();
			cout << testInstance << ", alpha: 0.2, kappa:" << kappa << ", time elapsed: " << end - begin << endl;
		}
		cout << endl;

		//fixed kappa=10
		for (auto alpha : alphas2) {
			begin = omp_get_wtime();
			KPathCentrality kpath(g, alpha, 10);
			kpath.run();
			end = omp_get_wtime();
			cout << testInstance << ", alpha:" << alpha << ", kappa: 10, time elapsed: " << end - begin << endl;
		}
		cout << endl;

		//fixed kappa=15
		for (auto alpha : alphas2) {
			begin = omp_get_wtime();
			KPathCentrality kpath(g, alpha, 15);
			kpath.run();
			end = omp_get_wtime();
			cout << testInstance << ", alpha:" << alpha << ", kappa: 15, time elapsed: " << end - begin << endl;
		}
		cout << endl;
	}

	return 0;
}

//TODO
double TestsNemes::correlation() {
	using namespace std;

	cout << "Testing Correlation:" << endl;

	METISGraphReader reader;
	Graph g = reader.read("input/power.graph");
	cout << calculateCorrelation(g, 0, 0);

	return 0;
}

double TestsNemes::calculateCorrelation(const Graph& g, count measure1, count measure2) {
	using namespace std;

	double mean1 = 0;
	double mean2 = 0;
	double variance1 = 0;
	double variance2 = 0;
	double covariance = 0;
	double correlation = 0;

	//TODO
	KPathCentrality centrality1(g,0.0,10);
	centrality1.run();
	Betweenness centrality2(g,false);
	centrality2.run();

	for (count i = 0; i < g.numberOfNodes(); i++) {
		mean1 += centrality1.scores()[i];
		mean2 += centrality2.scores()[i];
	}

	mean1 /= g.numberOfNodes();
	mean2 /= g.numberOfNodes();

	for (count i = 0; i < g.numberOfNodes(); i++) {
		covariance += (centrality2.scores()[i] - mean2)*(centrality1.scores()[i] - mean1);
	}

	covariance /= g.numberOfNodes();

	for (count i = 0; i < g.numberOfNodes(); i++) {
		variance1 += pow((centrality1.scores()[i]-mean1),2);
		variance2 += pow((centrality2.scores()[i]-mean2),2);
	}

	variance1 /= g.numberOfNodes();
	variance2 /= g.numberOfNodes();

	correlation = covariance / ((pow(variance2, 0.5)) * (pow(variance1,0.5)));

	return correlation;
}

}
