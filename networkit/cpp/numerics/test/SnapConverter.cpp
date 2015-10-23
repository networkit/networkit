/*
 * SnapConverter.cpp
 *
 *  Created on: Sep 29, 2015
 *      Author: Michael
 */

#include "SnapConverter.h"
#include "../../io/EdgeListReader.h"
#include "../../io/METISGraphWriter.h"

#include <string>


using namespace std;

namespace NetworKit {

TEST(SnapConverter, convert) {
	vector<string> filePaths = {"Brightkite_edges.snap", "CA-AstroPh.snap", "CA-CondMat.snap", "CA-GrQc.snap", "CA-HepPh.snap", "CA-HepTh.snap",
								"Email-Enron.snap", "Gowalla_edges.snap", "as-skitter.snap", "com-amazon.ungraph.snap", "com-orkut.ungraph.snap",
								"com-youtube.ungraph.snap", "roadNet-CA.snap", "roadNet-PA.snap", "roadNet-TX.snap"};

	METISGraphWriter metisWriter;
	for (string filePath : filePaths) {
		EdgeListReader snapReader('\t',0,"#",false);
		filePath = "/Users/Michael/Downloads/SNAP/" + filePath;
		INFO("Reading ", filePath);
		Graph G = snapReader.read(filePath);
		INFO("n = ", G.numberOfNodes(), " m = ", G.numberOfEdges());
		string newFilePath = filePath.substr(0, filePath.find_last_of("."));
		metisWriter.write(G, newFilePath + ".graph");
	}
}

} /* namespace NetworKit */
