/*
 * SnapConverter.cpp
 *
 *  Created on: Sep 29, 2015
 *      Author: Michael
 */

#include "SnapConverter.h"
#include "../../io/SNAPGraphReader.h"
#include "../../io/METISGraphWriter.h"

#include <string>


using namespace std;

namespace NetworKit {

TEST(SnapConverter, convert) {
	vector<string> filePaths = {"Brightkite_edges.snap", "CA-AstroPh.snap", "CA-CondMat.snap", "CA-GrQc.snap", "CA-HepPh.snap", "CA-HepTh.snap",
								"Email-Enron.snap", "Gowalla_edges.snap", "as-skitter.snap", "com-amazon.ungraph.snap", "com-orkut.ungraph.snap",
								"com-youtube.ungraph.snap", "roadNet-CA.snap", "roadNet-PA.snap", "roadNet-TX.snap"};
	SNAPGraphReader snapReader;
	METISGraphWriter metisWriter;
	for (string filePath : filePaths) {
		filePath = "/Users/Michael/Downloads/SNAP/" + filePath;
		INFO("Reading ", filePath);
		Graph G = snapReader.read(filePath);
		string newFilePath = filePath.substr(0, filePath.find_last_not_of("."));
		metisWriter.write(G, newFilePath + ".graph");
	}
}

} /* namespace NetworKit */
