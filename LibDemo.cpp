#include <NetworKit/io/METISGraphReader.h>
#include <NetworKit/auxiliary/Log.h>
#include <NetworKit/community/PLM.h>
#include <iostream>


int main() {
	std::cout << "simple demonstration of NetworKit as a library\n";
	Aux::Log::Settings::setLogLevel(Aux::Log::LogLevel::info);
	NetworKit::METISGraphReader reader;
	NetworKit::Graph jazz = reader.read("./input/jazz.graph");
	NetworKit::PLM plm(jazz,true);
	plm.run();
	NetworKit::Partition communities = plm.getPartition();
	std::cout << communities.numberOfSubsets() << " communities have been found\n";	
	return 0;

}

