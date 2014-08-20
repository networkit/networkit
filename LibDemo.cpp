#include <NetworKit/io/FastMETISGraphReader.h>
#include <NetworKit/auxiliary/Log.h>
#include <NetworKit/community/PLM.h>
#include <iostream>


int main() {
	std::cout << "simple demonstration of NetworKit as a library\n";
	Aux::Log::Settings::setLogLevel(Aux::Log::LogLevel::info);
	NetworKit::FastMETISGraphReader reader;
	NetworKit::Graph jazz = reader.read("./input/jazz.graph");
	NetworKit::PLM plm(true);
	NetworKit::Partition communities = plm.run(jazz);

	std::cout << communities.numberOfSubsets() << " communities have been found\n";	
	return 0;

}

