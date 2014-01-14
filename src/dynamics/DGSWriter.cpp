/*
 * DGSWriter.cpp
 *
 *  Created on: 14.01.2014
 *      Author: cls
 */

#include "DGSWriter.h"

#include <fstream>

namespace NetworKit {

void DGSWriter::write(std::vector<GraphEvent>& stream, std::string path) {
	std::ofstream out(path);

	for (GraphEvent ev : stream) {
		switch (ev.type) {
			case GraphEvent::NODE_ADDITION : {
				out << "an " << ev.u << std::endl; 
				break;
			}
			case GraphEvent::NODE_REMOVAL : {
				out << "dn " << ev.u << std::endl;
				break;
			}
			case GraphEvent::EDGE_ADDITION : {
				out << "ae " << ev.u << "-" << ev.v << " " << ev.u << " " << ev.v << " weight=" << ev.w << std::endl;
				break;
			}
			case GraphEvent::EDGE_REMOVAL : {
				out << "de " << ev.u << "-" << ev.v << std::endl;
				break;
			}
			case GraphEvent::EDGE_WEIGHT_UPDATE : {
				out << "ce " << ev.u << "-" << ev.v  << " weight=" << ev.w << std::endl;
				break;
			}
			case GraphEvent::TIME_STEP : {
				out << "st" << std::endl;
				break;
			}
			default: {
				throw std::runtime_error("unknown event type");
			}
		}
	}

}


} /* namespace NetworKit */


