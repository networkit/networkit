/*
 * Parameters.cpp
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#include "Parameters.h"

namespace NetworKit {

void Parameters::setInt(std::string key, std::int64_t value) {
	intMap.insert(std::make_pair(key, value));
}

void Parameters::setDouble(std::string key, double value) {
	doubleMap.insert(std::make_pair(key, value));
}

void Parameters::setString(std::string key, std::string value) {
	stringMap.insert(std::make_pair(key, value));
}

void Parameters::setBool(std::string key, bool value) {
	boolMap.insert(std::make_pair(key, value));
}

std::int64_t Parameters::getInt(std::string key) const {
	return intMap.find(key)->second;
}

double Parameters::getDouble(std::string key) const {
	return doubleMap.find(key)->second;
}

std::string Parameters::getString(std::string key) const {
	return stringMap.find(key)->second;
}

bool Parameters::getBool(std::string key) const {
	return boolMap.find(key)->second;
}

} /* namespace NetworKit */
