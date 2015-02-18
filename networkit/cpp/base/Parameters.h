/*
 * Parameters.h
 *
 *  Created on: 24.06.2013
 *      Author: cls
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <unordered_map>
#include <string>
#include <cstdint>

namespace NetworKit {

/**
 * A class for transporting parameters.
 */
class Parameters {
public:
	void setInt(std::string key, std::int64_t value);

	void setDouble(std::string key, double value);

	void setString(std::string key, std::string value);

	void setBool(std::string key, bool value);


	std::int64_t getInt(std::string key) const;

	double getDouble(std::string key) const;

	std::string getString(std::string key) const;

	bool getBool(std::string key) const;



protected:

	std::unordered_map<std::string, std::int64_t> intMap;
	std::unordered_map<std::string, double> doubleMap;
	std::unordered_map<std::string, std::string> stringMap;
	std::unordered_map<std::string, bool> boolMap;

};

} /* namespace NetworKit */
#endif /* PARAMETERS_H_ */
