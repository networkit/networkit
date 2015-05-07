#include <string>
#include <stdexcept>

#ifndef cython_call_operator
#define cython_call_operator operator()

void throw_runtime_error(std::string message) {
	throw std::runtime_error(message);
};
#endif