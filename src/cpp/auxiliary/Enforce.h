#ifndef ENFORCE_H_
#define ENFORCE_H_

#include <cassert>
#include <stdexcept>
#include <string>

namespace Aux {

/**
 * Enforces that b is true and throws an Exception otherwise.
 *
 * If provided, msg must not be null, otherwise the behavior of this
 * function is undefined.
 */
template<typename Exception = std::runtime_error>
inline void enforce(bool b, const char* msg = "") {
	assert(msg && "Message to enforce must not be nullptr");
	if (!b) {
		throw Exception{msg};
	}
}

/**
 * Overload that accepts a std::string. This is mainly for convenience
 * while keeping the default-version free of unneeded allocations.
 */
template<typename Exception = std::runtime_error>
inline void enforce(bool b, const std::string& msg) {
	enforce<Exception>(b, msg.c_str());
}

/**
 * Checks that the provided fstream is opened and throws an exception otherwise.
 */
template<typename Stream>
inline void enforceOpened(const Stream& stream) {
	enforce(stream.is_open());
}

} // namespace Aux



#endif // ENFORCE_H_
