#include "Loglevel.h"

namespace Aux {

	int LogLevel::current = 4;
	void LogLevel::setLevel(int level) {
		current = level;
	}
	int LogLevel::getLevel() {
		return current;
	}
}
