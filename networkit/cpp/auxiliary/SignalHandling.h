#ifndef SIGNALHANDLING_H
#define SIGNALHANDLING_H

#include <csignal>
//#include <atomic>

namespace Aux {
	class SignalHandling {
	private:
		//static std::atomic<bool> running;
		//static std::atomic<bool> handlerInitialized;
		static bool running;
		static bool handlerInitialized;

	public:
		static bool isRunning();

		static void setRunning(bool running);

		static void init();
	};
}

#endif /* SIGNALHANDLING_H */