/*
 * Timer.h
 *
 *  Created on: 28.11.2012
 *      Author: cls
 */

#ifndef TIMER_H_
#define TIMER_H_


//#include <time.h>
//#include <sys/time.h>
//
//#ifdef __MACH__
//#include <mach/clock.h>
//#include <mach/mach.h>
//#endif
//
//
//
//struct timespec ts;
//
//#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
//clock_serv_t cclock;
//mach_timespec_t mts;
//host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
//clock_get_time(cclock, &mts);
//mach_port_deallocate(mach_task_self(), cclock);
//ts.tv_sec = mts.tv_sec;
//ts.tv_nsec = mts.tv_nsec;
//
//#else
//clock_gettime(CLOCK_REALTIME, &ts);
//#endif

namespace EnsembleClustering {


/**
 * TODO: Platform-agnostic timer class
 */
class Timer {
public:
	Timer();
	virtual ~Timer();
};

} /* namespace EnsembleClustering */
#endif /* TIMER_H_ */
