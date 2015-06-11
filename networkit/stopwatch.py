# -*- coding: utf-8 -*-
#
# Copyright (C) 2008 John Paulett (john -at- 7oars.com)
# All rights reserved.
#
# This software is licensed as described in the file COPYING, which
# you should have received as part of this distribution.

import time

"""stopwatch is a very simple Python module for measuring time.
Great for finding out how long code takes to execute.

>>> import stopwatch
>>> t = stopwatch.Timer()
>>> t.elapsed
3.8274309635162354
>>> print t
15.9507198334 sec
>>> t.stop()
30.153270959854126
>>> print t
30.1532709599 sec

Decorator exists for printing out execution times:
>>> from stopwatch import clockit
>>> @clockit
    def mult(a, b):
        return a * b
>>> print mult(2, 6)
mult in 1.38282775879e-05 sec
6
       
"""

__version__ = '0.3.1'
__author__ = 'John Paulett <http://blog.7oars.com>'

class Timer(object):
    def __init__(self):
        self.__stopped = None
        self.__start = self.__time()

    def stop(self):
        """Stops the clock permanently for the instance of the Timer.
        Returns the time at which the instance was stopped.
        """
        self.__stopped = self.__last_time()
        return self.elapsed

    def elapsed(self):
        """The number of seconds since the current time that the Timer
        object was created.  If stop() was called, it is the number
        of seconds from the instance creation until stop() was called.
        """
        return self.__last_time() - self.__start
    elapsed = property(elapsed)
    
    def start_time(self):
        """The time at which the Timer instance was created.
        """
        return self.__start
    start_time = property(start_time)
    
    def stop_time(self):
        """The time at which stop() was called, or None if stop was 
        never called.
        """
        return self.__stopped 
    stop_time = property(stop_time)        
    
    def __last_time(self):
        """Return the current time or the time at which stop() was call,
        if called at all.
        """
        if self.__stopped is not None:
            return self.__stopped
        return self.__time()
    
    def __time(self):
        """Wrapper for time.time() to allow unit testing.
        """
        return time.time()
    
    def __str__(self):
        """Nicely format the elapsed time
        """
        return str(self.elapsed) + ' sec'

def clockit(func):
    """Function decorator that times the evaluation of *func* and prints the
    execution time.
    """
    def new(*args, **kw):
        t = Timer()
        retval = func(*args, **kw)
        t.stop()
        print('%s in %s') % (func.__name__, t)
        del t
        return retval
    return new
