#!/usr/bin/env python

import time
import subprocess
import sys
from watchdog.observers import Observer
from watchdog.events import PatternMatchingEventHandler

# get test name (or parts of it) from first command parameter, other leave it empty (and don't run tests)
testName = sys.argv[1] if len(sys.argv) > 1 else ''

# ERROR | WARN | INFO | DEBUG | TRACE
logLevel = 'ERROR'

sourceDirectory = 'src/'
workingDirectory = None

def buildAntTest():
    # compile
    result = subprocess.call(['scons', '--optimize=O', '--target=Tests'], cwd = workingDirectory)

    # run test to if compilation was successful and a test name was specified
    if result == 0 and len(testName) > 0:
        subprocess.call(['./NetworKit-Tests-O', '--tests', '--loglevel=' + logLevel, '--gtest_filter=*' + testName + '*'], cwd = workingDirectory)
    print("\nWaiting for code changes...")

if __name__ == "__main__":
    buildAntTest()

    class ContinousTesting(PatternMatchingEventHandler):
        def __init__(self):
            super(ContinousTesting, self).__init__(patterns = ['*.h', '*.cpp'])
        
        def on_any_event(self, event):
            buildAntTest()
    
    event_handler = ContinousTesting()

    observer = Observer()
    observer.schedule(event_handler, path = sourceDirectory, recursive = True)
    observer.start()

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()

    observer.join()
