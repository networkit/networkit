#!/usr/bin/env python3

# README
# This script should make your life as NetworKit developer easier. Once started it will recompile NetworKit each time you change a source file (*.h, *.cpp). You can give it also a test case name and it will run it after a successful compilation. Just run the script, start working and when you want to see the results of you work, safe and the script will trigger compilation (and testing).
#
# To only compile:
# 	./test.py
# To compile and run all test cases use:
# 	./test.py test
# To compile and run test cases in GraphGTest with optimize level set to Dbg and loglevel to DEBUG:
# 	./test.py -debug GraphGTest
# To compile and run testParallelPartitionCoarsening test case with loglevel set to TRACE use:
# 	./test.py -loglevel=TRACE testParallelPartitionCoarsening
# To see more options see the help:
# 	./test.py -h
#
# The script requires watch. You can install watchdog for you local user using 'pip3 install --user watchdog'

import sys
import time
import argparse
import subprocess
import multiprocessing
from watchdog.observers import Observer
from watchdog.events import PatternMatchingEventHandler

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--optimize', choices=['Dbg', 'Pro', 'Opt'], default='Opt', help='optimize level')
parser.add_argument('-t', '--target', choices=['Tests', 'Core', 'Lib'], default='Tests', help='target for compilation')
parser.add_argument('-l', '--loglevel', choices=['ERROR', 'WARN', 'INFO', 'DEBUG', 'TRACE'], default='ERROR', help='log level')
parser.add_argument('-j', '--jobs', type=int, help='number of jobs to use for compilation, default value is half of the number of CPUs')
parser.add_argument('-debug', action='store_true', help='alias for \'--optimize=Dbg --target=Tests --loglevel=DEBUG\'')
parser.add_argument('-once', action='store_true', help='compile and execute only once')
parser.add_argument('gtest_filter', nargs='?', help='filter passed to gtest, will be surrounded by *, to run all test cases enter \'test\'')
args = parser.parse_args()

# additional argument parsing
if args.debug:
	args.optimize = 'Dbg'
	args.target = 'Tests'
	args.loglevel = 'DEBUG'
if not args.jobs:
	args.jobs = multiprocessing.cpu_count() // 2
if args.gtest_filter and args.target != 'Tests':
	print('Target \'' + args.target + '\' does not support running test cases. Switch to target \'Test\' or remove argument \'' + args.gtest_filter + '\'')
	sys.exit()

# hard coded configuration
sourceDirectory = 'networkit/cpp'
workingDirectory = None

def buildAntTest():
	# compile
	result = subprocess.call(['scons', '--optimize=' + args.optimize, '--target=' + args.target, '-j', str(args.jobs)], cwd = workingDirectory)

	# run test to if compilation was successful and a test name was specified
	if result == 0 and args.target == 'Tests' and args.gtest_filter:
		subprocess.call(['./NetworKit-' + args.target + '-' + args.optimize, '--tests', '--loglevel=' + args.loglevel, '--gtest_filter=*' + args.gtest_filter + '*'], cwd = workingDirectory)

class ContinousTesting(PatternMatchingEventHandler):
	def __init__(self):
		super(ContinousTesting, self).__init__(patterns = ['*.h', '*.cpp'])
	
	def on_any_event(self, event):
		buildAntTest()
		print('\nWaiting for code changes...')

buildAntTest()
if args.once:
	sys.exit()

# listen for file changes
observer = Observer()
observer.schedule(ContinousTesting(), path = sourceDirectory, recursive = True)
observer.start()
print('\nWaiting for code changes...')

# run forever
try:
	while True:
		time.sleep(1)
except KeyboardInterrupt:
	observer.stop()

observer.join()
