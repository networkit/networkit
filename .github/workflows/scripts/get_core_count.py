import os
try:
	# For Linux
	print(len(os.sched_getaffinity(0)))
except AttributeError:
	# For macOS
	print(os.cpu_count())