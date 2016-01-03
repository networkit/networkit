import time


class Timer:
	""" Use the Python with-statement to time your code
	with this timer. """

	def __enter__(self):
		self.start = time.time()
		return self

	def __exit__(self, *args):
		self.end = time.time()
		self.elapsed = round(self.end - self.start, 6)


# Logging

def info(message):
	print(message)

def error(message):
	print(message)

def debug(message):
	pass
	# print(message)
