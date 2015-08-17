#
# file: threadpool.py
# author: Mark Erb
#

import multiprocessing


def numberOfProcessors():
	return multiprocessing.cpu_count()


class Worker(multiprocessing.Process):
	""" TODO: """

	def __init__(self, tasks, results):
		multiprocessing.Process.__init__(self)
		self.__tasks = tasks
		self.__results = results


	def run(self):
		while True:
			task = self.__tasks.get()
			if task is None:
				self.__tasks.task_done()
				break
			data = "Error"
			try:
				data = task.run()
			except Exception as e:
				print("Error: " + task.getType() + " - " + task.getName(), flush=True)
				print(str(e), flush=True)
			result = (task.getType(), task.getName(), data)
			self.__tasks.task_done()
			self.__results.put(result)


class ThreadPool():
	"""
		Suggested Syntax:

		n = multiprocessing.numberOfProcessors() * 2
		pool = multiprocessing.ThreadPool(n)

		class Job:
			def run(self):
				...
				return ...

			def getType(self):
				...

			def getName(self):
				...

		tasks = []
		tasks.append(Job())

		for task in tasks:
			pool.put(task)
		while pool.numberOfTasks() > 0:
			(type, name, data) = pool.get()
			try:
				... unpack data & post processing ...
			except:
				pass
		pool.join()
	"""

	def __init__(self, numberOfWorkers):
		self.__numberOfWorkers = numberOfWorkers
		self.__numberOfTasks = 0
		self.__tasks = multiprocessing.JoinableQueue()
		self.__results = multiprocessing.Queue()
		self.__workers = [Worker(self.__tasks, self.__results) for i in range(self.__numberOfWorkers)]
		for w in self.__workers:
			w.deamon = True
			w.start()


	def numberOfTasks(self):
		""" Current number of unfinished tasks """
		return self.__numberOfTasks


	def put(self, task):
		""" Assign a task """
		self.__tasks.put(task)
		self.__numberOfTasks += 1


	def get(self):
		result = self.__results.get()
		self.__numberOfTasks -= 1
		return result;


	def join(self):
		for i in range(self.__numberOfWorkers):
			self.__tasks.put(None)
		self.__tasks.join()
