#
# file: multiprocessing.py
# author: Mark Erb
#

import multiprocessing
from collections import deque


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
	def __init__(self, numberOfWorkers, isParallel=True):
		self.__numberOfWorkers = numberOfWorkers
		self.__numberOfTasks = 0
		self.__isParallel = isParallel
		if self.__isParallel:
			self.__tasks = multiprocessing.JoinableQueue()
			self.__results = multiprocessing.Queue()
			self.__workers = [Worker(self.__tasks, self.__results) for i in range(self.__numberOfWorkers)]
			count = 0
			for w in self.__workers:
				w.deamon = True
				try:
					w.start()
					count += 1
				except Exception as e:
					if count < 1:
						raise RuntimeError(e)
					break
		else:
			self.__tasks = deque()


	def numberOfTasks(self):
		""" Current number of unfinished tasks """
		return self.__numberOfTasks

		
	def numberOfWorkers(self):
		""" Initilized number of workers """
		return self.__numberOfWorkers
		

	def put(self, task):
		""" Assign a task """
		if self.__isParallel:
			self.__tasks.put(task)
		else:
			self.__tasks.append(task)
		self.__numberOfTasks += 1


	def get(self):
		if self.__isParallel:
			result = self.__results.get()
		else:
			task = self.__tasks.popleft()
			try:
				data = task.run()
			except Exception as e:
				print("Error: " + task.getType() + " - " + task.getName(), flush=True)
				print(str(e), flush=True)
			result = (task.getType(), task.getName(), data)
		self.__numberOfTasks -= 1
		return result;


	def join(self):
		if self.__isParallel:
			for i in range(self.__numberOfWorkers):
				self.__tasks.put(None)
			self.__tasks.join()
			self.__tasks.close()
			self.__results.close()