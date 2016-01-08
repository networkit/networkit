#
# file: multiprocessing_helper.py
# author: Mark Erb
#

import multiprocessing
from collections import deque

import sys, traceback


def numberOfProcessors():
	""" return the number of Processors """
	return multiprocessing.cpu_count()


class Worker(multiprocessing.Process):
	""" Worker Process """

	def __init__(self, tasks, results):
		""" constructor: init worker with task and result queue """
		multiprocessing.Process.__init__(self)
		self.__tasks = tasks
		self.__results = results


	def run(self):
		""" gain tasks and compute them until ThreadPool sends termination signal """
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
				print(traceback.format_exc(), flush=True)
			result = (task.getType(), task.getName(), data)
			self.__tasks.task_done()
			self.__results.put(result)


class ThreadPool():
	""" Threadpool

	suggested usage:

	from . import job

	n = multiprocessing.numberOfProcessors() * 2
	pool = multiprocessing.ThreadPool(n)

	class Foo(job.Job):
		...

	tasks = []
	tasks.append(Foo())

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
		""" constructor
		Args:
			numberOfWorkers: number of workers to create
			isParallel: parallel or sequential task computation
		"""
		self.__numberOfTasks = 0
		self.__numberOfWorkers = numberOfWorkers
		self.__isParallel = isParallel
		if self.__numberOfWorkers < 1:
			self.__numberOfWorkers = 0
			self.__isParallel = False
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
					self.__numberOfWorkers = count
					if count < 1:
						self.__tasks.close()
						self.__results.close()
						self.__isParallel = False
						self.__tasks = None
						self.__results = None
						self.__workers = None
					print(e)
					break
		if not self.__isParallel:
			self.__tasks = deque()


	def numberOfTasks(self):
		""" current number of unfinished tasks """
		return self.__numberOfTasks


	def numberOfWorkers(self):
		""" initilized number of workers """
		return self.__numberOfWorkers


	def isParallel(self):
		""" are the workers run parallel or sequencial  """
		return self.__isParallel


	def put(self, task):
		""" assign a task """
		if self.__isParallel:
			self.__tasks.put(task)
		else:
			self.__tasks.append(task)
		self.__numberOfTasks += 1


	def get(self):
		""" gains results from workers """
		if self.__isParallel:
			result = self.__results.get()
		else:
			task = self.__tasks.popleft()
			try:
				data = task.run()
			except Exception as e:
				print("Error: " + task.getType() + " - " + task.getName(), flush=True)
				print(str(e), flush=True)
				print(traceback.format_exc(), flush=True)
				data = None
			result = (task.getType(), task.getName(), data)
		self.__numberOfTasks -= 1
		return result;


	def join(self):
		""" wait for workers to finish and close ThreadPool """
		if self.__isParallel:
			for i in range(self.__numberOfWorkers):
				self.__tasks.put(None)
			self.__tasks.join()
			self.__tasks.close()
			self.__results.close()
