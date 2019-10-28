import os
import sys

class MissingDependencyError (RuntimeError):
	def __init__(self, package):
		self.package = package
		super().__init__("Optional dependency {} is not installed".format(package))

