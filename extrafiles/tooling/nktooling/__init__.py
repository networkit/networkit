import os
import glob
import sys
import tempfile

READONLY = True

def setup(args = None):
	"""
	Parses command-line arguments and sets module runtime flags.
	"""
	if args is None:
		args = sys.argv

	if "-h" in args or "--help" in args:
		print("-h | --help    Print this message\n"
			  "-w | --write   Apply fixes to source code file\n")
		sys.exit(0)

	global READONLY
	READONLY = not ("-w" in args)

def isReadonly():
	"""
	Returns whether module is in read-only mode (default: True).
	In this case, FileRewriter.commit is silently deactivated.
	"""
	return READONLY

def failIfReadonly(prgname):
	"""
	Helper function indented to give a summary after a test in read-only mode.
	In read-only mode, fixes could not be applied, hence, we fail with a negative
	return code. This can be used during CI to check whether repository is compliant.
	"""
	if not isReadonly():
		return

	print("Changes are necessary.")
	print("Use '%s -w' to apply changes." % prgname)
	sys.exit(1)

def getNetworKitRoot():
	"""
	Returns an absolute path to the project root (which is derived from the module's absolute path)
	"""
	return os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..", ".."))

def getCXXFiles(includeSources = True, includeHeaders = True):
	"""
	Returns a list of all NetworKit C++ source files.

	includeSource : bool
		controls whether files from the source directory (networkit/cpp) are included.
	includeHeaders : bool
		controls whether file from the include directory (include/networkit) are included.
	"""
	files = []
	if includeSources:
		newfiles = glob.glob(os.path.join(getNetworKitRoot(), "networkit", "cpp", "**", "*.[ch]pp"), recursive=True)
		assert(len(newfiles) > 0)
		files += newfiles

	if includeHeaders:
		newfiles = glob.glob(os.path.join(getNetworKitRoot(), "include", "networkit", "**", "*.[ch]pp"), recursive=True)
		assert(len(newfiles) > 0)
		files += newfiles

	return files

class FileRewriter:
	"""
	The FileRewriter acts as an wrapper to safely
	read and update a text file provided to the constructor.
	To read the file iterate for lines().
	To write a (un)modified line call write(line).
	Once your done, call commit() to replace the file atomically.
	Observe that commit will do nothing if the nktooling module
	is in readonly mode (which is the default).
	"""
	def __init__(self, filepath):
		self.path = filepath
		self.in_lines = []
		self.out_lines = []
		self.reading = False

	def lines(self):
		"""
		Iterates over all input lines.
		"""
		with open(self.path, "r") as fh:
			self.reading = True
			for line in fh:
				self.in_lines.append(line)
				yield line

			self.reading = False

	def write(self, line):
		"""
		Adds data to the output buffer.
		"""
		self.out_lines.append(line)

	def isIdentical(self):
		"""
		Checks whether the content written matches the content read (so far).
		"""
		return self.in_lines == self.out_lines

	def commit(self, force = False):
		"""
		Writes buffered lines into the specified in the constructor.
		Writes are carried out atomically, i.e. we first write into a temporary file
		and then replace the original.
		If the module is currently readonly (which is the default, see setup) no change
		is carried out.
		This function may only be called if not reads are taking place currently (e.g.,
		via lines())
		"""
		assert(not self.reading)
		if isReadonly() and not force:
			return

		with tempfile.TemporaryDirectory(dir=getNetworKitRoot()) as dir:
			filename = os.path.join(dir, "Writer")
			with open(filename, "w") as out:
				out.write("".join(self.out_lines))
			os.replace(filename, self.path)