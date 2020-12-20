import os
import glob
import sys
import tempfile
import subprocess

READONLY = True
VERBOSE = False
SHOW_DIFFS = True
MIN_CLANG_FORMAT_VERSION = 8
MAX_CLANG_FORMAT_VERSION = 11

def setup(args = None):
	"""
	Parses command-line arguments and sets module runtime flags.
	"""
	if args is None:
		args = sys.argv

	if "-h" in args or "--help" in args:
		print("-h | --help    Print this message\n"
			  "-v | --verbose Report more details\n"
			  "-q | --nodiff  Do not report diffs of changes\n"
			  "-w | --write   Apply fixes to source code file\n")
		sys.exit(0)

	global READONLY
	READONLY = not ("-w" in args)

	global VERBOSE
	VERBOSE = ("-v" in args or "--verbose" in args)

	global SHOW_DIFFS
	SHOW_DIFFS = not ("-q" in args or "--nodiff" in args)

def isReadonly():
	"""
	Returns whether module is in read-only mode (default: True).
	In this case, FileRewriter.commit is silently deactivated.
	"""
	return READONLY

def isVerbose():
	"""
	Returns whether user requested verbose output (default: False).
	"""
	return VERBOSE

def doReportDiff():
	"Returns whether user requested diffs (default: True)."
	return SHOW_DIFFS

def reportChange(msg):
	"""
	Informs the user (currently via STDOUT) about a change that was carried out.
	The message is prefixed with either [ERROR] (if we operate in readonly mode,
	since the change cannot be carried out) or [FIXED] (if write are allowed)
	"""
	msgType = "[ERROR]" if isReadonly() else "[FIXED]"
	print(msgType, msg)

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

def computeAndReportDiff(originalFilename, formatedFilename):
	"""Compute a colorful diff between the original file and the formatted one"""
	p = subprocess.Popen(["diff", "-a", "--color=always", originalFilename, formatedFilename],
						 stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=-1)
	output, _ = p.communicate()

	print("-" * 20, "Begin of diff", "-" * 20)
	print("Input file: ", originalFilename)
	sys.stdout.buffer.write(output)
	print("-" * 21, "End of diff", "-" * 21)

class FileRewriter:
	"""
	The FileRewriter acts as a wrapper to safely read and update a text file provided to the constructor.
	To read the file iterate for lines(). To write a (un)modified line call write(line).
	Once you are done, call commit() to replace the file atomically.
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
		return "".join( self.in_lines ) == "".join( self.out_lines )

	def commit(self, force = False, outputDiff = None):
		"""
		Writes buffered lines into the file specified in the constructor.
		Writes are carried out atomically, i.e. we first write into a temporary file
		and then replace the original.
		If the module is currently readonly (which is the default, see setup),
		the original file is not replaced.
		This function may only be called if not reads are taking place currently (e.g.,
		via lines())

		If outputDiff is None, its value is copied from doReportDiff(). If it is True,
		Diff is reported to the user
		"""
		assert(not self.reading)

		doReport = doReportDiff() if outputDiff is None else outputDiff
		doWrite = isReadonly() and not force

		if not (doWrite or doReport):
			return

		with tempfile.TemporaryDirectory(dir=getNetworKitRoot()) as dir:
			filename = os.path.join(dir, "Writer")
			with open(filename, "w") as out:
				out.write("".join(self.out_lines))

			if doReport:
				computeAndReportDiff(self.path, filename)

			if doWrite:
				os.replace(filename, self.path)

	def reportDiff(self, force = False):
		"""
		Prints a diff between the original file and the changes to STDOUT.
		If the module is configured not to show diffs and force is False (default),
		no output is generated
		"""
		if not doReportDiff() and not force:
			return

		with tempfile.TemporaryDirectory(dir=getNetworKitRoot()) as dir:
			filename = os.path.join(dir, "Diff")
			with open(filename, "w") as out:
				out.write("".join(self.out_lines))
			reportDiff(self.path, filename)

	def _writeToFile(self, fileHandle):
		"""INTERNAL USE ONLY. Writes commited lines into a previously opened file."""
		fileHandle.write("".join(self.out_lines))
