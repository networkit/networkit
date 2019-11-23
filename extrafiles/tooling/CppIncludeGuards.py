#!/usr/bin/env python3
"""
This tool iterates over all C++ files and replaces tab-based
indentation with spaces (tabwidth = 4). If used in read-only
mode it returns a negative exit code if a change is necessary.
"""
import nktooling as nkt
import re
import os, sys

nkt.setup()

numScannedFiled = 0
numNoncompliantFiles = 0
numIllegalFiles = 0

files = nkt.getCXXFiles(includeSources = False)

includePragmaOnce = False

rePragmaOnce = re.compile("#pragma\s+once$")

for file in files:
	if "/ext/" in file:
		continue

	numScannedFiled += 1

	# compute new include guard name
	newGuard = os.path.relpath(file, os.path.join(nkt.getNetworKitRoot(), "include"))
	newGuard = re.sub("include/", "", newGuard)
	newGuard = re.sub("[/\.]", "_", newGuard)
	newGuard = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', newGuard)
	newGuard = re.sub('([a-z0-9])([A-Z])', r'\1_\2', newGuard)
	newGuard = re.sub('__+', '_', newGuard + "_")
	newGuard = newGuard.upper()

	rw = nkt.FileRewriter(file)

	oldGuard = None
	requireDefine = False
	labelStack = []
	prevLineEmpty = True

	try:
		if includePragmaOnce:
			rw.write("#pragma once\n")

		for lineno, line in enumerate(rw.lines()):
			sline = line.strip()
			prevLineEmpty = not sline

			if rePragmaOnce.match(sline):
				continue

			# search for the first #ifndef (which we interpret as guard)
			if oldGuard is None:
				if sline.startswith("#ifndef"):
					oldGuard = re.split("\\s", sline)[-1]
					if not "_H" in oldGuard:
						raise Exception("%s:%d expected #define of include guard. Label needs to contain '_H'. Found: %s"% (file, lineno+1, oldGuard))
					requireDefine = True
					labelStack = [newGuard]
					continue

			if requireDefine:
				if not sline.startswith("#define"):
					raise Exception("%s:%d expected #define of include guard." % (file, lineno+1))
				if not sline.endswith(oldGuard):
					raise Exception("%s:%d line assumed to be #define of include guard. label mismatch. expected %s" % (file, lineno+1, oldGuard))
				requireDefine = False
				rw.write("#ifndef %s\n#define %s\n" % (newGuard, newGuard))
				continue

			if sline.startswith("#ifdef") or sline.startswith("#ifndef"):
				label = re.split("\\s", sline)[-1]
				labelStack.append(label)

			elif sline.startswith("#if "):
				labelStack.append(None)

			elif sline.startswith("#endif"):
				label = labelStack.pop()
				if not label is None:
					indent = line[:line.find("#endif")]
					rw.write("%s#endif // %s\n" % (indent, label))
					continue

			rw.write(line)

		if oldGuard is None:
			raise Exception("%s no include guard found" % file)

		if not rw.isIdentical():
			nkt.reportChange(file + " changed.")
			rw.commit()
			numNoncompliantFiles += 1

	except Exception as e:
		if str(e):
			print("[PARSE_ERROR] ", e)
		else:
			print("[PARSE_ERROR] ", file, "unknown error")

		numIllegalFiles += 1

print("Scanned %d files. Non-compliant files: %d. Illegal files: %d" % (len(files), numNoncompliantFiles, numIllegalFiles))

if numNoncompliantFiles > 0:
	nkt.failIfReadonly(__file__)

if numIllegalFiles > 0:
	sys.exit(-1)

