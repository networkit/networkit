#!/usr/bin/env python3
"""
This tool iterates over all C++ files and replaces tab-based
indentation with spaces (tabwidth = 4). If used in read-only
mode it returns a negative exit code if a change is necessary.
"""
import nktooling as nkt
import sys

nkt.setup()

TabWidth = 4
noncompliantFiles = 0

files = nkt.getCXXFiles()

for file in files:
	rw = nkt.FileRewriter(file)
	for line in rw.lines():
		x = 0
		skip = 0
		for (i, c) in enumerate(line):
			if c == '\t':
				x += TabWidth - (x % TabWidth)
			elif c == ' ':
				x += 1
			else:
				skip = i
				break

		rw.write((" " * int(x)) + line[i:])

	if not rw.isIdentical():
		nkt.reportChange(file + " is (partially) indented with tabs")
		rw.commit()
		noncompliantFiles += 1

print("Scanned %d files. Non-compliant files: %d." % (len(files), noncompliantFiles))

if noncompliantFiles > 0:
	nkt.failIfReadonly(__file__)