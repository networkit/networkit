#!/usr/bin/env python3
"""
This tool iterates over all C++ files and applies clang-format
to them. Source files need to opt-in by containing the term
"networkit-format" (preferably in a comment near the header).

If used in read-only mode it returns a negative exit code if
a change is necessary.
"""
import filecmp
import nktooling as nkt
import tempfile
import sys
import shutil
import subprocess
import os

def findClangFormat():
    """Tries to find clang-format-XXX variants within the path"""
    allowed = ["clang-format" + x for x in ["-9", "-8", "-7", ""]]
    for candidate in allowed:
        if not shutil.which(candidate) is None:
            return candidate

    raise "clang-format binary not found. We searched for:\n " + "\n ".join(allowed)

def subscribedToFormat(filename, pattern = "networkit-format"):
    """If pattern is present within the file, this file subscribed to auto formatting."""
    with open(filename, 'r') as file:
        return any( ((pattern in line) for line in file) )

def runClangFormat(inputFilename, outputFilename, clangFormatBinary = 'clang-format-8'):
    """Execute clang-format onto inputFilename and stores the result in outputFilename"""
    with open(outputFilename, "w") as outfile:
        subprocess.call([clangFormatBinary, '-style=file', inputFilename], stdout=outfile)


nkt.setup()

# we need to set pwd in order for clang-format to find the .clang-format file!
os.chdir(nkt.getNetworKitRoot())

numberNonCompliant = 0
numberFileSkipped = 0

clangFormatCommand = findClangFormat()
with tempfile.TemporaryDirectory(dir=nkt.getNetworKitRoot()) as tempDir:
    files = nkt.getCXXFiles()
    for file in files:
        if not subscribedToFormat(file):
            numberFileSkipped += 1
            continue

        tempFile = os.path.join(tempDir, 'cfOutput')
        runClangFormat(file, tempFile, clangFormatCommand)

        if not filecmp.cmp(file, tempFile, shallow=False):
            numberNonCompliant += 1
            if not nkt.isReadonly():
                os.replace(tempFile, file)

            nkt.reportChange(file + " is non-compliant")

print("Scanned %d files (skipped %d files without subscription). Non-compliant files: %d." %
      (len(files), numberFileSkipped, numberNonCompliant))

if numberNonCompliant > 0:
    nkt.failIfReadonly(__file__)