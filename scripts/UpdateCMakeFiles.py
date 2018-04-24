#!/usr/bin/env python3
import os
import fnmatch

def stripCppExtension(fn):
    assert(fn[-4:] == ".cpp")
    return fn[:-4]


srcDir = os.path.join("networkit", "cpp")
for (dirpath, dirnames, filenames) in os.walk(srcDir):
    dirname = os.path.basename(dirpath)
    isDirTest = (dirname == "test")

    with open(os.path.join(dirpath, "CMakeLists.txt"), "w") as f:
        f.write("# This file was automatically generated with UpdateCMakeFiles.py\n" +
                "# Changes will be overridden; do not edit.\n\n")

        f.write("".join( ("add_subdirectory(\"%s\")\n" % dn for dn in sorted(dirnames)) ))
        if (dirnames): f.write("\n")

        if (dirpath == srcDir):
            # We are at root level and do not include local file
            continue

        cppFiles = sorted(map(stripCppExtension, fnmatch.filter(filenames, "*.cpp")))

        # check Test name
        for fn in cppFiles:
            isFileTest = fn.endswith("Test")
            if (isDirTest and not isFileTest):
                print("Warning: File in test directory is missing suffix Test.cpp: ", fn)
            elif  (not isDirTest and isFileTest):
                print("Warning: File with Test.cpp suffix not in test directory: ", fn)

        cmakeCmd = "networkit_add_test" if isDirTest else "networkit_add_source"

        f.write("".join(("%s(%s)\n" % (cmakeCmd, sym) for sym in cppFiles)))

