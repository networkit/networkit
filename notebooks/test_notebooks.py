#!/usr/bin/env python3
import nbformat
import os
import sys
from nbconvert.preprocessors import ExecutePreprocessor

def run_notebook(path):
    nb_name, _ = os.path.splitext(os.path.basename(path))
    dirname = os.path.dirname(path)

    with open(path) as f:
        nb = nbformat.read(f, as_version=4)

    print("Start ", path)
    sys.stdout.flush()

    proc = ExecutePreprocessor(timeout=600, kernel_name='python3')
    proc.allow_errors = True

    proc.preprocess(nb, {'metadata': {'path': dirname}})
    errors = []
    for cell in nb.cells:
        if 'outputs' in cell:
            for output in cell['outputs']:
                if output.output_type == 'error':
                    errors.append(output)
    if errors == []:
        print(" " + path + " test successfully completed.")
        return 0
    else:
        print(" " + path + " test exited with errors.")
        return 1

if __name__ == '__main__':
    if len(sys.argv) == 1:
        raise Exception("Expected path to notebook")
    else:
        path = sys.argv[1]
    status = 0
    for file in os.listdir(path):
        if file.endswith(".ipynb"):
            status += run_notebook(os.path.join(path, file))
    sys.exit(status)
