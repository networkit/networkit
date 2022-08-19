#!/usr/bin/env python3
import nbformat
import os
import sys
from nbconvert.preprocessors import ExecutePreprocessor
import asyncio

def run_notebook(path):
    nb_name, _ = os.path.splitext(os.path.basename(path))
    dirname = os.path.dirname(path)

    with open(path) as f:
        nb = nbformat.read(f, as_version=4)

    print("Start ", path)
    sys.stdout.flush()

    proc = ExecutePreprocessor(timeout=600, kernel_name='python3')
    proc.allow_errors = True

    # For Python >= 3.8 on Windows, the Preprocessor throws a NotImplementedError
    # This check is inserted as a workaround
    # See: https://github.com/jupyter/nbconvert/issues/1372
    if sys.version_info[0] == 3 and sys.version_info[1] >= 8 and sys.platform.startswith('win'):
        asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())

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
        print(errors)
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
