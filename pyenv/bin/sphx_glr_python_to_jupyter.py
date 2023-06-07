#!/home/petholza/Documents/networkit/networkit/pyenv/bin/python3
r"""
Sphinx Gallery Notebook converter
=================================

Exposes the Sphinx-Gallery Notebook renderer to directly convert Python
scripts into Jupyter Notebooks.

"""
# Author: Óscar Nájera
# License: 3-clause BSD


from sphinx_gallery.notebook import python_to_jupyter_cli


if __name__ == '__main__':
    python_to_jupyter_cli()
