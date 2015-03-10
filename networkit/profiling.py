from networkit import *

import matplotlib.pyplot as plt
from matplotlib._pylab_helpers import Gcf
from IPython.core.pylabtools import print_figure
from base64 import b64encode
from IPython.core.display import HTML


def asImage(plotFunction, *args, **kwargs):
    """
    Call any plot function with the given argument and return the image in an HTML <img> tag.
    """
    plotFunction(*args, **kwargs)
    # Get a handle for the plot that was just generated
    fig = Gcf.get_all_fig_managers()[-1].canvas.figure
    # Generate a data URL for the image
    imageData = "data:image/png;base64,{0}".format(b64encode(print_figure(fig)).decode("utf-8"))
    # Remove the plot from the list of plots for the current cell
    Gcf.destroy_fig(fig)
    # generate img tag
    image = "<img src='{0}'\>".format(imageData)
    return image


profileTemplate = """
	<h1>Profile of Network</h1>

	<h2>Degree Distribution</h2>

	{ddImage}
"""

def profile(G):
	"""
	Output profile page of network as HTML
	"""
	ddImage  = asImage(plot.degreeDistribution, G)
	page = HTML(profileTemplate.format(**locals()))
	return page
