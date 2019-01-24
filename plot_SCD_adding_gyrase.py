"""
Plot the supercoiling density before and after adding Gyrase (Chong experiment)
"""

import sys
import numpy as np
import TCDS.simulation as sim
import TCDS.plotting as plotting
import matplotlib
matplotlib.rcParams.update({'font.size': 13})

# parameters file
INI_file=sys.argv[1]
# location of the npz files before adding the Gyrase (TopoI only)
output_dir_pre=sys.argv[2]
# location of the npz files after adding the Gyrase
output_dir_post=sys.argv[3]

plotting.plot_topoI_gyrase(INI_file, output_dir_pre, output_dir_post)