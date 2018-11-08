#! /usr/bin/python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import warnings
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# Define 2d-gaussian, centered at (0,0), scale = 1
def gaussian(x,y, sigmax, sigmay):
    return np.exp( -( (x**2 / (2*sigmax)) + (y**2 / (2*sigmay) ) ) )

# Parse command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawTextHelpFormatter,
    description='''
SSAGES 2d Metadynamics Post-Analysis
====================================

This script generates a colored plot from the
SSAGES metadynamics output.

''')

parser.add_argument('-i', '--input', type=str,
                    help='SSAGES Metadynamics output file.',
                    default='hills.out')

parser.add_argument('-o', '--output', type=str,
                    help='Output file for free energy landscape.')

parser.add_argument('--original', type=str,
                    help='Original free energy landscape for comparison.',
                    default='fes.dat')

args = parser.parse_args()

# Import original data if it exists
try:
    originalData = np.loadtxt(args.original)
except IOError:
    print("File {} does not exist, ignoring original result.".format(args.original))

# Import SSAGES metadynamics data
# xcenter, ycenter, sigma, sigma, height
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    try:
        SsagesData = np.loadtxt(args.input)
    except IOError:
        print("File {} does not exist. Exiting.".format(args.input))

    if not SsagesData.any():
        print("No data in {}! Exiting.".format(args.input))
        exit(1)

# Sum hills using gaussian function
X =  np.arange(-1.5, 2.0, 0.05)
Y =  np.arange(-1.5, 2.0, 0.05)
X, Y = np.meshgrid(X, Y)
Z = sum( [height*gaussian(X - xc, Y - yc, sigmax, sigmay) for it, xc, yc, sigmax, sigmay, height in SsagesData ] )

# Plot the inverse free energy landscape
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=True)
cbar = fig.colorbar(surf, shrink=0.7, aspect=10)

ax.set_title('Bias energy SSAGES')
ax.set_xlabel('x position')
ax.set_ylabel('y position')
cbar.set_label(r'$K_{b}T$')

# Plot the free energy landscape
fig2 = plt.figure()
ax2 = fig2.gca(projection='3d')
surf2 = ax2.plot_surface(X, Y, -Z, cmap=cm.coolwarm, linewidth=0, antialiased=True)
cbar2 = fig2.colorbar(surf2, shrink=0.7, aspect=10)

ax2.set_title('Free energy SSAGES')
ax2.set_xlabel('x position')
ax2.set_ylabel('y position')
cbar2.set_label(r'$K_{b}T$')

# Generate the free energy landscape from original data
# -- This part is missing yet

plt.show()
