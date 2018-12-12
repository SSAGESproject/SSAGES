#! /usr/bin/python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import warnings
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# Define 2d-gaussian, centered at (0,0), scale = 1
def gaussian(x,y, sigmax, sigmay):
    return np.exp( -( (x**2 / (2*sigmax**2)) + (y**2 / (2*sigmay**2) ) ) )

def original(X, Y):
    heights = [-18.4631620773, -39.8889505555, 9.8889505555]
    centers = [ -0.9790631338,   0.9790631338, 0           ]
    sigmas  = [  0.3412727171,   0.2420191493, 0.5         ]

    Z = sum( [heights[i]*gaussian(X - centers[i], Y - centers[i], sigmas[i], sigmas[i])
                for i in range(3)] )
    Z -= np.max(Z)

    return Z

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

args = parser.parse_args()

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
Ztrue = original(X, Y)
fig3 = plt.figure()
ax3 = fig3.gca(projection='3d')
surf3 = ax3.plot_surface(X, Y, Ztrue, cmap=cm.coolwarm, linewidth=0, antialiased=True)
cbar3 = fig3.colorbar(surf3, shrink=0.7, aspect=10)

ax3.set_title('Original free energy')
ax3.set_xlabel('x position')
ax3.set_ylabel('y position')
cbar3.set_label(r'$K_{b}T$')

# Show difference plot
Zdiff = Ztrue + Z # Z is the bias potential Ztrue the free energy
fig4 = plt.figure()
ax4 = fig4.gca(projection='3d')
surf4 = ax4.plot_surface(X, Y, Zdiff, cmap=cm.coolwarm, linewidth=0, antialiased=True)
cbar4 = fig4.colorbar(surf4, shrink=0.7, aspect=10)

ax4.set_title('Difference plot:\nOriginal free energy - SSAGES free energy')
ax4.set_xlabel('x position')
ax4.set_ylabel('y position')
cbar4.set_label(r'$K_{b}T$')

plt.show()
