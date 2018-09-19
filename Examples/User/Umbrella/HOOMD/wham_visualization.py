"""
Author: Bradley Dice (bdice@umich.edu)
Description: The SSAGES Umbrella Sampling example for HOOMD-blue runs a
             simulation of a butane molecule. This script is a visualization
             tool that reads output files from SSAGES and wham to plot the free
             energy as a function of the dihedral angle of the carbon backbone.
"""

import numpy as np
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

starting_row = 0
rad_to_deg = 180/np.pi
angle_ticks = np.linspace(-180, 180, 7)

with open('multiwalker_umbrella_input.json', 'r') as multi_input:
    nwalkers = json.load(multi_input)['walkers']

# Plot collective variables for all trajectories
for i in range(nwalkers):
    data = np.genfromtxt('{:03d}_umbrella.dat'.format(i))
    plt.plot(data[starting_row:, 0], rad_to_deg*data[starting_row:, 1])

plt.title('Collective Variable')
plt.xlabel('Timesteps')
plt.ylabel('Dihedral C-C-C-C angle (deg)')
plt.xlim((data[starting_row, 0], data[-1, 0]))
plt.ylim((-180, 180))
plt.yticks(angle_ticks)
plt.gcf().savefig('cv_vs_time.png', facecolor='white')
plt.close()

# Plot individual trajectory histograms
bin_totals = None
for i in range(nwalkers):
    sample_data = np.genfromtxt(
        '{:03d}_umbrella.dat'.format(i))[starting_row:, 1]
    binned_data, bins = np.histogram(
        sample_data, bins=100, range=(-np.pi, np.pi))
    bins = rad_to_deg * (bins[1:] + bins[:-1]) / 2
    if bin_totals is None:
        bin_totals = binned_data
    else:
        bin_totals += binned_data
    plt.plot(bins, binned_data, label=i)

plt.title('Individual trajectory histograms')
plt.xlabel('Dihedral C-C-C-C angle (deg)')
plt.ylabel('Bin Counts')
plt.xlim((-180, 180))
plt.xticks(angle_ticks)
plt.gcf().savefig('histogram_trajectories.png', facecolor='white')
plt.close()

# Plot combined histogram
plt.plot(bins, bin_totals)
plt.title('Combined histogram')
plt.xlabel('Dihedral C-C-C-C angle (deg)')
plt.ylabel('Bin Counts')
plt.xlim((-180, 180))
plt.ylim((0, np.max(bin_totals)))
plt.xticks(angle_ticks)
plt.gcf().savefig('histogram_combined.png', facecolor='white')
plt.close()

# Make free energy plot from output of ./wham_analysis.sh
data = np.genfromtxt('wham_freefile')
plt.plot(rad_to_deg * data[:, 0], data[:, 1])
plt.title('Butane PMF')
plt.xlabel('Dihedral C-C-C-C angle (deg)')
plt.ylabel('Free energy (kcal/mol)')
plt.xlim((-180, 180))
plt.xticks(angle_ticks)
plt.savefig('wham_free_energy.png', facecolor='white')
plt.close()
