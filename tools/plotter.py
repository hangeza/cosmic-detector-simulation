import matplotlib.pyplot as plt
import numpy as np
import sys

from matplotlib import colors
from matplotlib.ticker import PercentFormatter

# plotter.py
#
# plot data from a file path specified in command line arguments
# INPUT FILETYPE MUST BE .dat or .hist !!
##############################################################

# get these command line arguments
file = open(sys.argv[1], 'r')
filename = sys.argv[1]

# determine the file name and the file type
temp = filename.split('.')
filename = temp[0]
filetype = temp[1]

# this function is for .dat files and generates a dot plot of the given values
def data_plots(file, filename):
	# create a distribution to hold all of the data points
	distr = np.empty(shape = (181, 3))	# 181 rows (-90 deg. to +90) and 3 columns (angle, value, error)

	line_no = 0
	max_rate = -1.0

	# iterate through every line in the file and add each data point to the distribution
	for line in file:
		if line.startswith('#'):	# ignore headers
			continue

		cols = line.split()		# split lines and add relevant values to distr
		distr[line_no][0] = cols[0]
		distr[line_no][1] = cols[2]
		distr[line_no][2] = cols[3]
		line_no += 1

		if float(cols[2]) > max_rate:	# determine the max rate for fitting the cos^2 curve later
			max_rate = float(cols[2])


	fig, axs = plt.subplots(1, 1, tight_layout=True)

	# create cos^2 curve to be plotted on the same graph
	x = np.linspace(-2.0, 2.0, 100)
	y = max_rate * ( np.cos(x) ** 2)

	# plot each point from the distribution
	for point in distr:
		axs.errorbar(point[0], point[1], point[2], fmt='o', linewidth=2, capsize=6)

	# plot cos^2
	axs.plot(x, y, linewidth=2)

	# save figure as pdf w/ appropriate name and display it
	plt.savefig("simulated_" + filename.split('_')[-1] + "_vs_cos2.pdf", dpi=400)
	plt.show()

# this program is for .hist files and generates a histogram of given bins and entries
def histogram(file, filename):
	num_of_lines = 0
	bin_center = []
	entries = []

	# iterate over every line and add bin center and entry values to respective lists
	for line in file:
		if line.startswith('#'):		# ignore headers
			continue

		num_of_lines += 1
		cols = line.split()

		bin_center.append( float(cols[1]) )	# add bin center value
		entries.append( int(cols[2]) )		# add entry value

	# by default, number of bins is given by the .hist file
	# this can be changed as a command line argument after the file name
	n_bins = int(num_of_lines)
	if len(sys.argv) > 2:
		if sys.argv[2].isdigit():
			n_bins = int(sys.argv[2])
		else:
			print("Invalid argument [3] -- number of bins should be integer")

	# plot histogram using bin centers as values and weighting them by number of entries
	plt.hist(bin_center, n_bins, weights=entries)

	# save figure as pdf w/ appropriate name and display it
	plt.savefig("simulated_hist_" + filename + ".pdf", dpi=400)
	plt.show()


# input files must be either .dat or .hist
# program will do nothing otherwise (there is no functionality for other filetypes)
if filetype == "dat":
	data_plots(file, filename)

if filetype == "hist":
	histogram(file, filename)
