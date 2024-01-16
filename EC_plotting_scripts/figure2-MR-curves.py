#import stability_curve as scc # import custom python file which computes the stability curve
import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'pyfbs'))
import data	#import file to read in data

# load in data for the EOS comparison plot:
number_of_files = 4
filenames = [None]*number_of_files

filenames[0] = "../output/ECstar_curve_EOS_DD2_beta_0.0000000000_gamma_2.0000000000.txt"
filenames[1] = "../output/ECstar_curve_EOS_DD2_beta_10.0000000000_gamma_2.0000000000.txt"
filenames[2] = "../output/ECstar_curve_EOS_DD2_beta_20.0000000000_gamma_2.0000000000.txt"
filenames[3] = "../output/ECstar_curve_EOS_DD2_beta_100.0000000000_gamma_2.0000000000.txt"


#label_names = ["APR","DD2","FSG","KDE0v1","LNS"]
label_names = ["Pure DD2 NS","$\\beta = 10$","$\\beta = 20$","$\\beta = 100$"]

# read in data:
df = [None]*number_of_files
indices = [None]*number_of_files
masses = [None]*number_of_files
radii = [None]*number_of_files

for i in range(len(filenames)):
	df[i], indices[i] = data.load_file(filenames[i]) # load the files

for i in range(len(filenames)):
	masses[i] = df[i][:,indices[i]['M_T']]
	radii[i] = df[i][:,indices[i]['R_NS']]


tickFontSize = 18

fig = plt.figure()
plt.xlim([7., 19.0])
plt.ylim([0., 3.5])
plt.yscale("linear")

plt.xlabel(r"Radius [km]", fontsize = tickFontSize)
plt.ylabel(r"Total Gravitational Mass [M$_\odot$]", fontsize = tickFontSize)
#plt.tick_params(axis='both', which='major', labelsize=tickFontSize)
plt.grid(alpha=0.2, linestyle="--")

#plt.scatter(13.27,1.7, label = "", marker = "*", s=200, c = "orange",zorder=2.5)

linecolors = ["black","b","b","#DF5327","#DF5327"]
linestyles = ["-","-.","-.","--","--"]
for i in range(len(filenames)):
	plt.plot(radii[i], masses[i], label = label_names[i], c=linecolors[i],ls=linestyles[i], lw=1.6)

#plt.text(13.7, 1.7, r"1.7M$_\odot$, 13.27 km", fontsize=14)
plt.legend(loc="upper right", fontsize = 13)

# add DM fractions:
#plt.text(11.9, 2.2, r"$0\%$", fontsize=12, color=linecolors[0],rotation=-20)
#plt.text(11.1, 1.8, r"$10\%$", fontsize=12, color=linecolors[1],rotation=-20)
#plt.text(10.53, 1.43, r"$20\%$", fontsize=12, color=linecolors[1],rotation=-20)
#plt.text(11.9, 2.6, r"$60\%$", fontsize=12, color=linecolors[3],rotation=-20)
#plt.text(11.9, 3.0, r"$75\%$", fontsize=12, color=linecolors[3],rotation=-20)

# add observational constraints
xarr = [0.0, 40.]
yarr1 = [2.18, 2.18]
yarr2 = [2.52, 2.52]
plt.fill_between(xarr, yarr1, yarr2, alpha = 0.55, color="orange")
plt.text(7.4, 2.4, r"PSR J0952−0607", fontsize=8)


myfigname = "Figure2-MR_plot_test.pdf"
plt.savefig(myfigname, dpi=400, bbox_inches='tight')
#plt.show()
