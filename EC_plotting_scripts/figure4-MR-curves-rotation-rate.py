#import stability_curve as scc # import custom python file which computes the stability curve
import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'pyfbs'))
import data	#import file to read in data

# load in data for the EOS comparison plot:
number_of_files = 6
filenames = [None]*number_of_files

filenames[0] = "../output/ECstar_curve_rotation_rate_EOS_DD2_frequencyHz_0.0000000000.txt"
filenames[1] = "../output/ECstar_curve_rotation_rate_EOS_DD2_frequencyHz_714.0000000000.txt"
filenames[2] = "../output/ECstar_curve_rotation_rate_EOS_DD2_frequencyHz_1.0000000000.txt" # Keplerian rotation rate
filenames[3] = "../output/ECstar_curve_rotation_rate_EOS_APR_frequencyHz_0.0000000000.txt"
filenames[4] = "../output/ECstar_curve_rotation_rate_EOS_APR_frequencyHz_714.0000000000.txt"
filenames[5] = "../output/ECstar_curve_rotation_rate_EOS_APR_frequencyHz_1.0000000000.txt" # Keplerian rotation rate


#label_names = ["APR","DD2","FSG","KDE0v1","LNS"]
label_names = ["Pure DD2","$\Omega=714\,Hz$","$\Omega = \Omega_K$","Pure APR","$\Omega=714\,Hz$","$\Omega = \Omega_K$"]

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
plt.xlim([9., 14.0])
plt.ylim([0., 2.5])
plt.yscale("linear")

plt.xlabel("Radius [km]", fontsize = tickFontSize)
plt.ylabel("Gravitational Mass [M$_\odot$]", fontsize = tickFontSize)
#plt.tick_params(axis='both', which='major', labelsize=tickFontSize)
plt.grid(alpha=0.2, linestyle="--")

#plt.scatter(13.27,1.7, label = "", marker = "*", s=200, c = "orange",zorder=2.5)

linecolors = ["black","b","b", "black","r","r"] #"#DF5327","#DF5327"]
linestyles = ["-",":","-.","-",":","-."]
for i in range(len(filenames)):
	plt.plot(radii[i], masses[i], label = label_names[i], c=linecolors[i],ls=linestyles[i], lw=1.6)

plt.legend(loc="lower left", fontsize = 13)

# add observational constraints
xarr = [0.0, 40.]
yarr1 = [2.18, 2.18]
yarr2 = [2.52, 2.52]
plt.fill_between(xarr, yarr1, yarr2, alpha = 0.55, color="orange")
plt.text(9.1, 2.4, r"PSR J0952−0607", fontsize=8)

plt.text(10.7, 1.7, "Preliminary", fontsize = 24, rotation='horizontal', c="gray", alpha=0.5)


myfigname = "Figure4-MR_for_different_rotation_rates.pdf"
plt.savefig(myfigname, dpi=400, bbox_inches='tight')
#plt.show()
