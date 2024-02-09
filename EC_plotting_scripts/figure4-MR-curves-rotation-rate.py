#import stability_curve as scc # import custom python file which computes the stability curve
import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'pyfbs'))
import data	#import file to read in data

def strip_mr_curves(masses, radii):
	# remove the unstable configurations after the maximum mass was reached
	mmax = 0.0
	maxindex = 0
	for i in range(len(masses)):
		maxindex = i
		if masses[i] > mmax:
			mmax = masses[i]
		else:
			break
	masses_out = np.zeros(maxindex)
	radii_out = np.zeros(maxindex)
	for j in range(maxindex):
		masses_out[j] = masses[j]
		radii_out[j] = radii[j]
	return masses_out, radii_out


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
	masses[i], radii[i] = strip_mr_curves(masses[i], radii[i])


tickFontSize = 18

fig = plt.figure()
plt.xlim([10., 14.0])
plt.ylim([0., 2.5])
plt.yscale("linear")

plt.xlabel("Radius [km]", fontsize = tickFontSize)
plt.ylabel("Gravitational Mass [M$_\odot$]", fontsize = tickFontSize)
#plt.tick_params(axis='both', which='major', labelsize=tickFontSize)
plt.grid(alpha=0.2, linestyle="--")

#plt.scatter(13.27,1.7, label = "", marker = "*", s=200, c = "orange",zorder=2.5)

c1 = "forestgreen"#"green"
c2 = "mediumblue"#"mediumslateblue"#"blue"
linecolors = [c1,c1,c1,c2,c2,c2]#"#DF5327","#DF5327"]
linestyles = ["-",":","-.","-",":","-."]
for i in range(len(filenames)):
	plt.plot(radii[i], masses[i], label = label_names[i], c=linecolors[i],ls=linestyles[i], lw=1.6)

plt.legend(loc="lower left", fontsize = 13)

# add observational constraints
xarr = [0.0, 40.]
yarr1 = [2.18, 2.18]
yarr2 = [2.52, 2.52]
plt.fill_between(xarr, yarr1, yarr2, alpha = 0.3, color="sandybrown")
plt.text(10.1, 2.4, r"PSR J0952−0607", fontsize=8,color="sandybrown")

# add NICER results here in form of 1-2-sigma clouds
number_of_NICER_files = 5
filenames_NICER = [None]*number_of_NICER_files
df_NICER = [None]*number_of_NICER_files
indices_NICER = [None]*number_of_NICER_files

filenames_NICER[0] = "Posteriors_NS-measurments/A_NICER_VIEW_OF_PSR_J0740+6620/NICER_x_XMM_J0740__XPSI_STU_NSX_FIH__contour_radius_mass_68.txt" # 1 sigma
filenames_NICER[1] = "Posteriors_NS-measurments/A_NICER_VIEW_OF_PSR_J0740+6620/NICER_x_XMM_J0740__XPSI_STU_NSX_FIH__contour_radius_mass_95.txt" # 2 sigma
filenames_NICER[2] = "Posteriors_NS-measurments/updated_analyses_PSRJ0030_up_to_2018_NICER_data/J0030_PDTU_NxXMM_68.txt" # 1 sigma
filenames_NICER[3] = "Posteriors_NS-measurments/updated_analyses_PSRJ0030_up_to_2018_NICER_data/J0030_PDTU_NxXMM_95.txt" # 2 sigma
filenames_NICER[4] = "Posteriors_NS-measurments/updated_analyses_PSRJ0030_up_to_2018_NICER_data/J0030_PDTU_NxXMM_99.txt" # 2 sigma

for j in range(len(filenames_NICER)):
	df_NICER[j], indices_NICER[j] = data.load_file(filenames_NICER[j]) # load the files

col_nicer1 = "orangered"
col_nicer2 = "darkorange"
plt.fill(df_NICER[0][:,indices_NICER[0]['R_NS']], df_NICER[0][:,indices_NICER[0]['M_T']], alpha=0.15, c=col_nicer1, ls =":")
plt.fill(df_NICER[1][:,indices_NICER[1]['R_NS']], df_NICER[1][:,indices_NICER[1]['M_T']], alpha=0.15, c=col_nicer1, ls ="-.")
plt.fill(df_NICER[2][:,indices_NICER[2]['R_NS']], df_NICER[2][:,indices_NICER[2]['M_T']], alpha=0.10, c=col_nicer2, ls =":")
plt.fill(df_NICER[3][:,indices_NICER[3]['R_NS']], df_NICER[3][:,indices_NICER[3]['M_T']], alpha=0.10, c=col_nicer2, ls ="-.")
plt.fill(df_NICER[4][:,indices_NICER[4]['R_NS']], df_NICER[4][:,indices_NICER[4]['M_T']], alpha=0.10, c=col_nicer2, ls ="-")
plt.text(11.7, 2.05, "PSR J0740+6620", fontsize = 8, rotation='horizontal', c=col_nicer1)
plt.text(12.2, 1.4, "PSR J0030+0451 ", fontsize = 8, rotation='horizontal', c=col_nicer2)


plt.text(10.7, 1.7, "Preliminary", fontsize = 26, rotation='horizontal', c="gray", alpha=0.5)


myfigname = "Figure4-MR_for_different_rotation_rates-2.pdf"
plt.savefig(myfigname, dpi=400, bbox_inches='tight')
#plt.show()
