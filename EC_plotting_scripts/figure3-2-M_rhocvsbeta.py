import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'pyfbs'))
import data	#import file to read in data

def strip_curves(masses, radii, rhos, betas):
	# remove the unstable configurations after the maximum density was reached
	rhomax = 0.0
	maxindex = 0
	for i in range(len(rhos)):
		maxindex = i
		if rhos[i] > rhomax:
			rhomax = rhos[i]
		else:
			break
	masses_out = np.zeros(maxindex)
	radii_out = np.zeros(maxindex)
	rhos_out = np.zeros(maxindex)
	betas_out = np.zeros(maxindex)
	for j in range(maxindex):
		masses_out[j] = masses[j]
		radii_out[j] = radii[j]
		rhos_out[j] = rhos[j]
		betas_out[j] = betas[j]
	return masses_out, radii_out, rhos_out, betas_out


# load in data:
number_of_files = 8 #
filenames = [None]*number_of_files
df = [None]*number_of_files
indices = [None]*number_of_files
masses = [None]*number_of_files
radii = [None]*number_of_files
rhos = [None]*number_of_files
betas = [None]*number_of_files

filenames[0] = "../output/ECstar_curve_const_mass_M_rest_with_different_beta_EOS_DD2_nsmass_0.8000000000_gamma_2.0000000000.txt"
filenames[1] = "../output/ECstar_curve_const_mass_M_rest_with_different_beta_EOS_DD2_nsmass_1.0000000000_gamma_2.0000000000.txt"
filenames[2] = "../output/ECstar_curve_const_mass_M_rest_with_different_beta_EOS_DD2_nsmass_1.4000000000_gamma_2.0000000000.txt"
filenames[3] = "../output/ECstar_curve_const_mass_M_rest_with_different_beta_EOS_DD2_nsmass_2.0000000000_gamma_2.0000000000.txt"
#ECstar_curve_const_mass_M_rest_with_different_beta_EOS_DD2_nsmass_2.0000000000_gamma_2.0000000000
filenames[4] = "../output/ECstar_curve_const_mass_M_rest_with_different_beta_EOS_DD2_nsmass_0.8000000000_gamma_3.0000000000.txt"
filenames[5] = "../output/ECstar_curve_const_mass_M_rest_with_different_beta_EOS_DD2_nsmass_1.0000000000_gamma_3.0000000000.txt"
filenames[6] = "../output/ECstar_curve_const_mass_M_rest_with_different_beta_EOS_DD2_nsmass_1.4000000000_gamma_3.0000000000.txt"
filenames[7] = "../output/ECstar_curve_const_mass_M_rest_with_different_beta_EOS_DD2_nsmass_2.0000000000_gamma_3.0000000000.txt"

for i in range(len(filenames)):
	df[i], indices[i] = data.load_file(filenames[i])
	#print(indices[i])

indices_star = indices[0]
print(indices_star)

for i in range(len(filenames)):
	masses[i] = df[i][:,indices[i]['M_T']]
	radii[i]  = df[i][:,indices[i]['R_NS']]
	rhos[i]   = df[i][:,indices[i]['rho_0']]
	betas[i]  = df[i][:,indices[i]['beta']]
	masses[i], radii[i], rhos[i], betas[i] = strip_curves(masses[i], radii[i], rhos[i], betas[i])
	

linew = 1.6
# unit conversion:
code_to_km = 1.476625061
code_to_MeV_fm3 = 346455.096775659
code_to_GeV_fm3 = code_to_MeV_fm3/1.0e3
code_to_rho_nuc = 1. / (0.16 * 2.886376934e-6 * 939.565379)
kappa = 8*np.pi


nrows = 6
ncols = 1
fig = plt.figure(figsize=(6.5,16.5))
gs = fig.add_gridspec(nrows, ncols, hspace=0.15, wspace=0)
axs = gs.subplots()#sharex='col', sharey='row')

###################################################################################################################################################
c1 = "forestgreen"#"green"
c2 = "mediumblue"#"mediumslateblue"#"blue"
c3 = "black"
c4 = "crimson"

# upper plots with beta=2:
axs[0].plot(betas[0], radii[0], lw=linew, ls = "-", color =c1, label="$M_{rest}=0.8 M_\odot$")
axs[0].plot(betas[1], radii[1], lw=linew, ls = "-", color =c2, label="$1.0 M_\odot$")
axs[0].plot(betas[2], radii[2], lw=linew, ls = "-", color =c3, label="$1.4 M_\odot$")
axs[0].plot(betas[3], radii[3], lw=linew, ls = "-", color =c4, label="$2.0 M_\odot$")

axs[1].plot(betas[0], rhos[0]*code_to_rho_nuc, lw=linew, ls = "-", color =c1, label="$M_{rest}=0.8 M_\odot$")
axs[1].plot(betas[1], rhos[1]*code_to_rho_nuc, lw=linew, ls = "-", color =c2, label="$1.0 M_\odot$")
axs[1].plot(betas[2], rhos[2]*code_to_rho_nuc, lw=linew, ls = "-", color =c3, label="$1.4 M_\odot$")
axs[1].plot(betas[3], rhos[3]*code_to_rho_nuc, lw=linew, ls = "-", color =c4, label="$2.0 M_\odot$")

axs[2].plot(betas[0], masses[0]/masses[0][0]*100 -100, lw=linew, ls = "-", color =c1, label="$M_{rest}=0.8 M_\odot$")
axs[2].plot(betas[1], masses[1]/masses[1][0]*100 -100, lw=linew, ls = "-", color =c2, label="$1.0 M_\odot$")
axs[2].plot(betas[2], masses[2]/masses[2][0]*100 -100, lw=linew, ls = "-", color =c3, label="$1.4 M_\odot$")
axs[2].plot(betas[3], masses[3]/masses[3][0]*100 -100, lw=linew, ls = "-", color =c4, label="$2.0 M_\odot$")

# plots with beta=3:
axs[3].plot(betas[4], radii[4], lw=linew, ls = "-", color =c1, label="$M_{rest}=0.8 M_\odot$")
axs[3].plot(betas[5], radii[5], lw=linew, ls = "-", color =c2, label="$1.0 M_\odot$")
axs[3].plot(betas[6], radii[6], lw=linew, ls = "-", color =c3, label="$1.4 M_\odot$")
axs[3].plot(betas[7], radii[7], lw=linew, ls = "-", color =c4, label="$2.0 M_\odot$")

axs[4].plot(betas[4], rhos[4]*code_to_rho_nuc, lw=linew, ls = "-", color =c1, label="$M_{rest}=0.8 M_\odot$")
axs[4].plot(betas[5], rhos[5]*code_to_rho_nuc, lw=linew, ls = "-", color =c2, label="$1.0 M_\odot$")
axs[4].plot(betas[6], rhos[6]*code_to_rho_nuc, lw=linew, ls = "-", color =c3, label="$1.4 M_\odot$")
axs[4].plot(betas[7], rhos[7]*code_to_rho_nuc, lw=linew, ls = "-", color =c4, label="$2.0 M_\odot$")

axs[5].plot(betas[4], masses[4]/masses[4][0]*100 -100, lw=linew, ls = "-", color =c1, label="$M_{rest}=0.8 M_\odot$")
axs[5].plot(betas[5], masses[5]/masses[5][0]*100 -100, lw=linew, ls = "-", color =c2, label="$1.0 M_\odot$")
axs[5].plot(betas[6], masses[6]/masses[6][0]*100 -100, lw=linew, ls = "-", color =c3, label="$1.4 M_\odot$")
axs[5].plot(betas[7], masses[7]/masses[7][0]*100 -100, lw=linew, ls = "-", color =c4, label="$2.0 M_\odot$")

#######
# add small circle at the pointof critical density/radius/mass for the relevant configurations:
markerstyle = "*" # "o"
markersize = 30
axs[0].scatter(betas[2][-1], radii[2][-1], s=markersize, c=c3, marker=markerstyle)
axs[0].scatter(betas[3][-1], radii[3][-1], s=markersize, c=c4, marker=markerstyle)
axs[1].scatter(betas[2][-1], rhos[2][-1]*code_to_rho_nuc, s=markersize, c=c3, marker=markerstyle)
axs[1].scatter(betas[3][-1], rhos[3][-1]*code_to_rho_nuc, s=markersize, c=c4, marker=markerstyle)
axs[2].scatter(betas[2][-1], masses[2][-1]/masses[2][0]*100 -100, s=markersize, c=c3, marker=markerstyle)
axs[2].scatter(betas[3][-1], masses[3][-1]/masses[3][0]*100 -100, s=markersize, c=c4, marker=markerstyle)

axs[3].scatter(betas[6][-1], radii[6][-1], s=markersize, c=c3, marker=markerstyle)
axs[3].scatter(betas[7][-1], radii[7][-1], s=markersize, c=c4, marker=markerstyle)
axs[4].scatter(betas[6][-1], rhos[6][-1]*code_to_rho_nuc, s=markersize, c=c3, marker=markerstyle)
axs[4].scatter(betas[7][-1], rhos[7][-1]*code_to_rho_nuc, s=markersize, c=c4, marker=markerstyle)
axs[5].scatter(betas[6][-1], masses[6][-1]/masses[6][0]*100 -100, s=markersize, c=c3, marker=markerstyle)
axs[5].scatter(betas[7][-1], masses[7][-1]/masses[7][0]*100 -100, s=markersize, c=c4, marker=markerstyle)

###################################################################################################################################################


# set scale, plotlim, legend ...
axs[0].set_yscale("linear") # "log" or "linear"
axs[1].set_yscale("linear") # "log" or "linear"
axs[2].set_yscale("linear") # "log" or "linear"
axs[3].set_yscale("linear") # "log" or "linear"
axs[4].set_yscale("linear") # "log" or "linear"
axs[5].set_yscale("linear") # "log" or "linear"

axs[0].set_ylim((12.8,13.4))
axs[1].set_ylim((1,3.5))
axs[2].set_ylim((-1.5,1.0))
axs[3].set_ylim((12.8,13.4))
axs[4].set_ylim((1,3.5))
axs[5].set_ylim((-1.5,1.0))


axs[0].set_xlim((0.,100.))
axs[1].set_xlim((0.,100.))
axs[2].set_xlim((0.,100.))
axs[3].set_xlim((0.,0.6e6))
axs[4].set_xlim((0.,0.6e6))
axs[5].set_xlim((0.,0.6e6))
#axs[2].ticklabel_format(axis="x",scilimits=(0,0),useOffset=True)
#axs[3].ticklabel_format(axis="x",scilimits=(0,0),useMathText=True)
#axs[1].set_yticks(np.geomspace(ylim2[0],ylim2[1], 51))

axs[0].grid(alpha=0.15, linestyle="--", c="black")
axs[1].grid(alpha=0.15, linestyle="--", c="black")
axs[2].grid(alpha=0.15, linestyle="--", c="black")
axs[3].grid(alpha=0.15, linestyle="--", c="black")
axs[4].grid(alpha=0.15, linestyle="--", c="black")
axs[5].grid(alpha=0.15, linestyle="--", c="black")

# ----- add legends:
legend_ax0 = axs[0].legend(loc="upper right", fontsize = 9,ncol=2)
axs[0].add_artist(legend_ax0)
legend_line0, = axs[0].plot([-2.], [-2.], label ="$\gamma = 2$", c="white")
axs[0].legend(handles=[legend_line0], handlelength=0, handletextpad=0, loc="upper left", fontsize = 12)

legend_ax1 = axs[1].legend(loc="upper right", fontsize = 9,ncol=2)
axs[1].add_artist(legend_ax1)
legend_line1, = axs[1].plot([-2.], [-2.], label ="$\gamma = 2$", c="white")
axs[1].legend(handles=[legend_line1], handlelength=0, handletextpad=0, loc="upper left", fontsize = 12)

legend_ax2 = axs[2].legend(loc="upper right", fontsize = 9,ncol=2)
axs[2].add_artist(legend_ax2)
legend_line2, = axs[2].plot([-2.], [-2.], label ="$\gamma = 2$", c="white")
axs[2].legend(handles=[legend_line2], handlelength=0, handletextpad=0, loc="upper left", fontsize = 12)

legend_ax3 = axs[3].legend(loc="upper right", fontsize = 9,ncol=2)
axs[3].add_artist(legend_ax3)
legend_line3, = axs[3].plot([-2.], [-2.], label ="$\gamma = 3$", c="white")
axs[3].legend(handles=[legend_line3], handlelength=0, handletextpad=0, loc="upper left", fontsize = 12)

legend_ax4 = axs[4].legend(loc="upper right", fontsize = 9,ncol=2)
axs[4].add_artist(legend_ax4)
legend_line4, = axs[4].plot([-2.], [-2.], label ="$\gamma = 3$", c="white")
axs[4].legend(handles=[legend_line4], handlelength=0, handletextpad=0, loc="upper left", fontsize = 12)

legend_ax5 = axs[5].legend(loc="upper right", fontsize = 9,ncol=2)
axs[5].add_artist(legend_ax5)
legend_line5, = axs[5].plot([-2.], [-2.], label ="$\gamma = 3$", c="white")
axs[5].legend(handles=[legend_line5], handlelength=0, handletextpad=0, loc="upper left", fontsize = 12)


# shared x-axis label of the plot
plt.xlabel("$\\beta$ [$(M_\odot / R_{g\odot}^3)^{1-\\gamma}] / [\kappa]$", fontsize = 18)
# shared y-axis label
axs[0].set_ylabel("$R_{NS}$ [km]", fontsize = 18)
axs[1].set_ylabel("$\\rho_{c}$ [$\\rho_{sat}$]", fontsize = 18)
axs[2].set_ylabel("$\Delta M_{grav,rel}$ [$\\%$]", fontsize = 18)
axs[3].set_ylabel("$R_{NS}$ [km]", fontsize = 18)
axs[4].set_ylabel("$\\rho_{c}$ [$\\rho_{sat}$]", fontsize = 18)
axs[5].set_ylabel("$\Delta M_{grav,rel}$ [$\\%$]", fontsize = 18)

axs[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axs[4].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# prevent labels on the "inside" between the subplots:
#for ax in axs:
#	ax.label_outer()
		
# plot labeling etc.:
# text in the plot:
axs[1].text(35., 2.25, "Preliminary", fontsize = 24, rotation='horizontal', c="gray", alpha=0.5)
#axs[1].text(35., 13.25, "maybe try other eos or $\gamma$", fontsize = 24, rotation='horizontal', c="gray", alpha=0.5)
#axs[0].text(6.6, 0.000055, "APR EOS", fontsize = 12, rotation='horizontal', c="blue")


figname = "Figure3-2-Mrhocvsbeta-mt.pdf"
plt.savefig(figname, dpi=400, bbox_inches='tight')
print("Saved figure as: " + figname)
#plt.show()
