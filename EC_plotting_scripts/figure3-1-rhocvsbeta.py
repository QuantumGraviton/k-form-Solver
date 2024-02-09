import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'pyfbs'))
import data	#import file to read in data

# load in data:
number_of_files = 8 #
filenames = [None]*number_of_files
df = [None]*number_of_files
indices = [None]*number_of_files

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


linew = 1.6
# unit conversion:
code_to_km = 1.476625061
code_to_MeV_fm3 = 346455.096775659
code_to_GeV_fm3 = code_to_MeV_fm3/1.0e3
code_to_rho_nuc = 1. / (0.16 * 2.886376934e-6 * 939.565379)
kappa = 8*np.pi


nrows = 2
ncols = 1
fig = plt.figure(figsize=(6.5,5))
gs = fig.add_gridspec(nrows, ncols, hspace=0.15, wspace=0)
axs = gs.subplots()#sharex='col', sharey='row')

###################################################################################################################################################


axs[0].plot(df[0][:,indices_star['beta']], df[0][:,indices_star['rho_0']]*code_to_rho_nuc, lw=linew, ls = "-", color ="green", label="$M_{rest}=0.8 M_\odot$")
axs[0].plot(df[1][:,indices_star['beta']], df[1][:,indices_star['rho_0']]*code_to_rho_nuc, lw=linew, ls = "-", color ="blue", label="$1.0 M_\odot$")
axs[0].plot(df[2][:,indices_star['beta']], df[2][:,indices_star['rho_0']]*code_to_rho_nuc, lw=linew, ls = "-", color ="black", label="$1.4 M_\odot$")
axs[0].plot(df[3][:,indices_star['beta']], df[3][:,indices_star['rho_0']]*code_to_rho_nuc, lw=linew, ls = "-", color ="red", label="$2.0 M_\odot$")

axs[1].plot(df[4][:,indices_star['beta']], df[4][:,indices_star['rho_0']]*code_to_rho_nuc, lw=linew, ls = "-", color ="green", label="$0.8 M_\odot$")
axs[1].plot(df[5][:,indices_star['beta']], df[5][:,indices_star['rho_0']]*code_to_rho_nuc, lw=linew, ls = "-", color ="blue", label="$1.0 M_\odot$")
axs[1].plot(df[6][:,indices_star['beta']], df[6][:,indices_star['rho_0']]*code_to_rho_nuc, lw=linew, ls = "-", color ="black", label="$1.4 M_\odot$")
axs[1].plot(df[7][:,indices_star['beta']], df[7][:,indices_star['rho_0']]*code_to_rho_nuc, lw=linew, ls = "-", color ="red", label="$2.0 M_\odot$")
#axs[1].plot(df[0][:,indices_star['beta']], df[0][:,indices_star['R_NS']], lw=linew, ls = "-", color ="green", label="$\gamma = 2$")
#axs[1].plot(df[0][:,indices_star['r']]*code_to_km, kappas2 / df[0][:,indices_star['e']] * 100, lw=linew, ls = "-", color ="green", label="$\Omega = \Omega_{Kep}$")


###################################################################################################################################################


# set scale, plotlim, legend ...
axs[0].set_yscale("linear") # "log" or "linear"
axs[1].set_yscale("linear") # "log" or "linear"

axs[0].set_ylim((1,3.5))
axs[1].set_ylim((1,3.5))

axs[0].set_xlim((0.,100.))
axs[1].set_xlim((0.,1e6))
#axs[1].set_yticks(np.geomspace(ylim2[0],ylim2[1], 51))

axs[0].grid(alpha=0.15, linestyle="--", c="black")
axs[1].grid(alpha=0.15, linestyle="--", c="black")

# ----- add legends:
legend_ax0 = axs[0].legend(loc="upper right", fontsize = 9,ncol=2)
axs[0].add_artist(legend_ax0)
legend_line0, = axs[0].plot([-2.], [-2.], label ="$\gamma = 2$")
axs[0].legend(handles=[legend_line0], handlelength=0, handletextpad=0, loc="upper left", fontsize = 12)

legend_ax1 = axs[1].legend(loc="upper right", fontsize = 9)
axs[1].add_artist(legend_ax1)
legend_line1, = axs[1].plot([-2.], [-2.], label ="$\gamma = 3$")
axs[1].legend(handles=[legend_line1], handlelength=0, handletextpad=0, loc="upper left", fontsize = 12)

# shared x-axis label of the plot
plt.xlabel("$\\beta$ [$(M_\odot / R_{g\odot}^3)^{1-\\gamma}$]", fontsize = 18)
# shared y-axis label
axs[0].set_ylabel("$\\rho_{c}$ [$\\rho_{sat}$]", fontsize = 18)
axs[1].set_ylabel("$\\rho_{c}$ [$\\rho_{sat}$]", fontsize = 18)

# prevent labels on the "inside" between the subplots:
#for ax in axs:
#	ax.label_outer()
		
# plot labeling etc.:
# text in the plot:
axs[0].text(35., 2.25, "Preliminary", fontsize = 24, rotation='horizontal', c="gray", alpha=0.5)
#axs[1].text(35., 13.25, "maybe try other eos or $\gamma$", fontsize = 24, rotation='horizontal', c="gray", alpha=0.5)
#axs[0].text(6.6, 0.000055, "APR EOS", fontsize = 12, rotation='horizontal', c="blue")

#plt.title("relative error for $M_{tot}$ and $N_{B}$ inside stab region", fontsize = 18)
#plt.xlabel("Radius $r$ [km]", fontsize = 18)
#plt.ylabel("$e(r)$ [$M_\odot / R_g^3$]", fontsize = 18)
#plt.ylabel("Field $\phi(r)$ [$M_p$]", fontsize = 18)
#plt.grid(alpha=0.2, linestyle="--")#, c="black")

#plt.yscale("linear") # "log" or "linear"
#plt.ylim((1e-8,0.0025))
#plt.xlim((0.,15.))
#plt.legend(loc="upper right", fontsize = 11) #loc="lower left" or "upper right"

figname = "Figure3-1-Mvsbeta-Mrest.pdf"
plt.savefig(figname, dpi=400, bbox_inches='tight')
print("Saved figure as: " + figname)
#plt.show()
