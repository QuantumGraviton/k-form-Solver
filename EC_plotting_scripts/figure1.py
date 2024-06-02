import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'pyfbs'))
import data	#import file to read in data

# load in data:
number_of_files = 2 #
filenames = [None]*number_of_files
df = [None]*number_of_files
indices = [None]*number_of_files

filenames[0] = "../output/ECstar_profile_EOS_DD2_rho0_4.0000000000_beta_0.0000000000_gamma_2.0000000000.txt"
filenames[1] = "../output/ECstar_profile_EOS_APR_rho0_4.0000000000_beta_0.0000000000_gamma_2.0000000000.txt"
#
#filenames[0] = "../output/ECstar_profile_EOS_APR_rho0_4.0000000000_beta_0.0000000000_gamma_55.0000000000.txt"
#filenames[1] = "../output/ECstar_profile_EOS_APR_rho0_4.0000000000_beta_0.0300000000_gamma_55.0000000000.txt"

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
kappa = 8*np.pi


nrows = 2
ncols = 1
fig = plt.figure(figsize=(6.5,5.5)) # or (6.5,5)
gs = fig.add_gridspec(nrows, ncols, hspace=0.1, wspace=0)
axs = gs.subplots(sharex='col')#, sharey='row')

c1 = "forestgreen"#"green"
c2 = "mediumblue"#"mediumslateblue"#"blue"

##########################################################################################################

#Compute Kelparian rotation rate for APR and DD2 EOS:
Radius_DD2 = df[0][:,indices_star['r']][-1]
Mass_DD2 = 0.5*Radius_DD2 * (1. - 1. / pow(df[0][:,indices_star['a']][-1], 2.) )   # M = R/2 * (1 - 1/ a(R)^2)
Omega_kep_DD2 = np.sqrt( Mass_DD2 / pow(Radius_DD2,3) )
print(Radius_DD2, Mass_DD2, Omega_kep_DD2)

Radius_APR = df[1][:,indices_star['r']][-1]
Mass_APR = 0.5*Radius_APR * (1. - 1. / pow(df[1][:,indices_star['a']][-1], 2.) )   # M = R/2 * (1 - 1/ a(R)^2)
Omega_kep_APR = np.sqrt( Mass_APR / pow(Radius_APR,3) )
print(Radius_APR, Mass_APR, Omega_kep_APR)

Omega_kep = Omega_kep_DD2 #np.sqrt(2.2964173616 / pow(8.6584375001,3)) # Keplerian rotation DD2 model
print(Omega_kep)

s2 = [0.0]*len(df[0][:,indices_star['r']])
kappas2 = [0.0]*len(df[0][:,indices_star['r']])
Gamma = [0.0]*len(df[0][:,indices_star['r']])
r = df[0][:,indices_star['r']]
e = df[0][:,indices_star['e']]
P = df[0][:,indices_star['P']]

Omega = Omega_kep
for i in range(len(df[0][:,indices_star['r']])):
	Gamma[i] = 1. / np.sqrt(1. + Omega*Omega*r[i]*r[i])
	s2[i] = pow(r[i], 4) * Omega*Omega* pow(Gamma[i], 4) * pow(e[i] + P[i], 2)
	kappas2[i] = kappa * s2[i]

axs[0].plot(df[0][:,indices_star['r']]*code_to_km, kappas2*df[0][:,indices_star['e']]/df[0][:,indices_star['e']]*1, lw=linew, ls = "-", color =c1, label="$\Omega = \Omega_{Kep}$")
axs[1].plot(df[0][:,indices_star['r']]*code_to_km, kappas2 / df[0][:,indices_star['e']] * 100, lw=linew, ls = "-", color =c1, label="")
max_relative_torsion_strength_Omegakep_DD2 = max(kappas2 / df[0][:,indices_star['e']])
print( max_relative_torsion_strength_Omegakep_DD2 )


Omega = 0.4*Omega_kep # Keplerian rotation
for i in range(len(df[0][:,indices_star['r']])):
	Gamma[i] = 1. / np.sqrt(1. + Omega*Omega*r[i]*r[i])
	s2[i] = pow(r[i], 4) * Omega*Omega* pow(Gamma[i], 4) * pow(e[i] + P[i], 2)
	kappas2[i] = kappa * s2[i]
axs[0].plot(df[0][:,indices_star['r']]*code_to_km, kappas2*df[0][:,indices_star['e']]/df[0][:,indices_star['e']]*1, lw=linew, ls = "--", color =c1, label="$\Omega = 0.4\Omega_{Kep}$")
axs[1].plot(df[0][:,indices_star['r']]*code_to_km, kappas2 / df[0][:,indices_star['e']] * 100, lw=linew, ls = "--", color =c1, label="")
print( max(kappas2 / df[0][:,indices_star['e']]) )


Omega = 0.1*Omega_kep # Keplerian rotation
for i in range(len(df[0][:,indices_star['r']])):
	Gamma[i] = 1. / np.sqrt(1. + Omega*Omega*r[i]*r[i])
	s2[i] = pow(r[i], 4) * Omega*Omega* pow(Gamma[i], 4) * pow(e[i] + P[i], 2)
	kappas2[i] = kappa * s2[i]
axs[0].plot(df[0][:,indices_star['r']]*code_to_km, kappas2*df[0][:,indices_star['e']]/df[0][:,indices_star['e']]*1, lw=linew, ls = ":", color =c1, label="$\Omega = 0.1\Omega_{Kep}$")
axs[1].plot(df[0][:,indices_star['r']]*code_to_km, kappas2 / df[0][:,indices_star['e']] * 100, lw=linew, ls = ":", color =c1, label="")
print( max(kappas2 / df[0][:,indices_star['e']]) )

###################################################################################################################################################
Omega_kep = Omega_kep_APR #np.sqrt(1.5983200402 / pow(7.7735937501,3)) # Keplerian rotation APR model
print(Omega_kep)


s2 = [0.0]*len(df[1][:,indices_star['r']])
kappas2 = [0.0]*len(df[1][:,indices_star['r']])
Gamma = [0.0]*len(df[1][:,indices_star['r']])
r = df[1][:,indices_star['r']]
e = df[1][:,indices_star['e']]
P = df[1][:,indices_star['P']]

Omega = Omega_kep
for i in range(len(df[1][:,indices_star['r']])):
	Gamma[i] = 1. / np.sqrt(1. + Omega*Omega*r[i]*r[i])
	s2[i] = pow(r[i], 4) * Omega*Omega* pow(Gamma[i], 4) * pow(e[i] + P[i], 2)
	kappas2[i] = kappa * s2[i]

axs[0].plot(df[1][:,indices_star['r']]*code_to_km, kappas2*df[1][:,indices_star['e']]/df[1][:,indices_star['e']]*1, lw=linew, ls = "-", color =c2, label="$\Omega = \Omega_{Kep}$")
axs[1].plot(df[1][:,indices_star['r']]*code_to_km, kappas2 / df[1][:,indices_star['e']] * 100, lw=linew, ls = "-", color =c2, label="")
max_relative_torsion_strength_Omegakep_APR = max(kappas2 / df[1][:,indices_star['e']])
print( max_relative_torsion_strength_Omegakep_APR )

Omega = 0.4*Omega_kep # Keplerian rotation
for i in range(len(df[1][:,indices_star['r']])):
	Gamma[i] = 1. / np.sqrt(1. + Omega*Omega*r[i]*r[i])
	s2[i] = pow(r[i], 4) * Omega*Omega* pow(Gamma[i], 4) * pow(e[i] + P[i], 2)
	kappas2[i] = kappa * s2[i]
axs[0].plot(df[1][:,indices_star['r']]*code_to_km, kappas2*df[1][:,indices_star['e']]/df[1][:,indices_star['e']]*1, lw=linew, ls = "--", color =c2, label="$\Omega = 0.4\Omega_{Kep}$")
axs[1].plot(df[1][:,indices_star['r']]*code_to_km, kappas2 / df[1][:,indices_star['e']] * 100, lw=linew, ls = "--", color =c2, label="")
print( max(kappas2 / df[1][:,indices_star['e']]) )


Omega = 0.1*Omega_kep # Keplerian rotation
for i in range(len(df[1][:,indices_star['r']])):
	Gamma[i] = 1. / np.sqrt(1. + Omega*Omega*r[i]*r[i])
	s2[i] = pow(r[i], 4) * Omega*Omega* pow(Gamma[i], 4) * pow(e[i] + P[i], 2)
	kappas2[i] = kappa * s2[i]
axs[0].plot(df[1][:,indices_star['r']]*code_to_km, kappas2*df[1][:,indices_star['e']]/df[1][:,indices_star['e']]*1, lw=linew, ls = ":", color =c2, label="$\Omega = 0.1\Omega_{Kep}$")
axs[1].plot(df[1][:,indices_star['r']]*code_to_km, kappas2 / df[1][:,indices_star['e']] * 100, lw=linew, ls = ":", color =c2, label="")
print( max(kappas2 / df[1][:,indices_star['e']]) )

###################################################################################################################################################


# set scale, plotlim, legend ...
axs[0].set_yscale("linear") # "log" or "linear"
axs[1].set_yscale("log") # "log" or "linear"

#axs[0].set_ylim((1e-7,0.002))
axs[0].set_ylim((0,1.75e-4))
axs[0].ticklabel_format(axis="y",scilimits=(0,0))
axs[1].set_ylim((1e-3,30))

axs[0].set_xlim((0.,14.))
axs[1].set_xlim((0.,14.))
#axs[1].set_yticks(np.geomspace(ylim2[0],ylim2[1], 51))

axs[0].grid(alpha=0.15, linestyle="--", c="black")
axs[1].grid(alpha=0.15, linestyle="--", c="black")

axs[0].legend(loc="upper left", fontsize = 10)
#axs[1].legend(loc="upper left", fontsize = 9)

# shared x-axis label of the plot
plt.xlabel("Radius $r$ [km]", fontsize = 18)
# shared y-axis label
axs[0].set_ylabel("$\kappa s^2$ [$M_\odot / R_{g\odot}^3$]", fontsize = 18)
axs[1].set_ylabel("$\kappa s^2/e$ [$\%$]", fontsize = 18)

# prevent labels on the "inside" between the subplots:
for ax in axs:
	ax.label_outer()
		
# plot labeling etc.:
# text in the plot:
axs[0].text(11.2, 0.000095, "DD2 EOS", fontsize = 12, rotation='horizontal', c=c1)
axs[0].text(6.2, 0.000059, "APR EOS", fontsize = 12, rotation='horizontal', c=c2)
#axs[0].text(4.8, 0.000125, "Preliminary", fontsize = 24, rotation='horizontal', c="gray", alpha=0.5)

# plot horizontal line to indicate maximal relative torsion contribution:
axs[1].hlines(y=max_relative_torsion_strength_Omegakep_DD2*100, xmin=0.5, xmax=11.8, linewidth=1, color=c1, alpha=0.5)
percent_value = "{:.2f}".format(max_relative_torsion_strength_Omegakep_DD2*100) + "$\,\%$"
axs[1].text(0.5, max_relative_torsion_strength_Omegakep_DD2*100*0.5, percent_value, fontsize = 10, rotation='horizontal', c=c1)

axs[1].hlines(y=max_relative_torsion_strength_Omegakep_APR*100, xmin=2.5, xmax=10.3, linewidth=1, color=c2, alpha=0.5)
percent_value = "{:.2f}".format(max_relative_torsion_strength_Omegakep_APR*100) + "$\,\%$"
axs[1].text(2.5, max_relative_torsion_strength_Omegakep_APR*100*0.5, percent_value, fontsize = 10, rotation='horizontal', c=c2)

#############################
figname = "Figure1.pdf"
plt.savefig(figname, dpi=400, bbox_inches='tight')
print("Saved figure as: " + figname)
#plt.show()
