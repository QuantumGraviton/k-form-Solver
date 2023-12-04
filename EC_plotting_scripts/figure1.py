import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

# loads the custum FBS data values from a file written by the cpp code
def load_file(filename):
    f = open(filename)
    l0 = next(f)
    labels = l0.replace('#', '').strip().split('\t')
    indices = dict([(labels[i].strip(), i) for i in range(len(labels))])
    f.close()
    data = np.loadtxt(filename) #, delimeter = ' ')
    return data, indices


# load in data:
number_of_files = 2 #
filenames = [None]*number_of_files
df = [None]*number_of_files
indices = [None]*number_of_files

filenames[0] = "../output/ECstar_profile_EOS_DD2_rho0_4.0000000000_beta_0.0000000000_gamma_2.0000000000.txt"
filenames[1] = "../output/ECstar_profile_EOS_APR_rho0_4.0000000000_beta_0.0000000000_gamma_2.0000000000.txt"
#
filenames[0] = "../output/ECstar_profile_EOS_APR_rho0_4.0000000000_beta_0.0000000000_gamma_55.0000000000.txt"
filenames[1] = "../output/ECstar_profile_EOS_APR_rho0_4.0000000000_beta_0.0300000000_gamma_55.0000000000.txt"

for i in range(len(filenames)):
	df[i], indices[i] = load_file(filenames[i])
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
fig = plt.figure(figsize=(6.5,5))
gs = fig.add_gridspec(nrows, ncols, hspace=0.1, wspace=0)
axs = gs.subplots(sharex='col')#, sharey='row')



axs[0].plot(df[0][:,indices_star['r']]*code_to_km, df[0][:,indices_star['e']], lw=linew, ls = "-", color ="#DF5327", label="$e(r)$, DD2")

Omega_kep = np.sqrt(2.2964173616 / pow(8.6584375001,3)) # Keplerian rotation DD2 model
#Omega_kep = np.sqrt(1.5983200402 / pow(7.7735937501,3)) # Keplerian rotation APR model
print(Omega_kep)
#Omega = 713. / 71780.
#print(Omega)

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

axs[0].plot(df[0][:,indices_star['r']]*code_to_km, kappas2*df[0][:,indices_star['e']]/df[0][:,indices_star['e']]*1, lw=linew, ls = "-", color ="green", label="$\kappa s^2 (r)$, $\Omega_{Kep}$")
axs[1].plot(df[0][:,indices_star['r']]*code_to_km, kappas2 / df[0][:,indices_star['e']] * 100, lw=linew, ls = "-", color ="green", label="$\Omega = \Omega_{Kep}$")
print( max(kappas2 / df[0][:,indices_star['e']]) )

Omega = 0.4*Omega_kep # Keplerian rotation
for i in range(len(df[0][:,indices_star['r']])):
	Gamma[i] = 1. / np.sqrt(1. + Omega*Omega*r[i]*r[i])
	s2[i] = pow(r[i], 4) * Omega*Omega* pow(Gamma[i], 4) * pow(e[i] + P[i], 2)
	kappas2[i] = kappa * s2[i]
axs[0].plot(df[0][:,indices_star['r']]*code_to_km, kappas2*df[0][:,indices_star['e']]/df[0][:,indices_star['e']]*1, lw=linew, ls = "--", color ="green", label="$\kappa s^2 (r)$, $0.4\Omega_{Kep}$")
axs[1].plot(df[0][:,indices_star['r']]*code_to_km, kappas2 / df[0][:,indices_star['e']] * 100, lw=linew, ls = "--", color ="green", label="$\Omega = 0.4\Omega_{Kep}$")
print( max(kappas2 / df[0][:,indices_star['e']]) )


Omega = 0.1*Omega_kep # Keplerian rotation
for i in range(len(df[0][:,indices_star['r']])):
	Gamma[i] = 1. / np.sqrt(1. + Omega*Omega*r[i]*r[i])
	s2[i] = pow(r[i], 4) * Omega*Omega* pow(Gamma[i], 4) * pow(e[i] + P[i], 2)
	kappas2[i] = kappa * s2[i]
axs[0].plot(df[0][:,indices_star['r']]*code_to_km, kappas2*df[0][:,indices_star['e']]/df[0][:,indices_star['e']]*1, lw=linew, ls = ":", color ="green", label="$\kappa s^2 (r)$, $0.1\Omega_{Kep}$")
axs[1].plot(df[0][:,indices_star['r']]*code_to_km, kappas2 / df[0][:,indices_star['e']] * 100, lw=linew, ls = ":", color ="green", label="$\Omega = 0.1\Omega_{Kep}$")
print( max(kappas2 / df[0][:,indices_star['e']]) )

###################################################################################################################################################
#Omega_kep = np.sqrt(2.2964173616 / pow(8.6584375001,3)) # Keplerian rotation DD2 model
Omega_kep = np.sqrt(1.5983200402 / pow(7.7735937501,3)) # Keplerian rotation APR model
print(Omega_kep)
#Omega = 713. / 71780.
#print(Omega)

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

axs[0].plot(df[1][:,indices_star['r']]*code_to_km, kappas2*df[1][:,indices_star['e']]/df[1][:,indices_star['e']]*1, lw=linew, ls = "-", color ="blue", label="$\kappa s^2 (r)$, $\Omega_{Kep}$")
axs[1].plot(df[1][:,indices_star['r']]*code_to_km, kappas2 / df[1][:,indices_star['e']] * 100, lw=linew, ls = "-", color ="blue", label="$\Omega = \Omega_{Kep}$")
print( max(kappas2 / df[1][:,indices_star['e']]) )

Omega = 0.4*Omega_kep # Keplerian rotation
for i in range(len(df[1][:,indices_star['r']])):
	Gamma[i] = 1. / np.sqrt(1. + Omega*Omega*r[i]*r[i])
	s2[i] = pow(r[i], 4) * Omega*Omega* pow(Gamma[i], 4) * pow(e[i] + P[i], 2)
	kappas2[i] = kappa * s2[i]
axs[0].plot(df[1][:,indices_star['r']]*code_to_km, kappas2*df[1][:,indices_star['e']]/df[1][:,indices_star['e']]*1, lw=linew, ls = "--", color ="blue", label="$\kappa s^2 (r)$, $0.4\Omega_{Kep}$")
axs[1].plot(df[1][:,indices_star['r']]*code_to_km, kappas2 / df[1][:,indices_star['e']] * 100, lw=linew, ls = "--", color ="blue", label="$\Omega = 0.4\Omega_{Kep}$")
print( max(kappas2 / df[1][:,indices_star['e']]) )


Omega = 0.1*Omega_kep # Keplerian rotation
for i in range(len(df[1][:,indices_star['r']])):
	Gamma[i] = 1. / np.sqrt(1. + Omega*Omega*r[i]*r[i])
	s2[i] = pow(r[i], 4) * Omega*Omega* pow(Gamma[i], 4) * pow(e[i] + P[i], 2)
	kappas2[i] = kappa * s2[i]
axs[0].plot(df[1][:,indices_star['r']]*code_to_km, kappas2*df[1][:,indices_star['e']]/df[1][:,indices_star['e']]*1, lw=linew, ls = ":", color ="blue", label="$\kappa s^2 (r)$, $0.1\Omega_{Kep}$")
axs[1].plot(df[1][:,indices_star['r']]*code_to_km, kappas2 / df[1][:,indices_star['e']] * 100, lw=linew, ls = ":", color ="blue", label="$\Omega = 0.1\Omega_{Kep}$")
print( max(kappas2 / df[1][:,indices_star['e']]) )

###################################################################################################################################################


# set scale, plotlim, legend ...
axs[0].set_yscale("linear") # "log" or "linear"
axs[1].set_yscale("log") # "log" or "linear"

axs[0].set_ylim((1e-7,0.002))
axs[0].set_ylim((0,1.5e-4))
axs[1].set_ylim((1e-3,20))

axs[0].set_xlim((0.,14.))
axs[1].set_xlim((0.,14.))
#axs[1].set_yticks(np.geomspace(ylim2[0],ylim2[1], 51))

axs[0].grid(alpha=0.15, linestyle="--", c="black")
axs[1].grid(alpha=0.15, linestyle="--", c="black")

axs[0].legend(loc="upper left", fontsize = 9)
axs[1].legend(loc="upper left", fontsize = 9)

# shared x-axis label of the plot
plt.xlabel("Radius $r$ [km]", fontsize = 18)
axs[0].set_ylabel("$\kappa s^2$ [$M_\odot / R_{g\odot}^3$]", fontsize = 18)
axs[1].set_ylabel("$\kappa s^2/e$ [$\%$]", fontsize = 18)
# shared y-axis label
#fig.text(0.02, 0.5, ylabeltot, fontsize = 14, va='center', rotation='vertical')

# prevent labels on the "inside" between the subplots:
for ax in axs:
	ax.label_outer()
		
# plot labeling etc.:

#plt.title("relative error for $M_{tot}$ and $N_{B}$ inside stab region", fontsize = 18)
#plt.xlabel("Radius $r$ [km]", fontsize = 18)
#plt.ylabel("$e(r)$ [$M_\odot / R_g^3$]", fontsize = 18)
#plt.ylabel("Field $\phi(r)$ [$M_p$]", fontsize = 18)
#plt.grid(alpha=0.2, linestyle="--")#, c="black")

#plt.yscale("linear") # "log" or "linear"
#plt.ylim((1e-8,0.0025))
#plt.xlim((0.,15.))
#plt.legend(loc="upper right", fontsize = 11) #loc="lower left" or "upper right"

figname = "test4_rotationtest.pdf"
plt.savefig(figname, dpi=400, bbox_inches='tight')
print("Saved figure as: " + figname)
#plt.show()
