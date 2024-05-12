########################################################
#
#
#             Q2 Differential Cross Section Plot
#
#                LHC neutrino settings
#               
#            panel 1: cross section distribution
#            panel 2: ratio to NLO
#
#
#
########################################################

import matplotlib.pyplot as plt
import numpy as np


#                   NEUTRINO F I L E S

file_names_lo = [
                '../data/q2/pwg-lo.dat'
                  ]

file_names_nlo = [
                  '../data/q2/pwg-nlo.dat',
                  ]
 
file_names_min = [
                  '../data/q2/pwg-nlo-min.dat',
                  ]

file_names_max = [
                  '../data/q2/pwg-nlo-max.dat',
                  ] 

pythia_file = '../data/q2/pythia.dat'


#               EXTRACT DATA


def pwg_data(file_names):
    x_values_combined = []    # x and y values correspond to x and y axis, NOT dis kinematic variables
    x_errors_combined= []
    y_values_combined = []
    y_errors_combined = []
    
    pwg_events = 5000000 * 64
    pwg_sigtot = 4.29

    for file_name in file_names:
        x_values = []
        x_errors = []
        y_values = []
        y_errors = []

        with open(file_name, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue

                columns = line.split()
                if len(columns) >= 4:
                    x_values.append((float(columns[0]) + float(columns[1])) / 2)
                    x_errors.append((float(columns[1]) - float(columns[0])) / 2)
                    y_values.append(float(columns[2]))
                    y_errors.append((float(columns[3])**2) * (float(columns[2]) * pwg_events / pwg_sigtot))
                    
        x_values_combined.append(x_values)
        x_errors_combined.append(x_errors)
        y_values_combined.append(y_values)
        y_errors_combined.append(y_errors)
        

    # Convert lists to NumPy arrays
    x_values_combined = np.array(x_values_combined)
    x_errors_combined = np.array(x_errors_combined)
    y_values_combined = np.array(y_values_combined)
    y_errors_combined = np.array(y_errors_combined)

    # Calculate the mean and standard error of the mean
    x_values_mean = np.mean(x_values_combined, axis=0)
    x_errors_mean = np.mean(x_errors_combined, axis=0)
    y_values_mean = np.mean(y_values_combined, axis=0)
    y_errors_mean = np.mean(y_errors_combined, axis=0)

    return x_values_mean, x_errors_mean, y_values_mean, y_errors_mean


sigtot = 4.071      # pythia's total cross section, for conversion from number of events to pb
nevents = 1000000      # pythia's total number of events, also for conversion
binwidth = 49.8        # constant binwidth in GeV^2, also for conversion


def pythia_data(file_name):
    q2_pyth, pyth_sigma = [], []

    with open(file_name, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue

            columns = line.split()
            if len(columns) == 2:
                q2_pyth.append(float(columns[0]))
                pyth_sigma.append((float(columns[1])) * sigtot / nevents / binwidth)
    return np.array(q2_pyth), np.array(pyth_sigma)
            
q2_pyth, pyth_sigma = pythia_data(pythia_file)


# Extract LO and NLO data
x_values_lo, x_errors_lo, y_values_lo, y_errors_lo = pwg_data(file_names_lo)
x_values_nlo, x_errors_nlo, y_values_nlo, y_errors_nlo = pwg_data(file_names_nlo)

# Extract min and max data
nothing, nothing, y_values_min, nothing = pwg_data(file_names_min)
nothing, nothing, y_values_max, nothing = pwg_data(file_names_max)


# Calculate LO to NLO ratio
ratio = y_values_lo / y_values_nlo


# Calculate 7 point scale ratios
min_ratios = y_values_min / y_values_nlo
max_ratios = y_values_max / y_values_nlo

ratio_pyth_nlo = pyth_sigma / y_values_nlo

# Create a figure and two subplots (2 rows, 1 column)
fig, axs = plt.subplots(2, 1, figsize=(7.5, 11)) 

# fig, axs = plt.subplots(3, 1, figsize=(7,10))

# Get the default color cycle
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
default_blue = default_colors[0] 
            
# Plotting commands

#       TOP PANEL

# LO

axs[0].errorbar(x_values_lo, y_values_lo, xerr=x_errors_lo, yerr=None, fmt='o', markersize=5, 
   color='sandybrown', ecolor='sandybrown', capsize=4, label='LO pwg+pythia8')

# NLO

axs[0].errorbar(x_values_nlo, y_values_nlo, xerr=x_errors_nlo, yerr=None, fmt='o', markersize=5, 
    color=default_blue, ecolor=default_blue,
    capsize=4, label='NLO pwg+pythia8')

# PYTHIA

axs[0].errorbar(q2_pyth, pyth_sigma, xerr=(binwidth/2), fmt='x', capsize=5, 
             markersize=5, color='red', ecolor='red', label='pythia8')


axs[0].set_yscale('log')
axs[0].set_xlim(53.8, 1000)

axs[0].set_xticklabels([])
axs[0].set_ylabel(r'$\mathrm{d\sigma / dQ^2 \ [pb/GeV^2]}$')
axs[0].grid(linestyle=':', color='black', linewidth=1, alpha=0.8)

#   PRINT NEUTRINO PROCESS

axs[0].text(0.8, 0.67, r'${\nu}_e  p \rightarrow X  e^-$',        
            transform=axs[0].transAxes, fontsize=12, 
            horizontalalignment='center', verticalalignment='bottom')


#axs[0].set_title(r'Differential cross section w.r.t  $\mathrm{Q^2 \ [GeV^2]}$')


#       BOTTOM PANEL

axs[1].errorbar(x_values_lo, ratio, xerr=x_errors_lo, fmt='o', capsize=4,
                color='sandybrown')

# axs[1].errorbar(x_values, ratio_pyth_nlo, xerr=x_errors, fmt='o', capsize=4,
#                 color='red')

# # Plot the 7 point scale bars on the bottom panel

width = x_errors_nlo[0]*2
for x, min_ratio, max_ratio in zip(x_values_nlo, min_ratios, max_ratios):
    axs[1].bar(x, max_ratio - min_ratio, bottom=min_ratio, width=width, 
                color='gray', edgecolor='black', linewidth=1, align='center', 
                label='7-point scale variation')
    
axs[1].set_xlim(53.8, 1000)          
axs[1].set_ylim(0.75, 1.1)     
axs[1].axhline(y=1, color='black', linestyle='-')
axs[1].set_xlabel(r'$\mathrm{Q^2 \ [GeV^2]}$')
axs[1].set_ylabel('Ratio to NLO')
axs[1].grid(linestyle=':', color='black', linewidth=1, alpha=0.7)


#       EXTRA PANEL - RATIO PYTHIA / LO

#pyth_sigma_short = pyth_sigma[:len(y_values_lo)]    # to equalise the number of data points in pwg vs pythia
#q2_pyth_short = q2_pyth[:len(pyth_sigma_short)] 
#ratio_pyth = pyth_sigma_short / y_values_lo 

# axs[2].errorbar(q2_pyth_short, ratio_pyth, xerr=(binwidth/2), fmt='x', capsize=4,
#                color='red')

# axs[2].set_xlim(53.8, 1000)
# axs[2].set_yticks(np.arange(0.9, 1.1, 0.1))

# axs[2].axhline(y=1, color='black', linestyle='-')
# axs[2].grid(linestyle=':', color='black', linewidth=1, alpha=0.7)
# axs[2].set_ylabel('Ratio to LO')
# axs[2].set_xlabel(r'$\mathrm{Q^2} \ (GeV^2)$')


#  legend
handles, labels = [], []
for ax in axs:
    for h, l in zip(*ax.get_legend_handles_labels()):
        if l not in labels:  
            handles.append(h)
            labels.append(l)


fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1, 1), 
           bbox_transform=axs[0].transAxes, fontsize='large')

# Adjust tick parameters 
for ax in axs:
    ax.tick_params(axis='both', which='both', length=6, width=1, labelsize='medium', direction='in')
    ax.tick_params(axis='x', which='both', top=True)  
    ax.tick_params(axis='y', which='both', right=True)  


#adjust layout
fig.subplots_adjust(top=0.95, bottom=0.05, left=0.125, right=0.9, hspace=0.0)

# Show plot
plt.show()

