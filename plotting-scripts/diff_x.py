########################################################
#
#
#             x_B Differential Cross Section
#
#                LHC neutrino settings
#               
#             panel 1: differential distribution
#             panel 2: ratio to NLO
#
#
#
########################################################

import matplotlib.pyplot as plt
import numpy as np


#                   NEUTRINO F I L E S


file_names_lo = [
                '../DATA/x/LO/10cut.dat',
                  ]

file_names_nlo = [
                  '../DATA/x/NLO/10cut.dat',                  
                  ]
 
file_names_min = [
                  '../DATA/x/NLO/10cutmin.dat',
                  ]

file_names_max = [
                  '../DATA/x/NLO/10cutmax.dat',
                  ] 

pythia_file = '../DATA/x/Pythia/10cut.dat'


#               EXTRACT DATA

def extract_y_data(file_names):
    x_values_combined = []
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
                    y_errors.append(float(columns[3]) * (float(columns[2]) * pwg_events / pwg_sigtot)**0.5)

        x_values_combined.append(x_values)
        x_errors_combined.append(x_errors)                    
        y_values_combined.append(y_values)
        y_errors_combined.append(y_errors)
        

    # Convert lists to NumPy arrays
    x_values_combined = np.array(x_values_combined)
    x_errors_combined = np.array(x_errors_combined)
    y_values_combined = np.array(y_values_combined)
    y_erorrs_combined = np.array(y_errors_combined)

    # Calculate the mean and standard error of the mean
    x_values_mean = np.mean(x_values_combined, axis=0)
    x_errors_mean = np.mean(x_errors_combined, axis=0)
    y_values_mean = np.mean(y_values_combined, axis=0)
    y_errors_mean = np.mean(y_erorrs_combined, axis=0)

    return x_values_mean, x_errors_mean, y_values_mean, y_errors_mean


sigtot = 4.073     
nevents = 1000000      # high stats
binwidth = 0.05


def pythia_data(file_name):
    x_pyth, pyth_sigma = [], []

    with open(file_name, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue

            columns = line.split()
            if len(columns) == 2:
                x_pyth.append(float(columns[0]))
                pyth_sigma.append((float(columns[1])) * sigtot / nevents / binwidth)
    return np.array(x_pyth), np.array(pyth_sigma)
            
x_pyth, pyth_sigma = pythia_data(pythia_file)



# Extract LO and NLO data
x_values_lo, x_errors_lo, y_values_lo, y_errors_lo = extract_y_data(file_names_lo)
x_values_nlo, x_errors_nlo, y_values_nlo, y_errors_nlo = extract_y_data(file_names_nlo)


# Extract min and max data
nothing, nothing, y_values_min, nothing = extract_y_data(file_names_min)
nothing, nothing, y_values_max, nothing = extract_y_data(file_names_max)


# Calculate LO to NLO ratio
ratio = y_values_lo / y_values_nlo


# Calculate 7 point scale ratios
min_ratios = y_values_min / y_values_nlo
max_ratios = y_values_max / y_values_nlo

#line_widths = abs(y_values_max - y_values_min)


#   PYTHIA TO LO RATIO

pyth_sigma_short = pyth_sigma[:len(y_values_lo)] 
x_pyth_short = x_pyth[:len(pyth_sigma_short)] 
ratio_pyth = pyth_sigma_short / y_values_lo  
#print("Calculated ratio:", ratio_pyth)
#print("Length of the calculated ratio:", len(ratio_pyth))



#                           CREATE FIGURE            


# Create a figure and two subplots (2 rows, 1 column)
# fig, axs = plt.subplots(3, 1, figsize=(7,10))

            #   two panel settings
fig, axs = plt.subplots(2, 1, figsize=(7.5, 11))



# Get the default color cycle
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
default_blue = default_colors[0] 
            

# LO

axs[0].errorbar(x_values_lo, y_values_lo, xerr=x_errors_lo, yerr=None, fmt='o', markersize=5, 
   color='sandybrown', ecolor='sandybrown', capsize=4, label='LO pwg+pythia8 default')

# NLO

axs[0].errorbar(x_values_nlo, y_values_nlo, xerr=x_errors_nlo, yerr=None, fmt='o', markersize=5, 
    color=default_blue, ecolor=default_blue,
    capsize=4, label='NLO pwg+pythia8 default')

# PYTHIA

axs[0].errorbar(x_pyth, pyth_sigma, xerr=(binwidth/2), fmt='x', capsize=5, 
             markersize=5, color='red', ecolor='red', label='pythia8')


axs[0].set_yscale('log')
axs[0].set_ylim(10**-2.5,10**2)
axs[0].set_xlim(0,0.8)
axs[0].set_xticklabels([]) 
axs[0].set_ylabel(r'$\mathrm{d\sigma / dx_B \ [pb]}$')
axs[0].grid(linestyle=':', color='black', linewidth=1, alpha=0.8)

#       PRINT NEUTRINO PROCESS
axs[0].text(0.8, 0.67, r'${\nu}_e  p \rightarrow X  e^-$', 
            transform=axs[0].transAxes, fontsize=12, horizontalalignment='center', 
            verticalalignment='bottom')
#axs[0].set_title(r'Differential cross section w.r.t  $\mathrm{x_B}$')


#           RATIO TO NLO PANEL (MIDDLE PANEL)


axs[1].errorbar(x_values_lo, ratio, xerr=x_errors_lo, fmt='o', capsize=4,
                color='sandybrown')

#  7 point scale variation

width = binwidth
for x, min_ratio, max_ratio in zip(x_values_nlo, min_ratios, max_ratios):
    axs[1].bar(x, max_ratio - min_ratio, bottom=min_ratio, width=width, 
                color='gray', edgecolor='black', linewidth=1, align='center', 
                label='7-point scale variation')
 
axs[1].set_xlim(0,0.8)
axs[1].set_ylim(0.7,1.2)
#axs[1].set_xticklabels([])   # uncomment for third panel
#axs[1].set_yticklabels([])    # uncomment for third panel
axs[1].axhline(y=1, color='black', linestyle='-')
axs[1].set_xlabel(r'$\mathrm{x_B}$')
axs[1].set_ylabel('Ratio to NLO')
axs[1].grid(linestyle=':', color='black', linewidth=1, alpha=0.7)




#       EXTRA PANEL - RATIO PYTHIA / LO

# axs[2].errorbar(x_pyth_short, ratio_pyth, xerr=(binwidth/2), fmt='x', capsize=4,
#                color='red')

# #axs[2].set_xlim(0, 0.7)
# #axs[2].set_xlim(0.1,0.5)
# #axs[2].set_ylim(0.8,1.1)
# #axs[2].set_ylim(0.88,0.96)
# #axs[2].set_yticks(np.arange(0.8, 1.1, 0.1))

# axs[2].axhline(y=1, color='black', linestyle='-')
# axs[2].grid(linestyle=':', color='black', linewidth=1, alpha=0.7)
# axs[2].set_ylabel('Ratio to LO')
# axs[2].set_xlabel(r'$\mathrm{x_B}$')


# legend
handles, labels = [], []
for ax in axs:
    for h, l in zip(*ax.get_legend_handles_labels()):
        if l not in labels:  
            handles.append(h)
            labels.append(l)

fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(1, 1), 
           bbox_transform=axs[0].transAxes, fontsize='large')


for ax in axs:
    ax.tick_params(axis='both', which='both', length=6, width=1, labelsize='medium', direction='in')
    ax.tick_params(axis='x', which='both', top=True)  
    ax.tick_params(axis='y', which='both', right=True)  


#adjust layout
fig.subplots_adjust(top=0.95, bottom=0.05, left=0.125, right=0.9, hspace=0.0)

# Show plot
plt.show()
