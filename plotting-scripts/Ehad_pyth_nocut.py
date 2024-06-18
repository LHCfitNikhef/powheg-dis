########################################################
#
#
#                Outgoing Hadron Energy 
#                
#              PWG+PY8 LO, NLO - no cuts
#
#                    Pythia8 LO
#
#
#
########################################################

import matplotlib.pyplot as plt
import numpy as np


#                  NEUTRINO F I L E S

#           NO CUTS

file_names_lo = [
                '../DATA/Energy-had/LO/fixed_scales.dat',
                 '../DATA/Energy-had/LO/fixed_scales-1.dat',
                 ]

files_lo_min = [
                '../DATA/Energy-had/LO/fixed_scales_min.dat',
                '../DATA/Energy-had/LO/fixed_scales_min-1.dat',
               ]

files_lo_max = [
                '../DATA/Energy-had/LO/fixed_scales_max.dat',
                '../DATA/Energy-had/LO/fixed_scales_max-1.dat',
               ]


file_names_nlo = [
                '../DATA/Energy-had/NLO/fixed_scales.dat',
                '../DATA/Energy-had/NLO/fixed_scales-1.dat',
                  
                  ]
 
file_names_min = [
                '../DATA/Energy-had/NLO/fixed_scales_min.dat',
                '../DATA/Energy-had/NLO/fixed_scales_min-1.dat',
                  ]

file_names_max = [
                '../DATA/Energy-had/NLO/fixed_scales_max.dat',
                '../DATA/Energy-had/NLO/fixed_scales_min-1.dat',
                  ] 


pythia_file = '../DATA/Energy-had/Pythia/nocut.dat'





#               EXTRACT DATA

def extract_y_data(file_names):
    x_values_combined = []
    x_errors_combined= []
    y_values_combined = []
    

    for file_name in file_names:
        x_values = []
        x_errors = []
        y_values = []

        with open(file_name, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue

                columns = line.split()
                if len(columns) >= 4:
                    x_values.append((float(columns[0]) + float(columns[1])) / 2)
                    x_errors.append((float(columns[1]) - float(columns[0])) / 2)
                    y_values.append(float(columns[2]))

        x_values_combined.append(x_values)
        x_errors_combined.append(x_errors)                    
        y_values_combined.append(y_values)
        

    # Convert lists to NumPy arrays
    x_values_combined = np.array(x_values_combined)
    x_errors_combined = np.array(x_errors_combined)
    y_values_combined = np.array(y_values_combined)

    # Calculate the mean and standard error of the mean
    x_values_mean = np.mean(x_values_combined, axis=0)
    x_errors_mean = np.mean(x_errors_combined, axis=0)
    y_values_mean = np.mean(y_values_combined, axis=0)

    return x_values_mean, x_errors_mean, y_values_mean


nevents = 19999981
sigtot = 4.072
binwidth = 45


def pythia_data(file_name):
    Ehad_pyth, pyth_sigma = [], []

    with open(file_name, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue

            columns = line.split()
            if len(columns) == 2:
                Ehad_pyth.append(float(columns[0]))
                pyth_sigma.append((float(columns[1])) * sigtot / nevents / binwidth)
    return np.array(Ehad_pyth), np.array(pyth_sigma)
            
Ehad_pyth, pyth_sigma = pythia_data(pythia_file)



# Extract data
x_values_lo, x_errors_lo, y_values_lo = extract_y_data(file_names_lo)
x_values_nlo, x_errors_nlo, y_values_nlo = extract_y_data(file_names_nlo)

nothing, nothing, y_lo_min = extract_y_data(files_lo_min)
nothing, nothing, y_lo_max = extract_y_data(files_lo_max)

nothing, nothing, y_values_min = extract_y_data(file_names_min)
nothing, nothing, y_values_max = extract_y_data(file_names_max)


# Calculate ratios
ratio = y_values_nlo / y_values_lo
ratio_lo = y_values_lo / y_values_lo

min_lo = y_lo_min / y_values_lo
max_lo = y_lo_max / y_values_lo


min_ratios = y_values_min / y_values_lo
max_ratios = y_values_max / y_values_lo


# Create a figure and two subplots (2 rows, 1 column)
fig, axs =plt.subplots(2,1, figsize=(3.4,3.4), gridspec_kw={'height_ratios': [3, 1]}, dpi=300)

plt.rcParams['font.size'] = 9
plt.tight_layout()
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Computer Modern'],
    'mathtext.fontset': 'cm',
    'mathtext.rm': 'serif',
    'mathtext.it': 'serif:italic',
    'mathtext.bf': 'serif:bold',
    'axes.titlesize': 'x-large',      # Title font size
    'axes.labelsize': 'large',        # X and Y labels font size
    'xtick.labelsize': 'medium',      # X tick labels font size
    'ytick.labelsize': 'medium',      # Y tick labels font size
    'legend.fontsize': 'medium',      # Legend font size
    'figure.titlesize': 'x-large'     # Figure title font size
})



#       TOP PANEL

# LO

lines=[]
p1, = axs[0].plot(x_values_lo, y_values_lo, drawstyle="steps-post",
         label=r"$\rm POWHEG+PY8\; LO$",color="blue", linewidth=0.5)

axs[0].fill_between(x_values_lo, y_lo_min, y_lo_max,step="post",
          alpha=0.35,color="blue",linewidth=0.4,edgecolor="blue")

p2, = [plt.Rectangle((0, 0), 0, 0, facecolor="blue",linewidth=0.4, edgecolor="blue", alpha=0.35)]
lines+=((p1,p2),)

# NLO


p1, = axs[0].plot(x_values_nlo, y_values_nlo, drawstyle="steps-post",
         label=r"$\rm POWHEG+PY8\; NLO$",color="orange",linewidth=0.5)

axs[0].fill_between(x_values_nlo, y_values_min, y_values_max,step="post",
          alpha=0.35,color="orange", linewidth=0.4, edgecolor="orange")

p2, = [plt.Rectangle((0, 0), 0, 0, facecolor="orange",linewidth=0.4, edgecolor="orange", alpha=0.35)]
lines+=((p1,p2),)

#  PYTHIA

p1, = axs[0].plot(Ehad_pyth, pyth_sigma, drawstyle='steps-post',
            label=r'$\rm PY8\ STANDALONE$', color='red', linewidth=1)
lines+=(p1,)

ratio_pyth = pyth_sigma / y_values_lo
axs[1].plot(x_values_lo, ratio_pyth, drawstyle='steps-post',
            color='red', linewidth=1)

labels = [r"$\rm POWHEG+PY8\; LO$", r"$\rm POWHEG+PY8\; NLO$", r'$\rm PY8\ Stand\text{-}alone$']

axs[0].legend(lines, labels, loc='upper right', numpoints=1, fontsize='medium')



# neutrinos
axs[0].text(0.8, 0.1, r'${\nu}_e  p \rightarrow X  e^-$', 
            transform=axs[0].transAxes, horizontalalignment='center', 
            verticalalignment='bottom')




#   BOTTOM PANEL

axs[1].plot(x_values_nlo, ratio, drawstyle='steps-post',
            color='orange', linewidth=0.5)
axs[1].fill_between(x_values_nlo, min_ratios, max_ratios, step='post',
                    alpha=0.35, color='orange', linewidth=0.4, edgecolor='orange')



axs[1].plot(x_values_lo, ratio_lo, drawstyle='steps-post',
            color='blue', linewidth=0.5)
axs[1].fill_between(x_values_lo, min_lo, max_lo, step='post',
                    alpha=0.35, color='blue', linewidth=0.4, edgecolor='blue')

ratio_pyth = pyth_sigma / y_values_lo
axs[1].plot(x_values_lo, ratio_pyth, drawstyle='steps-post',
            color='red', linewidth=0.5)


#   PLOT SETTINGS 

axs[0].set_xlim(190, 910)
#axs[0].set_ylim(0.0036,0.0049)
axs[0].set_ylim(0.0005,0.0065)
axs[0].set_xticklabels([])
axs[0].set_ylabel(r'$d\sigma / dE_h \ \mathrm{[pb/GeV]}$')
axs[0].grid(linestyle='-', color='grey', linewidth=0.25)
axs[0].set_title(r'$E_\nu = \rm 1\ TeV$')



axs[1].set_xlim(190, 910)
axs[1].set_ylim(0.75, 1.25)
axs[1].set_xlabel(r'$E_h \ \mathrm{[GeV]}$')
axs[1].set_ylabel(r'$\rm Ratio \ to \ PW \ LO$', fontsize=8, labelpad=12)
axs[1].grid(linestyle='-', color='grey', linewidth=0.25)




for ax in axs:
    ax.tick_params(axis='both', which='both', width=1, labelsize='medium', direction='in')
    ax.tick_params(axis='x', which='both', top=True)  # Show ticks on the top
    ax.tick_params(axis='y', which='both', right=True)  # Show ticks on the right

for ax in axs:
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(9.8)


#adjust layout
fig.subplots_adjust(top=0.9, left=0.2, right=0.95, bottom=0.12, hspace=0.12)

# Show plot
#plt.show()

plt.savefig('Ehad_42_nocut.pdf', format='pdf')
