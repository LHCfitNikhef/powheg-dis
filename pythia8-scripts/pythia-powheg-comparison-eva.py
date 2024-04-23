# Python script to plot pythia and powheg in the same figure

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from numpy import *
import lhapdf

# Parameters

alfa = 1/137.035999084
xs = 0.003498
mb2gev = 2.56819
E1 = 1000
E2 = 0.938
G_F = 1.166 * 10**(-5)
M_W = 80.3

def y(x,q2):
        E = 2*E1*E2
        ys = q2 / (x*E)
        return ys	

# Pythia output

pp   = PdfPages('El-pythia-ycuts.pdf')
tmp1 = plt.figure(1)
tmp1.set_size_inches(8.00,6.00)
plot = open('pythia-El-data.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]

sigtot = 4.223 #pb
qbin = np.array(valx)[3]-np.array(valx)[2]
Ntot = 999905

# Get the total pythia (differential) cross section in pb / GeV
xsec = np.array(valy) / Ntot * sigtot / qbin

# POWHEG+PYTHIA output

Q2pwg1 = np.loadtxt('./powheg-LHEF-Elout-LO.dat', usecols=[0])
Q2pwg2 = np.loadtxt('./powheg-LHEF-Elout-LO.dat', usecols=[1])
outputpwg = np.loadtxt('./powheg-LHEF-Elout-LO.dat', usecols=[2])

Q2pwg12 = np.loadtxt('./powheg+pythia-Elout-NLO.dat', usecols=[0])
Q2pwg22 = np.loadtxt('./powheg+pythia-Elout-NLO.dat', usecols=[1])
outputpwg2 = np.loadtxt('./powheg+pythia-Elout-NLO.dat', usecols=[2])

Q2pwg = 0.5*(Q2pwg1 + Q2pwg2)
qbinpwg = Q2pwg2[2] - Q2pwg1[2]

Q2pwg2 = 0.5*(Q2pwg12 + Q2pwg22)
qbinpwg2 = Q2pwg22[2] - Q2pwg12[2]

xsecpwg = outputpwg			# pb / GeV
xsecpwg2 = outputpwg2			# pb/ GeV

# PLOT THE FIGURE

plt.errorbar(np.array(valx), np.array(xsec), fmt='.',yerr=None, xerr=qbin/2, label='Pythia')
plt.errorbar(Q2pwg, xsecpwg, fmt='.', yerr=None, xerr=qbinpwg/2, label='POWHEG LO')
plt.errorbar(Q2pwg2, xsecpwg2, fmt='.', yerr=None, xerr=qbinpwg/2, label='POWHEG NLO')
plt.ylim( 0.0030, 0.0050)
plt.grid()
plt.legend(frameon=False,loc='best')
plt.title('CC DIS with a 1 TeV neutrino')
plt.xlabel(r'$E_l$ (GeV)')
plt.ylabel(r'$\frac{d\sigma}{dE_l}$ (pb / GeV)')
pp.savefig(tmp1,bbox_inches='tight')
plt.clf()
pp.close()


