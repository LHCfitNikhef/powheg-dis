import numpy as np
import matplotlib.pyplot as plt
import lhapdf
from numpy import *

G_F = 1.1663787 * 10**(-5) 
M_W = 80.398
alfa = 0.007818E0
x = 0.056
xbin = 0.002

E1 = 0.938
E2 = 1000


def y(x,q2):
        E = 2*E1*E2
        ys = q2 / (x*E)
        return ys

# OUR POWHEG RESULTS

powhegpath = './pwgLHEF_analysis-11-Q2056.dat'
powhegpathmin = './pwgLHEF_analysis-min-Q2056.dat'
powhegpathmax = './pwgLHEF_analysis-max-Q2056.dat'

Qj1 = np.loadtxt(powhegpath, usecols=[0])
Qj2 = np.loadtxt(powhegpath, usecols=[1])
Q2 = 0.5*(Qj1 + Qj2)
qbin = Qj2[5] - Qj1[5]

sigpwg = np.loadtxt(powhegpath, usecols=[2])
sigmin = 1/ xbin * 2.56819e-09 * 4*np.pi*x / G_F**2 * (np.array(Q2) + M_W**2)**2 / M_W**4 * np.loadtxt(powhegpathmin, usecols=[2])
sigmax = 1 / xbin * 2.56819e-09 * 4*np.pi*x / G_F**2 * (np.array(Q2) + M_W**2)**2 / M_W**4 * np.loadtxt(powhegpathmax, usecols=[2])

sigpwgnorm = np.array(sigpwg) / xbin * 2.56819e-09
sigr_pwg = 4*np.pi*x / G_F**2 * (np.array(Q2) + M_W**2)**2 / M_W**4 * sigpwgnorm



# PYTHIA

plot = open('./Q2-bin-tests-3004-0.dat')
plot = [line.split() for line in plot]
valx = [float(x[0]) for x in plot]
valy = [float(x[1]) for x in plot]

qbinpyt = np.array(valx)[3]-np.array(valx)[2]

Q2pythia = []
for i in range(len(valx) - 1):
        Q2pythia.append(0.5*(valx[i]+valx[i+1]))
Q2pythia.append(valx[19]+qbinpyt/2)

sigtot = 4.072*2.56819e-09
Ntot = 5e7

# Get the total pythia (differential) cross section in pb / GeV
xsec = np.array(valy) / Ntot * sigtot / qbinpyt / xbin
sigr_pythia = 4*np.pi*x / G_F**2 * (np.array(Q2pythia) + M_W**2)**2 / M_W**4 * xsec 



# Analytical calculation

# PDF set
p = lhapdf.mkPDF("NNPDF40_nnlo_as_01180", 0)

# CKM matrix
Vud = 0.97383
Vus = 0.2272
Vub = 0.00396
Vdc = 0.2271
Vsc = 0.97296
Vcb = 0.04221
Vtd = 0.00814
Vts = 0.04161
Vtb = 0.9991

Q2s = []
F2 = []
xF3 = []
dsigmadxdQ2 = []
F2nu = []
xF3nu = []
sigr_nu = []
for k in range(len(np.array(Q2))):
        pdfd = p.xfxQ(1, x, np.sqrt(Q2[k]))
        pdfdbar = p.xfxQ(-1,x,np.sqrt(Q2[k]))
        pdfu = p.xfxQ(2, x, np.sqrt(Q2[k]))
        pdfubar = p.xfxQ(-2,x,np.sqrt(Q2[k]))
        pdfs = p.xfxQ(3, x, np.sqrt(Q2[k]))
        pdfsbar = p.xfxQ(-3,x,np.sqrt(Q2[k]))
        pdfc = p.xfxQ(4, x, np.sqrt(Q2[k]))
        pdfcbar = p.xfxQ(-4,x,np.sqrt(Q2[k]))
        pdfb = p.xfxQ(5, x, np.sqrt(Q2[k]))
        pdfbbar = p.xfxQ(-5,x,np.sqrt(Q2[k]))
        pdfgluon = p.xfxQ(21,x,np.sqrt(Q2[k]))	

        pdfdplus = pdfd + pdfdbar
        pdfuplus = pdfu + pdfubar
        pdfsplus = pdfs + pdfsbar
        pdfcplus = pdfc + pdfcbar
        pdfbplus = pdfb + pdfbbar
        #print('pdfdplus:',pdfdplus)
        #print('pdfuplus:',pdfuplus)
        #print('pdfsplus:',pdfsplus)
        #print('pdfcplus:',pdfcplus)
        #print('pdfbplus:',pdfbplus)
        Q2s.append(Q2[k])
        ys = y(x, Q2[k])
        
        # NEUTRINOS CC
        F2nuj = 2*x*( (Vud**2 + Vus**2 + Vub**2)*pdfubar + (Vud**2 + Vdc**2)*pdfd + (Vsc**2 + Vus**2)*pdfs + (Vsc**2 + Vdc**2 + Vcb**2)*pdfcbar + (Vub**2 + Vcb**2)*pdfb)
        xF3nuj = 2*x*( -(Vud**2 + Vus**2 + Vub**2)*pdfubar + (Vud**2 + Vdc**2)*pdfd + (Vsc**2 + Vus**2)*pdfs - (Vsc**2 + Vdc**2 + Vcb**2)*pdfcbar + (Vub**2 + Vcb**2)*pdfb)
        sigr_nuj = (1 + (1-ys)**2)*F2nuj + (1 - (1-ys)**2)*xF3nuj        
        sigr_nu.append(sigr_nuj)

        # ELECTRONS NC
        #F2j = (4/9*(pdfuplus + pdfcplus) + 1/9*(pdfdplus + pdfsplus + pdfbplus))
        #xF3j = 2*x*(-pdfubar + pdfd + pdfs - pdfcbar)                                  # electron structure fct F2
        #F2.append(F2j)
        #dsigmadxdQ2j = 4*np.pi*alfa**2 * ((0.5*(1 + (1 - ys1)**2)*F2j)) / (x * Q2[k]**2)        # LO cross section
        #dsigmadxdQ2_R = x*Q2[k]**2/(2*np.pi*alfa**2*(1+(1-ys1)**2))*dsigmadxdQ2j        # LO reduced cross section
        #dsigmadxdQ2.append(F2j)




# Calculate ratio
ratio = np.array(sigr_pwg) / np.array(sigr_pythia)
ratio2 = np.array(sigr_nu) / np.array(sigr_pythia)
print('pwg / pythia',ratio)
print('pythia / LO calc', 1/ratio2)


# Plotting
bin_edges = np.append(Qj1,Qj1[-1])
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8,10), gridspec_kw={'height_ratios': [3, 1]})

#fig.set_size_inches(12.00,12.00)
ax1.errorbar(Q2pythia, np.array(sigr_pythia), fmt='.', yerr=None, xerr=qbinpyt/2, color='orange', solid_capstyle='projecting', capsize=5, label='Pythia8')
ax1.errorbar(Q2, sigr_pwg, fmt='.', yerr=None, xerr=qbin/2, color='steelblue', solid_capstyle='projecting', capsize=5,label='POWHEG LO')
ax1.errorbar(Q2, sigr_nu, fmt='.', yerr=None, xerr=qbin/2, color='red', solid_capstyle='projecting', capsize=5,label='LO calc')
ax1.step(bin_edges, np.append(sigmin, sigmin[-1]), 'lightsteelblue', where='post')
ax1.step(bin_edges, np.append(sigmax, sigmax[-1]),'lightsteelblue', where='post')
#ax1.set_ylim(1.1, 1.6)
ax1.fill_between(bin_edges, np.append(sigmin, sigmin[-1]), np.append(sigmax, sigmax[-1]), step='post', color='lightsteelblue', alpha=0.5)

fig.text(0.2, 0.5, r'$\nu_e + p \rightarrow e_{+} + X$', fontsize='large')
ax1.grid('whitesmoke')
ax1.legend(frameon=False, loc='best', fontsize='large')
ax1.set_title(r'Reduced cross section for $E_{\nu} = 1$ TeV at $x = 0.056$')
ax1.set_ylabel(r'$\sigma_R^{CC}$', fontsize='large')

# Plot the ratio in the bottom panel
ax2.errorbar(np.array(Q2), ratio, fmt='.', yerr=None, xerr=qbin/2, color='steelblue', solid_capstyle='projecting', capsize=5)
ax2.errorbar(np.array(Q2pythia), ratio2, fmt='.', yerr=None, xerr=qbin/2, color='red', solid_capstyle='projecting', capsize=5)
ax2.axhline(y=1, color='black', linestyle='--')  # Add a horizontal line at y=1 for reference
ax2.set_xlabel(r'$Q^2 (GeV^2)$', fontsize='large')
ax2.set_ylabel('Ratio to pythia')
#ax2.set_ylim(0.9,1.2)
ax2.set_xlim(4,100)
ax2.grid()

#fig.tight_layout()
fig.subplots_adjust(hspace=0.05)
fig.savefig("sigr-comp-0305.png")



