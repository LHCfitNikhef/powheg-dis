# Script to calculate the LO structure fct F2 as a fct of Q2 for different x's from neutrino PDFs--> fig. 32.2 in QFT & the SM
# Second part makes double differential cross section for neutrino DIS  as a fct of Q2 for different x's --> benchmark for powheg

import numpy as np
import matplotlib.pyplot as plt
from numpy import *
from mpl_toolkits import mplot3d
import lhapdf

# get the LHAPDF set
pset = lhapdf.getPDFSet('CT10nlo')
pdfs1 = pset.mkPDFs()

flavours = [1,2,3,4,5]

# Q2 Q2 Q2 Q2 Q2 --> pdf's as a fct of Q2 for a few values of x

Q2set = np.linspace(10,1000,53)
#xs = [0.0013, 0.005, 0.013, 0.05, 0.13, 0.5, 0.85]
xs = [0.0001, 0.0005, 0.001, 0.01, 0.1, 0.3, 0.4]

# import other parameters
alfa = 1/137							# electron coupling
M = 1 								# proton mass
s = 1000							# Mandelstam s
y = 0.5								


# Tables for dd cross section and structure fct, columns are x-values, rows are Q2-values
dsigmadxdQ2 = [[],[],[],[],[],[],[]]
F2 = [[],[],[],[],[],[],[]]

for j in range(len(np.array(xs))):
	# make an empty table for the pdf's per x-value, columns are x-values, rows are Q2-values
	pdfdQ = [[],[],[],[],[],[],[]]
	pdfuQ = [[],[],[],[],[],[],[]]
	pdfsQ = [[],[],[],[],[],[],[]]
	pdfcQ = [[],[],[],[],[],[],[]]
	pdfbQ = [[],[],[],[],[],[],[]]
	
	# we will append cross sections and structure fct-s per value of x to this
	dsigmadxdQ2[j] = []
	F2[j] = []

	# import the pdf-s from lhapdf and calculate the f_q^+'s (sum of quark and antiquark contributions)
	for i in range(len(Q2set)):
		p_d = pdfs1[i].xfxQ2(1,xs[j],Q2set[i])
		p_antid = pdfs1[i].xfxQ2(-1,xs[j],Q2set[i])
		p_dplus = p_d + p_antid
		pdfdQ[j].append(p_dplus)

		p_u = pdfs1[i].xfxQ2(2,xs[j],Q2set[i])
		p_antiu = pdfs1[i].xfxQ2(-2,xs[j],Q2set[i])
		p_uplus = p_u + p_antiu
		pdfuQ[j].append(p_uplus)

		p_s = pdfs1[i].xfxQ2(3,xs[j],Q2set[i])
		p_antis = pdfs1[i].xfxQ2(-3,xs[j],Q2set[i])
		p_splus = p_s + p_antis
		pdfsQ[j].append(p_splus)

		p_c = pdfs1[i].xfxQ2(4,xs[j],Q2set[i])
		p_antic = pdfs1[i].xfxQ2(-4,xs[j],Q2set[i])
		p_cplus = p_c + p_antic
		pdfcQ[j].append(p_cplus)

		p_b = pdfs1[i].xfxQ2(5,xs[j],Q2set[i])
		p_antib = pdfs1[i].xfxQ2(-5,xs[j],Q2set[i])
		p_bplus = p_b + p_antib
		pdfbQ[j].append(p_bplus)
	
	# calculate the structure fct F_2 (neglect others in large Q2 regime) and double diff cross section and append these to the lists per x-value
	for k in range(len(Q2set)):
		F2k = xs[j]*(4/9*(pdfuQ[j][k] + pdfcQ[j][k]) + 1/9*(pdfdQ[j][k] + pdfsQ[j][k]))		# F_2 from the LHC as a Neutrino-Ion colllider paper for leptop-proton DIS
		F2[j].append(F2k)	# for filling the structure fct array --> e.g. for plotting a fig like 32.2 in QFT & the SM
		dsigmadxdy_k = 2 * np.pi * alfa**2 * s * (1 + (1 - y)**2) *  Q2set[k] * F2k	# double differential cross section from p.676 of QFT and the SM
		dsigmadxdQ2_k = 1/(xs[j] * (s - M**2)) * dsigmadxdy_k	# relation taken from Particle Data Book review ch. 18
		dsigmadxdQ2_R = dsigmadxdQ2_k * xs[j]*Q2set[k]**2 / (2*np.pi*alfa**2*(1+(1-y)**2))	# reduced dd cross section, defined in 3.3 of POWHEG paper	
		dsigmadxdQ2[j].append(dsigmadxdQ2_k)	# fill up cross section list per x-value


# Make the structure fct figure
fig = plt.figure(figsize =(14, 9))
for i in range(len(xs)):
	plt.plot(Q2set, 2**(i+1)*np.array(F2[i]), '.', label='x = %s'%xs[i])
	print(i)
	print('F2 for x=%s'%xs[i], F2[i])
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.title('Structure function F2(x,Q2) for different values of x')
plt.xlabel('Q2 (GeV^2)')
plt.ylabel('F2(x,Q2) * 2^i')
fig.savefig('20112023fig32_2_lowx.png')


# Make the cross section figure
fig = plt.figure(figsize =(14, 9))
for i in range(len(xs)):
        plt.plot(np.log(Q2set), 2**(25-i)* np.array(dsigmadxdQ2[i]), '.', label='x = %s'%xs[i])
plt.legend()
#plt.xscale('log')
#plt.yscale('log')
plt.title('Reduced double differential cross section for different values of x')
plt.xlabel('log Q2 (GeV^2)')
plt.ylabel('2^(25-i_x)dsigma/dxdQ2')
fig.savefig('20112023ddcrosssection_lowx.png')

