#!/usr/bin/env python
# coding: utf-8

import os,sys
import lhapdf
import numpy as np
import matplotlib.pyplot as py
from matplotlib import gridspec
from  matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)
from pylab import *
import scipy
from scipy.integrate import dblquad

print("\n *********************************************************")
print("      Plotting Neutrino PDFs                               ")
print(" ***********************************************************\n")


nfl=3
# Set x grid
xmin = 3e-4
xmax= 0.99
nx=1000
X = np.logspace(log(xmin),log(xmax),nx)
lhapdf.setVerbosity(0)
# number of experiments
nexp=3
# number of neutrino flavours
nfl=6
fl=[12,14,16,-12,-14,-16]

pdfset=["faserv",\
        "snd",\
        "faserv2"]
pdfsetlab=[r"${\rm FASER}\nu$",\
           r"${\rm SND@LHC}$",\
           r"${\rm FASER}\nu 2$"]

for iset in range(nexp):

    if(iset==0):
        fit1 = np.zeros((nfl,nx))
    if(iset==1):
        fit2 = np.zeros((nfl,nx))
    if(iset==2):
        fit3 = np.zeros((nfl,nx))
    if(iset==3):
        fit4 = np.zeros((nfl,nx))
 
    # Run over replicas
    p=lhapdf.mkPDF(pdfset[iset],0)
    print(pdfset[iset])
            
    # Run over x arrat
    for k in range(nx):
            
        x = X[k]
        q = 2 # Random number

        # run over flavours
        for ifl in range(nfl):

            if(iset==0):
                fit1[ifl][k] = p.xfxQ(fl[ifl],x,q)
            if(iset==1):
                fit2[ifl][k] = p.xfxQ(fl[ifl],x,q)
            if(iset==2):
                fit3[ifl][k] = p.xfxQ(fl[ifl],x,q)
            if(iset==3):
                fit4[ifl][k] = p.xfxQ(fl[ifl],x,q)
 

print("\n\n")

py.clf()
ncols,nrows=3,2
py.figure(figsize=(ncols*5,nrows*3.5))
gs = gridspec.GridSpec(nrows,ncols)
rescolors = py.rcParams['axes.prop_cycle'].by_key()['color']

# pdflabels
labelpdf=[r"$f_{\nu_e}(x_{\nu})$",\
    r"$f_{\bar{\nu}_e}(x_{\nu})$",\
    r"$f_{\nu_\mu}(x_{\nu})$",\
    r"$f_{\bar{\nu}_\mu}(x_{\nu})$",\
    r"$f_{\nu_\tau}(x_{\nu})$",\
    r"$f_{\bar{\nu}_\tau}(x_{\nu})$"]


for ifl in range(nfl):
   
    ax = py.subplot(gs[ifl])

    p1=ax.plot(X,
               fit1[ifl],
               ls="dashed",
               color="C1"
    )
    
    p2=ax.plot(X,
               fit2[ifl],
               ls="dotted",
               color="C2"
    )

    p3=ax.plot(X,
               fit3[ifl],
               ls="solid",
               color="C3"
    )

    ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(1e-1,5e5)
         
    ax.tick_params(which='both',direction='in',labelsize=12,right=True)
    ax.tick_params(which='major',length=7)
    ax.tick_params(which='minor',length=4)
    ax.set_ylabel(labelpdf[ifl],fontsize=16)
    if(ifl>2):
        ax.set_xlabel(r'$x_\nu$',fontsize=18)
    
    # Add the legend
    if(ifl==0):
        ax.legend([p1[0],p2[0],p3[0]],\
                  [pdfsetlab[0],pdfsetlab[1],pdfsetlab[2]], \
                  frameon=True,loc=2,prop={'size':11})
        

py.tight_layout(pad=1.4, w_pad=1.0, h_pad=1.0)
py.savefig('nupdfs-nue-exps.pdf')
print('output plot: nupdfs-nue-exps.pdf')
             
exit()




