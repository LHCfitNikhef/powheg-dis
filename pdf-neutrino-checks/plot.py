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
xmin = 1e-5
xmax= 0.99
nx=1000
X = np.logspace(log(xmin),log(xmax),nx)
lhapdf.setVerbosity(0)
# number of experiments
nexp=4
# number of neutrino flavours
nfl=1
fl=[12,-12,14,-14,16,-15]

pdfset=["faserv",\
        "snd",\
        "faserv2",\
        "flare10"]
pdfsetlab=[r"${\rm FASER}\nu$",\
           r"${\rm SND@LHC}$",\
           r"${\rm FASER}\nu 2$",\
           r"${\rm FLArE}$"]

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
            
    # Run over x arrat
    for k in range(nx):
            
        x = X[k]
        q = 10 # Random number

        # run over flavours
        for ifl in range(nfl):

            if(ifl==0):
                # Electron neutrino
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
ncols,nrows=1,1
py.figure(figsize=(ncols*5,nrows*3.5))
gs = gridspec.GridSpec(nrows,ncols)
rescolors = py.rcParams['axes.prop_cycle'].by_key()['color']

# pdflabels
labelpdf=[r"$f_{\nu_e}(x_{\nu})$"]

icount=0
for ifl in range(nfl):

    ax = py.subplot(gs[icount])

    
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

    p4=ax.plot(X,
               fit4[ifl],
               ls="dashdot",
               color="C4"
    )
     

    ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(1e5,1e14)
         
    ax.tick_params(which='both',direction='in',labelsize=12,right=True)
    ax.tick_params(which='major',length=7)
    ax.tick_params(which='minor',length=4)
    ax.set_ylabel(labelpdf[ifl],fontsize=16)
    ax.set_xlabel(r'$x_\nu$',fontsize=18)
    
    # Add the legend
    if(ifl==0):
        ax.legend([p1[0],p2[0],p3[0],p4[0]],\
                  [pdfsetlab[0],pdfsetlab[1],pdfsetlab[2],pdfsetlab[3]], \
                  frameon=True,loc=3,prop={'size':11})
        
    icount = icount + 1

py.tight_layout(pad=1.4, w_pad=1.0, h_pad=1.0)
py.savefig('neutrino-pdfs-comp.pdf')
print('output plot: neutrino-pdfs-comp.pdf')
             
exit()




