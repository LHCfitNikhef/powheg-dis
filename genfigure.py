import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

plt.rcParams['text.usetex'] = True

def read_hist(filename):
   low,high,val,err = np.loadtxt(filename, unpack=True)
   low = np.insert(low,len(low),high[-1],axis=0)
   val=np.insert(val,len(val),val[-1],axis=0)
   return low,high,val,err


def main():
   # set the names here
   LO_filename_n = "LO_nu_en_dipol_Ehcen.dat"
   LO_filename_min_n = "LO_nu_en_dipol_Ehmin.dat"
   LO_filename_max_n = "LO_nu_en_dipol_Ehmax.dat"
   LO_filename_p = "LO_nu_ep_dipol_Ehcen.dat"
   LO_filename_min_p = "LO_nu_ep_dipol_Ehmin.dat"
   LO_filename_max_p = "LO_nu_ep_dipol_Ehmax.dat"
   LO_low,LO_high,LO_nval,LO_err = read_hist(LO_filename_n)
   LO_low_min,LO_high_min,LO_nval_min,LO_err = read_hist(LO_filename_min_n)
   LO_low_max,LO_high_max,LO_nval_max,LO_err = read_hist(LO_filename_max_n)
   LO_low,LO_high,LO_pval,LO_err = read_hist(LO_filename_p)
   LO_low_min,LO_high_min,LO_pval_min,LO_err = read_hist(LO_filename_min_p)
   LO_low_max,LO_high_max,LO_pval_max,LO_err = read_hist(LO_filename_max_p)
   LO_val=74/183*LO_pval+(183-74)/183*LO_nval
   LO_val_min=74/183*LO_pval_min+(183-74)/183*LO_nval_min
   LO_val_max=74/183*LO_pval_max+(183-74)/183*LO_nval_max
   # names for NLO files
   NLO_filename_n = "NLO_nu_en_dipol_Ehcen.dat"
   NLO_filename_min_n = "NLO_nu_en_dipol_Ehmin.dat"
   NLO_filename_max_n = "NLO_nu_en_dipol_Ehmax.dat"
   NLO_filename_p = "NLO_nu_ep_dipol_Ehcen.dat"
   NLO_filename_min_p = "NLO_nu_ep_dipol_Ehmin.dat"
   NLO_filename_max_p = "NLO_nu_ep_dipol_Ehmax.dat"
   NLO_low,NLO_high,NLO_nval,NLO_err = read_hist(NLO_filename_n)
   NLO_low_min,NLO_high_min,NLO_nval_min,NLO_err = read_hist(NLO_filename_min_n)
   NLO_low_max,NLO_high_max,NLO_nval_max,NLO_err = read_hist(NLO_filename_max_n)
   NLO_low,NLO_high,NLO_pval,NLO_err = read_hist(NLO_filename_p)
   NLO_low_min,NLO_high_min,NLO_pval_min,NLO_err = read_hist(NLO_filename_min_p)
   NLO_low_max,NLO_high_max,NLO_pval_max,NLO_err = read_hist(NLO_filename_max_p)
   NLO_val=74/183*NLO_pval+(183-74)/183*NLO_nval
   NLO_val_min=74/183*NLO_pval_min+(183-74)/183*NLO_nval_min
   NLO_val_max=74/183*NLO_pval_max+(183-74)/183*NLO_nval_max

   # set size of figure
   fig=plt.figure(figsize=(3.4,3.4))
   gs=gridspec.GridSpec(2,1,height_ratios=[3,1])
   gs.update(left=0.16,right=0.95,top=0.95,hspace=0.12)
   ax=fig.add_subplot(gs[0])
   ax.grid(color="grey", linestyle='-', linewidth=0.25)
   ax.set_axisbelow(True)
   plt.setp(ax.get_xticklabels(), visible=False)
   ax.tick_params(which='both',left=True,right=True,bottom=True,top=True, direction='in')

   lines=[]
   p1, = ax.plot(LO_low,LO_val,drawstyle="steps-post",label=r"$\rm POWHEG+PYTHIA8\; LO$",color="blue",linewidth=0.5)
   ax.fill_between(LO_low,LO_val_min,LO_val_max,step="post",alpha=0.35,color="blue",linewidth=0.4,edgecolor="blue")
   p2, = [plt.Rectangle((0, 0), 0, 0, facecolor="blue",linewidth=0.4, edgecolor="blue",alpha=0.35)]
   lines+=((p1,p2),)

   p1, = ax.plot(NLO_low,NLO_val,drawstyle="steps-post",label=r"$\rm POWHEG+PYTHIA8\; NLO$",color="orange",linewidth=0.5)
   ax.fill_between(NLO_low,NLO_val_min,NLO_val_max,step="post",alpha=0.35,color="orange",linewidth=0.4,edgecolor="orange")
   p2, = [plt.Rectangle((0, 0), 0, 0, facecolor="orange",linewidth=0.4, edgecolor="orange",alpha=0.35)]
   lines+=((p1,p2),)
   labels=[r"$\rm LO\, dipole$",r"$\rm NLO \; dipole$"]
   ax.legend(lines, labels, loc='best', numpoints=1)
   ax.tick_params(axis="x", which="both", direction="in")
   ax.tick_params(axis="y", which="both", direction="in")

   # Adapt scales here ... set labels etc.
   #ax.set_yscale("log")
   #ax.set_ylim(6E2,2E5)
   #ax.set_ylim(0.175,0.2750)
   ax.set_ylim(0.15,0.35)
   #ax.set_xlim(0,0.011) # theta

   axr=fig.add_subplot(gs[1])
   axr.grid(color="grey", linestyle='-', linewidth=0.25)
   axr.tick_params(axis="x", which="both", direction="in")
   axr.tick_params(axis="y", which="both", direction="in")
   norm=LO_val
   axr.plot(LO_low,LO_val/norm,drawstyle="steps-post",color="blue",linewidth=0.5)
   axr.fill_between(LO_low,LO_val_min/norm,LO_val_max/norm,step="post",alpha=0.35,color="blue",linewidth=0.4)
   axr.plot(NLO_low,NLO_val/norm,drawstyle="steps-post",color="orange",linewidth=0.5)
   axr.fill_between(NLO_low,NLO_val_min/norm,NLO_val_max/norm,step="post",alpha=0.35,color="orange",linewidth=0.4)
   axr.tick_params(which='both',left=True,right=True,bottom=True,top=True, direction='in', zorder=200)

   # LOWER PANEL labels, scaling...
   # SCALE
   #axr.set_xlim(0,0.011) # theta
   #axr.set_ylim([0.75,1.25]) #Q2
   #axr.set_ylim([0.90,1.10]) #El
   axr.set_ylim([0.80,1.20]) #Eh
   #axr.set_ylim([0.50,1.50]) #theta
   # LABELS
   ax.set_ylabel(r"$d N_{ev} / d E_h$")
   #ax.set_ylabel(r"$d N_{ev} / d Q^2$") #Q2
   axr.set_ylabel(r"$\rm Ratio\; to\; LO$")
   #axr.set_xlabel(r"$Q^2$") #Q2
   #axr.set_xlabel(r"$\vartheta _\ell$") #theta
   axr.set_xlabel(r"$E_h$") #Eh

   plt.savefig("faser_Eh.pdf")

if __name__ == "__main__":
    main()

