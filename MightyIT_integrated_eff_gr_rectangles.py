#!/usr/bin/env python
# coding: utf-8

# In[15]:


import ROOT
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import root_pandas as rp
import root_numpy as rn
import itertools
from __future__ import division


# In[16]:


filenames_bfix = [ "nocuts_patchecker.root","530_100_patchecker.root", "530_200_patchecker.root",
                  "530_300_patchecker.root", "530_400_patchecker.root",
                  "530_500_patchecker.root", "1060_300_patchecker.root"]
xaxis_bfix = [0,200*1060,400*1060,600*1060,800*1060,1000*1060, 600*2120]
xaxis_bfix_labels = ["no cuts", "200*1060","400*1060","600*1060","800*1060","1000*1060", "600*2120"]


filenamelists = [filenames_bfix]
xaxislists = [xaxis_bfix]
xlabellists = [xaxis_bfix_labels]
xlabels = ["Area in IT cuts (mm^2)"]


# In[63]:


def get_ghostrate(title1, title2):
 
    #nPV_Ghosts
    nums1 = rn.hist2array(title1, copy=True)
    nghosts = np.sum(nums1)
    rel_errghosts = np.divide(np.sqrt(nghosts),nghosts)
    
    #nPV_Total
    nums2, edges2 = rn.hist2array(title2, return_edges=True, copy=True)
    ntot = np.sum(nums2)
    rel_errtot = np.divide(np.sqrt(ntot),ntot)
    
    ghostrate = nghosts/ntot
    ghostrate_err = ghostrate*(np.sqrt(np.add(np.square(rel_errghosts),np.square(rel_errtot))))
    print(ghostrate_err/ghostrate*100)
    #print nghosts, ntot
    
    return ntot, ghostrate, ghostrate_err


# In[64]:


#Ghost rate

fig, axes = plt.subplots()
fig.set_figheight(5)
fig.set_figwidth(10)
fig.suptitle("Ghost rate in ellipse IT cut")


path1 = "Track/PrChecker/TTrack/Eta_Ghosts"
path2 = "Track/PrChecker/TTrack/Eta_Total"
paths_list1 = path1.split("/")
paths_list2 = path2.split("/")

for filenamelist in filenamelists:
    totnumber = []
    ghostrate = []
    labelxbar = []
    for filename in filenamelist:
        xaxis = xaxislists[filenamelists.index(filenamelist)]
        f = ROOT.TFile.Open(filename, 'read')
        first = True
        for ipathlet in paths_list1:
            if first:
                current_item1 = f.Get(ipathlet)
                first = False
            else:
                current_item1 = current_item1.Get(ipathlet)
            if str(type(current_item1))=="<class '__main__.TH1D'>":
                break

        first = True
        for ipathlet in paths_list2:
            if first:
                current_item2 = f.Get(ipathlet)
                first = False
            else:
                current_item2 = current_item2.Get(ipathlet)
            if str(type(current_item2))=="<class '__main__.TH1D'>":
                break
        
        totnum, ghostratenum, ghostratenum_err = get_ghostrate(current_item1, current_item2)
        ghostrate.append(ghostratenum)
        #print(filename, (ghostratenum*100))
        totnumber.append(totnum)
        labelxbar.append((filename.strip(".root")).strip("_patchecker.root"))
    axes.set_title("Ghost rate")
    #axes.scatter(xaxis, ghostrate, marker='x', linewidth=1)
    axes.errorbar(xaxis, ghostrate, yerr=ghostratenum_err, fmt='x', ecolor='r', linewidth=1)
    axes.set_xticks(xaxis, minor=False)
    axes.set_xticklabels(xlabellists[filenamelists.index(filenamelist)])
    
    axes.set_ylabel("Ghost rate")


    fig.savefig("Ghost rate", bbox_inches='tight')


plt.show()


# In[65]:


def get_efficiency(title1, title2):
         
    #nPV_reconstructed
    nums1, edges1 = rn.hist2array(title1, return_edges=True, copy=True)
    n_ed = np.sum(nums1)
    rel_errned = np.divide(np.sqrt(n_ed),n_ed)
    
    #nPV_reconstructible
    nums2, edges2 = rn.hist2array(title2, return_edges=True, copy=True)
    n_ible = np.sum(nums2)
    rel_errnible = np.divide(np.sqrt(n_ed),n_ed)

    
    efficiency = n_ed/n_ible
    efficiency_err = efficiency*(np.sqrt(np.add(np.square(rel_errned),np.square(rel_errnible))))
    
    print(efficiency_err/efficiency*100)
    return efficiency, efficiency_err


# In[66]:


#efficiencies

commonpath = "Track/PrChecker/TTrack/"
types = ["hasT_", "long_", "long>5GeV_"]
end_path1 = "Eta_reconstructed"
end_path2 = "Eta_reconstructible"



fig, axes = plt.subplots(len(types),len(filenamelists))
fig.set_figheight(16)
fig.set_figwidth(10)
fig.suptitle("Reconstruction efficiency")


for typetrack in types:
    
    
    path1 = commonpath + typetrack + end_path1
    path2 = commonpath + typetrack + end_path2
    paths_list1 = path1.split("/")
    paths_list2 = path2.split("/")
    
    for filenamelist in filenamelists:
        xaxis = xaxislists[filenamelists.index(filenamelist)]
        efficiencies = []
        efficiencies_err = []
        labelxbar = []
        for filename in filenamelist:
            f = ROOT.TFile.Open(filename, 'read')
            first = True
            for ipathlet in paths_list1:
                if first:
                    current_item1 = f.Get(ipathlet)
                    first = False
                else:
                    current_item1 = current_item1.Get(ipathlet)
                if str(type(current_item1))=="<class '__main__.TH1D'>":
                    break

            first = True
            for ipathlet in paths_list2:
                if first:
                    current_item2 = f.Get(ipathlet)
                    first = False
                else:
                    current_item2 = current_item2.Get(ipathlet)
                if str(type(current_item2))=="<class '__main__.TH1D'>":
                    break

            efficiency, efficiencyerror = get_efficiency(current_item1, current_item2)
            efficiencies.append(efficiency)
            efficiencies_err.append(efficiencyerror)

        
        axes[types.index(typetrack)].set_title(typetrack.strip("_"))
        #axes[types.index(typetrack)].scatter(xaxis, efficiencies, marker='x', linewidth=1)
        axes[types.index(typetrack)].errorbar(xaxis, efficiencies, yerr=efficiencyerror, fmt='x', ecolor='r', linewidth=1)


        axes[types.index(typetrack)].set_xlabel(xlabels[filenamelists.index(filenamelist)])
        axes[types.index(typetrack)].set_ylabel("Efficiency")
        axes[types.index(typetrack)].set_xticks(xaxis, minor=False)
        axes[types.index(typetrack)].set_xticklabels(xlabellists[filenamelists.index(filenamelist)])
  
        
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    fig.savefig("Efficiencies", bbox_inches='tight')


# In[ ]:





# In[ ]:




