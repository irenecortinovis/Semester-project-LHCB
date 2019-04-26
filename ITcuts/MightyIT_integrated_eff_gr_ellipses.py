#!/usr/bin/env python
# coding: utf-8

# In[2]:


import ROOT
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import root_pandas as rp
import root_numpy as rn
import itertools
from __future__ import division


# In[19]:


#filenames = ["mag1_nocuts.root", "mag1_250_490.root", "mag1_300_540.root", "mag1_350_590.root"]
#filenames = ["Reco_magdown_1b_250_450.root", "Reco_magdown_1b_300_540.root", "Reco_magdown_1b_350_630.root", "Reco_magdown_1b_275_495.root", "Reco_magdown_1b_325_585.root", "Reco_magdown_1b_nocuts.root"]

filenames_bfix = ["Reco_magdown1b_ellipse_250_125.root", "Reco_magdown1b_ellipse_300_125.root", 
                  "Reco_magdown1b_ellipse_350_125.root", "Reco_magdown1b_ellipse_400_125.root", 
                  "Reco_magdown1b_ellipse_450_125.root", "Reco_magdown1b_ellipse_500_125.root",
                  "Reco_magdown1b_ellipse_550_125.root", "Reco_magdown1b_ellipse_600_125.root",
                  "Reco_magdown1b_ellipse_650_125.root"]
xaxis_bfix = [250,300,350,400,450,500,550,600,650]

filenames_afix = ["Reco_magdown1b_ellipse_300_100.root", "Reco_magdown1b_ellipse_300_125.root", 
                  "Reco_magdown1b_ellipse_300_150.root", "Reco_magdown1b_ellipse_300_200.root",
                  "Reco_magdown1b_ellipse_300_250.root", "Reco_magdown1b_ellipse_300_300.root",
                  "Reco_magdown1b_ellipse_300_350.root", "Reco_magdown1b_ellipse_300_400.root"]
xaxis_afix = [100,125,150,200,250,300,350,400]

filenamelists = [filenames_bfix, filenames_afix]
xaxislists = [xaxis_bfix, xaxis_afix]
xlabels = ["a in ellipse IT cuts (mm), b fixed to 125 mm", "b in ellipse IT cuts (mm), a fixed to 300 mm"]


# In[20]:


def get_ghostrate(title1, title2):
 
    #nPV_Ghosts
    nums1 = rn.hist2array(title1, copy=True)
    nghosts = np.sum(nums1)
    
    #nPV_Total
    nums2, edges2 = rn.hist2array(title2, return_edges=True, copy=True)
    ntot = np.sum(nums2)
    
    ghostrate = nghosts/ntot
    
    #print nghosts, ntot
    #check why they are slightly wrong
    
    return ntot, ghostrate


# In[21]:


#Ghost rate


fig, axes = plt.subplots(1,len(filenamelists))
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

        totnum, ghostratenum = get_ghostrate(current_item1, current_item2)
        ghostrate.append(ghostratenum)
        #print(filename, (ghostratenum*100))
        totnumber.append(totnum)
        #labelxbar.append((filename.strip(".root")).strip("mag1_"))
        #labelxbar.append((filename.strip(".root")).strip("Reco_magdown_1b_"))
        labelxbar.append((filename.strip(".root")).strip(" Reco_magdown1b_ellipse_"))
    axes[filenamelists.index(filenamelist)].set_title("Ghost rate")
    #axes[filenamelists.index(filenamelist)].plot(xaxis, ghostrate, marker='x', linewidth=1, markersize=10)
    axes[filenamelists.index(filenamelist)].scatter(xaxis, ghostrate, marker='x', linewidth=1)
    axes[filenamelists.index(filenamelist)].set_xlabel(xlabels[filenamelists.index(filenamelist)])
    axes[filenamelists.index(filenamelist)].set_ylabel("Ghost rate")


    fig.savefig("Ghost rate", bbox_inches='tight')


plt.show()


# In[22]:


def get_efficiency(title1, title2):
         
    #nPV_reconstructed
    nums1, edges1 = rn.hist2array(title1, return_edges=True, copy=True)
    n_ed = np.sum(nums1)
    
    #nPV_reconstructible
    nums2, edges2 = rn.hist2array(title2, return_edges=True, copy=True)
    n_ible = np.sum(nums2)
    
    efficiency = n_ed/n_ible
    
    return efficiency   


# In[24]:


#efficiencies



commonpath = "Track/PrChecker/TTrack/"
types = ["hasT_", "long_", "long>5GeV_"]
end_path1 = "Eta_reconstructed"
end_path2 = "Eta_reconstructible"





fig, axes = plt.subplots(len(filenamelists),len(types))
fig.set_figheight(8)
fig.set_figwidth(15)
fig.suptitle("Reconstruction efficiency")


for typetrack in types:
    
    
    path1 = commonpath + typetrack + end_path1
    path2 = commonpath + typetrack + end_path2
    paths_list1 = path1.split("/")
    paths_list2 = path2.split("/")
    
    for filenamelist in filenamelists:
        xaxis = xaxislists[filenamelists.index(filenamelist)]
        efficiencies = []
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

            efficiencies.append(get_efficiency(current_item1, current_item2))
            #print(filename, paths_list1, get_efficiency(current_item1, current_item2)*100)

        
        axes[filenamelists.index(filenamelist), types.index(typetrack)].set_title(typetrack.strip("_"))
        #axes[filenamelists.index(filenamelist), types.index(typetrack)].plot(xaxis, efficiencies, marker='x', linewidth=1, markersize=10)
        axes[filenamelists.index(filenamelist), types.index(typetrack)].scatter(xaxis, efficiencies, marker='x', linewidth=1)

        axes[filenamelists.index(filenamelist), types.index(typetrack)].set_xlabel(xlabels[filenamelists.index(filenamelist)])
        axes[filenamelists.index(filenamelist), types.index(typetrack)].set_ylabel("Efficiency")
        
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.4)

    fig.savefig("Efficiencies", bbox_inches='tight')


# In[ ]:




