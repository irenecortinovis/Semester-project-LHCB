#!/usr/bin/env python
# coding: utf-8

# In[92]:


from ROOT import *
import numpy as np
import math
from shapely.geometry import LineString


# In[93]:


ch1=TChain("MCParticleNTuple/Tracks")
ch1.Add("~/MightyIT/MCtracks_oldfile_200ev.root")
ch2=TChain("MCParticleNTuple/Tracks")
ch2.Add("~/MightyIT/MCtracks_oldfile_200ev.root")

gROOT.ProcessLine(".x ~/lhcbStyle.C")
gStyle.SetPaintTextFormat("1.3f")


# In[94]:


#function that finds intersection of velo and t track for a longtrack, on the bending plane
#outputs the z coordinate (zmag)

def findZmag(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1): #to ask: velo p?
    tyV=(yv1 - yv0)/(zv1 - zv0)
    txV=(xv1 - xv0)/(zv1 - zv0)
    tyT=(yt1 - yt0)/(zt1 - zt0)
    txT=(xt1 - xt0)/(zt1 - zt0)
    
    #extrapolate velo track from v1 to t0 with velo slope
    xV_zt0 = xv1 + txV*(zt0-zv1)
    
    #extrapolate t track from t0 to v1 with t track slope
    xT_zv1 = xt0 + txT*(zv1-zt0)
    
    velo_extr = LineString([(xv1,zv1), (xV_zt0,zt0)])
    t_extr = LineString([(xT_zv1,zv1), (xt0,zt0)])
    
    #todo: sometimes no intersection, but this is not optimized?
    if velo_extr.intersection(t_extr):
        zmag = velo_extr.intersection(t_extr).y #it is actually z
    else:
        zmag = -1
   
    return zmag


# In[97]:



Nbins_p=100
Nbins_zmag=100
zmagvsp=TH2F("zmagvsp","",Nbins_p,Nbins_zmag,0,5000,4000,6000)
list_zmag = []


# In[98]:


#find zmag for many long tracks

nEvent = 0
for event in ch1:
    
    #selecting some long tracks for event in ch1 and cut for !=0 else division/0
    if (event.HitVeloZpos[1]-event.HitVeloZpos[0] != 0 and event.HitZpos[1]-event.HitZpos[0] != 0
       and  0.<event.HitVeloZpos[0]<event.HitVeloZpos[1]<800. and 7825.<event.HitZpos[0]<7875. and 7900.<event.HitZpos[1]<7950. and event.p > 5000):
        nEvent = nEvent+1
        
        xv0 = event.HitVeloXpos[0]; yv0 = event.HitVeloYpos[0]; zv0 = event.HitVeloZpos[0];
        xv1 = event.HitVeloXpos[1]; yv1 = event.HitVeloYpos[1]; zv1 = event.HitVeloZpos[1];
        xt0 = event.HitXpos[0]; yt0 = event.HitYpos[0]; zt0 = event.HitZpos[0];
        xt1 = event.HitXpos[1]; yt1 = event.HitYpos[1]; zt1 = event.HitZpos[1];
        vp = event.p
        
      
        zmag = findZmag(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1) 
        #zmagvsp.Fill(zmag[0], zmag[1])
        zmagvsp.Fill(vp, zmag)
        list_zmag.append(zmag)
        #cut at n good events (computational time issue)
        if nEvent > 100000:
            break


# In[99]:


#plot zmag as function of momentum

c5=TCanvas("c5","",1400,1000)
zmagvsp.Fit("pol0")
zmagvsp.Draw()
c5.SaveAs("dist_testing_zmag.png")



# In[102]:


#mean and sigma for the distribution (since it is quite constant)
print(np.mean(list_zmag))
print(np.std(list_zmag))


# In[ ]:




