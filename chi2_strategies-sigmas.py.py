#!/usr/bin/env python
# coding: utf-8

# In[270]:


from ROOT import *
import numpy as np
import math


# In[282]:


ch1=TChain("MCParticleNTuple/Tracks")
ch1.Add("~/MightyIT/MCtracks_oldfile_200ev.root")
ch2=TChain("MCParticleNTuple/Tracks")
ch2.Add("~/MightyIT/MCtracks_oldfile_200ev.root")

#gROOT.ProcessLine(".x ~/lhcbStyle.C")
#gStyle.SetPaintTextFormat("1.3f")

#first parameter is zmag, to be optimized (fine tuning from 5107 +- 77 found from zmag program)
m_zMagParams=[5287.6, -7.98878, 317.683, 0.0119379, -1418.42]

Par=[4943170,6314610]
ParF=[-0.0001666775,-1304.68]


#tested from scatterx and y
rangex=40
rangey=40

zS=10000.
z0=0.


# In[283]:


def ChiSquareVeloT(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1 , zS, m_zMagParams, Par, rangex, rangey, vp): #to ask: velo p?
    tyV=(yv1 - yv0)/(zv1 - zv0)
    txV=(xv1 - xv0)/(zv1 - zv0)
    tyT=(yt1 - yt0)/(zt1 - zt0)
    txT=(xt1 - xt0)/(zt1 - zt0)

    #zmag calculations
    xTzS=xt1+(zS-zt1)*txT
    dSlopex=txV-txT
    zmag=m_zMagParams[0]+m_zMagParams[1]*abs(dSlopex)+m_zMagParams[2]*dSlopex*dSlopex+m_zMagParams[3]*abs(xTzS)+m_zMagParams[4]*txV*txV
    xVzmag=xv1+(zmag-zv1)*txV
    yVzmag=yv1+(zmag-zv1)*tyV
       
    
    #chi2 calculations
    txPre=(xt0-xVzmag)/(zt0-zmag)
    tyPre=(yt0-yVzmag)/(zt0-zmag)
    sigma2tx=Par[0]/pow(vp,3.)
    sigma2ty=Par[1]/pow(vp,3.)
    
    #chi2 = (txPre-txT)*(txPre-txT)/sigma2tx+(tyPre-tyT)*(tyPre-tyT)/sigma2ty

    #ŧesting: chi2 only in y
    chi2 = (tyPre-tyT)*(tyPre-tyT)/sigma2ty

    #extrapolated point from zmag to t station, with the same slope of the t track
    xT_ex=xVzmag+(zt0-zmag)*txT
    #yT_ex=yVzmag+(zt0-zmag)*tyT

    #from velo to t track with velo slope
    yT_ex=yv1+(zt0-zv1)*tyT

    return chi2, (xT_ex-xt0), (yT_ex-yt0), (txPre-txT), (tyPre-tyT)


# In[302]:


Nbins=1000
chihisto=TH1F("chihisto","",Nbins,0, 3)



Nbins=1000
histolimx = 50
histolimy = rangey
scatterx=TH1F("scatterx","",Nbins,-histolimx, +histolimx)
scattery=TH1F("scattery","",Nbins,-histolimy, +histolimy)



Nbins=1000
histolimtx = 0.01
histolimty = 0.01
scattertx=TH1F("scattertx","",Nbins,-histolimtx, +histolimtx)
scatterty=TH1F("scatterty","",Nbins,-histolimty, +histolimty)


# In[303]:


list_chi = []
list_scx = []
list_scy = []
list_sctx = []
list_scty = []


# In[304]:



#TESTING SEARCH WINDOWS
#after noting an asymmetry in distribution depending on the sign of the charge, the zmag has been found from the data. around 5100
#fine tuning: looking for the best zmag coordinate to have a symmetric distribution of deltax and deltatx with charge>0

m_zMagParams[0] = 5092.5
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

        #testing asymmetry
        #if event.qop>0:
        chi2real = ChiSquareVeloT(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1 , zS, m_zMagParams, Par, rangex, rangey, vp) #to ask: velo p?
        chihisto.Fill(chi2real[0])
        scatterx.Fill(chi2real[1])
        scattery.Fill(chi2real[2])
        scattertx.Fill(chi2real[3])
        scatterty.Fill(chi2real[4])
        
        #finding sigmas (not as a function of momentum)        
        chi2real = ChiSquareVeloT(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1 , zS, m_zMagParams, Par, rangex, rangey, vp) #to ask: velo p?
        list_chi.append(chi2real[0])
        list_scx.append(chi2real[1])
        list_scy.append(chi2real[2])
        list_sctx.append(chi2real[3])
        list_scty.append(chi2real[4])
            
    if(nEvent > 10000):
        break


#found 5092.5 with 100000 events


# In[305]:


print(np.mean(list_scx))
print(np.mean(list_scy))
print(np.mean(list_sctx))
print(np.mean(list_scty))

print(np.std(list_scx))
print(np.std(list_scy))
print(np.std(list_sctx))
print(np.std(list_scty))


# In[316]:


#looking for a threshold in the chi2 to include most real chi2 in the final algorithm

list_chi.sort()
percentage = 95
threshold = list_chi[len(list_chi)/100*percentage]
print(threshold)


# In[307]:


#plot chi distribution

c4=TCanvas("c4","",1300,900)
chihisto.Draw()
c4.SaveAs("dist_testing_chi2.png")


# In[276]:


#plot deltax, y tx, ty distributions

c5=TCanvas("c5","",1300,900)
funcx= TF1("funcx", "gaus", -histolimx, +histolimx)
funcx.SetParameters(500,0,0)
fitresult_x = scatterx.Fit(funcx)
print(funcx.GetParameter(2))

scatterx.Draw()
c5.SaveAs("dist_testing_windowx.png")





c6=TCanvas("c6","",1300,900)

funcy = TF1("funcy", "gaus", -histolimy, +histolimy)
funcy.SetParameters(500,0,0)
fitresult_y = scattery.Fit(funcy)
print(funcy.GetParameter(2))

scattery.Draw()
c6.SaveAs("dist_testing_windowy.png")





c7=TCanvas("c7","",1300,900)
funct = TF1("functx", "gaus", -histolimtx, +histolimtx)
functx.SetParameters(500,0,0)
fitresult_tx = scattertx.Fit(functx)
print(functx.GetParameter(2))
scattertx.Draw()
c7.SaveAs("dist_testing_windowtx.png")




c8=TCanvas("c8","",1300,900)
scatterty.Draw()

functy = TF1("functy", "gaus", -histolimty, +histolimty)
functy.SetParameters(500,0,0)
fitresult_ty = scatterty.Fit(functy)
print(functy.GetParameter(2))

c8.SaveAs("dist_testing_windowty.png")


# In[ ]:




