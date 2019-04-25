#!/usr/bin/env python
# coding: utf-8

# In[19]:


from ROOT import *
import numpy as np
import math

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


# In[199]:


ch1=TChain("MCParticleNTuple/Tracks")
ch1.Add("~/MightyIT/MCtracks_oldfile_200ev.root")
ch2=TChain("MCParticleNTuple/Tracks")
ch2.Add("~/MightyIT/MCtracks_oldfile_200ev.root")


zmag = 5092.5

#sigma for chi2 calculation (found in chi2_strategies-sigmas.py)

sigmax_nosm=9.2
sigmay_nosm=5.4
sigmatx_nosm=0.0033
sigmaty_nosm=0.0017

sigmax_wsm=9.5
sigmay_wsm=22.3
sigmatx_wsm=0.0034
sigmaty_wsm=0.0034



# pixel smearing in mightyIT
resx_mit = 0.1/TMath.Sqrt(12)
resy_mit = 0.4/TMath.Sqrt(12)
# strips smearing outside mightyIT
resx_sf = 0.25/TMath.Sqrt(12)
resy_sf = 0.25/TMath.Cos(5.*TMath.Pi()/180.) #takes into account 5° tilting






        


# In[200]:


#function that checks if a point (format: (x,y)) is inside the MightyIT region

def MightyITregion():
    dx = 540
    dy = 200
    
    polygon = Polygon([(-4*dx,1*dy),(-3*dx,1*dy),(-3*dx,2*dy),(-2*dx,2*dy),(-2*dx,3*dy),(-1*dx,3*dy),(-1*dx,4*dy),
                      (1*dx,4*dy),(1*dx,3*dy),(2*dx,3*dy),(2*dx,2*dy),(3*dx,2*dy),(3*dx,1*dy),(4*dx,1*dy),
                      (4*dx,-1*dy),(3*dx,-1*dy),(3*dx,-2*dy),(2*dx,-2*dy),(2*dx,-3*dy),(1*dx,-3*dy),(1*dx,-4*dy),
                      (-1*dx,-4*dy),(-1*dx,-3*dy),(-2*dx,-3*dy),(-2*dx,-2*dy),(-3*dx,-2*dy),(-3*dx,-1*dy),(-4*dx,-1*dy)])
    return polygon
    
def is_inregion(point):
    point = Point(point[0], point[1])
           
    region = MightyITregion()
    if region.contains(point):
        return True
    else:
        return False
    


# In[201]:


def ChiSquareVeloT(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1, zmag, smearing=False):
    
    #gaussian smearing - pixel resolutions x and y according to if in mightyit region or outside
    if smearing == True:
        if (is_inregion([xt0,yt0]) == True):
            xt0 = gRandom.Gaus(xt0,resx_mit)
            yt0 = gRandom.Gaus(yt0,resy_mit)
        else:
            xt0 = gRandom.Gaus(xt0,resx_sf)
            yt0 = gRandom.Gaus(yt0,resy_sf)
            
        if (is_inregion([xt1,yt1]) == True):
            xt1 = gRandom.Gaus(xt1,resx_mit)
            yt1 = gRandom.Gaus(yt1,resy_mit)
        else:
            xt1 = gRandom.Gaus(xt1,resx_sf)
            yt1 = gRandom.Gaus(yt1,resy_sf)
            
 
        
    tyV=(yv1 - yv0)/(zv1 - zv0)
    txV=(xv1 - xv0)/(zv1 - zv0)
    tyT=(yt1 - yt0)/(zt1 - zt0)
    txT=(xt1 - xt0)/(zt1 - zt0)

    #zmag calculations
    xVzmag=xv1+(zmag-zv1)*txV
    yVzmag=yv1+(zmag-zv1)*tyV
       
    
    #chi2 calculations
    txPre=(xt0-xVzmag)/(zt0-zmag)
    tyPre=(yt0-yVzmag)/(zt0-zmag)

    #extrapolated point from zmag to t station, with the same slope of the t track
    xT_ex=xVzmag+(zt0-zmag)*txT
    #from velo to t track with velo slope
    yT_ex=yv1+(zt0-zv1)*tyT
    
    chi2 = (txPre-txT)*(txPre-txT)/sigma2tx+(tyPre-tyT)*(tyPre-tyT)/sigma2ty + (yT_ex-yt0)*(yT_ex-yt0)/sigma2y + (xT_ex-xt0)*(xT_ex-xt0)/sigma2x 


    return chi2, (xT_ex-xt0), (yT_ex-yt0), (txPre-txT), (tyPre-tyT)


# In[202]:


Nbins=1000; nrange=5
chihisto_nosm=TH1F("chihisto_nosm","",Nbins,0, nrange)
chihisto_wsm=TH1F("chihisto_wsm","",Nbins,0, nrange)



Nbins=1000
histolimx = 50
histolimy = 50
scatterx_nosm=TH1F("scatterx_nosm","",Nbins,-histolimx, +histolimx)
scattery_nosm=TH1F("scattery_nosm","",Nbins,-histolimy, +histolimy)



Nbins=1000
histolimtx = 0.01
histolimty = 0.01
scattertx_nosm=TH1F("scattertx_nosm","",Nbins,-histolimtx, +histolimtx)
scatterty_nosm=TH1F("scatterty_nosm","",Nbins,-histolimty, +histolimty)



Nbins=1000
histolimx = 50
histolimy = 50
scatterx_wsm=TH1F("scatterx_wsm","",Nbins,-2*histolimx, +2*histolimx)
scattery_wsm=TH1F("scattery_wsm","",Nbins,-2*histolimy, +2*histolimy)



Nbins=1000
histolimtx = 0.01
histolimty = 0.01
scattertx_wsm=TH1F("scattertx_wsm","",Nbins,-2*histolimtx, +2*histolimtx)
scatterty_wsm=TH1F("scatterty_wsm","",Nbins,-2*histolimty, +2*histolimty)


# In[172]:


def reject_outliers(data, nsigma=3):
    data = np.array(data)
    return data[abs(data - np.mean(data)) < nsigma * np.std(data)]


# In[205]:


stdsx_nosm = []
stdsy_nosm = []
stdstx_nosm = []
stdsty_nosm = []

stdsx_wsm = []
stdsy_wsm = []
stdstx_wsm = []
stdsty_wsm = []

#rangex =range(100,2100,100)+range(2000,20000,1000)
rangex =range(8000,8001,1)
for maxNevent in rangex:
    
    print('\n', maxNevent)

    list_chi_wsm = []
    list_chi_nosm = []
    list_scx_wsm = []
    list_scx_nosm = []
    list_scy_wsm = []
    list_scy_nosm = []
    list_sctx_wsm = []
    list_sctx_nosm = []
    list_scty_wsm = []
    list_scty_nosm = []



    #TESTING SEARCH WINDOWS
    #after noting an asymmetry in distribution depending on the sign of the charge, the zmag has been found from the data. around 5100
    #fine tuning: looking for the best zmag coordinate to have a symmetric distribution of deltax and deltatx with charge>0

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

            #testing smearing
            for smearing in [False,True]:

                #change sigmas and windows depending on smearing on/off
                if smearing == False:
                    #sigma squared
                    sigma2x=pow(sigmax_nosm,2)
                    sigma2y=pow(sigmay_nosm,2)
                    sigma2tx=pow(sigmatx_nosm,2)
                    sigma2ty=pow(sigmaty_nosm,2)

                else:
                    #sigma squared
                    sigma2x=pow(sigmax_wsm,2)
                    sigma2y=pow(sigmay_wsm,2)
                    sigma2tx=pow(sigmatx_wsm,2)
                    sigma2ty=pow(sigmaty_wsm,2)

                chi2_smtest = ChiSquareVeloT(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1 , zmag, smearing=smearing)
                if smearing == True:
                    chihisto_wsm.Fill(chi2_smtest[0])
                    list_chi_wsm.append(chi2_smtest[0])
                    scatterx_wsm.Fill(chi2_smtest[1])
                    list_scx_wsm.append(chi2_smtest[1])
                    scattery_wsm.Fill(chi2_smtest[2])
                    list_scy_wsm.append(chi2_smtest[2])
                    scattertx_wsm.Fill(chi2_smtest[3])
                    list_sctx_wsm.append(chi2_smtest[3])
                    scatterty_wsm.Fill(chi2_smtest[4])
                    list_scty_wsm.append(chi2_smtest[4])
                else:
                    chihisto_nosm.Fill(chi2_smtest[0])
                    list_chi_nosm.append(chi2_smtest[0])
                    scatterx_nosm.Fill(chi2_smtest[1])
                    list_scx_nosm.append(chi2_smtest[1])
                    scattery_nosm.Fill(chi2_smtest[2])
                    list_scy_nosm.append(chi2_smtest[2])
                    scattertx_nosm.Fill(chi2_smtest[3])
                    list_sctx_nosm.append(chi2_smtest[3])
                    scatterty_nosm.Fill(chi2_smtest[4])
                    list_scty_nosm.append(chi2_smtest[4])
                    
                    #testing asymmetry: add if event.qop>0



        if(nEvent > maxNevent):
            break
            
    stdsx_nosm.append(np.std(reject_outliers(list_scx_nosm)))
    stdsy_nosm.append(np.std(reject_outliers(list_scy_nosm)))
    stdstx_nosm.append(np.std(reject_outliers(list_sctx_nosm)))
    stdsty_nosm.append(np.std(reject_outliers(list_scty_nosm)))
    
    stdsx_wsm.append(np.std(reject_outliers(list_scx_wsm)))
    stdsy_wsm.append(np.std(reject_outliers(list_scy_wsm)))
    stdstx_wsm.append(np.std(reject_outliers(list_sctx_wsm)))
    stdsty_wsm.append(np.std(reject_outliers(list_scty_wsm)))

    



# In[197]:


plt.scatter(rangex,stdsx_nosm )
plt.show()
plt.scatter(rangex,stdsy_nosm )
plt.show()
plt.scatter(rangex,stdstx_nosm )
plt.show()
plt.scatter(rangex,stdsty_nosm )
plt.show()


# In[195]:


#mean of stds after it becomes stable (around 8000)
index_thr = rangex.index(8000)
print 'no smearing'
print(np.mean(reject_outliers(stdsx_nosm[index_thr:])))
print(np.mean(reject_outliers(stdsy_nosm[index_thr:])))
print(np.mean(reject_outliers(stdstx_nosm[index_thr:])))
print(np.mean(reject_outliers(stdsty_nosm[index_thr:])))

print '\n'
print 'with smearing'
print(np.mean(reject_outliers(stdsx_wsm[index_thr:])))
print(np.mean(reject_outliers(stdsy_wsm[index_thr:])))
print(np.mean(reject_outliers(stdstx_wsm[index_thr:])))
print(np.mean(reject_outliers(stdsty_wsm[index_thr:])))


# In[ ]:





# In[206]:


#print the threshold for the chi2 to include a percentage of long track chi2s when cutting below threshold

perc_list = []
thr_wsm_list = []
thr_nosm_list = []

list_chi_wsm.sort()
list_chi_nosm.sort()

for percentage in range(80,100):
    threshold_wsm = list_chi_wsm[len(list_chi_wsm)/100*percentage]
    threshold_nosm = list_chi_nosm[len(list_chi_nosm)/100*percentage]
    perc_list.append(percentage)
    thr_wsm_list.append(round(threshold_wsm,2))
    thr_nosm_list.append(round(threshold_nosm,2))
    
    
print(perc_list, thr_nosm_list, thr_wsm_list)


# In[ ]:


c9=TCanvas("c9","",1300,900)
chihisto_nosm.SetLineColor(kBlue)
chihisto_nosm.Draw()
chihisto_wsm.SetLineColor(kRed)
chihisto_wsm.Draw("SAME")
c9.SaveAs("dist_chi2_wwo_smearing.png")


# In[14]:


#looking for a threshold in the chi2 to include most real chi2 in the final algorithm

list_chi.sort()
percentage = 95
threshold = list_chi[len(list_chi)/100*percentage]
print(threshold)

#3.3698
#79.91


# In[307]:


#plot chi distribution

c4=TCanvas("c4","",1300,900)
chihisto.Draw()
c4.SaveAs("dist_testing_chi2.png")


# In[98]:


#plot deltax, y tx, ty distributions



c5=TCanvas("c5","",1300,900)

c5.Divide(2,2)

c5.cd(1)
funcx= TF1("funcx", "gaus", -histolimx, +histolimx)
funcx.SetParameters(500,0,0)
fitresult_x = scatterx_nosm.Fit(funcx)
print(funcx.GetParameter(2))
#scatterx_nosm.Draw()




c5.cd(2)
funcy = TF1("funcy", "gaus", -histolimy, +histolimy)
funcy.SetParameters(500,0,0)
fitresult_y = scattery_nosm.Fit(funcy)
print(funcy.GetParameter(2))
#scattery_nosm.Draw()



c5.cd(3)
funct = TF1("functx", "gaus", -histolimtx, +histolimtx)
functx.SetParameters(500,0,0)
fitresult_tx = scattertx_nosm.Fit(functx)
print(functx.GetParameter(2))
#scattertx_nosm.Draw()

c5.cd(4)
functy = TF1("functy", "gaus", -histolimty, +histolimty)
functy.SetParameters(500,0,0)
fitresult_ty = scatterty_nosm.Fit(functy)
print(functy.GetParameter(2))
scatterty_nosm.Draw()


c5.SaveAs("dist_testing_deltas_nosmearing.png")


# In[94]:


#plot deltax, y tx, ty distributions



c1=TCanvas("c5","",1300,900)

c1.Divide(2,2)

c1.cd(1)
funcx= TF1("funcx", "gaus", -histolimx, +histolimx)
funcx.SetParameters(500,0,0)
fitresult_x = scatterx_wsm.Fit(funcx)
print(funcx.GetParameter(2))
scatterx_wsm.Draw()




c1.cd(2)
funcy = TF1("funcy", "gaus", -histolimy, +histolimy)
funcy.SetParameters(500,0,0)
fitresult_y = scattery_wsm.Fit(funcy)
print(funcy.GetParameter(2))
scattery_wsm.Draw()



c1.cd(3)
funct = TF1("functx", "gaus", -histolimtx, +histolimtx)
functx.SetParameters(500,0,0)
fitresult_tx = scattertx_wsm.Fit(functx)
print(functx.GetParameter(2))
scattertx_wsm.Draw()

c1.cd(4)
functy = TF1("functy", "gaus", -histolimty, +histolimty)
functy.SetParameters(500,0,0)
fitresult_ty = scatterty_wsm.Fit(functy)
print(functy.GetParameter(2))
scatterty_wsm.Draw()


c1.SaveAs("dist_testing_deltas_smearing.png")


# In[ ]:




