#!/usr/bin/env python
# coding: utf-8

# In[64]:


from ROOT import *
import numpy as np
import math
import time

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


# In[65]:


ch1=TChain("MCParticleNTuple/Tracks")
ch1.Add("~/MightyIT/MCtracks_oldfile_20ev.root")
ch2=TChain("MCParticleNTuple/Tracks")
ch2.Add("~/MightyIT/MCtracks_oldfile_20ev.root")

'''ch1=TChain("MCParticleNTuple/Tracks")
ch1.Add("~/MightyIT/MCtracks_oldfile_nocuts_20ev.root")
ch2=TChain("MCParticleNTuple/Tracks")
ch2.Add("~/MightyIT/MCtracks_oldfile_nocuts_20ev.root")'''

#77981
#817971

nentries = ch1.GetEntries()
print(nentries)


# In[66]:




gROOT.ProcessLine(".x ~/lhcbStyle.C")
gStyle.SetPaintTextFormat("1.3f")

#zmag optimised (from chi2_strategies-zmag_p.py and fine tuning chi2_strategies-sigmas.py )
zmag=5092.5

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
resy_sf = 0.25/TMath.Cos(5.*np.pi/180.) #takes into account 5° tilting


#chi2 threshold to include percentage long tracks when cutting below chi2 (from chi2_strategies-sigmas.py)
perc_list = [80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]
thr_nosm_list = [4.02, 4.33, 4.61, 4.99, 5.44, 5.93, 6.41, 7.13, 7.96, 9.07, 10.3, 11.7, 13.18, 16.25, 20.38, 27.76, 36.46, 54.25, 87.58, 220.49]
thr_wsm_list = [4.67, 4.89, 5.15, 5.4, 5.69, 6.05, 6.38, 6.79, 7.18, 7.77, 8.35, 9.01, 10.0, 11.55, 13.78, 16.73, 21.7, 29.91, 44.23, 116.12]



# In[68]:


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
    


# In[69]:


#function which calculates the chi2, also checks if t track is in window if flag is set True

def ChiSquareVeloT(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1 , zmag, rangex, rangey, sigma2x, sigma2y, sigma2tx, sigma2ty, check_if_inwindow=False, smearing=False): #to ask: velo p?
    
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
    #yT_ex=yVzmag+(zt0-zmag)*tyT
    #from velo to t track with velo slope
    yT_ex=yv1+(zt0-zv1)*tyT

    #complete chi2 (slopes and coordinates)
    chi2 = (txPre-txT)*(txPre-txT)/sigma2tx+(tyPre-tyT)*(tyPre-tyT)/sigma2ty + (yT_ex-yt0)*(yT_ex-yt0)/sigma2y + (xT_ex-xt0)*(xT_ex-xt0)/sigma2x 

    #check if in window
    if (check_if_inwindow == True):
        

        #if in window: return chi2, else placeholder -1
        if(abs(xT_ex-xt0)<rangex and abs(yT_ex-yt0)<rangey):
            return chi2
        else:
            return (-1)
            
    #no need to check 
    else:
        return chi2
    
 


# In[70]:


smearing = False

for percentage in [94,95,96,97,98,99]:
    print(percentage)
    nLongTracks = 0
    real=0.;fake=0.
    realbelow=0.;fakebelow=0.

    timestart = time.time()
    nentries = ch1.GetEntries()

    #percentage = 97
    index_perc = perc_list.index(percentage)

    #change chi threshold, sigmas and windows depending on smearing on/off
    if smearing == False:
        threshold_chi2 = thr_nosm_list[index_perc]
        #squared
        sigma2x=pow(sigmax_nosm,2)
        sigma2y=pow(sigmay_nosm,2)
        sigma2tx=pow(sigmatx_nosm,2)
        sigma2ty=pow(sigmaty_nosm,2)
        #2 sigma acceptance
        rangex=2*sigmax_nosm
        rangey=2*sigmay_nosm
    else:
        threshold_chi2 = thr_wsm_list[index_perc]
        sigma2x=pow(sigmax_wsm,2)
        sigma2y=pow(sigmay_wsm,2)
        sigma2tx=pow(sigmatx_wsm,2)
        sigma2ty=pow(sigmaty_wsm,2)
        #2 sigma acceptance
        rangex=2*sigmax_wsm
        rangey=2*sigmay_wsm

    '''Nbins=100
    #plots for strategy min chi2 (need threshold defined)
    chi2R=TH1F("chi2R","",Nbins,0,threshold_chi2)
    chi2F=TH1F("chi2F","",Nbins,0,threshold_chi2)
    #plots for strategy cut below threshold
    percbelow_th=TH1F("chi2R_th","",Nbins,0,1)
    percF_th=TH1F("chi2F_th","",Nbins,0,1)'''

    print('percentage: ', percentage)
    print('smearing: ', smearing)
    print('threshold chi2: : ', threshold_chi2)



    #loop over long tracks
    for event in ch1:

        ChisquareTry=[]
        Chi2Try_entryNumbers=[]

        #selecting some long tracks for event in ch1 and cut for !=0 else division/0
        if (event.HitVeloZpos[1]-event.HitVeloZpos[0] != 0 and event.HitZpos[1]-event.HitZpos[0] != 0
           and  0.<event.HitVeloZpos[0]<event.HitVeloZpos[1]<800. and event.HitZpos[0]>7000. and event.HitZpos[1]>7000 and event.p > 5000):
            nLongTracks = nLongTracks+1

            xv0 = event.HitVeloXpos[0]; yv0 = event.HitVeloYpos[0]; zv0 = event.HitVeloZpos[0];
            xv1 = event.HitVeloXpos[1]; yv1 = event.HitVeloYpos[1]; zv1 = event.HitVeloZpos[1];
            xt0 = event.HitXpos[0]; yt0 = event.HitYpos[0]; zt0 = event.HitZpos[0];
            xt1 = event.HitXpos[1]; yt1 = event.HitYpos[1]; zt1 = event.HitZpos[1];
            vp = event.p

            chi2real = ChiSquareVeloT(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1, zmag, rangex, rangey, sigma2x, sigma2y, sigma2tx, sigma2ty)
            entryreal = event.GetReadEntry()
            eventreal = event.eventNumber


            #loop over all events, cut for momentum and dz !=0 else division/0 (later cut in window)
            for event2 in ch2:
                #select T tracks in same event
                if (event2.eventNumber == eventreal and event2.HitZpos[0] > 7000 and event2.HitZpos[1] > 7000  and event2.p > 5000 and event2.HitZpos[1]-event2.HitZpos[0] != 0):

                    #overwrite the t track coordinates (velo coordinates and momentum (?) remain the same)
                    xt0 = event2.HitXpos[0]; yt0 = event2.HitYpos[0]; zt0 = event2.HitZpos[0];
                    xt1 = event2.HitXpos[1]; yt1 = event2.HitYpos[1]; zt1 = event2.HitZpos[1];

                    chi2try = ChiSquareVeloT(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1 , zmag, rangex, rangey, sigma2x, sigma2y, sigma2tx, sigma2ty, check_if_inwindow=True, smearing=smearing)
                    if chi2try != -1: #cut events in search window (see def of chi2 function)
                        ChisquareTry.append(chi2try)
                        Chi2Try_entryNumbers.append(event2.GetReadEntry())


            #STRATEGY MINSQUARE: find min chi2

            if len(ChisquareTry) != 0:
                ChisquareSM=min(ChisquareTry)
                indexMin= ChisquareTry.index(ChisquareSM)

                #check if it real or fake match (compare number of entry)
                if (Chi2Try_entryNumbers[indexMin]==entryreal):
                    real=real+1.
                    #chi2R.Fill(ChisquareSM)
                else:
                    fake=fake+1.
                    #chi2F.Fill(ChisquareSM)


            #STRATEGY CUT BELOW THRESHOLD

            chi_below_thr = 0
            fakebelow = 0
            count = 0
            for chitry in ChisquareTry:
                count +=1
                indexchitry = ChisquareTry.index(chitry)
                if (chitry <= threshold_chi2):
                    #find how many tracks you find below threshold
                    chi_below_thr=chi_below_thr+1.
                    #find how many times the real one is amongst those below threshold
                    if (Chi2Try_entryNumbers[indexchitry]==entryreal):
                        realbelow=realbelow+1
                    #find how many wrong ones are below the threshold:
                    else:
                        fakebelow=fakebelow+1

            #plot some kind of (un)efficiency 
            #if chi_below_thr != 0:
                #percF_th.Fill(fakebelow/chi_below_thr)
            #plot % of t tracks below threshold 
            #if len(ChisquareTry) != 0:
                #percbelow_th.Fill(chi_below_thr/len(ChisquareTry))

        #cut at n good events (computational time issue)
            if (nLongTracks == 2000 or entryreal == nentries):
                break

    deltat = time.time() - timestart
    print(deltat, nLongTracks)

    fraction_min=real/(real+fake) #real+fake is == nevents since only one outcome
    print ('MINIMUM CHI: real, fake, fraction')
    print real,fake,fraction_min

    fraction_thr = realbelow / nLongTracks #here real+fake != nevents since more than one wrong outcome (and not always right outcome)
    print ('CUT BELOW THRESHOLD: realbelow, n tot events, fraction')
    print realbelow, nLongTracks, fraction_thr


# In[ ]:





# In[285]:



#draw plots for minimum chi2 strategy

c4=TCanvas("c4","",1300,900)
chi2R.Draw()
chi2R.GetXaxis().SetTitle("#chi^{2}_{smallest}")
chi2R.GetYaxis().SetTitle("Events")
chi2R.GetYaxis().SetTitleOffset(0.8)
chi2R.SetLineColor(kBlue)
chi2R.SetMinimum(0.)
chi2F.Draw("SAME")
chi2F.SetLineColor(kRed)

legendi=TLegend(0.6,0.75,0.8,0.9)
legendi.AddEntry("chi2R","Real match","f")
legendi.AddEntry("chi2F","Fake match","f")
legendi.Draw()

la=TLatex()
la.SetTextSize(0.04)
la.DrawLatex(60,80,"N_{real}/N_{tot}="+"{0:.2f}%".format(fraction_min*100))
c4.SaveAs("distinguish_range_chi2fr.png")




#draw plots for cut below threshold strategy

c5=TCanvas("c5","",1300,900)
percF_th.Draw()
percF_th.GetXaxis().SetTitle("#chi^{2}_{smallest}")
percF_th.GetYaxis().SetTitle("Events")
percF_th.GetYaxis().SetTitleOffset(0.8)
percF_th.SetLineColor(kBlue)
percF_th.SetMinimum(0.)
c5.SaveAs("distinguish_range_fakebelow.png")


c6=TCanvas("c6","",1300,900)
percbelow_th.Draw()
percbelow_th.GetXaxis().SetTitle("#chi^{2}_{smallest}")
percbelow_th.GetYaxis().SetTitle("Events")
percbelow_th.GetYaxis().SetTitleOffset(0.8)
percbelow_th.SetLineColor(kRed)
percbelow_th.SetMinimum(0.)
c6.SaveAs("distinguish_range_totbelow.png")


# In[ ]:




