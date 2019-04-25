#!/usr/bin/env python
# coding: utf-8

# In[40]:


from ROOT import *
import numpy as np
import math
import time

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


# In[41]:


ch1=TChain("MCParticleNTuple/Tracks")
ch1.Add("~/MightyIT/MCtracks_oldfile_200ev.root")
ch2=TChain("MCParticleNTuple/Tracks")
ch2.Add("~/MightyIT/MCtracks_oldfile_200ev.root")


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



# In[43]:


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
    


# In[44]:


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

    #ŧesting: complete chi2 (slopes and coordinates)
    chi2 = (txPre-txT)*(txPre-txT)/sigma2tx+(tyPre-tyT)*(tyPre-tyT)/sigma2ty + (yT_ex-yt0)*(yT_ex-yt0)/sigma2y + (xT_ex-xt0)*(xT_ex-xt0)/sigma2x 

    #check if in window
    if (check_if_inwindow == True):
        

        #if in window: return chi2, else placeholder -1
        if(abs(xT_ex-xt0)<rangex and abs(yT_ex-yt0)<rangey):
            return chi2, (xT_ex-xt0), (yT_ex-yt0), (txPre-txT), (tyPre-tyT)
        else:
            return (-1, -1, -1, -1, -1)
            
    #no need to check 
    else:
        return chi2, (xT_ex-xt0), (yT_ex-yt0), (txPre-txT), (tyPre-tyT) 


# In[45]:


smearing = False
percentage = 96

nLongTracks = 0
real=0.;fake=0.
realbelow=0.;fakebelow=0.

nentries = ch1.GetEntries()

index_perc = perc_list.index(percentage)

threshold_chi2 = 0
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


# In[46]:



Nbins=60

#plots for strategy min chi2
chi2R_true=TH1F("chi2R_true","",Nbins,0,threshold_chi2)
chi2R=TH1F("chi2R","",Nbins,0,threshold_chi2)
chi2F=TH1F("chi2F","",Nbins,0,threshold_chi2)

pR_true=TH1F("pF_true","",Nbins,5000,50000)
pR=TH1F("pR","",Nbins,5000,50000)
pF=TH1F("pF","",Nbins,5000,50000)

xyzt0R_true=TH2F("xyzt0R_true","",Nbins,-3000,+3000,Nbins,-3000,+3000)
xyzt0R=TH2F("xyzt0R","",Nbins,-3000,+3000,Nbins,-3000,+3000)
xyzt0F=TH2F("xyzt0F","",Nbins,-3000,+3000,Nbins,-3000,+3000)

dxR_true=TH1F("dxR_true","",Nbins,-rangex,rangex)
dxR=TH1F("dxR","",Nbins,-rangex,rangex)
dxF=TH1F("dxF","",Nbins,-rangex,rangex)

dyR_true=TH1F("dyR_true","",Nbins,-rangey,rangey)
dyR=TH1F("dyR","",Nbins,-rangey,rangey)
dyF=TH1F("dyR","",Nbins,-rangey,rangey)

sigmatx = np.sqrt(sigma2tx)
dtxR_true=TH1F("dtxR_true","",Nbins,-2*sigmatx,2*sigmatx)
dtxR=TH1F("dtxR","",Nbins,-2*sigmatx,2*sigmatx)
dtxF=TH1F("dtxF","",Nbins,-2*sigmatx,2*sigmatx)

sigmaty = np.sqrt(sigma2ty)
dtyR_true=TH1F("dtyR_true","",Nbins,-2*sigmaty,2*sigmaty)
dtyR=TH1F("dtyR","",Nbins,-2*sigmaty,2*sigmaty)
dtyF=TH1F("dtyF","",Nbins,-2*sigmaty,2*sigmaty)


#plots for strategy cut below threshold
percbelow_th=TH1F("chi2R_th","",Nbins,0,1)
percF_th=TH1F("chi2F_th","",Nbins,0,1)

print('percentage: ', percentage)
print('smearing: ', smearing)
print('threshold chi2: : ', threshold_chi2)


# In[47]:


timestart = time.time()
#loop over long tracks
for event in ch1:
    
    ChisquareTry=[]
    Chi2Try_entryNumbers=[]
    pTry=[]
    xzt0Try=[]
    yzt0Try=[]
    dxTry=[]
    dyTry=[]
    dtxTry=[]
    dtyTry=[]

    #selecting some long tracks for event in ch1 and cut for !=0 else division/0
    if (event.HitVeloZpos[1]-event.HitVeloZpos[0] != 0 and event.HitZpos[1]-event.HitZpos[0] != 0
       and  0.<event.HitVeloZpos[0]<event.HitVeloZpos[1]<800. and event.HitZpos[0]>7000. and event.HitZpos[1]>7000 and event.p > 5000):
        nLongTracks = nLongTracks+1
        print(nLongTracks)

        xv0 = event.HitVeloXpos[0]; yv0 = event.HitVeloYpos[0]; zv0 = event.HitVeloZpos[0];
        xv1 = event.HitVeloXpos[1]; yv1 = event.HitVeloYpos[1]; zv1 = event.HitVeloZpos[1];
        xt0 = event.HitXpos[0]; yt0 = event.HitYpos[0]; zt0 = event.HitZpos[0];
        xt1 = event.HitXpos[1]; yt1 = event.HitYpos[1]; zt1 = event.HitZpos[1];
        vp = event.p

        chi2real_l = ChiSquareVeloT(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1 , zmag, rangex, rangey, sigma2x, sigma2y, sigma2tx, sigma2ty, check_if_inwindow=False,smearing=smearing)
        chi2real = chi2real_l[0]
        entryreal = event.GetReadEntry()
        eventreal = event.eventNumber
        p = event.p
        xzt0 = xt0 
        yzt0 = yt0
        dx = chi2real_l[1]
        dy = chi2real_l[2]
        dtx = chi2real_l[3]
        dty = chi2real_l[4]


        #loop over all events, cut for momentum and dz !=0 else division/0 (later cut in window)
        for event2 in ch2:
            #select T tracks in same event
            if (event2.eventNumber == eventreal and event2.HitZpos[0] > 7000 and event2.HitZpos[1] > 7000  and event2.p > 5000 and event2.HitZpos[1]-event2.HitZpos[0] != 0):

                #overwrite the t track coordinates (velo coordinates and momentum (?) remain the same)
                xt0 = event2.HitXpos[0]; yt0 = event2.HitYpos[0]; zt0 = event2.HitZpos[0];
                xt1 = event2.HitXpos[1]; yt1 = event2.HitYpos[1]; zt1 = event2.HitZpos[1];

                chi2try_l = ChiSquareVeloT(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1 , zmag, rangex, rangey, sigma2x, sigma2y, sigma2tx, sigma2ty, check_if_inwindow=True, smearing=smearing)
                if chi2try_l[0] != -1: #cut events in search window (see def of chi2 function)
                    ChisquareTry.append(chi2try_l[0])
                    Chi2Try_entryNumbers.append(event2.GetReadEntry())
                    pTry.append(event2.p)
                    xzt0Try.append(xt0)
                    yzt0Try.append(yt0)
                    dxTry.append(chi2try_l[1])
                    dyTry.append(chi2try_l[2])
                    dtxTry.append(chi2try_l[3])
                    dtyTry.append(chi2try_l[4])



        if len(ChisquareTry) != 0:
            
            #STRATEGY MINSQUARE: find min chi2
            
            ChisquareSM=min(ChisquareTry)
            indexMin= ChisquareTry.index(ChisquareSM)

            #check if it real or fake match (compare number of entry)
            if (Chi2Try_entryNumbers[indexMin]==entryreal):
                real=real+1.
                #fill all the real plots
                chi2R_true.Fill(chi2real)
                pR_true.Fill(p)
                xyzt0R_true.Fill(xzt0,yzt0)
                dxR_true.Fill(dx)
                dyR_true.Fill(dy)
                dtxR_true.Fill(dtx)
                dtyR_true.Fill(dty)
                
            else:
                fake=fake+1.
                #fill all the real plots
                chi2R.Fill(chi2real)
                pR.Fill(p)
                xyzt0R.Fill(xzt0,yzt0)
                dxR.Fill(dx)
                dyR.Fill(dy)
                dtxR.Fill(dtx)
                dtyR.Fill(dty)
                #fill all the fake plots
                chi2F.Fill(ChisquareSM)
                pF.Fill(pTry[indexMin])
                xyzt0F.Fill(xzt0Try[indexMin],yzt0Try[indexMin])
                dxF.Fill(dxTry[indexMin])
                dyF.Fill(dyTry[indexMin])
                dtxF.Fill(dtxTry[indexMin])
                dtyF.Fill(dtyTry[indexMin])


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
            if chi_below_thr != 0:
                percF_th.Fill(fakebelow/chi_below_thr)
            #plot % of t tracks below threshold 
            if len(ChisquareTry) != 0:
                percbelow_th.Fill(chi_below_thr/len(ChisquareTry))


        #cut at n good events (computational time issue)
        if (nLongTracks == 2500 or entryreal == nentries):
            break

deltat = time.time() - timestart
print(deltat)

fraction_min=real/(real+fake) #real+fake is == nevents since only one outcome
print ('MINIMUM CHI: real, fake, fraction')
print real,fake,fraction_min

fraction_thr = realbelow / nLongTracks #here real+fake != nevents since more than one wrong outcome (and not always right outcome)
print ('CUT BELOW THRESHOLD: realbelow, n tot events, fraction')
print realbelow, nLongTracks, fraction_thr


# In[48]:



#draw plots for minimum chi2 strategy

c4=TCanvas("c4","",2000,2000)

c4.Divide(3,3)

c4.cd(1)
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

c4.cd(2)
pR.Draw()
pR.GetXaxis().SetTitle("p")
pR.GetYaxis().SetTitle("Events")
pR.SetLineColor(kBlue)
pF.Draw("SAME")
pF.SetLineColor(kRed)

c4.cd(3)
xyzt0R.Draw()
xyzt0R.GetXaxis().SetTitle("x")
xyzt0R.GetYaxis().SetTitle("y")
xyzt0R.SetMarkerColor(kBlue)
xyzt0F.Draw("SAME")
xyzt0F.SetMarkerColor(kRed)


c4.cd(4)
dxR.Draw()
dxR.GetXaxis().SetTitle("dx")
dxR.GetYaxis().SetTitle("Events")
dxR.SetLineColor(kBlue)
dxF.Draw("SAME")
dxF.SetLineColor(kRed)


c4.cd(5)

dyR.Draw()
dyR.GetXaxis().SetTitle("dy")
dyR.GetYaxis().SetTitle("Events")
dyR.SetLineColor(kBlue)
dyF.Draw("SAME")
dyF.SetLineColor(kRed)

c4.cd(6)
dtxR.Draw()
dtxR.GetXaxis().SetTitle("dtx")
dtxR.GetYaxis().SetTitle("Events")
dtxR.SetLineColor(kBlue)
dtxF.Draw("SAME")
dtxF.SetLineColor(kRed)

c4.cd(7)
dtyR.Draw()
dtyR.GetXaxis().SetTitle("dty")
dtyR.GetYaxis().SetTitle("Events")
dtyR.SetLineColor(kBlue)
dtyF.Draw("SAME")
dtyF.SetLineColor(kRed)

c4.SaveAs("distinguish_range_fr_fakemin.png")


# In[49]:



#draw plots for minimum chi2 strategy

c3=TCanvas("c3","",2000,2000)

c3.Divide(3,3)

c3.cd(1)
chi2R_true.Draw()
chi2R_true.GetXaxis().SetTitle("#chi^{2}_{smallest}")
chi2R_true.GetYaxis().SetTitle("Events")
chi2R_true.GetYaxis().SetTitleOffset(0.8)
chi2R_true.SetLineColor(kBlue)
chi2R_true.SetMinimum(0.)

c3.cd(2)
pR_true.Draw()
pR_true.GetXaxis().SetTitle("p")
pR_true.GetYaxis().SetTitle("Events")
pR_true.SetLineColor(kBlue)

c3.cd(3)
xyzt0R_true.Draw()
xyzt0R_true.GetXaxis().SetTitle("x")
xyzt0R_true.GetYaxis().SetTitle("y")
xyzt0R.SetMarkerColor(kBlue)


c3.cd(4)
dxR_true.Draw()
dxR_true.GetXaxis().SetTitle("dx")
dxR_true.GetYaxis().SetTitle("Events")
dxR_true.SetLineColor(kBlue)

c3.cd(5)
dyR_true.Draw()
dyR_true.GetXaxis().SetTitle("dy")
dyR_true.GetYaxis().SetTitle("Events")
dyR_true.SetLineColor(kBlue)

c3.cd(6)
dtxR_true.Draw()
dtxR_true.GetXaxis().SetTitle("dtx")
dtxR_true.GetYaxis().SetTitle("Events")
dtxR_true.SetLineColor(kBlue)

c3.cd(7)
dtyR_true.Draw()
dtyR_true.GetXaxis().SetTitle("dty")
dtyR_true.GetYaxis().SetTitle("Events")
dtyR_true.SetLineColor(kBlue)

c3.SaveAs("distinguish_range_fr_truemin.png")


# In[50]:



#draw plots for cut below threshold strategy

c5=TCanvas("c5","",1300,900)
c5.Divide(2,1)
c5.cd(1)
percF_th.Draw()
percF_th.GetXaxis().SetTitle("fakebelow/chibelowthr")
percF_th.GetYaxis().SetTitle("Events")
percF_th.GetYaxis().SetTitleOffset(0.8)
percF_th.SetLineColor(kBlue)
percF_th.SetMinimum(0.)

c5.cd(2)
percbelow_th.Draw()
percbelow_th.GetXaxis().SetTitle("totbelow/chi2try")
percbelow_th.GetYaxis().SetTitle("Events")
percbelow_th.GetYaxis().SetTitleOffset(0.8)
percbelow_th.SetLineColor(kRed)
percbelow_th.SetMinimum(0.)
c5.SaveAs("distinguish_range_belowthr.png")


# In[ ]:




