#!/usr/bin/env python
# coding: utf-8

# In[195]:


from ROOT import *
import numpy as np
import math


# In[216]:


zS=10000.
z0=0.

ch1=TChain("MCParticleNTuple/Tracks")
ch1.Add("~/MightyIT/MCtracks_oldfile_200ev.root")
ch2=TChain("MCParticleNTuple/Tracks")
ch2.Add("~/MightyIT/MCtracks_oldfile_200ev.root")

gROOT.ProcessLine(".x ~/lhcbStyle.C")
gStyle.SetPaintTextFormat("1.3f")

#optimized zmag
m_zMagParams=[5092.5, -7.98878, 317.683, 0.0119379, -1418.42]

Par=[4943170,6314610]

#tested from scatterx, y, tx, ty
#2 sigma acceptance
rangex=2*22.63
rangey=2*9.982
#sigma squared for chi2 calculation
sigma2x=pow(22.63,2)
sigma2y=pow(9.982,2)
sigma2tx=pow(0.009824,2)
sigma2ty=pow(0.002917,2)


Nbins=80;Nrange=80.
#plots for strategy min chi2
chi2R=TH1F("chi2R","",Nbins,0,Nrange)
chi2F=TH1F("chi2F","",Nbins,0,Nrange)
#plots for strategy cut below threshold
percbelow_th=TH1F("chi2R_th","",Nbins,0,1)
percF_th=TH1F("chi2F_th","",Nbins,0,1)

#threshold_chi2 = 30 #coming from sigmas program. 99%
threshold_chi2 = 3.4 #coming from sigmas program: 95%


# In[217]:


#function which calculates the chi2, also checks if t track is in window if flag is set True
def ChiSquareVeloT(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1 , zS, m_zMagParams, Par, rangex, rangey, vp, sigma2x, sigma2y, sigma2tx, sigma2ty, check_if_inwindow): #to ask: velo p?
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
    #sigma2tx=Par[0]/pow(vp,3.)
    #sigma2ty=Par[1]/pow(vp,3.)
    
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
            return chi2
        else:
            return (-1)
            
    #no need to check 
    else:
        return chi2
    
 


# In[218]:


print(ch1.GetEntries())


# In[223]:


nEvent = 0
real=0.;fake=0.
realbelow=0.;fakebelow=0.


#loop over long tracks
for event in ch1:

    ChisquareTry=[]
    Chi2Try_entryNumbers=[]
    
    #selecting some long tracks for event in ch1 and cut for !=0 else division/0
    if (event.HitVeloZpos[1]-event.HitVeloZpos[0] != 0 and event.HitZpos[1]-event.HitZpos[0] != 0
       and  0.<event.HitVeloZpos[0]<event.HitVeloZpos[1]<800. and 7825.<event.HitZpos[0]<7875. and 7900.<event.HitZpos[1]<7950. and event.p > 5000):
        nEvent = nEvent+1
        
        xv0 = event.HitVeloXpos[0]; yv0 = event.HitVeloYpos[0]; zv0 = event.HitVeloZpos[0];
        xv1 = event.HitVeloXpos[1]; yv1 = event.HitVeloYpos[1]; zv1 = event.HitVeloZpos[1];
        xt0 = event.HitXpos[0]; yt0 = event.HitYpos[0]; zt0 = event.HitZpos[0];
        xt1 = event.HitXpos[1]; yt1 = event.HitYpos[1]; zt1 = event.HitZpos[1];
        vp = event.p
       
        chi2real = ChiSquareVeloT(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1 , zS, m_zMagParams, Par, rangex, rangey, vp, sigma2x, sigma2y, sigma2tx, sigma2ty, check_if_inwindow=False) #to ask: velo p?
        entryreal = event.GetReadEntry()
        print('chi2real ', chi2real, 'entryreal ', entryreal)

        
        #loop over all events, cut for momentum and dz !=0 else division/0 (later cut in window)
        for event2 in ch2:
            #select T tracks 
            if (event2.HitZpos[0] > 7000 and event2.HitZpos[1] > 7000  and event2.p > 5000 and event2.HitZpos[1]-event2.HitZpos[0] != 0):

                #overwrite the t track coordinates (velo coordinates and momentum (?) remain the same)
                xt0 = event2.HitXpos[0]; yt0 = event2.HitYpos[0]; zt0 = event2.HitZpos[0];
                xt1 = event2.HitXpos[1]; yt1 = event2.HitYpos[1]; zt1 = event2.HitZpos[1];

                chi2try = ChiSquareVeloT(xv0, yv0, zv0, xv1, yv1, zv1, xt0, yt0, zt0, xt1, yt1, zt1 , zS, m_zMagParams, Par, rangex, rangey, vp, sigma2x, sigma2y, sigma2tx, sigma2ty, check_if_inwindow=True)
                if chi2try != -1: #cut events in search window (see def of chi2 function)
                    ChisquareTry.append(chi2try)
                    Chi2Try_entryNumbers.append(event2.GetReadEntry())
                    
                                         
        #STRATEGY MINSQUARE: find min chi2
        ChisquareSM=min(ChisquareTry)
        indexMin= ChisquareTry.index(ChisquareSM)
        
        print len(ChisquareTry), ChisquareSM, indexMin
        
        #check if it real or fake match (compare number of entry)
        if (Chi2Try_entryNumbers[indexMin]==entryreal):
            real=real+1.
            chi2R.Fill(ChisquareSM)
        else:
            fake=fake+1.
            chi2F.Fill(ChisquareSM)
        
       
        #STRATEGY CUT BELOW THRESHOLD
        
        #testing
        #threshold_chi2 = chi2real #should be fixed
        
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
        print('fake below thr: ', fakebelow)
        #plot some kind of (un)efficiency <- will tune threshold to minimise this
        if chi_below_thr != 0:
            percF_th.Fill(fakebelow/chi_below_thr)
        #plot % of t tracks below threshold (separate study? complementary)
        if len(ChisquareTry) != 0:
            percbelow_th.Fill(chi_below_thr/len(ChisquareTry))
            
        
        #cut at n good events (computational time issue)
        if nEvent == 3:
            break


# In[224]:


fraction_min=real/(real+fake) #real+fake is == nevents since only one outcome
print ('MINIMUM CHI: real, fake, fraction')
print real,fake,fraction_min



fraction_thr = realbelow / nEvent #here real+fake != nevents since more than one wrong outcome (and not always right outcome)
print ('CUT BELOW THRESHOLD: realbelow, n tot events, fraction')
print realbelow, nEvent, fraction_thr
print chi_below_thr, count


# In[60]:



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
c4.SaveAs("distinguish_range_tcut_500ev.png")




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




