from ROOT import *
import numpy as np
import math

zS=10000.
z0=0.

ch1=TChain("MCParticleNTuple/Tracks")
ch1.Add("~/MightyIT/MCtracks_oldfile_200ev.root")
ch2=TChain("MCParticleNTuple/Tracks")
ch2.Add("~/MightyIT/MCtracks_oldfile_200ev.root")

gROOT.ProcessLine(".x ~/lhcbStyle.C")
gStyle.SetPaintTextFormat("1.3f")

m_zMagParams=[5287.6, -7.98878, 317.683, 0.0119379, -1418.42]

Par=[4943170,6314610]
ParF=[-0.0001666775,-1304.68]
rangex=160.
rangey=50.
real=0.;fake=0.

Nbins=80;Nrange=80.
chi2R=TH1F("chi2R","",Nbins,0,Nrange)
chi2F=TH1F("chi2F","",Nbins,0,Nrange)

nEvent = 0

for event in ch1:
    i=1;j=0
    Chisquare=[]
    #selecting some long tracks for event in ch1 and cut for !=0 else division/0
    if (event.HitVeloZpos[i]-event.HitVeloZpos[j] != 0 and event.HitZpos[1]-event.HitZpos[0] != 0 
       and  0.<event.HitVeloZpos[j]<event.HitVeloZpos[i]<800. and 7825.<event.HitZpos[0]<7875. and 7900.<event.HitZpos[1]<7950.):
        nEvent = nEvent+1
        
        #calculate the slopes
        tyV=(event.HitVeloYpos[i]-event.HitVeloYpos[j])/(event.HitVeloZpos[i]-event.HitVeloZpos[j])
        txV=(event.HitVeloXpos[i]-event.HitVeloXpos[j])/(event.HitVeloZpos[i]-event.HitVeloZpos[j])
        tyT=(event.HitYpos[1]-event.HitYpos[0])/(event.HitZpos[1]-event.HitZpos[0])
        txT=(event.HitXpos[1]-event.HitXpos[0])/(event.HitZpos[1]-event.HitZpos[0])
        
        #zmag calculations
        xTzS=event.HitXpos[1]+(zS-event.HitZpos[1])*txT
        dSlopex=txV-txT
        zmag=m_zMagParams[0]+m_zMagParams[1]*abs(dSlopex)+m_zMagParams[2]*dSlopex*dSlopex+m_zMagParams[3]*abs(xTzS)+m_zMagParams[4]*txV*txV
        xVzmag=event.HitVeloXpos[i]+(zmag-event.HitVeloZpos[i])*txV
        yVzmag=event.HitVeloYpos[i]+(zmag-event.HitVeloZpos[i])*tyV

        #chi2 calculations
        txPre=(event.HitXpos[0]-xVzmag)/(event.HitZpos[0]-zmag)
        tyPre=(event.HitYpos[0]-yVzmag)/(event.HitZpos[0]-zmag)
        sigma2tx=Par[0]/pow(event.p,3.)
        sigma2ty=Par[1]/pow(event.p,3.)
        Chisquare.append((txPre-txT)*(txPre-txT)/sigma2tx+(tyPre-tyT)*(tyPre-tyT)/sigma2ty)
        
        #loop over all events in the window, !!!cut for !=0 else division/0
        for k in range(0,ch2.GetEntries()):
            #start from the same event so i can check later if the min chi2 is the same
            ch2.GetEntry(event.GetReadEntry()+k)
            #cut in window
            if 0<abs(ch2.GetLeaf("HitXpos").GetValue(0)-event.HitXpos[0])<rangex and 0<abs(ch2.GetLeaf("HitYpos").GetValue(0)-event.HitYpos[0])<rangey and ch2.GetLeaf("p").GetValue()>5000.:
                #cut for !=0 else division/0
                if (ch2.GetLeaf("HitZpos").GetValue(1)-ch2.GetLeaf("HitZpos").GetValue(0) != 0 ):
                    #calculate slopes
                    tyT2=(ch2.GetLeaf("HitYpos").GetValue(1)-ch2.GetLeaf("HitYpos").GetValue(0))/(ch2.GetLeaf("HitZpos").GetValue(1)-ch2.GetLeaf("HitZpos").GetValue(0))
                    txT2=(ch2.GetLeaf("HitXpos").GetValue(1)-ch2.GetLeaf("HitXpos").GetValue(0))/(ch2.GetLeaf("HitZpos").GetValue(1)-ch2.GetLeaf("HitZpos").GetValue(0))
                    #zmag calculations
                    xTzS=ch2.GetLeaf("HitXpos").GetValue(1)+(zS-ch2.GetLeaf("HitZpos").GetValue(1))*txT2
                    dSlopex=txV-txT2
                    zmag=m_zMagParams[0]+m_zMagParams[1]*abs(dSlopex)+m_zMagParams[2]*dSlopex*dSlopex+m_zMagParams[3]*abs(xTzS)+m_zMagParams[4]*txV*txV
                    xVzmag=event.HitVeloXpos[i]+(zmag-event.HitVeloZpos[i])*txV
                    yVzmag=event.HitVeloYpos[i]+(zmag-event.HitVeloZpos[i])*tyV
                    #chisquare calculations
                    txPre=(ch2.GetLeaf("HitXpos").GetValue(0)-xVzmag)/(ch2.GetLeaf("HitZpos").GetValue(0)-zmag)
                    tyPre=(ch2.GetLeaf("HitYpos").GetValue(0)-yVzmag)/(ch2.GetLeaf("HitZpos").GetValue(0)-zmag)
                    sigma2tx=Par[0]/pow(event.p,3.)
                    sigma2ty=Par[1]/pow(event.p,3.)
                    Chisquare.append((txPre-txT2)*(txPre-txT2)/sigma2tx+(tyPre-tyT2)*(tyPre-tyT2)/sigma2ty)
        #find min chi2
        ChisquareSM=min(Chisquare)
        print len(Chisquare),ChisquareSM, Chisquare.index(ChisquareSM)
        #check if it real or fake match
        if Chisquare.index(ChisquareSM)==0 and ChisquareSM!=0.:
            real=real+1.
            chi2R.Fill(ChisquareSM)
        elif Chisquare.index(ChisquareSM)!=0 and ChisquareSM!=0.:
            fake=fake+1.
            chi2F.Fill(ChisquareSM)
        if event.GetReadEntry()>=20000:
            break
        #cut at 50 good events, just a computational time issue
        if nEvent > 1:
            break

#draw plots
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
fraction=real/(real+fake)
print real,fake,fraction
la.DrawLatex(60,80,"N_{real}/N_{tot}="+"{0:.2f}%".format(fraction*100))
c4.SaveAs("distinguish_range.png")
