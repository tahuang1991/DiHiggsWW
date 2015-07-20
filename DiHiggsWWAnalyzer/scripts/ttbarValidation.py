import random
import ROOT
import os

ROOT.gROOT.SetBatch(1)
#gStyle from TStyle
ROOT.gStyle.SetStatW(0.17)
ROOT.gStyle.SetStatH(0.15)

ROOT.gStyle.SetOptStat(111110)

ROOT.gStyle.SetTitleStyle(0)
ROOT.gStyle.SetTitleAlign(13) ## coord in top left
ROOT.gStyle.SetTitleX(0.)
ROOT.gStyle.SetTitleY(1.)
ROOT.gStyle.SetTitleW(1)
#ROOT.gStyle.SetTitleTextColor(4)
ROOT.gStyle.SetTitleXSize(0.05)
ROOT.gStyle.SetTitleYSize(0.05)
ROOT.gStyle.SetTitleH(0.058)
ROOT.gStyle.SetTitleBorderSize(0)

ROOT.gStyle.SetPadLeftMargin(0.126)
ROOT.gStyle.SetPadRightMargin(0.14)
ROOT.gStyle.SetPadTopMargin(0.06)
ROOT.gStyle.SetPadBottomMargin(0.13)



#___________________________________________
def draw1D(file,dir,todraw,x_bins,x_title,cut,pic_name):
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    f = ROOT.TFile(file)
    t = f.Get(dir)
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    
    b1 = ROOT.TH1F("b1","b1",xBins,xminBin,xmaxBin)
    b1.SetTitle("h2#rightarrow hh#rightarrow BBWW, B3"+" "*12 + "CMS Simulation Preliminary")
    b1.GetYaxis().SetTitle("Events")
    b1.GetXaxis().SetTitle("%s"%x_title)
    #b1.SetStats(0)

    b1.Sumw2() 
    t.Draw(todraw+">>b1",cut)
   # b1.Draw("colz")
#    b2 = ROOT.TF2("b2","x^2+y^2",xminBin,xmaxBin,yminBin,ymaxBin)
 #   contour = [1,2,3,4]
    #list = [None]*3
    #list.append(1)
    #list.append(2)
    #list.append(3)
    #list = list[-3:]    
    #b2.SetContour(4, contour)
   # b2.SetContourLevel(0,1)
    #b2.SetContourLevel(1,2)
    #b2.SetContourLevel(2,3)
    #b2.SetContourLevel(3,4)
    #b2.Draw("CONT3 same")
    b1.Fit("gaus")
    legend = ROOT.TLegend(0.15,0.46,0.45,0.64)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetFillStyle(0)
#    legend.SetHeader("PU140,simTrack Pt(%s"%pt_min+",%s)"%pt_max)
#legend.AddEntry(e1,"","l")
#    legend.Draw("same")
 #   line1 = "PU140,simTrack Pt(%s"%pt_min+",%s)"%pt_max
  #  line2 = "98 dphicut %f"%dphi_cut
   # tex = ROOT.TLatex(0.15,0.45,line1)
    #tex.SetTextSize(0.05)
    #tex.SetNDC()
    #tex.Draw("same")
    #tex2 = ROOT.TLatex(0.15,0.35,line2)
    #tex2.SetTextSize(0.05)
    #tex2.SetNDC()
    #tex2.Draw("same")
	
    c1.SaveAs("Dihiggs_%s"%pic_name+"_B3.pdf")
    c1.SaveAs("Dihiggs_%s"%pic_name+"_B3.png")


#____________________________________________________________________
def draw2D(file,dir,num,xaxis,yaxis,x_bins,y_bins):
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    f = ROOT.TFile(file)
    t = f.Get(dir)
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    
    b1 = ROOT.TH2F("b1","b1",xBins,xminBin,xmaxBin,yBins,yminBin,ymaxBin)
    b1.GetYaxis().SetTitle("%s"%yaxis)
    b1.GetXaxis().SetTitle("%s"%xaxis)
    b1.SetTitle("h1#rightarrow BB or WW, B3"+" "*12 + "CMS Simulation Preliminary")
   # b1.SetStats(1)
    todraw = "(%s)"%yaxis+":"+"(%s)>>b1"%xaxis
    t.Draw(todraw,num,"colz")
#    b1.SetMaximum(150)
    b1.Draw("colz")
    ROOT.gPad.SetLogz() 
    legend = ROOT.TLegend(0.15,0.56,0.40,0.64)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetHeader("PU140")
#legend.AddEntry(e1,"","l")
# legend.Draw("same")
    Num = b1.Integral(1,xBins,1,yBins)
    print "number of event ", Num

    tex = ROOT.TLatex(0.15,0.30,"#splitline{p_{T}^{sim}>20GeV,#frac{(pt-trackpt)}{pt}<-0.5}{stubs in TF matcehd to simtracks, Entries %s}"%Num)
#    tex = ROOT.TLatex(0.15,0.30,"p_{T}^{sim}>20GeV, #frac{abs(pt-trackpt)}{pt}<0.5, Entries %s"%Num)
#    tex = ROOT.TLatex(0.20,0.30,"#frac{(pt-trackpt)}{pt}>0.5, Entries %s"%Num)
#    tex = ROOT.TLatex(.2,0.3,"all stubs in TF matched to simtrack ")
    tex.SetTextSize(0.05)
    tex.SetTextFont(42)
    tex.SetNDC()
    tex.Draw("same")
	
    c1.SaveAs("Dihiggs_%s"%xaxis+"_VS_%s.pdf"%yaxis)
    c1.SaveAs("Dihiggs_%s"%xaxis+"_VS_%s.png"%yaxis)



#___________________________________________
def draw1D_combined(file,dir,todraw1,todraw2,x_bins,x_title,cut1,cut2, pic_name):
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    f = ROOT.TFile(file)
    t = f.Get(dir)
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    
    b2 = ROOT.TH1F("b2","b2",xBins,xminBin,xmaxBin)
    b2.SetTitle("h2#rightarrow hh#rightarrow BBWW, B3"+" "*12 + "CMS Simulation Preliminary")
    b2.GetYaxis().SetTitle("Events")
    b2.GetXaxis().SetTitle("%s"%x_title)
    t.Draw(todraw2+">>b2",cut2)
    #b1.SetStats(0)
    #wmassout = ROOT.TFile("onshellwmassout.root","recreate")
    #wmassout.cd()
    b1 = ROOT.TH1F("onshellWmasspdf","b1",xBins,xminBin,xmaxBin)
    b1.SetTitle("h2#rightarrow hh#rightarrow BBWW, B3"+" "*12 + "CMS Simulation Preliminary")
    b1.GetYaxis().SetTitle("Events")
    b1.GetXaxis().SetTitle("%s"%x_title)
    #b1.SetStats(0)


    #ROOT.gStyle.SetOptFit(1)
    t.Draw(todraw1+">>onshellWmasspdf",cut1)
    mean = 80.1
    signma = 1.5
    myfun = ROOT.TF1("myfun","exp(x*[0]+[1])+[2]*exp(-0.5*((x-80.1)/2.00)**2)",xminBin,xmaxBin)
    b1.Add(b2)
    b1.Sumw2()
    integral = b1.Integral("width")
    
    print "b1 integral ", integral," max",b1.GetBinContent(b1.GetMaximumBin())
    b1.Scale(1/integral)
    b1.Fit("pol6")
    #st =ROOT.TPaveStats(b1.FindObject("stats"))
    #st.SetX1NDC(0.1)
    #st.SetX2NCD(0.4)
    #print " -1 bin: ",b1.FindBin(-1)," 1000, bin: ",b1.FindBin(1000)," GetnBin",b1.GetNbinsX()
    """
    i = 0
    while i<10:
	r = random.uniform(50,90)
	bin = b1.FindBin(r)
	print "r ",r ," bin ",bin," bincenter1 ",b1.GetBinCenter(bin)," center2 ",b1.GetBinCenter(bin+1)
	i = i+1
    legend = ROOT.TLegend(0.15,0.46,0.45,0.64)
    legend.SetFillColor(ROOT.kWhite)
    #b1.Draw()
    """
    #b2.Draw("same")
    b1.Draw()
    #wmassout.Write()
    #wmassout.Close()
    
    c1.SaveAs("Dihiggs_%s"%pic_name+"combined_B3.pdf")
    c1.SaveAs("Dihiggs_%s"%pic_name+"combined_B3.png")
    c1.SaveAs("Dihiggs_%s"%pic_name+"combined_B3.C")
#    c1.SaveAs("Dihiggs_%s"%pic_name+"combined_B3.ROOT")


#___________________________________________
def draw2D_combined(file,dir,todraw1,todraw2,x_bins,y_bins,x_title,y_title,cut1,cut2, pic_name):
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    f = ROOT.TFile(file)
    t = f.Get(dir)
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    
    b2 = ROOT.TH2F("b2","b2",xBins,xminBin,xmaxBin,yBins,yminBin,ymaxBin)
    b2.SetTitle("h2#rightarrow hh#rightarrow BBWW, B3"+" "*12 + "CMS Simulation Preliminary")
    b2.GetYaxis().SetTitle("Events")
    b2.GetXaxis().SetTitle("%s"%x_title)
    t.Draw(todraw2+">>b2",cut2)
    #b1.SetStats(0)
#    wmassout = ROOT.TFile("onshellwmassout.root","recreate")
#    wmassout.cd()
    b1 = ROOT.TH2F("b1","b1",xBins,xminBin,xmaxBin,yBins,yminBin,ymaxBin)
    b1.SetTitle("h2#rightarrow hh#rightarrow BBWW, B3"+" "*12 + "CMS Simulation Preliminary")
    b1.GetYaxis().SetTitle("%s"%y_title)
    b1.GetXaxis().SetTitle("%s"%x_title)
    #b1.SetStats(0)


    #ROOT.gStyle.SetOptFit(1)
    t.Draw(todraw1+">>b1",cut1)
    mean = 80.1
    signma = 1.5
    b1.Add(b2)
    #b1.Sumw2()
   # b1.Fit("myfun")
    #st =ROOT.TPaveStats(b1.FindObject("stats"))
    #st.SetX1NDC(0.1)
    #st.SetX2NCD(0.4)
    """
    i = 0
    while i<10:
	r = random.uniform(50,90)
	bin = b1.FindBin(r)
	print "r ",r ," bin ",bin," bincenter1 ",b1.GetBinCenter(bin)," center2 ",b1.GetBinCenter(bin+1)
	i = i+1
    legend = ROOT.TLegend(0.15,0.46,0.45,0.64)
    legend.SetFillColor(ROOT.kWhite)
    #b1.Draw()
    """
    #b2.Draw("same")
    b1.Draw("colz")
 #   wmassout.Write()
 #   wmassout.Close()
    
    #c1.SaveAs("Dihiggs_%s"%pic_name+"combined_B3.pdf")
    #c1.SaveAs("Dihiggs_%s"%pic_name+"combined_B3.png")
    #c1.SaveAs("Dihiggs_%s"%pic_name+"combined_B3.C")
#    c1.SaveAs("Dihiggs_%s"%pic_name+"combined_B3.ROOT")




#deltaeta, deltaphi distribution
#___________________________________________
def deltaR1(file,dir,x_bins,y_bins,cut,pic_name):
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    f = ROOT.TFile(file)
    t = f.Get(dir)
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    
    b1 = ROOT.TH2F("b1","b1",xBins,xminBin,xmaxBin,yBins,yminBin,ymaxBin)
    b1.SetTitle("h2#rightarrow h1h1#rightarrow BBWW, B3"+" "*12 + "CMS Simulation Preliminary")
    b1.GetYaxis().SetTitle("#Delta#phi(#mu,#nu)")
    b1.GetXaxis().SetTitle("#Delta#eta(#mu,#nu)")
    b1.SetStats(0)
    
    todraw = "TVector2::Phi_mpi_pi(mu1_phi-nu1_phi):(mu1_eta-nu1_eta)>>b1"
    t.Draw(todraw,cut)
    b1.Draw("colz")
#    b2 = ROOT.TF2("b2","x^2+y^2",xminBin,xmaxBin,yminBin,ymaxBin)
 #   contour = [1,2,3,4]
    #list = [None]*3
    #list.append(1)
    #list.append(2)
    #list.append(3)
    #list = list[-3:]    
    #b2.SetContour(4, contour)
   # b2.SetContourLevel(0,1)
    #b2.SetContourLevel(1,2)
    #b2.SetContourLevel(2,3)
    #b2.SetContourLevel(3,4)
    #b2.Draw("CONT3 same")
    legend = ROOT.TLegend(0.15,0.46,0.45,0.64)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetFillStyle(0)
#    legend.SetHeader("PU140,simTrack Pt(%s"%pt_min+",%s)"%pt_max)
#legend.AddEntry(e1,"","l")
#    legend.Draw("same")
 #   line1 = "PU140,simTrack Pt(%s"%pt_min+",%s)"%pt_max
  #  line2 = "98 dphicut %f"%dphi_cut
   # tex = ROOT.TLatex(0.15,0.45,line1)
    #tex.SetTextSize(0.05)
    #tex.SetNDC()
    #tex.Draw("same")
    #tex2 = ROOT.TLatex(0.15,0.35,line2)
    #tex2.SetTextSize(0.05)
    #tex2.SetNDC()
    #tex2.Draw("same")
	
    c1.SaveAs("Dihiggs_deltaR1_%s"%pic_name+"_B3.pdf")
    c1.SaveAs("Dihiggs_deltaR1_%s"%pic_name+"_B3.png")


#___________________________________________
def deltaR2(file,dir,x_bins,y_bins,cut,pic_name):
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    f = ROOT.TFile(file)
    t = f.Get(dir)
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    
    b1 = ROOT.TH2F("b1","b1",xBins,xminBin,xmaxBin,yBins,yminBin,ymaxBin)
    b1.SetTitle("h2#rightarrow h1h1#rightarrow BBWW, B3"+" "*12 + "CMS Simulation Preliminary")
    b1.GetYaxis().SetTitle("#Delta#phi(#mu,#nu)")
    b1.GetXaxis().SetTitle("#Delta#eta(#mu,#nu)")
    b1.SetStats(0)
    
    todraw = "TVector2::Phi_mpi_pi(mu2_phi-nu2_phi):(mu2_eta-nu2_eta)>>b1"
    t.Draw(todraw,cut)
    b1.Draw("colz")
#    b2 = ROOT.TF2("b2","x^2+y^2",xminBin,xmaxBin,yminBin,ymaxBin)
 #   contour = [1,2,3,4]
    #list = [None]*3
    #list.append(1)
    #list.append(2)
    #list.append(3)
    #list = list[-3:]    
    #b2.SetContour(4, contour)
   # b2.SetContourLevel(0,1)
    #b2.SetContourLevel(1,2)
    #b2.SetContourLevel(2,3)
    #b2.SetContourLevel(3,4)
    #b2.Draw("CONT3 same")
    legend = ROOT.TLegend(0.15,0.46,0.45,0.64)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetFillStyle(0)
#    legend.SetHeader("PU140,simTrack Pt(%s"%pt_min+",%s)"%pt_max)
#legend.AddEntry(e1,"","l")
#    legend.Draw("same")
 #   line1 = "PU140,simTrack Pt(%s"%pt_min+",%s)"%pt_max
  #  line2 = "98 dphicut %f"%dphi_cut
   # tex = ROOT.TLatex(0.15,0.45,line1)
    #tex.SetTextSize(0.05)
    #tex.SetNDC()
    #tex.Draw("same")
    #tex2 = ROOT.TLatex(0.15,0.35,line2)
    #tex2.SetTextSize(0.05)
    #tex2.SetNDC()
    #tex2.Draw("same")
	
    c1.SaveAs("Dihiggs_deltaR2_%s"%pic_name+"_B3.pdf")
    c1.SaveAs("Dihiggs_deltaR2_%s"%pic_name+"_B3.png")



#_______________________________________________________
def getPurity(file,dir,den,num,xaxis,h_bins):
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    f = ROOT.TFile(file)
    t = f.Get(dir)
    nBins = int(h_bins[1:-1].split(',')[0])
    minBin = float(h_bins[1:-1].split(',')[1])
    maxBin = float(h_bins[1:-1].split(',')[2])
    
    b1 = ROOT.TH1F("b1","b1",nBins,minBin,maxBin)
    b1.GetYaxis().SetRangeUser(0.60,1.02)
    b1.GetYaxis().SetTitleOffset(1.2)
    b1.GetYaxis().SetNdivisions(520)
    b1.GetYaxis().SetTitle("TFTrack GE11dphicut Efficiency")
    b1.GetXaxis().SetTitle("%s of Simtrack"%xaxis)
    b1.SetTitle("Track reco in TrackFinder"+" "*16 + "CMS Simulation Preliminary")
    b1.SetStats(0)

    h1 = ROOT.TH1F("h1","h1",nBins,minBin,maxBin)
    t.Draw("abs(%s) >> h1"%xaxis,den)
    h2 = ROOT.TH1F("h2","h2",nBins,minBin,maxBin)
    t.Draw("abs(%s) >> h2"%xaxis,num)
    e = ROOT.TEfficiency(h2,h1)
    
    b1.Draw()
    e.SetLineColor(ROOT.kRed)
    e.Draw("same")
    legend = ROOT.TLegend(0.23,0.60,0.80,0.82)
#legend.SetFillColor(ROOT.kWhite)
    legend.SetHeader("PU140,simPt > 10Gev")
#legend.AddEntry(e1,"","l")
#legend.Draw("same")
    tex = ROOT.TLatex(0.30,0.50,"PU140, p_{T}^{sim} > 15Gev")
    tex.SetTextSize(0.05)
    tex.SetNDC()
    tex.Draw("same")
	
    c1.SaveAs("TFTrack_reco_eff_%s_Pt10_GE11dphi_PU140.pdf"%xaxis)
    c1.SaveAs("TFTrack_reco_eff_%s_Pt10_GE11dphi_PU140.png"%xaxis)

#_____________________________________________________________________________
def drawAll_1D(dir, treename, todraw,x_bins,x_title,cut,pic_name,text):
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    if not os.path.isdir(dir):
          print "ERROR: This is not a valid directory: ", dir
    ls = os.listdir(dir)
    tot = len(ls)
    rootfile = dir[:]+ls[0]
    tfile0 = ROOT.TFile(rootfile)
    t = tfile0.Get(treename)
    m = 0
    chain = ROOT.TChain(treename)
    e = ROOT.TH1F("e","e",xBins,xminBin,xmaxBin)
    b1 = ROOT.TH1F("b1","b1",xBins,xminBin,xmaxBin)
    b1.SetTitle("ttbar#rightarrow bbbarWW, B3"+" "*12 + "CMS Simulation Preliminary")
    b1.GetYaxis().SetTitle("Events")
    b1.GetXaxis().SetTitle("%s"%x_title)
    for x in ls:
	x = dir[:]+x
	chain.Add(x)
    chain.Draw(todraw+">>b1", cut)
    b1.Draw() 
    #print "chain ",chain, " b1 entries ",b1.GetEntries()
    tex = ROOT.TLatex(0.15,0.45,text)
    tex.SetTextSize(0.05)
    tex.SetNDC()
    tex.Draw("same")

    c1.SaveAs("ttbar_%s"%pic_name+"_All_B3.pdf")
    c1.SaveAs("ttbar_%s"%pic_name+"_All_B3.png")



#_____________________________________________________________________________
def drawAll_2D(dir, treename, todraw,x_bins,y_bins, x_title, y_title,cut,pic_name):

    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    
    b1 = ROOT.TH2F("b1","b1",xBins,xminBin,xmaxBin, yBins,yminBin,ymaxBin)
    b1.SetTitle("ttbar#rightarrow bbbarWW, B3"+" "*12 + "CMS Simulation Preliminary")
    b1.GetYaxis().SetTitle("%s"%y_title)
    b1.GetXaxis().SetTitle("%s"%x_title)
    if not os.path.isdir(dir):
          print "ERROR: This is not a valid directory: ", dir
    ls = os.listdir(dir)
    tot = len(ls)
    
    chain = ROOT.TChain(treename)
    for x in ls:
	x = dir[:]+x
	chain.Add(x)
    chain.Draw(todraw+">>b1",cut,"colz")
    b1.Draw("colz") 
    ROOT.gPad.SetLogz()

    
    c1.SaveAs("ttbar_%s"%pic_name+"_All_B3.pdf")
    c1.SaveAs("ttbar_%s"%pic_name+"_All_B3.png")
   

#______________________________________________________________________________
def buildTChain(dir,treename, rootfilename):

    if not os.path.isdir(dir):
          print "ERROR: This is not a valid directory: ", dir
    ls = os.listdir(dir)
    tot = len(ls)
    
    chain = ROOT.TChain(treename)
    for x in ls:
	x = dir[:]+x
	chain.Add(x)
   
    file = ROOT.TFile(rootfilename,"recreate") 
    chain.CloneTree(-1,"fast")
    file.Write()
    file.Close()

#_______________________________________________________________________________
if __name__ == "__main__":
     
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs-1M-B3-1071409.root"
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_100k_correctnu_0324_B3.root"
    file = "/eos/uscms/store/user/tahuang/DiHiggs/ttbar_100k_0414_B3.root"
    filedir = "/fdata/hepx/store/user/taohuang/ttbar_run2_PU20_bbana_cuts_1M_finalStates_V1/"

    dir = "ttbarAna/evtree"
    buildTChain(filedir,dir,"ttbar_TChain.root") 
    #htoWW_mass = "sqrt((mu1_energy+mu2_energy+nu1_energy+nu2_energy)**2-(mu1_px+mu2_px+nu1_px+nu2_px)**2-(mu1_py+mu2_py+nu1_py+nu2_py)**2-(mu1_pz+mu2_pz+nu1_pz+nu2_pz)**2)"
    
    
    W1_mass = "sqrt((mu1_energy+nu1_energy)**2-(mu1_px+nu1_px)**2-(mu1_py+nu1_py)**2-(mu1_pz+nu1_pz)**2)"
    W2_mass = "sqrt((mu2_energy+nu2_energy)**2-(mu2_px+nu2_px)**2-(mu2_py+nu2_py)**2-(mu2_pz+nu2_pz)**2)"
    Wmass_bins = "(40,60,100)"
    
    #drawAll_1D(filedir,dir,W1_mass,Wmass_bins,"reconstructed mass of W_{1}#rightarrow #mu_{1}#nu_{1} (W_{1} from #bar{t})", "1","W1_mass_1M_0605_finalStates")
    #drawAll_1D(filedir,dir,W2_mass,Wmass_bins,"reconstructed mass of W_{2}#rightarrow #mu_{2}#nu_{2} (W_{2} from t)", "1","W2_mass_1M_0605_finalStates")
    W1cand_mass = "mu1_mother_mass"
    W2cand_mass = "mu2_mother_mass"
    #draw1D(file,dir,W1cand_mass,Wmass_bins,"mass of W_{1} candidates ","1","W1cand_mass_100k_0413")
    #draw1D(file,dir,W2cand_mass,Wmass_bins,"mass of W_{2} candidates ","1","W2cand_mass_100k_0413")


    drawAll_1D(filedir,dir,"bjet_decendant_energy/bjet_energy","(101,-.005,1.005)","#frac{E_{b decendants}}{E_{bjet}}", "dR_bjet<0.1","decendantenergyratio_dR1_0605_V1", "#DeltaR<0.1, p_T>30, |#eta|<2.5")
    drawAll_1D(filedir,dir,"bjet_decendant_energy/bjet_energy","(101,-.005,1.005)","#frac{E_{b decendants}}{E_{bjet}}", "dR_bjet>0.4","decendantenergyratio_dR2_0605_V1", "#DeltaR>0.4, p_T>30, |#eta|<2.5")
    drawAll_1D(filedir,dir,"bbarjet_decendant_energy/bbarjet_energy","(101,-.005,1.005)","#frac{E_{b decendants}}{E_{bbarjet}}", "dR_bbarjet<0.1","decendantenergyratio_bjetdR1_0605_V1", "#DeltaR<0.1, p_T>30, |#eta|<2.5")
    drawAll_1D(filedir,dir,"bbarjet_decendant_energy/bbarjet_energy","(101,-.005,1.005)","#frac{E_{b decendant}}{E_{bbarjet}}", "dR_bbarjet>0.4","decendantenergyratio_bbarjetdR2_0605_V1", "#DeltaR>0.4, p_T>30, |#eta|<2.5")

    #drawAll_2D(filedir,dir,"(bjet_decendant_energy/bjet_energy):dR_bjet","(50,0,2)","(55,0,1.1)", "dR(b genParticle, b GenJet)", "#frac{E_{b decendants}}{E_{bjet}}","1","energyratioVsdR_bjet_0605_V1")

    tbar_mass = "sqrt((mu1_energy+nu1_energy+b2_energy)**2-(mu1_px+nu1_px+b2_px)**2-(mu1_py+nu1_py+b2_py)**2-(mu1_pz+nu1_pz+b2_pz)**2)"
    t_mass = "sqrt((mu2_energy+nu2_energy+b1_energy)**2-(mu2_px+nu2_px+b1_px)**2-(mu2_py+nu2_py+b1_py)**2-(mu2_pz+nu2_pz+b1_pz)**2)"
    tmass_bins = "(60,140,200)"
    #draw1D(file,dir,t_mass,tmass_bins,"reconstructed mass of t#rightarrowWb", "1","tmass_100k_0413")
    #draw1D(file,dir,tbar_mass,tmass_bins,"reconstructed mass of #bar{t}#rightarrowW#bar{b}", "1","tbarmass_100k_0413")
    tbarcand_mass = "tbar_mass"
    tcand_mass = "t_mass" 
    #draw1D(file,dir,tcand_mass,tmass_bins,"mass of t cand", "1","tcandmass_100k_0413")
    #draw1D(file,dir,tbarcand_mass,tmass_bins,"mass of #bar{t} cand", "1","tbarcandmass_100k_0413")
   # h2toh1h1_mass = "sqrt((mu1_energy+mu2_energy+nu1_energy+nu2_energy+b1_energy+b2_energy)**2-(mu1_px+mu2_px+nu1_px+nu2_px+b1_px+b2_px)**2-(mu1_py+mu2_py+nu1_py+nu2_py+b1_py+b2_py)**2-(mu1_pz+mu2_pz+nu1_pz+nu2_pz+b1_pz+b2_pz)**2)"
    #h2toh1h1_cut = "h2tohh"
    #h2mass_bins_6 = "(200,400,600)"
    #h2mass_bins_3 = "(200,250,450)"
    #draw1D(file,dir,h2toh1h1_mass,h2mass_bins_3,"reconstructed mass of h2#rightarrow BB#mu#nu#mu#nu", h2toh1h1_cut,"h2toh1h1_mass_1M_0413")
#draw jets mass and draw h2 reconstruction mass from jets and muons neutrinos
    jets_mass_bins = "(200, 100,400)"
    #draw1D(file,dir,"jets_mass",jets_mass_bins,"invariant mass of all decendants from b+#bar{b}", htoWW_cut,"jets_mass_1M_0413")
    h2_mass_jets = "sqrt((mu1_energy+mu2_energy+nu1_energy+nu2_energy+jets_energy)**2-(mu1_px+mu2_px+nu1_px+nu2_px+jets_px)**2-(mu1_py+mu2_py+nu1_py+nu2_py+jets_py)**2-(mu1_pz+mu2_pz+nu1_pz+nu2_pz+jets_pz)**2)"
    #draw1D(file,dir,h2_mass_jets,h2mass_bins_3,"invariant mass of all decendants from b+#bar{b}", htoWW_cut,"h2_mass_jets_1M_0413")
    met_bins = "(150,0,150)"
    met = "met"
    #draw1D(file,dir,met,met_bins,"Simulated #slash{E}_{T}", "1","MET_1M_0413")
    #drawAll_1D(filedir,dir,"dR_bjet","(50,0,2)","deltaR(bjet, b genParticle)", "1","dR_bjet_cuts_1M_0605_finalStates")
    #drawAll_1D(filedir,dir,"dR_bbarjet","(50,0,2)","deltaR(bbarjet, bbar genParticle)", "1","dR_bbarjet_cuts_1M_0605_finalStates")
    
### as a reference for monitoring plots in MMC
    wmass_offshell_bins = "(60,0.0,60.0)" 
    wmass_onshell_bins = "(50,40.0,90.0)" 
    eta_bins = "(30,-6,6)"
    nu1_eta = "nu1_eta"
    nu2_eta = "nu2_eta"

    #offshell_nupt_bins = "(25,0,100)"
    offshell_nupt_bins = "(25,0,100)"
    onshell_nupt_bins = "(25,0,125)"
    nu1_pt = "sqrt(nu1_px**2+nu1_py**2)"
    nu2_pt = "sqrt(nu2_px**2+nu2_py**2)"
    onshellW_1_cut = "mu1_mother_mass>mu2_mother_mass"
    offshellW_1_cut = "mu2_mother_mass>mu1_mother_mass"
    delta_phi = "(25,-3.1415,3.1415)"
    delta_eta = "(50,-5.0,5.0)"
    #deltaR1(file,dir,delta_eta,delta_phi,h2toh1h1_cut,"h2toh1h1_0223")
    #deltaR2(file,dir,delta_eta,delta_phi,h2toh1h1_cut,"h2toh1h1_0223")

