
import ROOT
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
    b1.SetTitle("h2#rightarrow hh#rightarrow BBWW, B3"+" "*8 + "CMS Simulation Preliminary")
    b1.GetYaxis().SetTitle("Events")
    b1.GetXaxis().SetTitle("%s"%x_title)
    #b1.SetStats(0)
    
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



#_______________________________________________________________________________
if __name__ == "__main__":
     
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs-1M-B3-1071409.root"
    file = "DiHiggs_test.root"
    dir = "DiHiggsWWAna/evtree"
  
    #htoWW_mass = "sqrt((mu1_energy+mu2_energy+nu1_energy+nu2_energy)**2-(mu1_px+mu2_px+nu1_px+nu2_px)**2-(mu1_py+mu2_py+nu1_py+nu2_py)**2-(mu1_pz+mu2_pz+nu1_pz+nu2_pz)**2)"
    
    
    htoWW_mass = "sqrt((mu1_energy+mu2_energy+nu1_energy+nu2_energy)**2-(mu1_px+mu2_px+nu1_px+nu2_px)**2-(mu1_py+mu2_py+nu1_py+nu2_py)**2-(mu1_pz+mu2_pz+nu1_pz+nu2_pz)**2)"
    hmass_bins = "(50,100,150)"
    hmass_bins1 = "(10,124.5,125.5)" 
    htoWW_cut = "1" 
    
    draw1D(file,dir,htoWW_mass,hmass_bins1,"mass of h#rightarrow WW, reconstruction from(#mu#mu#nu#nu)", htoWW_cut,"htoWW_final_mass_100k_cmssw-0322")
    
    htoBB_mass = "sqrt((b1_energy+b2_energy)**2-(b1_px+b2_px)**2-(b1_py+b2_py)**2-(b1_pz+b2_pz)**2)"
    htoBB_cut = "1"
    draw1D(file,dir,htoBB_mass,hmass_bins1,"mass of h#rightarrow BB", htoBB_cut,"htoBB_mass_100k_cmssw-0322")
   
    h2toh1h1_mass = "sqrt((mu1_energy+mu2_energy+nu1_energy+nu2_energy+b1_energy+b2_energy)**2-(mu1_px+mu2_px+nu1_px+nu2_px+b1_px+b2_px)**2-(mu1_py+mu2_py+nu1_py+nu2_py+b1_py+b2_py)**2-(mu1_pz+mu2_pz+nu1_pz+nu2_pz+b1_pz+b2_pz)**2)"
    h2toh1h1_cut = "1"
    h2mass_bins_6 = "(100,400,600)"
    h2mass_bins_3 = "(100,250,450)"
    draw1D(file,dir,h2toh1h1_mass,h2mass_bins_3,"mass of h2#rightarrow BB#mu#nu#mu#nu", h2toh1h1_cut,"h2toh1h1_mass_100k_cmssw-0322")
    
    htoWW_mass1 = "sqrt((mu1_mother_energy+mu2_mother_energy)**2-(mu1_mother_px+mu2_mother_px)**2-(mu1_mother_py+mu2_mother_py)**2-(mu1_mother_pz+mu2_mother_pz)**2)"
    draw1D(file,dir,htoWW_mass1,hmass_bins1,"mass of h#rightarrow WW, reconstruction from (WW)", htoWW_cut,"htoWW_median_mass_100k_cmssw-0322")
#    draw1D(file,dir,"htoWW_mass",hmass_bins1,"mass of h in generation", htoWW_cut,"h_mass_100k_rivet-0321")

    phi_bins = "(25,-3.1415,3.1415)"
    nu1_phi = "nu1_phi";
    nu2_phi = "nu2_phi"
    #draw1D(file,dir,nu1_phi,phi_bins,"phi of #nu1", htoWW_cut,"htoWW_nu1_phi_0223")
    #draw1D(file,dir,nu2_phi,phi_bins,"phi of #nu2", htoWW_cut,"htoWW_nu2_phi_0223")
    delta_phi = "(25,-3.1415,3.1415)"
    delta_eta = "(50,-5.0,5.0)"
    #deltaR1(file,dir,delta_eta,delta_phi,h2toh1h1_cut,"h2toh1h1_0223")
    #deltaR2(file,dir,delta_eta,delta_phi,h2toh1h1_cut,"h2toh1h1_0223")
