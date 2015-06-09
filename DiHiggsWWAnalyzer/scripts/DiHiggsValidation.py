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

#_____________________________________________________________________________
def drawAll_1D(dir, treename, todraw,x_bins,x_title,cut,pic_name):
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    if not os.path.isdir(dir):
          print "ERROR: This is not a valid directory: ", inputDir
    ls = os.listdir(dir)
    tot = len(ls)
    rootfile = dir[:]+ls[0]
    tfile0 = ROOT.TFile(rootfile)
    t = tfile0.Get(treename)
    m = 0
    chain = ROOT.TChain(treename)
    e = ROOT.TH1F("e","e",xBins,xminBin,xmaxBin)
    b1 = ROOT.TH1F("b1","b1",xBins,xminBin,xmaxBin)
    b1.SetTitle("h2#rightarrow hh#rightarrow BBWW, B3"+" "*12 + "CMS Simulation Preliminary")
    b1.GetYaxis().SetTitle("Events")
    b1.GetXaxis().SetTitle("%s"%x_title)
    for x in ls:
	x = dir[:]+x
	chain.Add(x)
	"""
	if m>0:
		x = dir[:]+x
#	b1.Add(draw1D(x, treename, todraw,x_bins,x_title,cut,pic_name, m))
    		f = ROOT.TFile(x)
    		t.AddFriend(treename,f)
	#t.Draw(todraw+">>e",cut)
	#print "e entries ",e.GetEntries(), " b1 entries ", b1.GetEntries()
		f.Close()
	#b1.Add(e)
	m = m+1
	"""
    chain.Draw(todraw+">>b1", cut)
    b1.Draw() 
    print "chain ",chain, " b1 entries ",b1.GetEntries()

    c1.SaveAs("Dihiggs_%s"%pic_name+"_All_B3.pdf")
    c1.SaveAs("Dihiggs_%s"%pic_name+"_All_B3.png")

#_______________________________________________________________________________
if __name__ == "__main__":
     
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs-1M-B3-1071409.root"
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_100k_correctnu_0324_B3.root"
    filedir = "/fdata/hepx/store/user/taohuang/DiHiggs_run2_PU0_htobbana_cuts_750k_B3/"
    dir = "DiHiggsWWAna/evtree"
  
    #htoWW_mass = "sqrt((mu1_energy+mu2_energy+nu1_energy+nu2_energy)**2-(mu1_px+mu2_px+nu1_px+nu2_px)**2-(mu1_py+mu2_py+nu1_py+nu2_py)**2-(mu1_pz+mu2_pz+nu1_pz+nu2_pz)**2)"
    
    
    htoWW_mass = "sqrt((mu1_energy+mu2_energy+nu1_energy+nu2_energy)**2-(mu1_px+mu2_px+nu1_px+nu2_px)**2-(mu1_py+mu2_py+nu1_py+nu2_py)**2-(mu1_pz+mu2_pz+nu1_pz+nu2_pz)**2)"
    hmass_bins = "(50,100,300)"
    hmass_bins1 = "(50,0,200)" 
    htoWW_cut = "h2tohh" 
    
    #draw1D(file,dir,htoWW_mass,hmass_bins1,"reconstructed mass of h#rightarrow WW, reconstruction from(#mu#mu#nu#nu)", htoWW_cut,"htoWW_final_mass_1M_mediateStates_0325")
    bjet_pt = "sqrt(bjet_px**2+bjet_py**2)" 
    bbarjet_pt = "sqrt(bbarjet_px**2+bbarjet_py**2)" 
    htoBB_mass = "sqrt((b1_energy+b2_energy)**2-(b1_px+b2_px)**2-(b1_py+b2_py)**2-(b1_pz+b2_pz)**2)"
    htoBBJets_mass = "sqrt((bjet_energy+bbarjet_energy)**2-(bjet_px+bbarjet_px)**2-(bjet_py+bbarjet_py)**2-(bjet_pz+bbarjet_pz)**2)"
    htoBBJetstot_mass = "sqrt((bjet_energy_tot+bbarjet_energy_tot)**2-(bjet_px_tot+bbarjet_px_tot)**2-(bjet_py_tot+bbarjet_py_tot)**2-(bjet_pz_tot+bbarjet_pz_tot)**2)"
    htoBB_cut = "h2tohh"
    drawAll_1D(filedir,dir,"dR_bjet","(50,0,2)","deltaR(bjet, b genParticle)", htoBB_cut,"dR_bjet_cuts_1M_0605")
    drawAll_1D(filedir,dir,"dR_bbarjet","(50,0,2)","deltaR(bbarjet, bbar genParticle)", htoBB_cut,"dR_bbarjet_cuts_1M_0605")
    drawAll_1D(filedir,dir,bjet_pt, "(100,0,200)","p#_T(bjet)",htoBB_cut,"bjetPt_cuts_1M_0605")
    drawAll_1D(filedir,dir,bbarjet_pt, "(100,0,200)","p#_T(bbarjet)",htoBB_cut,"bbarjetPt_cuts_1M_0605")
    #draw1D(file,dir,htoBB_mass,hmass_bins1,"reconstructed mass of h#rightarrow BB", htoBB_cut,"htoBB_mass_1M_mediateStates_0325")
    drawAll_1D(filedir,dir,htoBBJets_mass,hmass_bins1,"reconstructed mass of h#rightarrow BB, from GenJet (closest)", htoBB_cut,"htoBBjets_cuts_mass_1M_0606")
    #drawAll_1D(filedir,dir,htoBBJetstot_mass,hmass_bins1,"reconstructed mass of h#rightarrow BB, from GenJet (sum)", htoBB_cut,"htoBBjetstot_mass_1M_0606")
   
    h2toh1h1_mass = "sqrt((mu1_energy+mu2_energy+nu1_energy+nu2_energy+b1_energy+b2_energy)**2-(mu1_px+mu2_px+nu1_px+nu2_px+b1_px+b2_px)**2-(mu1_py+mu2_py+nu1_py+nu2_py+b1_py+b2_py)**2-(mu1_pz+mu2_pz+nu1_pz+nu2_pz+b1_pz+b2_pz)**2)"
    h2toh1h1Jets_mass = "sqrt((mu1_energy+mu2_energy+nu1_energy+nu2_energy+bjet_energy+bbarjet_energy)**2-(mu1_px+mu2_px+nu1_px+nu2_px+bjet_px+bbarjet_px)**2-(mu1_py+mu2_py+nu1_py+nu2_py+bjet_py+bbarjet_py)**2-(mu1_pz+mu2_pz+nu1_pz+nu2_pz+bjet_pz+bbarjet_pz)**2)"
    h2toh1h1_cut = "h2tohh"
    h2mass_bins_6 = "(200,400,600)"
    h2mass_bins_3 = "(200,250,450)"
    #draw1D(file,dir,h2toh1h1_mass,h2mass_bins_3,"reconstructed mass of h2#rightarrow BB#mu#nu#mu#nu", h2toh1h1_cut,"h2toh1h1_mass_1M_mediateStates_0325")
    
    htoWW_mass1 = "sqrt((mu1_mother_energy+mu2_mother_energy)**2-(mu1_mother_px+mu2_mother_px)**2-(mu1_mother_py+mu2_mother_py)**2-(mu1_mother_pz+mu2_mother_pz)**2)"
    #draw1D(file,dir,htoWW_mass1,hmass_bins1," reconstructed mass of h#rightarrow WW, reconstruction from (WW)", htoWW_cut,"htoWW_median_mass_1M_mediateStates_0325")
    #draw1D(file,dir,"htoWW_mass",hmass_bins1,"true h mass from  h candidates in generation", htoWW_cut,"htoWW_mass_100k_Gen_0324")
    #draw1D(file,dir,"h2tohh_mass",h2mass_bins_3,"true h2 mass from h2 candidates in generation", htoWW_cut,"h2_mass_100k_Gen_0324")

#draw jets mass and draw h2 reconstruction mass from jets and muons neutrinos
    jets_mass_bins = "(200, 100,400)"
    #draw1D(file,dir,"jets_mass",jets_mass_bins,"invariant mass of all decendants from b+#bar{b}", htoWW_cut,"jets_mass_1M_mediateStates_0325")
    h2_mass_jets = "sqrt((mu1_energy+mu2_energy+nu1_energy+nu2_energy+jets_energy)**2-(mu1_px+mu2_px+nu1_px+nu2_px+jets_px)**2-(mu1_py+mu2_py+nu1_py+nu2_py+jets_py)**2-(mu1_pz+mu2_pz+nu1_pz+nu2_pz+jets_pz)**2)"
    #draw1D(file,dir,h2_mass_jets,h2mass_bins_3,"invariant mass of all decendants from b+#bar{b}", htoWW_cut,"h2_mass_jets_1M_mediateStates_0325")
    met_bins = "(150,0,150)"
    met = "met"
    #draw1D(file,dir,met,met_bins,"Simulated #slash{E}_{T}", "1","MET_1M_mediateStates_0325")
    
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
#    draw1D_combined(file,dir,"mu1_mother_mass","mu2_mother_mass", wmass_onshell_bins,"Simulated M_{W}^{onshell}",onshellW_1_cut,offshellW_1_cut,"onshellW_mass_1M_mediateStates_0325")
    
    #draw1D_combined(file,dir,nu1_pt,nu2_pt, onshell_nupt_bins,"Simulated p_{T#nu}^{onshellW}",onshellW_1_cut,offshellW_1_cut,"onshell_nupt_1M_mediateStates_0325")
    delta_phi = "(25,-3.1415,3.1415)"
    delta_eta = "(50,-5.0,5.0)"
    #deltaR1(file,dir,delta_eta,delta_phi,h2toh1h1_cut,"h2toh1h1_0223")
    #deltaR2(file,dir,delta_eta,delta_phi,h2toh1h1_cut,"h2toh1h1_0223")
   
    onoffshellWmass1 = "mu1_mother_mass:mu2_mother_mass"
    onoffshellWmass2 = "mu2_mother_mass:mu1_mother_mass"
     
    #draw2D_combined(file,dir, onoffshellWmass2, onoffshellWmass1, wmass_onshell_bins,wmass_offshell_bins,"Simulated M_{W}^{onshell}","Simulated M_{W}^{offshell}",onshellW_1_cut,offshellW_1_cut,"onshellVsoffshell_Wmass_1M_mediateStates_0325")

