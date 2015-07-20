from array import *
import ROOT
import math
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

ROOT.gStyle.SetPadLeftMargin(0.12)
#ROOT.gStyle.SetPadRightMargin(0.14)
ROOT.gStyle.SetPadTopMargin(0.15)
ROOT.gStyle.SetPadBottomMargin(0.13)


#used to scale hist
count = 100000.0

h2tohhmass = "h2tohh_Mass"
h2truemass = "mass_h2_true"
offshellWmass = "offshellW_Mass"
offshellWtruemass = "mass_offshellW_true"


#___________________________________________
def draw1D(file,dir,todraw,title,x_bins,x_title,tex,pic_name):
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    f = ROOT.TFile(file)
    tree = f.Get(dir)
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    
    cut = ["weight","weight*(control<2)","weight*(control>1)"]
    b1 = ROOT.TH1F("b1","b1",xBins,xminBin,xmaxBin)
    b1.SetTitle("%s"%title+" "*8 + "CMS Simulation Preliminary")
 #   b1.Draw()
    combined = 1
    if combined:
    #b1.SetStats(0)
        h2mass_bins = "(80,260,420)"
	todraw9 = "h2tohh_Mass"
        todrawtrue9 = "mass_h2_true"
	hs_title9 = "PDF of M_{H}"
	hs = gethiststack(tree, todraw9,todrawtrue9, h2mass_bins, hs_title9)
	#hs9.SetName("%s"%tree.GetTitle()+"_%s"%todraw9)
        """
        hs = ROOT.THStack("hs","hs");
        hs.SetTitle("%s"%title+" "*8 + "CMS Simulation Preliminary")
        hist1 = hist_1D(tree, h2tohhmass, x_bins, cut,0)
        hist1.SetLineColor(ROOT.kRed)
	hist1.SetLineWidth(4)
       # hist1.SetFillColor(ROOT.kRed)
      #  hist1.SetFillStyle(4050)
        max = hist1.GetMaximum()
        #hist1.Draw("E3") 
        
        hist2 = hist_1D(tree, h2tohhmass, x_bins, cut,1)
        #hist2.SetLineColor(ROOT.kBlue)
        #hist2.SetLineWidth(4)
	#hist2.SetMarkerStyle(21)
        hist2.SetFillColor(ROOT.kBlue)
        hist2.SetFillStyle(3445)
        #hist2.Draw("sameE3") 
   
        hist3 = hist_1D(tree, h2tohhmass, x_bins, cut,2)
        #hist3.SetLineColor(ROOT.kGreen)#3
        #hist3.SetLineWidth(4)#3
#	hist3.SetMarkerStyle(20)
#	hist3.SetMarkerColor(ROOT.kGreen)
        hist3.SetFillColor(ROOT.kGreen)
        hist3.SetFillStyle(3454)
        #hist3.Draw("sameE3") 

        hist4 = hist_1D(tree, h2truemass, x_bins, cut,0)
        hist4.SetFillColor(ROOT.kMagenta)
   #     hist4.SetFillStyle(4050)
        bin_max = hist4.GetMaximumBin()
	hist4.SetBinContent(bin_max, max/4.0)
        
        #hs.GetXaxis().SetTitle("%s"%x_title)
        #hs.GetYaxis().SetTitle("Possibility")
        ROOT.gStyle.SetHatchesLineWidth(2)
	hs.Add(hist1)
	hs.Add(hist2)
	hs.Add(hist3)
	hs.Add(hist4)
        """
     	hs.Draw("nostack")
     #   c1.SaveAs("contour_%d"%(n%9)+".png") 
        hs.GetHistogram().GetXaxis().SetTitle("%s"%x_title)
        hs.GetHistogram().GetYaxis().SetTitle("Possibility")
        c1.Update()
    tex2 = ROOT.TLatex(0.15,0.6,"%s "%tex)
    tex2.SetTextSize(0.04)
    tex2.SetTextFont(42)
    tex2.SetNDC()
    #tex2.Draw("same")
	
    c1.SaveAs("Dihiggs_MMC_style1_%s"%pic_name+"_B3.pdf")
    c1.SaveAs("Dihiggs_MMC_style1_%s"%pic_name+"_B3.png")
    c1.SaveAs("Dihiggs_MMC_style1_%s"%pic_name+"_B3.C")

#___________________________________________
def draw2D(file,dir,todraw,xaxis,yaxis,title,x_bins,x_title,y_bins,y_title,tex,pic_name):
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    f = ROOT.TFile(file)
    tree = f.Get(dir)
#    todraw = "%s"%yaxis+":%s"%xaxis
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    
    cut = ["weight","weight*(control<2)","weight*(control>1)"]
    b1 = ROOT.TH2F("b1","b1",xBins,xminBin,xmaxBin, yBins, yminBin, ymaxBin)
    b1.SetTitle("%s"%title+" "*8 + "CMS Simulation Preliminary")
 #   b1.Draw()
    combined = 1
    if combined:
    #b1.SetStats(0)
        hist1 = hist_2D(tree, todraw, x_bins, y_bins, cut,0)
	hist1.SetLineColor(ROOT.kRed)
        
        hist2 = hist_2D(tree, todraw, x_bins, y_bins, cut,1)
        hist2.SetLineColor(ROOT.kBlue)
        #hist2.SetFillColor(42)
        #hist2.SetFillStyle(4050)
        #hist2.Draw("sameE3")
	 
   
        hist3 = hist_2D(tree, todraw, x_bins, y_bins, cut,2)
        hist3.SetLineColor(ROOT.kGreen)#3

        #hist3.SetFillColor(46)
      #  hist3.SetFillStyle(4050)
        #hist3.Draw("sameE3") 
	var_x = getTrueValue(tree,xaxis,"(100,-6,6)")
	var_y = getTrueValue(tree,yaxis,"(100,-3.1415,3.1415)")
	array_x = array('d',[var_x])
	array_y = array('d',[var_y])
	graph = ROOT.TGraph(1, array_x,array_y)
        graph.SetMarkerColor(ROOT.kMagenta)
 	graph.SetMarkerStyle(15)
	graph.SetMarkerSize(2)

        hist2.GetXaxis().SetTitle("#eta")
        hist2.GetYaxis().SetTitle("#phi")
	hist2.Draw("CONT3")
	hist3.Draw("CONT3same")
	graph.Draw("psame")
	c1.Update()


    tex2 = ROOT.TLatex(0.15,0.6,"%s "%tex)
    tex2.SetTextSize(0.04)
    tex2.SetTextFont(42)
    tex2.SetNDC()
    #tex2.Draw("same")
	
    c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3.pdf")
    c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3.png")


#___________________________________________
def draw2DAll(file,todraw,title,x_bins,x_title,y_bins,y_title,tex,pic_name,cut = "weight"):
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    f = ROOT.TFile(file)
    list = f.GetListOfKeys()
    key = list.At(0)
    obj = key.ReadObj() #DiHiggsWWAna
    sub_list = obj.GetListOfKeys()
    sub_key = sub_list.At(1)

#    todraw = "%s"%yaxis+":%s"%xaxis
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    
    c1 = ROOT.TCanvas()
    b1 = ROOT.TH2F("b1","b1",xBins,xminBin,xmaxBin, yBins, yminBin, ymaxBin)
    b1.SetTitle("%s"%title+" "*8 + "CMS Simulation Preliminary")
    #cut = ["weight","weight*(control<2)","weight*(control>1)"]
    n = 1
    while sub_key: 
	c1.cd()
        c1.Clear()

	tree = sub_key.ReadObj()
        tree.Draw(todraw+">>b1",cut)
    	b1.Draw("colz")
        #b1.SetStats(0)
        c1.Update()
    	b1.GetXaxis().SetTitle("%s"%x_title)
    	b1.GetYaxis().SetTitle("%s"%y_title)
	stats1 = b1.GetListOfFunctions().FindObject("stats")#.Clone("stats1")
        stats1.SetX1NDC(.4)
        stats1.SetX2NDC(.55)
        stats1.SetY1NDC(.6)
        stats1.SetY2NDC(.9)
       # stats1.Draw("same")
	var_x1 = getTrueValue(tree,"mass_onshellW_true",x_bins)
	var_y1 = getTrueValue(tree,"pt_nuonshellW_true",y_bins)
	#array_x = array('d',[var_x1])
	#array_y = array('d',[var_y1])
    	b2 = ROOT.TH2F("b2_%s"%tree.GetTitle(),"b2",xBins,xminBin,xmaxBin, yBins, yminBin, ymaxBin)
	#graph2d = ROOT.TGraph(1, array_x,array_y)
        b2.SetMarkerColor(ROOT.kMagenta)
	bin = b2.FindBin(var_x1,var_y1)
	b2.SetBinContent(bin,1)
 	b2.SetMarkerStyle(15)
	b2.SetMarkerSize(1)
        b2.Draw("same")
        #print "array x ", array_x, " array y", array_y
        c1.Update()
    	c1.SaveAs("MMC_%s"%tree.GetTitle()+"_%s"%pic_name+"_0406_weight1.pdf")
    	c1.SaveAs("MMC_%s"%tree.GetTitle()+"_%s"%pic_name+"_0406_weight1.png")
    	c1.SaveAs("MMC_%s"%tree.GetTitle()+"_%s"%pic_name+"_0406_weight1.C")
        n = n+1
        sub_key = sub_list.At(n) 

#___________________________________________
def drawh2Mass_combined(file,x_bins,pic_name):
   
    
    f = ROOT.TFile(file)
    list = f.GetListOfKeys()
    key = list.At(0)
    obj = key.ReadObj() #DiHiggsWWAna
    sub_list = obj.GetListOfKeys()
    sub_key = sub_list.At(1)

    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
     
     
    cut = ["weight","weight*(control<2)","weight*(control>1)"]
    c1 = ROOT.TCanvas()
    #ROOT.gPad.SetFillStyle(4050)
    n = 0;
    tex = "Probability Density Function for M_{H} in MMC"
    tex2 = ROOT.TLatex(0.25,0.95,"%s "%tex)
    tex2.SetTextSize(0.04)
    tex2.SetTextFont(42)
    tex2.SetNDC()
    #contour = sub_key.ReadObj()
    #contourhist = hist_1D(contour, h2tohhmass, x_bins, cut,0)
    #contourhist.Draw()
    #c1.SaveAs("contour.png")
    hslist = []
    while sub_key: 
	print " n", n,"n%9", n%9 
    	if n%9 == 0:
                c1.Clear()
		c1.Divide(3,3)
		print  "c1 divided "
        c1.cd(n%9+1)
	tree = sub_key.ReadObj()
        #tree.Print()
        hs = ROOT.THStack("hs_%d"%(n%9)," ");
        hist1 = hist_1D(tree, h2tohhmass, x_bins, cut,0)
	hist1.SetMarkerStyle(21)
        hist1.SetFillColor(ROOT.kRed)
        #hist1.SetFillColor(16)
      #  hist1.SetFillStyle(4050)
        max = hist1.GetMaximum()
        #hist1.Draw("E3") 
	hs.Add(hist1)
        
        hist2 = hist_1D(tree, h2tohhmass, x_bins, cut,1)
        hist2.SetLineColor(ROOT.kBlue)
        hist2.SetLineWidth(4)
        #hist2.SetFillColor(42)
        #hist2.SetFillStyle(4050)
        #hist2.Draw("sameE3") 
	hs.Add(hist2)
   
        hist3 = hist_1D(tree, h2tohhmass, x_bins, cut,2)
        hist3.SetLineColor(ROOT.kGreen)#3
        hist3.SetLineWidth(4)
        #hist3.SetFillColor(46)
      #  hist3.SetFillStyle(4050)
        #hist3.Draw("sameE3") 
	hs.Add(hist3)

        hist4 = hist_1D(tree, h2truemass, x_bins, cut,0)
        hist4.SetFillColor(ROOT.kMagenta)
        hist4.SetFillStyle(4050)
        bin_max = hist4.GetMaximumBin()
	hist4.SetBinContent(bin_max, max/4.0)
        hs.Add(hist4)
        
        #hs.SetTitle("%s"%title+" "*8 + "CMS Simulation Preliminary")
        hslist.append(hs)
     	hs.Draw("nostack")
     #   c1.SaveAs("contour_%d"%(n%9)+".png") 
	if n%9 == 8:
		c1.cd()
        	tex2.Draw("same")
		c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".pdf")
    		c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".png")
    		c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".C")
        n = n+1
        sub_key = sub_list.At(n+1)
    if n%9 != 0:
        c1.cd()
        tex2.Draw("same")
	c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".pdf")
    	c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".png")
    	c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".C")



#___________________________________________
def draw_combined(file,todraw,todrawtrue,x_bins,tex,pic_name):
    
    f = ROOT.TFile(file)
    list = f.GetListOfKeys()
    key = list.At(0)
    obj = key.ReadObj() #DiHiggsWWAna
    sub_list = obj.GetListOfKeys()
    sub_key = sub_list.At(1)

    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    
    cut = ["weight","weight*(control<2)","weight*(control>1)"]
    c1 = ROOT.TCanvas()
    c1.Divide(3,3)
    n = 0;
    #tex = "contour"
    tex2 = ROOT.TLatex(0.25,.950,"%s "%tex)
    tex2.SetTextSize(0.04)
    tex2.SetTextFont(42)
    tex2.SetNDC()
    #contour = sub_key.ReadObj()
    #contourhist = hist_1D(contour, h2tohhmass, x_bins, cut,0)
    #contourhist.Draw()
    #c1.SaveAs("contour.png")
    hslist = []
    while sub_key: 
	print " n", n,"n%9", n%9 
    	if n%9 == 0:
                c1.Clear()
		c1.Divide(3,3)
		print  "c1 divided "
        c1.cd(n%9+1)
	tree = sub_key.ReadObj()
        #tree.Print()
        hs = ROOT.THStack("hs_%d"%(n%9)," ");
        hist1 = hist_1D(tree, todraw, x_bins, cut,0)
        hist1.SetFillColor(ROOT.kRed)
        #hist1.Draw("E3") 
	hs.Add(hist1)
        
        hist2 = hist_1D(tree, todraw, x_bins, cut,1)
        hist2.SetLineColor(ROOT.kBlue)
        hist2.SetLineWidth(6)
        #hist2.SetLineStyle(20)
        #hist2.Draw("sameE3") 
	hs.Add(hist2)
   
        hist3 = hist_1D(tree, todraw, x_bins, cut,2)
        hist3.SetLineColor(ROOT.kGreen)
        hist3.SetLineWidth(4)
        #hist2.SetLineStyle(21)
        #hist3.Draw("sameE3") 
	hs.Add(hist3)

        hist4 = hist_1D(tree, todrawtrue, x_bins, cut,0)
        hist4.SetFillColor(ROOT.kMagenta)
        max = hist1.GetMaximum()
        bin_max = hist4.GetMaximumBin()
	hist4.SetBinContent(bin_max, max/4.0)
	hist4.SetBinContent(bin_max-1, max/4.0)
	#hist4.SetBinContent(bin_max+1, max/4.0)
        hs.Add(hist4)
        
        hslist.append(hs)
     	hs.Draw("nostack")
     #   c1.SaveAs("contour_%d"%(n%9)+".png") 
	if n%9 == 8:
		c1.cd()
        	tex2.Draw("same")
		c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".pdf")
    		c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".png")
    		c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".C")
        n = n+1
        sub_key = sub_list.At(n+1)
    if n%9 != 0:
        c1.cd()
        tex2.Draw("same")
	c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".pdf")
    	c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".png")
    	c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".C")

        

#____________________________________________________________________________
def monitoringMMC(file,pic, cut = ["weight","weight*(control<2)","weight*(control>1)"]):

    num = gettreenum(file)
    #cut = ["weight*(nu_onshellW_pt>10)","weight*(control<2&&nu_onshellW_pt>10)","weight*(control>1&&nu_onshellW_pt>10)"]
    #cut = ["weight","weight*(control<2)","weight*(control>1)"]
   # cut_weight1 = ["weight1","weight1*(control<2)","weight1*(control>1)"]
    c1 = ROOT.TCanvas("c1","c1",1200,900)
    n = 0
    m = 1
    while m<num:
    	f = ROOT.TFile(file)
    	list = f.GetListOfKeys()
    	key = list.At(0)
    	obj = key.ReadObj() #DiHiggsWWAna
    	sub_list = obj.GetListOfKeys()
    	sub_key = sub_list.At(m)
	tree = sub_key.ReadObj()
	if tree.GetEntries()<100:
		m =m+1
		f.close()
		continue
        #h2_mass_true = getTrueValue(tree,"mass_h2_true","(400,200,600)")
#        h2_mass_reconstructed = geth2MassMostProba(tree,cut[0]) 
#        if (h2_mass_reconstructed < 450 and h2_mass_reconstructed > 250):
#		m = m+1
#		f.Close()
#		continue
    
	#eta-phi of nu_onshellW contour plot, nu_onshellW_ete/phi
        #mmctree_618->Draw("nu_onshellW_phi:nu_onshellW_eta>>th2","control<2")
        eta_bins = "(50,-6,6)"
        phi_bins = "(40,-3.1415,3.1415)"
	todraw1 = "nu_onshellW_phi:nu_onshellW_eta"
        
        todrawtrue1 = "phi_nuonshellW_true:eta_nuonshellW_true"
        hs_title1 = "#phi_{#nu}^{onshellW} Vs #eta_{#nu}^{onshellW}"
        hs_Xtitle1 = "#eta_{#nu}^{onshellW}"
        hs_Ytitle1 = "#phi_{#nu}^{onshellW}"
        hs1 = gethiststack_2d(tree, todraw1, todrawtrue1, eta_bins, phi_bins, hs_title1, cut)
        hist1_true = hist_2D(tree, todrawtrue1, eta_bins, phi_bins, cut,1)
        hist1_true.SetMarkerStyle(7) 
        hist1_true.SetMarkerSize(1)
        hist1_true.SetMarkerColor(ROOT.kMagenta)
	onshellmuon_eta = getTrueValue(tree,"mu_onshellW_eta","(200,-6,6)")
	onshellmuon_phi = getTrueValue(tree,"mu_onshellW_phi","(200,-3.1415,3.1415)")
	array_mu1_eta = array('d',[onshellmuon_eta])
	array_mu1_phi = array('d',[onshellmuon_phi])
	graph1_muon = ROOT.TGraph(1, array_mu1_eta,array_mu1_phi)
        graph1_muon.SetMarkerColor(ROOT.kBlack)
 	graph1_muon.SetMarkerStyle(16)
	graph1_muon.SetMarkerSize(1)


	#"nu_onshellW_phi:nu_onshellW_eta"
	#eta of nuonshellW	
        todraw10 = "nu_onshellW_eta"
	todrawtrue10 = "eta_nuonshellW_true"
        hs_title10 = "#eta_{#nu}^{onshellW}"
	hs10 = gethiststack(tree, todraw10, todrawtrue10, eta_bins, hs_title10, cut)

	tex10 = ROOT.TLatex(0.4,0.7,"#eta_{#mu}^{onshell}=%.2f"%onshellmuon_eta)
        tex10.SetTextSize(0.08)
        tex10.SetTextFont(42)
        tex10.SetNDC()
	
	#phi of nuonshellW	
        todraw11 = "nu_onshellW_phi"
	todrawtrue11 = "phi_nuonshellW_true"
        hs_title11 = "#phi_{#nu}^{onshellW}"
	hs11 = gethiststack(tree, todraw11, todrawtrue11, phi_bins, hs_title11, cut)
	tex11 = ROOT.TLatex(0.4,0.7,"#phi_{#mu}^{onshell}=%.2f"%onshellmuon_phi)
        tex11.SetTextSize(0.08)
        tex11.SetTextFont(42)
        tex11.SetNDC()
	

        #pt of nu_onshellW, nu_onshellW_pt,
	#c1.cd(2)
        todraw2 = "nu_onshellW_pt"
        todrawtrue2 = "pt_nuonshellW_true"
	nupt_onshell_bins = "(60,0,120)"
	hs_title2 = "p_{T#nu}^{onshellW}"
	hs2 = gethiststack(tree, todraw2, todrawtrue2, nupt_onshell_bins, hs_title2, cut)
	#hs2.SetName("%s"%tree.GetTitle()+"_%s"%todraw2)
     	#hs2.Draw("nostack")
        onshellmuon_pt =getTrueValue(tree,"mu_onshellW_pt","(400,0,200)")
	tex2 = ROOT.TLatex(0.4,0.7,"p_{T#mu}^{onshell}=%.2f"%onshellmuon_pt)
        tex2.SetTextSize(0.08)
        tex2.SetTextFont(42)
        tex2.SetNDC()
        

        #mass of onshellW, onshellW_Mass, 
	#c1.cd(3)
        todraw3 = "onshellW_Mass"
        todrawtrue3 = "mass_onshellW_true"
	wmass_onshell_bins = "(50,45,95)"
	hs_title3 = "M_{W}^{onshell}"
	hs3 = gethiststack(tree, todraw3, todrawtrue3, wmass_onshell_bins, hs_title3, cut)
	#hs3.SetName("%s"%tree.GetTitle()+"_%s"%todraw3)
     	#hs3.Draw("nostack")


        #eta-phi of nu_offshellW
	#c1.cd(4)
	todraw4 = "nu_offshellW_phi:nu_offshellW_eta"
        todrawtrue4 = "phi_nuoffshellW_true:eta_nuoffshellW_true"
        hs_title4 = "#phi_{#nu}^{offshellW} Vs #eta_{#nu}^{offshellW}"
        hs_Xtitle4 = "#eta_{#nu}^{offshellW}"
        hs_Ytitle4 = "#phi_{#nu}^{offshellW}"
        hs4 = gethiststack_2d(tree, todraw4, todrawtrue4, eta_bins, phi_bins, hs_title4, cut)
        hist4_true = hist_2D(tree, todrawtrue4, eta_bins, phi_bins, cut,1)
        hist4_true.SetMarkerStyle(7) 
        hist4_true.SetMarkerSize(1)
        hist4_true.SetMarkerColor(ROOT.kMagenta)
	offshellmuon_eta = getTrueValue(tree,"mu_offshellW_eta","(200,-6,6)")
	offshellmuon_phi = getTrueValue(tree,"mu_offshellW_phi","(200,-3.1415,3.1415)")
	array_mu2_eta = array('d',[offshellmuon_eta])
	array_mu2_phi = array('d',[offshellmuon_phi])
	graph4_muon = ROOT.TGraph(1, array_mu2_eta,array_mu2_phi)
        graph4_muon.SetMarkerColor(ROOT.kBlack)
 	graph4_muon.SetMarkerStyle(16)
	graph4_muon.SetMarkerSize(1)
	
	#eta of nuoffshellW	
        todraw12 = "nu_offshellW_eta"
	todrawtrue12 = "eta_nuoffshellW_true"
        hs_title12 = "#eta_{#nu}^{offshellW}"
	hs12 = gethiststack(tree, todraw12, todrawtrue12, eta_bins, hs_title12, cut)
	tex12 = ROOT.TLatex(0.4,0.7,"#eta_{#mu}^{offshell}=%.2f"%offshellmuon_eta)
        tex12.SetTextSize(0.08)
        tex12.SetTextFont(42)
        tex12.SetNDC()
	
	#phi of nuoffshellW	
        todraw13 = "nu_offshellW_phi"
	todrawtrue13 = "phi_nuoffshellW_true"
        hs_title13 = "#phi_{#nu}^{offshellW}"
	hs13 = gethiststack(tree, todraw13, todrawtrue13, phi_bins, hs_title13, cut)
	tex13 = ROOT.TLatex(0.4,0.7,"#phi_{#mu}^{offshell}=%.2f"%offshellmuon_phi)
        tex13.SetTextSize(0.08)
        tex13.SetTextFont(42)
        tex13.SetNDC()
        	

        #pt of nu_offshellW
	#c1.cd(5)
        todraw5 = "nu_offshellW_pt"
        todrawtrue5 = "pt_nuoffshellW_true"
	nupt_offshell_bins = "(40,0,80)"
	hs_title5 = "p_{T#nu}^{offshellW}"
	hs5 = gethiststack(tree, todraw5, todrawtrue5, nupt_offshell_bins, hs_title5, cut)
	#hs5.SetName("%s"%tree.GetTitle()+"_%s"%todraw5)
     	#hs5.Draw("nostack")
        offshellmuon_pt =getTrueValue(tree,"mu_offshellW_pt","(400,0,200)")
	tex5 = ROOT.TLatex(0.4,0.7,"p_{T#mu}^{offshell}=%.2f"%offshellmuon_pt)
        tex5.SetTextSize(0.08)
        tex5.SetTextFont(42)
        tex5.SetNDC()
        


	#mass of offshellW
	#c1.cd(6)
        todraw6 = "offshellW_Mass"
	todrawtrue6 = "mass_offshellW_true"
	wmass_offshell_bins = "(60,0,60)"
	hs_title6 = "M_{W}^{offshell}"
	hs6 = gethiststack(tree, todraw6, todrawtrue6, wmass_offshell_bins, hs_title6, cut)
	hs6.SetName("%s"%tree.GetTitle()+"_%s"%todraw6)
     	#hs6.Draw("nostack")


	#mass of h, htoWW_Mass
	#c1.cd(7)
	hmass_bins = "(20,124,126)"
        todraw7 = "htoWW_Mass"
        todrawtrue7 = "mass_htoWW_true"
	hs_title7 = "M_{h#rightarrow WW}"
	hs7 = gethiststack(tree, todraw7,todrawtrue7, hmass_bins, hs_title7, cut)
	#hs7.SetName("%s"%tree.GetTitle()+"_%s"%todraw7)
	#hs7.Draw("nostack")
      
        #pt of h2
        #c1.cd(8)
        pt_h2_bins = "(100,0,0.00001)"
	todraw8 = "h2tohh_Pt"
	todrawtrue8 = "pt_h2_true"
	hs_title8 = "p_{TH}"
	hs8 = gethiststack(tree, todraw8,todrawtrue8, pt_h2_bins, hs_title8, cut)
	#hs8.SetName("%s"%tree.GetTitle()+"_%s"%todraw8)
        #hs8.Draw("nostack")
 	#mass of h2,h2tohh_Mass
	#c1.cd(9)
        #todraw
        h2mass_bins = "(80,100,900)"
	todraw9 = "h2tohh_Mass"
        todrawtrue9 = "mass_h2_true"
	hs_title9 = "M_{H}"
	hs9 = gethiststack(tree, todraw9,todrawtrue9, h2mass_bins, hs_title9, cut)
#	hs14 = gethiststack(tree, todraw9,todrawtrue9, h2mass_bins, hs_title9,cut_weight1)
       # hs9.GetHistogram().GetXaxis().SetTitleSize(0.06)
	#hs9.SetName("%s"%tree.GetTitle()+"_%s"%todraw9)
        
        todraw14 = "onshellW_Mass:h2tohh_Mass"
        todrawtrue14 = "mass_onshellW_true:mass_h2_true"
        hs_title14 = "M_{W}^{onshell} Vs M_{H}"
        hs_Xtitle14 = "M_{H}"
        hs_Ytitle14 = "M_{W}^{onshell}"
        hs14 = gethiststack_2d(tree, todraw14, todrawtrue14, h2mass_bins, wmass_onshell_bins, hs_title14, cut)
        hist14_true = hist_2D(tree, todrawtrue14, h2mass_bins, wmass_onshell_bins, cut,1)
        hist14_true.SetMarkerStyle(7) 
        hist14_true.SetMarkerSize(1)
        hist14_true.SetMarkerColor(ROOT.kMagenta)
        
        
        todraw15 = "offshellW_Mass:h2tohh_Mass"
        todrawtrue15 = "mass_offshellW_true:mass_h2_true"
        hs_title15 = "M_{W}^{offshell} Vs M_{H}"
        hs_Xtitle15 = "M_{H}"
        hs_Ytitle15 = "M_{W}^{offshell}"
        hs15 = gethiststack_2d(tree, todraw15, todrawtrue15, h2mass_bins, wmass_offshell_bins, hs_title15, cut)
        hist15_true = hist_2D(tree, todrawtrue15, h2mass_bins, wmass_offshell_bins, cut,1)
        hist15_true.SetMarkerStyle(7) 
        hist15_true.SetMarkerSize(1)
        hist15_true.SetMarkerColor(ROOT.kMagenta)
        
        
        todraw16 = "nu_onshellW_pt:h2tohh_Mass"
        todrawtrue16 = "pt_nuonshellW_true:mass_h2_true"
        hs_title16 = "p_{T#nu}^{onshell} Vs M_{H}"
        hs_Xtitle16 = "M_{H}"
        hs_Ytitle16 = "p_{T#nu}^{onshell}"
        hs16 = gethiststack_2d(tree, todraw16, todrawtrue16, h2mass_bins, nupt_onshell_bins, hs_title16, cut)
        hist16_true = hist_2D(tree, todrawtrue16, h2mass_bins, nupt_onshell_bins, cut,1)
        hist16_true.SetMarkerStyle(7) 
        hist16_true.SetMarkerSize(1)
        hist16_true.SetMarkerColor(ROOT.kMagenta)
        
        
        todraw17 = "nu_offshellW_pt:h2tohh_Mass"
        todrawtrue17 = "pt_nuoffshellW_true:mass_h2_true"
        hs_title17 = "p_{T#nu}^{offshell} Vs M_{H}"
        hs_Xtitle17 = "M_{H}"
        hs_Ytitle17 = "p_{T#nu}^{offshell}"
        hs17 = gethiststack_2d(tree, todraw17, todrawtrue17, h2mass_bins, nupt_offshell_bins, hs_title17, cut)
        hist17_true = hist_2D(tree, todrawtrue17, h2mass_bins, nupt_offshell_bins, cut,1)
        hist17_true.SetMarkerStyle(7) 
        hist17_true.SetMarkerSize(1)
        hist17_true.SetMarkerColor(ROOT.kMagenta)
        

	#tot
	c1.Clear()
	c1.Divide(4,4)
	c1.cd(1)
     	hs1.Draw("nostackcolz")
        hs1.GetHistogram().GetXaxis().SetTitle("%s"%hs_Xtitle1)
        hs1.GetHistogram().GetYaxis().SetTitle("%s"%hs_Ytitle1)
        hs1.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs1.GetHistogram().GetXaxis().SetTitleSize(0.08)
        hist1_true.Draw("same")
	#hist1_2.Draw("CONT3")
	#hist1_3.Draw("CONT3same")
	#graph1.Draw("psame")
	graph1_muon.Draw("psame")

	c1.cd(2)
     	hs10.Draw("nostack")
        hs10.GetHistogram().GetXaxis().SetTitle("%s"%hs_title10)
        hs10.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs10.GetHistogram().GetXaxis().SetTitleSize(0.08)
        tex10.Draw("same")

	c1.cd(3)
     	hs11.Draw("nostack")
        hs11.GetHistogram().GetXaxis().SetTitle("%s"%hs_title11)
        hs11.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs11.GetHistogram().GetXaxis().SetTitleSize(0.08)
        tex11.Draw("same")

	c1.cd(4)
	hs2.Draw("nostack")
        hs2.GetHistogram().GetXaxis().SetTitle("%s"%hs_title2)
        hs2.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs2.GetHistogram().GetXaxis().SetTitleSize(0.08)
        tex2.Draw("same")

	c1.cd(5)
	#hist4_2.Draw("CONT3")
	#hist4_3.Draw("CONT3same")
	#graph4.Draw("psame")
     	hs4.Draw("nostackcolz")
        hs4.GetHistogram().GetXaxis().SetTitle("%s"%hs_Xtitle4)
        hs4.GetHistogram().GetYaxis().SetTitle("%s"%hs_Ytitle4)
        hs4.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs4.GetHistogram().GetXaxis().SetTitleSize(0.08)
        hist4_true.Draw("same")
	graph4_muon.Draw("psame")

	c1.cd(6)
     	hs12.Draw("nostack")
        hs12.GetHistogram().GetXaxis().SetTitle("%s"%hs_title12)
        hs12.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs12.GetHistogram().GetXaxis().SetTitleSize(0.08)
	tex12.Draw("same")

	c1.cd(7)
     	hs13.Draw("nostack")
        hs13.GetHistogram().GetXaxis().SetTitle("%s"%hs_title13)
        hs13.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs13.GetHistogram().GetXaxis().SetTitleSize(0.08)
	tex13.Draw("same")

	c1.cd(8)
     	hs5.Draw("nostack")
        hs5.GetHistogram().GetXaxis().SetTitle("%s"%hs_title5)
        hs5.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs5.GetHistogram().GetXaxis().SetTitleSize(0.08)
	tex5.Draw("same")

	c1.cd(9)
     	hs3.Draw("nostack")
        hs3.GetHistogram().GetXaxis().SetTitle("%s"%hs_title3)
        hs3.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs3.GetHistogram().GetXaxis().SetTitleSize(0.08)

	c1.cd(10)
     	hs6.Draw("nostack")
        hs6.GetHistogram().GetXaxis().SetTitle("%s"%hs_title6)
        hs6.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs6.GetHistogram().GetXaxis().SetTitleSize(0.08)

        c1.cd(11) 
     	hs7.Draw("nostack")
        hs7.GetHistogram().GetXaxis().SetTitle("%s"%hs_title7)
        hs7.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs7.GetHistogram().GetXaxis().SetTitleSize(0.08)

        c1.cd(12)
     	hs9.Draw("nostack")
        #hs10.Draw("nostacksame")
        hs9.GetHistogram().GetXaxis().SetTitle("%s"%hs_title9)
        hs9.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs9.GetHistogram().GetXaxis().SetTitleSize(0.08)
       
        c1.cd(13)
     	hs14.Draw("nostackcolz")
        hs14.GetHistogram().GetXaxis().SetTitle("%s"%hs_Xtitle14)
        hs14.GetHistogram().GetYaxis().SetTitle("%s"%hs_Ytitle14)
        hs14.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs14.GetHistogram().GetXaxis().SetTitleSize(0.08)
        hist14_true.Draw("same")

        c1.cd(14)
     	hs15.Draw("nostackcolz")
        hs15.GetHistogram().GetXaxis().SetTitle("%s"%hs_Xtitle15)
        hs15.GetHistogram().GetYaxis().SetTitle("%s"%hs_Ytitle15)
        hs15.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs15.GetHistogram().GetXaxis().SetTitleSize(0.08)
        hist15_true.Draw("same")

        c1.cd(15)
     	hs16.Draw("nostackcolz")
        hs16.GetHistogram().GetXaxis().SetTitle("%s"%hs_Xtitle16)
        hs16.GetHistogram().GetYaxis().SetTitle("%s"%hs_Ytitle16)
        hs16.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs16.GetHistogram().GetXaxis().SetTitleSize(0.08)
        hist16_true.Draw("same")

        c1.cd(16)
     	hs17.Draw("nostackcolz")
        hs17.GetHistogram().GetXaxis().SetTitle("%s"%hs_Xtitle17)
        hs17.GetHistogram().GetYaxis().SetTitle("%s"%hs_Ytitle17)
        hs17.GetHistogram().GetXaxis().SetTitleOffset(0.6)
        hs17.GetHistogram().GetXaxis().SetTitleSize(0.08)
        hist17_true.Draw("same")

	c1.Update()

        c1.cd()
        met = getTrueValue(tree,"met_true","(400,0,200)")
    #    print "met ", met
    	tex = ROOT.TLatex(0.15,.958,"Performance of MMC,#slash{E}_{T}=%.2f"%met+", %s"%(tree.GetTitle()))
    	#tex = ROOT.TLatex(0.25,.954,"Performance of MMC"+", %s"%(tree.GetTitle()))
    	tex.SetTextSize(0.04)
    	tex.SetTextFont(42)
    	tex.SetNDC()
	tex.Draw("same")
	c1.Update()
        #c1.SaveAs("%s"%(tree.GetTitle())+"_%s_%s_B3.pdf"%(pic, cut[0]))
	#c1.SaveAs("%s"%(tree.GetTitle())+"_%s_%s_B3.png"%(pic, cut[0]))
    #    c1.SaveAs("%s"%(tree.GetTitle())+"_%s_%s_B3.C"%(pic, cut[0]))
        ##########################################
        ##### fill hist_h2
        
	if m>200:
		break
	m = m+1
        f.Close()


#____________________________________________________________________________
def drawh2Mass_combined(file, cut="weight"):


    num = gettreenum(file)
    print "num of trees ", num
    hist_h2 = ROOT.TH1F("hist_h2"," ",150,200,500)
    hist_h2true = ROOT.TH1F("hist_h2true"," ",150,200,500)
    hist_h2Mass_RMS = ROOT.TH2F("hist_h2Mass_RMS"," ",150,200,500,50,0,50)
    hist_h2Mass_RMS2 = ROOT.TH2F("hist_h2Mass_RMS2"," ",150,200,500,50,0,0.5)
    hist_h2Mass_deltaR = ROOT.TH2F("hist_h2Mass_deltaR"," ",150,200,500,60,0,6)
    hist_h2Mass_mean = ROOT.TH2F("hist_h2Mass_mean"," ",150,200,500,150,200,500)
    hist_h2Mass_integral = ROOT.TH2F("hist_h2Mass_integral"," ",150,200,500,100,0,totintegral[cut]/count)
    hist_h2Mass_integral_peak5 = ROOT.TH2F("hist_h2Mass_integral_peak5"," ",150,200,500,75,0,totintegral[cut]/count*0.75)
    hist_h2Mass_integral_peak10 = ROOT.TH2F("hist_h2Mass_integral_peak10"," ",150,200,500,75,0,totintegral[cut]/count*0.75)
    hist_rms_integral_peak10 = ROOT.TH2F("hist_rms_integral_peak10"," ",50,0,50,75,0,totintegral[cut]/count*0.75)
    hist_h2Mass_met = ROOT.TH2F("hist_h2Mass_met"," ",150,200,500,50,0,200)
    hist_h2Mass_muonshellpt = ROOT.TH2F("hist_h2Mass_nuoffshellpt"," ",150,200,500,80,0,160)
    hist_h2Mass_muoffshellpt = ROOT.TH2F("hist_h2Mass_nuonshellpt"," ",150,200,500,80,0,160)
    hist_RMS_correct = ROOT.TH1F("hist_RMS_correct"," ",50,0,50)
    hist_RMS_incorrect = ROOT.TH1F("hist_RMS_incorrect"," ",50,0,50)
    hist_integralpeak10_correct = ROOT.TH1F("hist_integralpeak10_correct"," ",75,0,totintegral[cut]/count*0.75)
    hist_integralpeak10_incorrect = ROOT.TH1F("hist_integralpeak10_incorrect"," ",75,0,totintegral[cut]/count*0.75)
    hist_h2.SetXTitle("M_{H}")
    hist_h2Mass_RMS.SetXTitle("reconstruced M_{H}")
    hist_h2Mass_RMS.SetYTitle("MMC RMS")
    hist_h2Mass_RMS2.SetXTitle("reconstructed M_{H}")
    hist_h2Mass_RMS2.SetYTitle("MMC RMS/M_{H}")
    hist_h2Mass_deltaR.SetXTitle("reconstructed M_{H}")
    hist_h2Mass_deltaR.SetYTitle("#DeltaR_{#nu1,#nu2}")
    hist_h2Mass_mean.SetXTitle("reconstructed M_{H}")
    hist_h2Mass_mean.SetYTitle("average M_{H}")
    hist_h2Mass_integral.SetXTitle("reconstructed M_{H}")
    hist_h2Mass_integral.SetYTitle("integral over M_{H}")
    hist_h2Mass_integral_peak5.SetXTitle("reconstructed M_{H}")
    hist_h2Mass_integral_peak5.SetYTitle("integral over (maxbin-5, maxbin+5)")
    hist_h2Mass_integral_peak10.SetXTitle("reconstructed M_{H}")
    hist_h2Mass_integral_peak10.SetYTitle("integral over (maxbin-10, maxbin+10)")
    hist_rms_integral_peak10.SetXTitle("RMS")
    hist_rms_integral_peak10.SetYTitle("integral over (maxbin-10, maxb+10)")
    hist_h2Mass_met.SetXTitle("reconstructed M_{H}")
    hist_h2Mass_met.SetYTitle("#slash{E}_{T}")
    hist_h2Mass_muonshellpt.SetXTitle("reconstructed M_{H}")
    hist_h2Mass_muonshellpt.SetYTitle("p_{T#mu}^{onshell}")
    hist_h2Mass_muoffshellpt.SetXTitle("reconstructed M_{H}")
    hist_h2Mass_muoffshellpt.SetYTitle("p_{T#mu}^{offshell}")
    m = 1
    #test 
    #num = 10
    #print "totintegral ", totintegral[cut]
    while m<num:
    	f = ROOT.TFile(file)
    	list = f.GetListOfKeys()
    	key = list.At(0)
    	obj = key.ReadObj() #DiHiggsWWAna
    	sub_list = obj.GetListOfKeys()
    	sub_key = sub_list.At(m)
	tree = sub_key.ReadObj()
        h2_mass_true = getTrueValue(tree,"mass_h2_true","(800,200,1000)")
    	name = tree.GetTitle()+"_h2Mass"
   	b1 = ROOT.TH1F("%s"%name,"b1",800,200,1000)
    	tree.Draw("h2tohh_Mass>>%s"%name,cut)
        h2_mass_mean = b1.GetMean()
        h2_mass_integral = b1.Integral(0,800+1)/count
        maxbin = b1.GetMaximumBin() 
	integral_peak5 = b1.Integral(maxbin-5, maxbin+5)/count
	integral_peak10 = b1.Integral(maxbin-10, maxbin+10)/count
        h2_mass_reconstructed = geth2MassMostProba(tree,cut)
        h2_mass_reco_rms = getRMSh2Mass(tree,cut)
	nu1_eta = getTrueValue(tree,"eta_nuonshellW_true","(200,-6,6)")
	nu1_phi = getTrueValue(tree,"phi_nuonshellW_true","(200,-3.1415,3.1415)")
	nu2_eta = getTrueValue(tree,"eta_nuoffshellW_true","(200,-6,6)")
	nu2_phi = getTrueValue(tree,"phi_nuoffshellW_true","(200,-3.1415,3.1415)")
        met = getTrueValue(tree,"met_true","(400,0,200)")
        onshellmuon_pt =getTrueValue(tree,"mu_onshellW_pt","(400,0,200)")
        offshellmuon_pt =getTrueValue(tree,"mu_offshellW_pt","(400,0,200)")
        deltaphi = math.fabs(nu1_phi-nu2_phi)
        if deltaphi > math.pi :
		deltaphi = 2*math.pi - deltaphi
        deltaeta = nu1_eta-nu2_eta
        deltaR = math.sqrt(deltaphi*deltaphi + deltaeta*deltaeta)
#	print  " tree title ",tree.GetTitle()," entries ", tree.GetEntries()," rms ",h2_mass_reco_rms 
    #    if m>10:
#		break;    
	if (abs(h2_mass_reconstructed-h2_mass_true)<20):
	#	print "m ",m," ",tree.GetTitle()," reconstructed mass  ", h2_mass_reconstructed, " true mass ", h2_mass_true
		hist_integralpeak10_correct.Fill(integral_peak10)
		hist_RMS_correct.Fill(h2_mass_reco_rms)
	else:
		hist_integralpeak10_incorrect.Fill(integral_peak10)
		hist_RMS_incorrect.Fill(h2_mass_reco_rms)
		
 	hist_h2Mass_RMS.Fill(h2_mass_reconstructed, h2_mass_reco_rms)	
 	hist_h2Mass_RMS2.Fill(h2_mass_reconstructed, h2_mass_reco_rms/h2_mass_reconstructed)	
        hist_h2Mass_deltaR.Fill(h2_mass_reconstructed, deltaR)
        hist_h2Mass_mean.Fill(h2_mass_reconstructed, h2_mass_mean)
        hist_h2Mass_integral.Fill(h2_mass_reconstructed, h2_mass_integral)
        hist_h2Mass_integral_peak5.Fill(h2_mass_reconstructed, integral_peak5)
        hist_h2Mass_integral_peak10.Fill(h2_mass_reconstructed, integral_peak10)
        hist_rms_integral_peak10.Fill(h2_mass_reco_rms, integral_peak10)
        hist_h2Mass_met.Fill(h2_mass_reconstructed, met)
        hist_h2Mass_muonshellpt.Fill(h2_mass_reconstructed, onshellmuon_pt)
        hist_h2Mass_muoffshellpt.Fill(h2_mass_reconstructed, offshellmuon_pt)
       	hist_h2.Fill(h2_mass_reconstructed)
	hist_h2true.Fill(h2_mass_true)
        m = m+1
	#sub_key = sub_list.At(m)
	f.Close()
       
    hist_integralpeak10_correct.SetXTitle("intergral(M_{H}-20GeV, M_{H}+20GeV)")
    hist_RMS_correct.SetXTitle("RMS ")
    hist_integralpeak10_incorrect.SetXTitle("intergral(M_{H}-20GeV, M_{H}+20GeV)")
    hist_RMS_incorrect.SetXTitle("RMS ")
    hist_integralpeak10_correct.SetLineColor(ROOT.kRed)
    hist_RMS_correct.SetLineColor(ROOT.kRed)
    hist_integralpeak10_incorrect.SetLineColor(ROOT.kBlue)
    hist_RMS_incorrect.SetLineColor(ROOT.kBlue)
    leg_integralpeak10 = ROOT.TLegend(0.25,0.7,0.45,0.82)
    leg_integralpeak10.AddEntry(hist_integralpeak10_correct,"Good reconstruction'")
    leg_integralpeak10.AddEntry(hist_integralpeak10_incorrect,"bad reconstruction'")
    leg_RMS = ROOT.TLegend(0.25,0.7,0.45,0.82)
    leg_RMS.AddEntry(hist_RMS_correct,"Good reconstruction'")
    leg_RMS.AddEntry(hist_RMS_incorrect,"bad reconstruction'")

    h2massout = ROOT.TFile("h2mass_%s"%cut+"out.root","recreate")
    h2massout.cd()
    h2Mass_c = ROOT.TCanvas()
    hist_h2.SetLineColor(ROOT.kRed)
    hist_h2.GetXaxis().SetTitle("M_{H}")
    hist_h2.SetTitle("M_{H} reconstruction from MMC")

    hist_h2true.SetLineColor(ROOT.kBlue)
    hist_h2true.SetLineStyle(2)
    legend = ROOT.TLegend(0.25,0.7,0.45,0.82)
    legend.SetFillColor(ROOT.kWhite)
    #legend.SetFillStyle(0)
    legend.AddEntry(hist_h2,"Reconstructed ","l") 
    legend.AddEntry(hist_h2true,"True ","l") 
    hist_h2.Draw()
    hist_h2true.Draw("same")
    legend.Draw("same")
    h2Mass_c.SaveAs("MMC_h2Mass_0511_%s_1M_B3.pdf"%cut) 
    h2Mass_c.SaveAs("MMC_h2Mass_0511_%s_1M_B3.png"%cut) 
    h2Mass_c.SaveAs("MMC_h2Mass_0511_%s_1M_B3.C"%cut) 
    h2Mass_c.Write()
 
    h2_c2 = ROOT.TCanvas("h2_c2","h2_c2",1200,900)
    h2_c2.cd()
    h2_c2.Divide(4,3)
    h2_c2.cd(1)
    hist_h2.Draw()
    hist_h2true.Draw("same")
    legend.Draw("same")
    h2_c2.cd(2) 
    hist_h2Mass_RMS.Draw("colz")
    h2_c2.cd(3)
    hist_h2Mass_RMS2.Draw("colz")
    h2_c2.cd(4)
    hist_h2Mass_deltaR.Draw("colz")
    h2_c2.cd(5)
    hist_h2Mass_mean.Draw("colz")
    h2_c2.cd(6)
    hist_h2Mass_integral.Draw("colz")
    h2_c2.cd(7)
    hist_h2Mass_integral_peak5.Draw("colz")
    h2_c2.cd(8)
    hist_h2Mass_integral_peak10.Draw("colz")
    h2_c2.cd(9)
    hist_rms_integral_peak10.Draw("colz")
    h2_c2.cd(10)
    hist_RMS_correct.Draw()
    hist_RMS_incorrect.Draw("same")
    leg_RMS.Draw("same")
    
    h2_c2.cd(11)
    hist_integralpeak10_correct.Draw()
    hist_integralpeak10_incorrect.Draw("same")
    leg_integralpeak10.Draw("same")

    h2_c2.cd(12)
    hist_h2Mass_met.Draw("colz")
    #hist_h2Mass_muoffshellpt.Draw("colz")
    #hist_h2Mass_muonshellpt.Draw("colz")

    h2_c2.cd()
    h2_c2.SaveAs("MMC_h2Mass_combined_0511_%s_1M_B3.pdf"%cut)
    h2_c2.SaveAs("MMC_h2Mass_combined_0511_%s_1M_B3.png"%cut)
    h2_c2.SaveAs("MMC_h2Mass_combined_0511_%s_1M_B3.C"%cut)
    h2_c2.Write() 
    hist_h2.Write()
    hist_h2true.Write()
    hist_h2Mass_RMS.Write()
    hist_h2Mass_RMS2.Write()
    hist_h2Mass_deltaR.Write()
    hist_h2Mass_mean.Write()
    hist_h2Mass_integral.Write()
    hist_h2Mass_integral_peak5.Write()
    hist_h2Mass_integral_peak10.Write()
    hist_rms_integral_peak10.Write()
    hist_h2Mass_met.Write()
    hist_integralpeak10_correct.Write()
    hist_RMS_correct.Write()
    hist_integralpeak10_incorrect.Write()
    hist_RMS_incorrect.Write()
    hist_h2Mass_muonshellpt.Write()
    hist_h2Mass_muoffshellpt.Write()
    h2massout.Write()
    h2massout.Close()


#____________________________________________________________________________
def drawh2MassAll_combined(dir, pic, cut="weight"):


    #ROOT.gStyle.SetOptStat(0)

    hist_h2 = ROOT.TH1F("hist_h2"," ",150,200,500)
    hist_h2true = ROOT.TH1F("hist_h2true"," ",150,200,500)
    hist_h2ratio_metdiff = ROOT.TH2F("hist_h2ratio_metdiff"," ",40,0.6,1.4,50,0,30)
    hist_h2_Metdiff = ROOT.TH2F("hist_h2_Metdiff"," ",75,200,500,50,0,30)
    hist_h2_Metdiff.SetXTitle("reconstruced M_{H}")
    hist_h2_Metdiff.SetYTitle("|#slash{E}_{T}-#vec{nu1}-#vec{nu2}|")
    hist_Metdiffpx = ROOT.TH1F("hist_Metdiffpx", " ",100,-100,100)

    if not os.path.isdir(dir):
          print "ERROR: This is not a valid directory: ", dir
    ls = os.listdir(dir)
    tot = len(ls)
    met_diff_Vec = ROOT.TVector2()
    nu1_Vec = ROOT.TVector2()
    nu2_Vec = ROOT.TVector2()
    pi_root = ROOT.TMath.Pi()
    for x in ls:
	x = dir[:]+x
	print "rootfile X ",x
    	num = gettreenum(x)
	m = 1
    	while m<num:
    		f = ROOT.TFile(x)
    		list = f.GetListOfKeys()
    		key = list.At(0)
   	 	obj = key.ReadObj() #DiHiggsWWAna
    		sub_list = obj.GetListOfKeys()
    		sub_key = sub_list.At(m)
		tree = sub_key.ReadObj()
		print tree.GetTitle()
		if tree.GetEntries()<100:
			m = m+1
			f.Close()
			continue
			
        	h2_mass_true = getTrueValue(tree,"mass_h2_true","(900,100,1000)")
		met = getTrueValue(tree,"met_true","(100,0,300)")
		nu1_pt = getTrueValue(tree,"pt_nuoffshellW_true","(400,0,200)")
		nu2_pt = getTrueValue(tree,"pt_nuonshellW_true","(400,0,200)")
		nu1_phi = getTrueValue(tree,"phi_nuoffshellW_true","(100,%f,%f)"%(-1*pi_root,pi_root))
		nu2_phi = getTrueValue(tree,"phi_nuonshellW_true","(100,%f,%f)"%(-1*pi_root,pi_root))
		met_phi = getTrueValue(tree,"met_phi_true","(100,%f,%f)"%(0,2*pi_root))
		if met_phi>pi_root:
			met_phi = met_phi-2*pi_root
		#print "met ", met, " phi ", met_phi," nu1_pt ", nu1_pt, " nu1_phi ", nu1_phi," nu2_pt ",nu2_pt, " nu2_phi ",nu2_phi
			
		met_diff_Vec.SetMagPhi(met, met_phi)
		nu1_Vec.SetMagPhi(nu1_pt, nu1_phi)
		nu2_Vec.SetMagPhi(nu2_pt, nu2_phi)
		met_diff_Vec = met_diff_Vec-nu1_Vec-nu2_Vec
		#print "met ", met_diff_Vec.Print()
		#print "nu1 ", nu1_Vec.Print()
		#print "nu2 ", nu2_Vec.Print()
		#print "met diff ", met_diff_Vec.Print()
    		name = tree.GetTitle()+"_h2Mass"
   		b1 = ROOT.TH1F("%s"%name,"b1",800,200,1000)
    		tree.Draw("h2tohh_Mass>>%s"%name,cut)
		
        	h2_mass_mean = b1.GetMean()
        	maxbin = b1.GetMaximumBin() 
		h2_mass_reconstructed = b1.GetXaxis().GetBinCenter(maxbin)
			
		print name," ",h2_mass_reconstructed,"  true mass ", h2_mass_true	
       		hist_h2.Fill(h2_mass_reconstructed)
		hist_h2true.Fill(h2_mass_true)
		#hist_h2ratio_metdiff.Fill(h2_mass_reconstructed/h2_mass_true, met_diff_Vec.Mod())
		hist_h2_Metdiff.Fill(h2_mass_reconstructed, met_diff_Vec.Mod())
		hist_Metdiffpx.Fill(met_diff_Vec.Px())
        	m = m+1
		f.Close()
       
			
    h2Mass_c = ROOT.TCanvas()
    hist_h2.SetLineColor(ROOT.kRed)
    hist_h2.GetXaxis().SetTitle("M_{H}")
    hist_h2.SetTitle("M_{H} reconstruction from MMC")

    hist_h2true.SetLineColor(ROOT.kBlue)
    hist_h2true.SetLineStyle(2)
    legend = ROOT.TLegend(0.25,0.7,0.45,0.82)
    legend.SetFillColor(ROOT.kWhite)
    #legend.SetFillStyle(0)
    legend.AddEntry(hist_h2,"Reconstructed ","l") 
    legend.AddEntry(hist_h2true,"True ","l") 
    hist_h2.GetYaxis().SetRangeUser(0,hist_h2true.GetMaximum()+10)
    hist_h2.Draw()
    hist_h2true.Draw("same")
    legend.Draw("same")
    h2Mass_c.SaveAs("MMC_h2Mass_%s_B3_%s.pdf"%(cut,pic)) 
    h2Mass_c.SaveAs("MMC_h2Mass_%s_B3_%s.png"%(cut,pic)) 

    h2_c3 = ROOT.TCanvas()
    h2_c3.cd()
    hist_h2_Metdiff.Draw("colz")
    h2_c3.SaveAs("MMC_h2_Metdiff_%s_B3_%s.pdf"%(cut,pic))
    h2_c3.SaveAs("MMC_h2_Metdiff_%s_B3_%s.png"%(cut,pic))
    """
    h2_c2 = ROOT.TCanvas()
    h2_c2.cd()
    hist_h2ratio_metdiff.SetXTitle("#frac{reconstruced M_{H}}{true M_{H}}")
    hist_h2ratio_metdiff.SetYTitle("| #slash{E}_{T}-#vec{nu1}-#vec{nu2} |")
    hist_h2ratio_metdiff.Draw("colz")
    h2_c2.SaveAs("MMC_h2ratio_Metdiff_%s_B3_%s.pdf"%(cut,pic))
    h2_c2.SaveAs("MMC_h2ratio_Metdiff_%s_B3_%s.png"%(cut,pic))

    h2_c4 = ROOT.TCanvas()
    h2_c4.cd()
    hist_Metdiffpx.Draw()
    h2_c4.SaveAs("MMC_h2_MetdiffPx_%s_B3_%s.pdf"%(cut,pic))
    h2_c4.SaveAs("MMC_h2_MetdiffPx_%s_B3_%s.png"%(cut,pic))
    """
#___________________________________________________________________________
def geth2MassMostProba(t, cut = "weight"):
 
    name = t.GetTitle()
    b1 = ROOT.TH1F("%s"%name,"b1",800,200,1000)
    t.Draw("h2tohh_Mass>>%s"%name,cut)
    #print "maxbin ",b1.GetMaximumBin()

    mass = b1.GetXaxis().GetBinCenter(b1.GetMaximumBin())
    del b1
    return mass

#___________________________________________________________________________
def getRMSh2Mass(t, cut = "weight"):
 
    name = t.GetTitle()
    b1 = ROOT.TH1F("%s"%name,"b1",800,200,1000)
    t.Draw("h2tohh_Mass>>%s"%name,cut)
    #print "maxbin ",b1.GetMaximumBin()
    return b1.GetRMS()

#____________________________________________________________________________
def getTrueValue(t,var, x_bins, cut="control<2"):
   
    c = ROOT.TCanvas() 
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    name = "tmp_%s"%var
    b1 = ROOT.TH1F("%s"%name,"b1",xBins,xminBin,xmaxBin)
    t.Draw("%s"%var+">>%s"%name, cut)
    value = b1.GetMean()
    if value<xminBin or value>xmaxBin:
	print " getTrue value is not correct one ? ",var," value ", value
    	b1.Draw()
        c.SaveAs("getTrueVaule_%s"%var+"_test.png")
    #value = b1.GetXaxis().GetBinCenter(b1.GetMaximumBin())
    #print "xaxis", var, "b1 maximum", b1.GetMaximumBin()," value", value
    del b1
    return value
    

#_____________________________________________________________________________
def test(file, dir):
 #subdir
    f = ROOT.TFile(file)
    t = f.Get(dir)
    name = f.GetName()
    List = f.GetListOfKeys()
    print "at 0 ", List.At(0)
    key = List.At(0)
    c1 = ROOT.gROOT.GetClass(key.GetClassName())
    obj = key.ReadObj()
    sub_list = obj.GetListOfKeys()
    print "length of sub_list", len(sub_list)
    sub_key = sub_list.At(0)
    i = 0
    while sub_key:
	tree = sub_key.ReadObj()
#	tree.Print()
        i = i+1
        sub_key = sub_list.At(i)
    print "i ", i
    print "key getclassname ", key.GetClassName()
    print "c1 ",c1
    print "obj  ", obj
    print "sub_key ", sub_key
    print "sub_list 1", sub_list.At(1)
    print "dir ", dir, "  t ", t, " tree ", tree.GetTitle() 

#______________________________________________________________________________
def gettreenum(file):
	
    f = ROOT.TFile(file)
    List = f.GetListOfKeys()
    key = List.At(0)
    obj = key.ReadObj()
    sub_list = obj.GetListOfKeys()
    i = len(sub_list)
    f.Close()
    return i

#______________________________________________________________________________
def hist_1D(tree, todraw, x_bins, cut,i):
#constrol >1 , correct pairing 
#cut_weight should either "weight && (control<2)" or "weight" or "weight && (control>1)"
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    name = "hist1_%s"%tree.GetTitle()+"_%s"%todraw+"_%s"%cut[0]+"_%d"%i
    hist = ROOT.TH1F("%s"%name,"hist1",xBins,xminBin,xmaxBin)
    tree.Draw(todraw+">>%s"%name,cut[i])
    hist.Scale(1.0/count)
    hist.SetStats(0)

    return hist 

#______________________________________________________________________________
def hist_2D(tree, todraw, x_bins, y_bins, cut,i):
#constrol >1 , correct pairing 
#cut_weight should either "weight && (control<2)" or "weight" or "weight && (control>1)"
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    name = "hist2_%s"%tree.GetTitle()+"_%s"%todraw+"_%s"%cut[0]+"%d"%i
    hist = ROOT.TH2F("%s"%name,"%s"%todraw,xBins,xminBin,xmaxBin,yBins,yminBin,ymaxBin)
    tree.Draw(todraw+">>%s"%name,cut[i])
    #hist.Scale(1.0/count)
    hist.SetStats(0)

    return hist 




#______________________________________________________________________________
def h2tohhMass_hist(tree, x_bins, cut,i):
#constrol >1 , correct pairing 
#cut_weight should either "weight && (control<2)" or "weight" or "weight && (control>1)"
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    hist = ROOT.TH1F("hist","hist",xBins,xminBin,xmaxBin)
    tree.Draw("h2tohh_Mass>>hist",cut[i])
    hist.Scale(1.0/count)
    hist.SetStats(0)

    return hist 

#______________________________________________________________________________
def offshellWMass_hist(tree, x_bins, cut,i):
#constrol >1 , correct pairing 
#cut_weight should either "weight && (control<2)" or "weight" or "weight && (control>1)"
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    hist = ROOT.TH1F("hist","hist",xBins,xminBin,xmaxBin)
    tree.Draw("offshellW_Mass>>hist",cut[i])
    hist.Scale(1.0/count)
    hist.SetStats(0)


    return hist 
#______________________________________________________________________________________
def gethiststack(tree, todraw, todrawtrue, x_bins, hs_title, cut):
 

#        cut = ["weight","weight*(control<2)","weight*(control>1)"]
        #cut = ["weight*(nu_onshellW_pt>10)","weight*(control<2&&nu_onshellW_pt>10)","weight*(control>1&&nu_onshellW_pt>10)"]
        
	name = "hs_%s"%tree.GetTitle()+"_%s"%todraw+"_%s"%cut[0]
	print "name ",name
        #hs = ROOT.THStack("hs_%s"%tree.GetTitle()+"%s"%todraw,"%s"%hs_title);
        hs = ROOT.THStack("%s"%name,"  ");
        hist1 = hist_1D(tree, todraw, x_bins, cut,0)
        hist1.SetFillColor(ROOT.kRed)
        #hist1.SetLineColor(ROOT.kRed)
        #hist1.SetLineStyle(2)
        #hist1.SetLineWidth(4)
        #hist1.Draw("E3") 
        
        hist2 = hist_1D(tree, todraw, x_bins, cut,1)
        #hist2.SetLineColor(ROOT.kBlue)
        #ddhist2.SetLineWidth(2)
        hist2.SetFillColor(ROOT.kBlue)
        #hist2.Draw("sameE3") 
   
        hist3 = hist_1D(tree, todraw, x_bins, cut,2)
        #hist3.SetLineColor(ROOT.kGreen)
        #hist3.SetLineWidth(2)
        hist3.SetFillColor(ROOT.kGreen)
        #hist3.Draw("sameE3") 

        hist4 = hist_1D(tree, todrawtrue, x_bins, cut,0)
        max = hist1.GetMaximum()
        bin_max = hist4.GetMaximumBin()
	hist4.SetBinContent(bin_max, max/4.0)
        hist4.SetFillColor(ROOT.kMagenta)
#	hist4.SetBinContent(bin_max-1, max/4.0)
	#hist4.SetBinContent(bin_max+1, max/4.0)
        #add muon information
    #    mutodraw = todraw.replace("nu","mu",1)
     #   print "todraw ",todraw, " mutodraw ",mutodraw
	hs.Add(hist3)
	hs.Add(hist2)
	hs.Add(hist1)
        hs.Add(hist4)
         

    	ctemp = ROOT.TCanvas() 
	ctemp.cd()
    	hs.Draw()
    	ctemp.SaveAs("gethiststack_%s"%todraw+"_%s"%tree.GetTitle()+"_contour.png")
	#hs.SetTitleOffset(0.8)
        return hs   
     #   c1.SaveAs("contour_%d"%(n%9)+".png") 

#______________________________________________________________________________________
def gethiststack_2d(tree, todraw, todrawtrue, x_bins, y_bins, hs_title, cut):
 

#        cut = ["weight","weight*(control<2)","weight*(control>1)"]
        #cut = ["weight*(nu_onshellW_pt>10)","weight*(control<2&&nu_onshellW_pt>10)","weight*(control>1&&nu_onshellW_pt>10)"]
        
	name = "hs_%s"%tree.GetTitle()+"_%s"%todraw.replace(":","_")+"_%s"%cut[0]
        #hs = ROOT.THStack("hs_%s"%tree.GetTitle()+"%s"%todraw,"%s"%hs_title);
        hs = ROOT.THStack("%s"%name,"  ");

        hist1 = hist_2D(tree, todraw, x_bins, y_bins, cut,0)
        hist2 = hist_2D(tree, todraw, x_bins, y_bins, cut,1)
        hist2.SetLineColor(ROOT.kBlue)
        hist2.SetLineWidth(2)
        #hist2.Draw("CONT3")
        #truevars = todrawtrue.split(':')
        #truevars = todrawtrue.split(':') 
   
        hist3 = hist_2D(tree, todraw, x_bins, y_bins, cut,2)
        hist3.SetLineColor(ROOT.kGreen)
        hist3.SetLineWidth(2)
        #hist3.Draw("CONT3same") 
        hist3max = hist3.GetBinContent(hist3.GetMaximumBin())
        hist2max = hist2.GetBinContent(hist2.GetMaximumBin())
        hist1max = max(hist3max, hist2max)
        array_contour1 = array('d',[hist1max])
        i = 0.8
        while i > 0.0:
		array_contour1.append(hist1max*i) 
	 	i = i-0.1
        hist2.SetContour(len(array_contour1), array_contour1) 
        hist3.SetContour(len(array_contour1), array_contour1) 
       
        hist4 = hist_2D(tree, todrawtrue, x_bins, y_bins, cut,1)
 	hist4.SetMarkerStyle(15)
	hist4.SetMarkerSize(1)


	hs.Add(hist1)
	#hs.Add(hist3)
	#hs.Add(hist4)
         
	
  #  	c = ROOT.TCanvas() 
 #  	hs.Draw("nostackCONT3")
#    	c.SaveAs("gethiststack2d_%s"%todraw.replace(":","_")+"_%s"%tree.GetTitle()+"_contour.png")
	#hs.SetTitleOffset(0.8)
        return hs   
     #   c1.SaveAs("contour_%d"%(n%9)+".png") 



#_______________________________________________________________________________
if __name__ == "__main__":
    
    treename = "mmctree_4204" 
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs-100k-0406-mediateStates-B3-combined.root"
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_100k_iterations1M_0406_B3.root"
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_100k_weight_0408_B3.root"
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_100k_finalStates_checksolution_0511_B3.root"
    file = "/fdata/hepx/store/user/taohuang/Hhh/0504_h2tohh/DiHiggs-1M-0504-mediateStates-B3-combined.root"
    #filedir = "/fdata/hepx/store/user/taohuang/DiHiggs_run2_PU0_htobbana_cuts_1M_filter_B3_MMC/"
    filedir = "/fdata/hepx/store/user/taohuang/DiHiggs_run2_PU0_cuts_1M_filter_B3_MMC_useMET_METcorrection/"
    #filedir = "/fdata/hepx/store/user/taohuang/DiHiggs_run2_PU0_cuts_1M_filter_B3_MMC_PTconservation/"
    #filedir = "/fdata/hepx/store/user/taohuang/DiHiggs_run2_PU0_cuts_1M_filter_B3_MMC_useMET_METcorrection_V3MMC/"
    filedir = "/fdata/hepx/store/user/taohuang/DiHiggs_run2_PU0_cuts_1M_filter_B3_MMC_useMET_METcorrectionwithMMC/"
    #filedir = "/fdata/hepx/store/user/taohuang/DiHiggs_run2_PU0_cuts_1M_filter_B3_MMC_MetNubjets/"
    dir = "DiHiggsWWAna/%s"%treename
    #dir = "DiHiggsWWAna/"
    
    totintegral = {'weight':10000,'weight1':8000,'weight2':4000,'weight3':360} 
    cut_weight3 = ["weight3","weight3*(control<2)","weight3*(control>1)"]
    cut_weight2 = ["weight2","weight2*(control<2)","weight2*(control>1)"]
    cut_weight1 = ["weight1","weight1*(control<2)","weight1*(control>1)"]
    cut_weight = ["weight","weight*(control<2)","weight*(control>1)"]
    #test(file, dir)
    if not os.path.isdir(filedir):
          print "ERROR: This is not a valid directory: ", filedir
    ls = os.listdir(filedir)
    tot = len(ls)
    num=0
    for x in ls:
	x = filedir[:]+x
	print "gettreenum ",gettreenum(x)
    	num = num+gettreenum(x)
	if num >150:
		break
	m = 1
    #	while m<gettreenum(x):
    #		monitoringMMC(x,"0621_1M",cut_weight)  
    #monitoringMMC(file,cut_weight1)
    #monitoringMMC(file,cut_weight2)
    #monitoringMMC(file,cut_weight3)
    drawh2MassAll_combined(filedir,"0710_1M_useMET_METcorrectionwithMMC")
    #drawh2MassAll_combined(filedir,"0629_1M_PTconservation","weight1")
    #drawh2Mass_combined(file,"weight2")
    #drawh2Mass_combined(file,"weight3")
    title1 = "MMC PDF for M_{H}, Event 4204"
    h2massbins = "(80,260,420)"#for843 only
    
    pic_h2mass = "h2Mass_combined_0406_contour"
    #drawh2Mass_combined(file, h2massbins, pic_h2mass)

    h2mass = "M_{H}, [GeV]"
    #cut,0 = "control<2"
    tex1 = "correctly pair (muon candidates, muon lorentzVec in MMC)"
    pic1 = "evt_4204_h2Mass_0406_contour"
  #  draw1D(file, dir, h2tohhmass, title1, h2massbins, h2mass, tex1, pic1)

    h2massbins_2 = "(70,270,340)"#for843 only
    cut0 = "control>1"
    tex2 = "incorrectly pair (muon candidates, muon lorentzVec in MMC)"
    pic2 = "evt843_h2Mass_incorrect_0406_contour"
   # draw1D(file, dir, h2tohhmass, title1, h2massbins_2, h2mass, cut,0, tex2, pic2)

