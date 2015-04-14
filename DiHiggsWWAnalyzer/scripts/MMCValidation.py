from array import *
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
def monitoringMMC(file,cut = ["weight","weight*(control<2)","weight*(control>1)"]):

    f = ROOT.TFile(file)
    list = f.GetListOfKeys()
    key = list.At(0)
    obj = key.ReadObj() #DiHiggsWWAna
    sub_list = obj.GetListOfKeys()
    sub_key = sub_list.At(1)

    #cut = ["weight*(nu_onshellW_pt>10)","weight*(control<2&&nu_onshellW_pt>10)","weight*(control>1&&nu_onshellW_pt>10)"]
    #cut = ["weight","weight*(control<2)","weight*(control>1)"]
   # cut_weight1 = ["weight1","weight1*(control<2)","weight1*(control>1)"]
    c1 = ROOT.TCanvas("c1","c1",1200,900)
    n = 0
    m = 1
    while sub_key:
	tree = sub_key.ReadObj()
	#eta-phi of nu_onshellW contour plot, nu_onshellW_ete/phi
        #mmctree_618->Draw("nu_onshellW_phi:nu_onshellW_eta>>th2","control<2")
        eta_bins = "(50,-6,6)"
        phi_bins = "(40,-3.1415,3.1415)"
	todraw1 = "nu_onshellW_phi:nu_onshellW_eta"
        hist1_1 = hist_2D(tree, todraw1, eta_bins, phi_bins, cut,0)
	hist1_1.SetLineColor(ROOT.kRed)
        
        hist1_2 = hist_2D(tree, todraw1, eta_bins, phi_bins, cut,1)
        hist1_2.SetLineColor(ROOT.kBlue)
        
        hist1_3 = hist_2D(tree, todraw1, eta_bins, phi_bins, cut,2)
        hist1_3.SetLineColor(ROOT.kGreen)#3
        hist1_3max = hist1_3.GetBinContent(hist1_3.GetMaximumBin())
        hist1_2max = hist1_2.GetBinContent(hist1_2.GetMaximumBin())
        hist1max = max(hist1_3max, hist1_2max)
        array_contour1 = array('d',[hist1max])
        i = 0.8
        while i > 0.0:
		array_contour1.append(hist1max*i) 
	 	i = i-0.1
        """
        if (hist1max*i > hist1_3max):
		array_contour1.append(hist1_3max)
		array_contour1.append(hist1_3max*0.75)
		array_contour1.append(hist1_3max*0.5)
		array_contour1.append(hist1_3max*0.25)
        if (hist1max*i > hist1_2max):
		array_contour1.append(hist1_2max)
		array_contour1.append(hist1_2max*0.75)
		array_contour1.append(hist1_2max*0.5)
		array_contour1.append(hist1_2max*0.25)
        """
        hist1_2.SetContour(len(array_contour1), array_contour1) 
        hist1_3.SetContour(len(array_contour1), array_contour1) 

        
        #hist3.SetFillColor(46)
      #  hist3.SetFillStyle(4050)
        #hist3.Draw("sameE3") 
	var_x1 = getTrueValue(tree,"eta_nuonshellW_true","(200,-6,6)")
	var_y1 = getTrueValue(tree,"phi_nuonshellW_true","(200,-3.1415,3.1415)")
	array_x = array('d',[var_x1])
	array_y = array('d',[var_y1])
	graph1 = ROOT.TGraph(1, array_x,array_y)
        graph1.SetMarkerColor(ROOT.kMagenta)
 	graph1.SetMarkerStyle(15)
	graph1.SetMarkerSize(1)
	onshellmuon_eta = getTrueValue(tree,"mu_onshellW_eta","(200,-6,6)")
	onshellmuon_phi = getTrueValue(tree,"mu_onshellW_phi","(200,-3.1415,3.1415)")
	array_mu1_eta = array('d',[onshellmuon_eta])
	array_mu1_phi = array('d',[onshellmuon_phi])
	graph1_muon = ROOT.TGraph(1, array_mu1_eta,array_mu1_phi)
        graph1_muon.SetMarkerColor(ROOT.kBlack)
 	graph1_muon.SetMarkerStyle(16)
	graph1_muon.SetMarkerSize(1)
	
        hist1_2.GetXaxis().SetTitle("#eta_{#nu}^{onshellW}")
        hist1_2.GetXaxis().SetTitleSize(0.08)
        hist1_2.GetXaxis().SetTitleOffset(0.6)
        hist1_2.GetYaxis().SetTitle("#phi_{#nu}^{onshellW}")
        hist1_2.GetYaxis().SetTitleSize(0.08)
        hist1_2.GetYaxis().SetTitleOffset(0.6)
        hist1_2.SetTitle(" ")
        #hist1_2.SetTitleOffset(0.6)


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
        hist4_1 = hist_2D(tree, todraw4, eta_bins, phi_bins, cut,0)
	hist4_1.SetLineColor(ROOT.kRed)
         
        hist4_2 = hist_2D(tree, todraw4, eta_bins, phi_bins, cut,1)
        hist4_2.SetLineColor(ROOT.kBlue)
        hist4_3 = hist_2D(tree, todraw4, eta_bins, phi_bins, cut,2)
        hist4_3.SetLineColor(ROOT.kGreen)#3
        hist4_3max = hist4_3.GetBinContent(hist4_3.GetMaximumBin())
        hist4_2max = hist4_2.GetBinContent(hist4_2.GetMaximumBin())
        hist4max = max(hist4_2max, hist4_3max)
        array_contour4 = array('d',[hist4max])
        i = 0.8
        while i > 0.0:
		array_contour4.append(hist4max*i) 
	 	i = i-0.1
        """
        if (hist4max*i > hist4_3max):
		array_contour4.append(hist4_3max)
		array_contour4.append(hist4_3max*0.75)
		array_contour4.append(hist4_3max*0.5)
		array_contour4.append(hist4_3max*0.25)
        if (hist4max*i > hist4_2max):
		array_contour4.append(hist4_2max)
		array_contour4.append(hist4_2max*0.75)
		array_contour4.append(hist4_2max*0.5)
		array_contour4.append(hist4_2max*0.25)
        """
#        print array_contour4
        hist4_2.SetContour(len(array_contour4), array_contour4) 
        hist4_3.SetContour(len(array_contour4), array_contour4) 
        #hist3.SetFillColor(46)
      #  hist3.SetFillStyle(4050)
        #hist3.Draw("sameE3") 
	var_x4 = getTrueValue(tree,"eta_nuoffshellW_true","(200,-6,6)")
	var_y4 = getTrueValue(tree,"phi_nuoffshellW_true","(200,-3.1415,3.1415)")
	array_x4 = array('d',[var_x4])
	array_y4 = array('d',[var_y4])
	graph4 = ROOT.TGraph(1, array_x4,array_y4)
        graph4.SetMarkerColor(ROOT.kMagenta)
 	graph4.SetMarkerStyle(15)
	graph4.SetMarkerSize(1)
	offshellmuon_eta = getTrueValue(tree,"mu_offshellW_eta","(200,-6,6)")
	offshellmuon_phi = getTrueValue(tree,"mu_offshellW_phi","(200,-3.1415,3.1415)")
	array_mu2_eta = array('d',[offshellmuon_eta])
	array_mu2_phi = array('d',[offshellmuon_phi])
	graph4_muon = ROOT.TGraph(1, array_mu2_eta,array_mu2_phi)
        graph4_muon.SetMarkerColor(ROOT.kBlack)
 	graph4_muon.SetMarkerStyle(16)
	graph4_muon.SetMarkerSize(1)

        hist4_2.GetXaxis().SetTitle("#eta_{#nu}^{offshellW}")
        hist4_2.GetXaxis().SetTitleSize(0.08)
        hist4_2.GetXaxis().SetTitleOffset(0.6)
        hist4_2.GetYaxis().SetTitle("#phi_{#nu}^{offshellW}")
        hist4_2.GetYaxis().SetTitleSize(0.08)
        hist4_2.GetYaxis().SetTitleOffset(0.6)
        hist4_2.SetTitle("  ")
#	hist4_2.SetTitleOffset(0.6)

	
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
	#tot
	c1.Clear()
	c1.Divide(4,3)
	c1.cd(1)
	hist1_2.Draw("CONT3")
	hist1_3.Draw("CONT3same")
	graph1.Draw("psame")
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
	hist4_2.Draw("CONT3")
	hist4_3.Draw("CONT3same")
	graph4.Draw("psame")
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
	c1.Update()

        c1.cd()
        met = getTrueValue(tree,"met_true","(300,0,150)")
    	tex = ROOT.TLatex(0.15,.954,"Performance of MMC,#slash{E}_{T}=%.2f"%met+", %s"%(tree.GetTitle()))
    	#tex = ROOT.TLatex(0.25,.954,"Performance of MMC"+", %s"%(tree.GetTitle()))
    	tex.SetTextSize(0.04)
    	tex.SetTextFont(42)
    	tex.SetNDC()
	tex.Draw("same")
	c1.Update()
        c1.SaveAs("%s"%(tree.GetTitle())+"_0410_weight_newframe_B3.pdf")
        c1.SaveAs("%s"%(tree.GetTitle())+"_0410_weight_newframe_B3.png")
        c1.SaveAs("%s"%(tree.GetTitle())+"_0410_weight_newframe_B3.C")
        ##########################################
        ##### fill hist_h2
        
	#if m>1:
	#	break
	m = m+1
        sub_key = sub_list.At(m)


#____________________________________________________________________________
def drawh2Mass_combined(file, cut="weight"):

    f = ROOT.TFile(file)
    list = f.GetListOfKeys()
    key = list.At(0)
    obj = key.ReadObj() #DiHiggsWWAna
    sub_list = obj.GetListOfKeys()
    sub_key = sub_list.At(1)
    hist_h2 = ROOT.TH1F("hist_h2"," ",150,200,500)
    hist_h2true = ROOT.TH1F("hist_h2true"," ",150,200,500)
    m = 1
    while sub_key:
	tree = sub_key.ReadObj()
	print "tree title ",tree.GetTitle()," entries ", tree.GetEntries()
	hist_h2.Fill(geth2MassMostProba(tree,cut))
	hist_h2true.Fill(getTrueValue(tree,"mass_h2_true","(400,200,600)"))
        h2_mass_reconstructed = geth2MassMostProba(tree,cut)
        h2_mass_true = getTrueValue(tree,"mass_h2_true","(400,200,600)")
	if (abs(h2_mass_reconstructed-h2_mass_true)>5 or h2_mass_reconstructed > 450 or h2_mass_reconstructed < 250):
		print tree.GetTitle()," reconstructed mass  ", h2_mass_reconstructed, " true mass ", h2_mass_true
        m = m+1
	sub_key = sub_list.At(m)
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
    h2Mass_c.SaveAs("MMC_h2Mass_0410_weight_newframe_B3.pdf") 
    h2Mass_c.SaveAs("MMC_h2Mass_0410_weight_newframe_B3.png") 
    h2Mass_c.SaveAs("MMC_h2Mass_0410_weight_newframe_B3.C") 

#___________________________________________________________________________
def geth2MassMostProba(t, cut = "weight"):

    b1 = ROOT.TH1F("b1","b1",400,200,600)
    t.Draw("h2tohh_Mass>>b1",cut)
    return b1.GetXaxis().GetBinCenter(b1.GetMaximumBin())

#____________________________________________________________________________
def getTrueValue(t,var, x_bins):
   
    c = ROOT.TCanvas() 
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    name = "tmp_%s"%var
    b1 = ROOT.TH1F("%s"%name,"b1",xBins,xminBin,xmaxBin)
    t.Draw("%s"%var+">>%s"%name)
    b1.Draw()
   # c.SaveAs("getTrueVaule_%s"%var+"_contour.png")
    value = b1.GetXaxis().GetBinCenter(b1.GetMaximumBin())
    #print "xaxis", var, "b1 maximum", b1.GetMaximumBin()," value", value
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
    sub_key = sub_list.At(0)
    i = 0
    while sub_key:
	tree = sub_key.ReadObj()
#	tree.Print()
        i = i+1
        sub_key = sub_list.At(i)
    print "key getclassname ", key.GetClassName()
    print "c1 ",c1
    print "obj  ", obj
    print "sub_key ", sub_key
    print "sub_list 1", sub_list.At(1)
    print "dir ", dir, "  t ", t, " tree ", tree.GetTitle() 
    pt_h2_bins = "(100,0,0.00001)"
    todraw8 = "h2tohh_Pt"
    todrawtrue8 = "pt_h2_true"
    hs_title8 = "p_{TH}"
    hs8 = gethiststack(tree, todraw8,todrawtrue8, pt_h2_bins, hs_title8)
    c = ROOT.TCanvas() 
    hs8.Draw("nostack")
    c.SaveAs("gethiststack_contour.png")

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
def gethiststack(tree, todraw, todrawtrue, x_bins, hs_title, cut=["weight","weight*(control<2)","weight*(control>1)"]):
 

#        cut = ["weight","weight*(control<2)","weight*(control>1)"]
        #cut = ["weight*(nu_onshellW_pt>10)","weight*(control<2&&nu_onshellW_pt>10)","weight*(control>1&&nu_onshellW_pt>10)"]
        
	name = "hs_%s"%tree.GetTitle()+"_%s"%todraw+"_%s"%cut[0]
        #hs = ROOT.THStack("hs_%s"%tree.GetTitle()+"%s"%todraw,"%s"%hs_title);
        hs = ROOT.THStack("%s"%name,"  ");
        hist1 = hist_1D(tree, todraw, x_bins, cut,0)
        hist1.SetLineColor(ROOT.kRed)
        hist1.SetLineStyle(2)
        hist1.SetLineWidth(4)
        #hist1.Draw("E3") 
        
        hist2 = hist_1D(tree, todraw, x_bins, cut,1)
        hist2.SetLineColor(ROOT.kBlue)
        hist2.SetLineWidth(2)
        #hist2.Draw("sameE3") 
   
        hist3 = hist_1D(tree, todraw, x_bins, cut,2)
        hist3.SetLineColor(ROOT.kGreen)
        hist3.SetLineWidth(2)
        #hist3.Draw("sameE3") 

        hist4 = hist_1D(tree, todrawtrue, x_bins, cut,0)
        hist4.SetFillColor(ROOT.kMagenta)
        max = hist1.GetMaximum()
        bin_max = hist4.GetMaximumBin()
	hist4.SetBinContent(bin_max, max/4.0)
#	hist4.SetBinContent(bin_max-1, max/4.0)
	#hist4.SetBinContent(bin_max+1, max/4.0)
        #add muon information
    #    mutodraw = todraw.replace("nu","mu",1)
     #   print "todraw ",todraw, " mutodraw ",mutodraw
	hs.Add(hist1)
	hs.Add(hist2)
	hs.Add(hist3)
        hs.Add(hist4)
         
	
    	#c = ROOT.TCanvas() 
    #	hs.Draw("nostack")
    #	c.SaveAs("gethiststack_%s"%todraw+"_%s"%tree.GetTitle()+"_contour.png")
	#hs.SetTitleOffset(0.8)
        return hs   
     #   c1.SaveAs("contour_%d"%(n%9)+".png") 


#_______________________________________________________________________________
if __name__ == "__main__":
    
    treename = "mmctree_4204" 
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs-100k-0406-mediateStates-B3-combined.root"
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_100k_iterations1M_0406_B3.root"
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_100k_weight_0408_B3.root"
    file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_100k_0410_newframe_B3.root"
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_100k_0406_wmass_B3.root"
    dir = "DiHiggsWWAna/%s"%treename
    #dir = "DiHiggsWWAna/"
    
    cut_weight2 = ["weight2","weight2*(control<2)","weight2*(control>1)"]
    cut_weight1 = ["weight1","weight1*(control<2)","weight1*(control>1)"]
    cut_weight = ["weight","weight*(control<2)","weight*(control>1)"]
    #test(file, dir)
    monitoringMMC(file,cut_weight)  
    #monitoringMMC(file,cut_weight1)
    #monitoringMMC(file,cut_weight2)
    #drawh2Mass_combined(file)
    #drawh2Mass_combined(file,"weight2")
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

    h2massbins_3 = "(70,270,340)"#for843 only
    cut2 = " "
    tex3 = "inclduing both pairs (muon candidates, muon lorentzVec in MMC)"
    pic3 = "evt843_h2Mass_both_0406_contour"
  #  draw1D(file, dir, h2tohhmass, title1, h2massbins_3, h2mass, cut,2, tex3, pic3)
    #draw2D(file,dir,"nu_onshellW_phi:nu_onshellW_eta","eta_nuonshellW_true","phi_nuonshellW_true","#phi Vs #eta of nu_onshell","(50,-6,6)","#eta","(30,-3.14,3.14)","#phi","tex","etaVsphi_nuonshellW_0406_contour")
    #draw2D(file,dir,"nu_offshellW_phi:nu_offshellW_eta","eta_nuoffshellW_true","phi_nuoffshellW_true","#phi Vs #eta of nu_offshell","(50,-6,6)","#eta","(30,-3.14,3.14)","#phi","tex","etaVsphi_nuoffshellW_0406_contour")
    #draw2DAll(file,"nu_onshellW_pt:onshellW_Mass","onshell W#rightarrow #mu#nu ","(50,45,95)","M_{W}^{onshell}", "(125,0,125)","p_{T#nu}^{onshellW}","onshellnuptVsWmass","onshellnuptVsWmass")
    #draw2DAll(file,"nu_onshellW_pt:onshellW_Mass","onshell W#rightarrow #mu#nu ","(50,45,95)","M_{W}^{onshell}", "(125,0,125)","p_{T#nu}^{onshellW}","onshellnuptVsWmass","onshellnuptVsWmass","weight1")
    """ 
    offshellWmass = "offshellW_Mass"
    titleW1 = "off-shell W mass from MMC, Event843 "

    offshellWmassbins = "(50,0,50)"
    pic_offshellW = "offshellW_combined_0406_contour"
    drawoffshellWmass_combined(file, offshellWmassbins, pic_offshellW)
    
    offshellnu_eta = "nu_offshellW_eta"
    offshellnu_eta_true = "eta_nuoffshellW_true"
    nuetabins= "(100,-7,7)"
    pic_offshellnu_eta = "offshellW_nueta_0406_contour"
    tex_offshellnu_eta = "ensembles of allowable solutions for #eta_{#nu}^{offshellW}"
    draw_combined(file, offshellnu_eta, offshellnu_eta_true, nuetabins, tex_offshellnu_eta, pic_offshellnu_eta)

    onshellnu_eta = "nu_onshellW_eta"
    onshellnu_eta_true = "eta_nuonshellW_true"
    pic_onshellnu_eta = "onshellW_nueta_0406_contour"
    tex_onshellnu_eta = "ensembles of allowable solutions for #eta_{#nu}^{onshellW}"
    draw_combined(file, onshellnu_eta, onshellnu_eta_true, nuetabins, tex_onshellnu_eta, pic_onshellnu_eta)
 
    offshellnu_pt = "nu_offshellW_pt"
    offshellnu_pt_true = "pt_nuoffshellW_true"
    nuptbins = "(60,0,120)"
    pic_offshellnu_pt = "offshellW_nupt_0406_contour"
    tex_offshellnu_pt = "ensembles of allowable solutions for p_{T}_{#nu}^{offshellW}"
    draw_combined(file, offshellnu_pt, offshellnu_pt_true, nuptbins, tex_offshellnu_pt, pic_offshellnu_pt)
    
    onshellnu_pt = "nu_onshellW_pt"
    onshellnu_pt_true = "pt_nuonshellW_true"
    nuptbins = "(60,0,120)"
    pic_onshellnu_pt = "onshellW_nupt_0406_contour"
    tex_onshellnu_pt = "ensembles of allowable solutions for p_{T}_{#nu}^{onshellW}"
    draw_combined(file, onshellnu_pt, onshellnu_pt_true, nuptbins, tex_onshellnu_pt, pic_onshellnu_pt)
    

    offshellWmass_x = "off-shell W mass"
    cutW1 = "control<2"
    texW1 = "correctly pair (muon candidates, muon lorentzVec in MMC)"
    picW1 = "evt843_offshellWmass_correct_0406_contour"
  #  draw1D(file, dir, offshellWmass, titleW1, offshellWmassbins, offshellWmass_x, cutW1, texW1, picW1)

    cutW2 = "control>2"
    texW2 = "incorrectly pair (muon candidates, muon lorentzVec in MMC)"
    picW2 = "evt843_offshellWmass_incorrect_0406_contour"
  #  draw1D(file, dir, offshellWmass, titleW1, offshellWmassbins, offshellWmass_x, cutW2, texW2, picW2)

    cutW3 = ""
    texW3 = "incluing both pairs (muon candidates, muon lorentzVec in MMC)"
    picW3 = "evt843_offshellWmass_both_0406_contour"
   # draw1D(file, dir, offshellWmass, titleW1, offshellWmassbins, offshellWmass_x, cutW1, texW3, picW3)

    """
