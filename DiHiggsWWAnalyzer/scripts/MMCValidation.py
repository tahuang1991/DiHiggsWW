from array import *
import ROOT

ROOT.gROOT.SetBatch(1)
#gStyle from TStyle
ROOT.gStyle.SetStatW(0.17)
ROOT.gStyle.SetStatH(0.15)

#ROOT.gStyle.SetOptStat(111110)
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
ROOT.gStyle.SetPadTopMargin(0.18)
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

    cut1 = "weight"
    cut2 = "weight*(control<2)"
    cut3 = "weight*(control>1)"
    cut4 = "weight"
    f = ROOT.TFile(file)
    tree = f.Get(dir)
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    
    b1 = ROOT.TH1F("b1","b1",xBins,xminBin,xmaxBin)
    b1.SetTitle("%s"%title+" "*8 + "CMS Simulation Preliminary")
 #   b1.Draw()
    combined = 1
    if combined:
    #b1.SetStats(0)
        hs = ROOT.THStack("hs","hs");
        hs.SetTitle("%s"%title+" "*8 + "CMS Simulation Preliminary")
        hist1 = hist_1D(tree, h2tohhmass, x_bins, cut1)
        hist1.SetFillColor(ROOT.kRed)
        hist1.GetXaxis().SetTitle("%s"%x_title)
        hist1.GetYaxis().SetTitle("Possibility")
     	hs.Draw("nostack")
     #   c1.SaveAs("test_%d"%(n%9)+".png") 
        #hist1.SetFillColor(16)
      #  hist1.SetFillStyle(4050)
        max = hist1.GetMaximum()
        #hist1.Draw("E3") 
	hs.Add(hist1)
        
        hist2 = hist_1D(tree, h2tohhmass, x_bins, cut2)
        #hist2.SetLineColor(ROOT.kBlue)
        #hist2.SetLineWidth(4)
	#hist2.SetMarkerStyle(21)
        hist2.SetFillColor(ROOT.kBlue)
        hist2.SetFillStyle(3445)
        #hist2.Draw("sameE3") 
	hs.Add(hist2)
   
        hist3 = hist_1D(tree, h2tohhmass, x_bins, cut3)
        #hist3.SetLineColor(ROOT.kGreen)#3
        #hist3.SetLineWidth(4)#3
#	hist3.SetMarkerStyle(20)
#	hist3.SetMarkerColor(ROOT.kGreen)
        hist3.SetFillColor(ROOT.kGreen)
        hist3.SetFillStyle(3454)
        #hist3.Draw("sameE3") 
	hs.Add(hist3)

        hist4 = hist_1D(tree, h2truemass, x_bins, cut4)
        hist4.SetFillColor(ROOT.kMagenta)
   #     hist4.SetFillStyle(4050)
        bin_max = hist4.GetMaximumBin()
	hist4.SetBinContent(bin_max, max/4.0)
        hs.Add(hist4)
        
        #hs.GetXaxis().SetTitle("%s"%x_title)
        #hs.GetYaxis().SetTitle("Possibility")
        ROOT.gStyle.SetHatchesLineWidth(2)
     	hs.Draw("nostack")
     #   c1.SaveAs("test_%d"%(n%9)+".png") 
        hs.GetHistogram().GetXaxis().SetTitle("%s"%x_title)
        hs.GetHistogram().GetYaxis().SetTitle("Possibility")
        c1.Update()
    tex2 = ROOT.TLatex(0.15,0.6,"%s "%tex)
    tex2.SetTextSize(0.04)
    tex2.SetTextFont(42)
    tex2.SetNDC()
    #tex2.Draw("same")
	
    c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3.pdf")
    c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3.png")

#___________________________________________
def draw2D(file,dir,todraw,xaxis,yaxis,title,x_bins,x_title,y_bins,y_title,tex,pic_name):
    
    c1 = ROOT.TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()

    cut1 = "weight"
    cut2 = "weight*(control<2)"
    cut3 = "weight*(control>1)"
    cut4 = "weight"
    f = ROOT.TFile(file)
    tree = f.Get(dir)
#    todraw = "%s"%yaxis+":%s"%xaxis
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    
    b1 = ROOT.TH2F("b1","b1",xBins,xminBin,xmaxBin, yBins, yminBin, ymaxBin)
    b1.SetTitle("%s"%title+" "*8 + "CMS Simulation Preliminary")
 #   b1.Draw()
    combined = 1
    if combined:
    #b1.SetStats(0)
        hist1 = hist_2D(tree, todraw, x_bins, y_bins, cut1)
	hist1.SetLineColor(ROOT.kRed)
        
        hist2 = hist_2D(tree, todraw, x_bins, y_bins, cut2)
        hist2.SetLineColor(ROOT.kBlue)
        #hist2.SetFillColor(42)
        #hist2.SetFillStyle(4050)
        #hist2.Draw("sameE3")
	 
   
        hist3 = hist_2D(tree, todraw, x_bins, y_bins, cut3)
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
     
     
    cut1 = "weight"
    cut2 = "weight*(control<2)"
    cut3 = "weight*(control>1)"
    cut4 = "weight"
    c1 = ROOT.TCanvas()
    #ROOT.gPad.SetFillStyle(4050)
    n = 0;
    tex = "Probability Density Function for M_{H} in MMC"
    tex2 = ROOT.TLatex(0.25,0.95,"%s "%tex)
    tex2.SetTextSize(0.04)
    tex2.SetTextFont(42)
    tex2.SetNDC()
    #test = sub_key.ReadObj()
    #testhist = hist_1D(test, h2tohhmass, x_bins, cut1)
    #testhist.Draw()
    #c1.SaveAs("test.png")
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
        hist1 = hist_1D(tree, h2tohhmass, x_bins, cut1)
	hist1.SetMarkerStyle(21)
        hist1.SetFillColor(ROOT.kRed)
        #hist1.SetFillColor(16)
      #  hist1.SetFillStyle(4050)
        max = hist1.GetMaximum()
        #hist1.Draw("E3") 
	hs.Add(hist1)
        
        hist2 = hist_1D(tree, h2tohhmass, x_bins, cut2)
        hist2.SetLineColor(ROOT.kBlue)
        hist2.SetLineWidth(4)
        #hist2.SetFillColor(42)
        #hist2.SetFillStyle(4050)
        #hist2.Draw("sameE3") 
	hs.Add(hist2)
   
        hist3 = hist_1D(tree, h2tohhmass, x_bins, cut3)
        hist3.SetLineColor(ROOT.kGreen)#3
        hist3.SetLineWidth(4)
        #hist3.SetFillColor(46)
      #  hist3.SetFillStyle(4050)
        #hist3.Draw("sameE3") 
	hs.Add(hist3)

        hist4 = hist_1D(tree, h2truemass, x_bins, cut4)
        hist4.SetFillColor(ROOT.kMagenta)
        hist4.SetFillStyle(4050)
        bin_max = hist4.GetMaximumBin()
	hist4.SetBinContent(bin_max, max/4.0)
        hs.Add(hist4)
        
        #hs.SetTitle("%s"%title+" "*8 + "CMS Simulation Preliminary")
        hslist.append(hs)
     	hs.Draw("nostack")
     #   c1.SaveAs("test_%d"%(n%9)+".png") 
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
def drawoffshellWmass_combined(file,x_bins,pic_name):
    
    f = ROOT.TFile(file)
    list = f.GetListOfKeys()
    key = list.At(0)
    obj = key.ReadObj() #DiHiggsWWAna
    sub_list = obj.GetListOfKeys()
    sub_key = sub_list.At(1)

    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    
    cut1 = "weight"
    cut2 = "weight*(control<2)"
    cut3 = "weight*(control>1)"
    cut4 = "weight"
    c1 = ROOT.TCanvas()
    c1.Divide(3,3)
    n = 0;
    tex = "ensembles of allowable solutions for M_{W}^{offshell} in MMC"
    tex2 = ROOT.TLatex(0.25,.950,"%s "%tex)
    tex2.SetTextSize(0.04)
    tex2.SetTextFont(42)
    tex2.SetNDC()
    #test = sub_key.ReadObj()
    #testhist = hist_1D(test, h2tohhmass, x_bins, cut1)
    #testhist.Draw()
    #c1.SaveAs("test.png")
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
        hist1 = hist_1D(tree, offshellWmass, x_bins, cut1)
        hist1.SetFillColor(ROOT.kRed)
        #hist1.Draw("E3") 
	hs.Add(hist1)
        
        hist2 = hist_1D(tree, offshellWmass, x_bins, cut2)
        hist2.SetLineColor(ROOT.kBlue)
        hist2.SetLineWidth(4)
        #hist2.Draw("sameE3") 
	hs.Add(hist2)
   
        hist3 = hist_1D(tree, offshellWmass, x_bins, cut3)
        hist3.SetLineColor(ROOT.kGreen)
        hist3.SetLineWidth(4)
        #hist3.Draw("sameE3") 
	hs.Add(hist3)

        hist4 = hist_1D(tree, offshellWtruemass, x_bins, cut4)
        hist4.SetFillColor(ROOT.kMagenta)
        max = hist1.GetMaximum()
        bin_max = hist4.GetMaximumBin()
	hist4.SetBinContent(bin_max, max/4.0)
        hs.Add(hist4)
        
        hslist.append(hs)
     	hs.Draw("nostack")
     #   c1.SaveAs("test_%d"%(n%9)+".png") 
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
    
    cut1 = "weight"
    cut2 = "weight*(control<2)"
    cut3 = "weight*(control>1)"
    cut4 = "weight"
    c1 = ROOT.TCanvas()
    c1.Divide(3,3)
    n = 0;
    #tex = "test"
    tex2 = ROOT.TLatex(0.25,.950,"%s "%tex)
    tex2.SetTextSize(0.04)
    tex2.SetTextFont(42)
    tex2.SetNDC()
    #test = sub_key.ReadObj()
    #testhist = hist_1D(test, h2tohhmass, x_bins, cut1)
    #testhist.Draw()
    #c1.SaveAs("test.png")
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
        hist1 = hist_1D(tree, todraw, x_bins, cut1)
        hist1.SetFillColor(ROOT.kRed)
        #hist1.Draw("E3") 
	hs.Add(hist1)
        
        hist2 = hist_1D(tree, todraw, x_bins, cut2)
        hist2.SetLineColor(ROOT.kBlue)
        hist2.SetLineWidth(6)
        #hist2.SetLineStyle(20)
        #hist2.Draw("sameE3") 
	hs.Add(hist2)
   
        hist3 = hist_1D(tree, todraw, x_bins, cut3)
        hist3.SetLineColor(ROOT.kGreen)
        hist3.SetLineWidth(4)
        #hist2.SetLineStyle(21)
        #hist3.Draw("sameE3") 
	hs.Add(hist3)

        hist4 = hist_1D(tree, todrawtrue, x_bins, cut4)
        hist4.SetFillColor(ROOT.kMagenta)
        max = hist1.GetMaximum()
        bin_max = hist4.GetMaximumBin()
	hist4.SetBinContent(bin_max, max/4.0)
	hist4.SetBinContent(bin_max-1, max/4.0)
	#hist4.SetBinContent(bin_max+1, max/4.0)
        hs.Add(hist4)
        
        hslist.append(hs)
     	hs.Draw("nostack")
     #   c1.SaveAs("test_%d"%(n%9)+".png") 
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
def monitoringMMC(file):

    f = ROOT.TFile(file)
    list = f.GetListOfKeys()
    key = list.At(0)
    obj = key.ReadObj() #DiHiggsWWAna
    sub_list = obj.GetListOfKeys()
    sub_key = sub_list.At(1)

    cut1 = "weight"
    cut2 = "weight*(control<2)"
    cut3 = "weight*(control>1)"
    cut4 = "weight"

    c1 = ROOT.TCanvas("c1","c1",700,700)
    n = 0
    m = 1
    while sub_key:
	c1.Clear()
	c1.Divide(3,3)
	tree = sub_key.ReadObj()
	#eta-phi of nu_onshellW contour plot, nu_onshellW_ete/phi
        #mmctree_618->Draw("nu_onshellW_phi:nu_onshellW_eta>>th2","control<2")
	c1.cd(1)
        eta_bins = "(50,-6,6)"
        phi_bins = "(40,-3.1415,3.1415)"
	todraw1 = "nu_onshellW_phi:nu_onshellW_eta"
        hist1_1 = hist_2D(tree, todraw1, eta_bins, phi_bins, cut1)
	hist1_1.SetLineColor(ROOT.kRed)
        
        hist1_2 = hist_2D(tree, todraw1, eta_bins, phi_bins, cut2)
        hist1_2.SetLineColor(ROOT.kBlue)
        #hist2.SetFillColor(42)
        #hist2.SetFillStyle(4050)
        #hist2.Draw("sameE3")
	 
   
        hist1_3 = hist_2D(tree, todraw1, eta_bins, phi_bins, cut3)
        hist1_3.SetLineColor(ROOT.kGreen)#3

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
	graph1.SetMarkerSize(2)

        hist1_2.GetXaxis().SetTitle("#eta_{#nu}^{onshellW}")
        hist1_2.GetYaxis().SetTitle("#phi_{#nu}^{onshellW}")
	hist1_2.Draw("CONT3")
	hist1_3.Draw("CONT3same")
	graph1.Draw("psame")
	"""
        #pt of nu_onshellW, nu_onshellW_pt,
	c1.cd(2)
        todraw2 = "nu_onshellW_pt"
        todrawtrue2 = "pt_nuonshellW_true"
	nupt_onshell_bins = "(60,0,120)"
	hs_title2 = "p_{T#nu}^{onshellW}"
	hs2 = gethiststack(tree, todraw2, todrawtrue2, nupt_onshell_bins, hs_title2)
	#hs2.SetName("%s"%tree.GetTitle()+"_%s"%todraw2)
     	hs2.Draw("nostack")


        #mass of onshellW, onshellW_Mass, 
	c1.cd(3)
        todraw3 = "onshellW_Mass"
        todrawtrue3 = "mass_onshellW_true"
	wmass_onshell_bins = "(20,75,85)"
	hs_title3 = "M_{W}^{onshell}"
	hs3 = gethiststack(tree, todraw3, todrawtrue3, wmass_onshell_bins, hs_title3)
	#hs3.SetName("%s"%tree.GetTitle()+"_%s"%todraw3)
     	hs3.Draw("nostack")


        #eta-phi of nu_offshellW
	c1.cd(4)
	todraw4 = "nu_offshellW_phi:nu_offshellW_eta"
        hist4_1 = hist_2D(tree, todraw4, eta_bins, phi_bins, cut1)
	hist4_1.SetLineColor(ROOT.kRed)
        
        hist4_2 = hist_2D(tree, todraw4, eta_bins, phi_bins, cut2)
        hist4_2.SetLineColor(ROOT.kBlue)
        #hist2.SetFillColor(42)
        #hist2.SetFillStyle(4050)
        #hist2.Draw("sameE3")
	 
   
        hist4_3 = hist_2D(tree, todraw4, eta_bins, phi_bins, cut3)
        hist4_3.SetLineColor(ROOT.kGreen)#3

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
	graph4.SetMarkerSize(2)

        hist4_2.GetXaxis().SetTitle("#eta_{#nu}^{offshellW}")
        hist4_2.GetYaxis().SetTitle("#phi_{#nu}^{offshellW}")
	hist4_2.Draw("CONT3")
	hist4_3.Draw("CONT3same")
	graph4.Draw("psame")
	

        #pt of nu_offshellW
	c1.cd(5)
        todraw5 = "nu_offshellW_pt"
        todrawtrue5 = "pt_nuoffshellW_true"
	nupt_offshell_bins = "(40,0,80)"
	hs_title5 = "p_{T#nu}^{offshellW}"
	hs5 = gethiststack(tree, todraw5, todrawtrue5, nupt_offshell_bins, hs_title5)
	#hs5.SetName("%s"%tree.GetTitle()+"_%s"%todraw5)
     	hs5.Draw("nostack")


	#mass of offshellW
	c1.cd(6)
        todraw6 = "offshellW_Mass"
	todrawtrue6 = "mass_offshellW_true"
	wmass_offshell_bins = "(50,0,50)"
	hs_title6 = "M_{W}^{offshell}"
	hs6 = gethiststack(tree, todraw6, todrawtrue6, wmass_offshell_bins, hs_title6)
	hs6.SetName("%s"%tree.GetTitle()+"_%s"%todraw6)
     	hs6.Draw("nostack")


	#mass of h, htoWW_Mass
	c1.cd(7)
	hmass_bins = "(20,124,126)"
        todraw7 = "htoWW_Mass"
        todrawtrue7 = "mass_htoWW_true"
	hs_title7 = "M_{h#rightarrow WW}"
	hs7 = gethiststack(tree, todraw7,todrawtrue7, hmass_bins, hs_title7)
	#hs7.SetName("%s"%tree.GetTitle()+"_%s"%todraw7)
	hs7.Draw("nostack")
      
        #pt of h2
        c1.cd(8)
        pt_h2_bins = "(100,0,0.00001)"
	todraw8 = "h2tohh_Pt"
	todrawtrue8 = "pt_h2_true"
	hs_title8 = "p_{TH}"
	hs8 = gethiststack(tree, todraw8,todrawtrue8, pt_h2_bins, hs_title8)
	#hs8.SetName("%s"%tree.GetTitle()+"_%s"%todraw8)
        hs8.Draw("nostack")
 	#mass of h2,h2tohh_Mass
	c1.cd(9)
        #todraw
        h2mass_bins = "(30,280,370)"
	todraw9 = "h2tohh_Mass"
        todrawtrue9 = "mass_h2_true"
	hs_title9 = "PDF of M_{H}"
	hs9 = gethiststack(tree, todraw9,todrawtrue9, h2mass_bins, hs_title9)
	#hs9.SetName("%s"%tree.GetTitle()+"_%s"%todraw9)
     	hs9.Draw("nostack")
	"""
	#tot
	c1.cd()
    	tex = ROOT.TLatex(0.25,.950,"Monitoring Performance of MMC, %s"%(tree.GetTitle()))
    	tex.SetTextSize(0.04)
    	tex.SetTextFont(42)
    	tex.SetNDC()
	#tex.Draw("same")
	c1.Update()
        c1.SaveAs("%s"%(tree.GetTitle())+"_test_0325_B3.pdf")
        c1.SaveAs("%s"%(tree.GetTitle())+"_test_0325_B3.png")
        c1.SaveAs("%s"%(tree.GetTitle())+"_test_0325_B3.C")
	if m>1:
		break
	m = m+1
        sub_key = sub_list.At(m)

#____________________________________________________________________________
def getTrueValue(t,var, x_bins):
   
    c = ROOT.TCanvas() 
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    b1 = ROOT.TH1F("b1","b1",xBins,xminBin,xmaxBin)
    t.Draw("%s"%var+">>b1")
    b1.Draw()
   # c.SaveAs("getTrueVaule_%s"%var+"_test.png")
    value = b1.GetXaxis().GetBinCenter(b1.GetMaximumBin())
    print "xaxis", var, "b1 maximum", b1.GetMaximumBin()," value", value
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
    c.SaveAs("gethiststack_test.png")

#______________________________________________________________________________
def hist_1D(tree, todraw, x_bins, cut_weight):
#constrol >1 , correct pairing 
#cut_weight should either "weight && (control<2)" or "weight" or "weight && (control>1)"
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    name = "hist1_%s"%tree.GetTitle()+"_%s"%todraw+"_%s"%cut_weight
    hist = ROOT.TH1F("%s"%name,"hist1",xBins,xminBin,xmaxBin)
    tree.Draw(todraw+">>%s"%name,cut_weight)
    hist.Scale(1.0/count)
    hist.SetStats(0)

    return hist 

#______________________________________________________________________________
def hist_2D(tree, todraw, x_bins, y_bins, cut_weight):
#constrol >1 , correct pairing 
#cut_weight should either "weight && (control<2)" or "weight" or "weight && (control>1)"
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])
    name = "hist2_%s"%tree.GetTitle()+"_%s"%todraw
    hist = ROOT.TH2F("%s"%name,"hist2",xBins,xminBin,xmaxBin,yBins,yminBin,ymaxBin)
    tree.Draw(todraw+">>%s"%name,cut_weight)
    #hist.Scale(1.0/count)
    hist.SetStats(0)

    return hist 




#______________________________________________________________________________
def h2tohhMass_hist(tree, x_bins, cut_weight):
#constrol >1 , correct pairing 
#cut_weight should either "weight && (control<2)" or "weight" or "weight && (control>1)"
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    hist = ROOT.TH1F("hist","hist",xBins,xminBin,xmaxBin)
    tree.Draw("h2tohh_Mass>>hist",cut_weight)
    hist.Scale(1.0/count)
    hist.SetStats(0)

    return hist 

#______________________________________________________________________________
def offshellWMass_hist(tree, x_bins, cut_weight):
#constrol >1 , correct pairing 
#cut_weight should either "weight && (control<2)" or "weight" or "weight && (control>1)"
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    hist = ROOT.TH1F("hist","hist",xBins,xminBin,xmaxBin)
    tree.Draw("offshellW_Mass>>hist",cut_weight)
    hist.Scale(1.0/count)
    hist.SetStats(0)


    return hist 
#______________________________________________________________________________________
def gethiststack(tree, todraw, todrawtrue, x_bins, hs_title):


    	cut1 = "weight"
   	cut2 = "weight*(control<2)"
    	cut3 = "weight*(control>1)"
    	cut4 = "weight"
	name = "hs_%s"%tree.GetTitle()+"_%s"%todraw
        #hs = ROOT.THStack("hs_%s"%tree.GetTitle()+"%s"%todraw,"%s"%hs_title);
        hs = ROOT.THStack("%s"%name,"%s"%hs_title);
        hist1 = hist_1D(tree, todraw, x_bins, cut1)
        hist1.SetFillColor(ROOT.kRed)
        #hist1.Draw("E3") 
	hs.Add(hist1)
        
        hist2 = hist_1D(tree, todraw, x_bins, cut2)
        hist2.SetLineColor(ROOT.kBlue)
        hist2.SetLineWidth(4)
        #hist2.Draw("sameE3") 
	hs.Add(hist2)
   
        hist3 = hist_1D(tree, todraw, x_bins, cut3)
        hist3.SetLineColor(ROOT.kGreen)
        hist3.SetLineWidth(4)
        #hist3.Draw("sameE3") 
	hs.Add(hist3)

        hist4 = hist_1D(tree, todrawtrue, x_bins, cut4)
        hist4.SetFillColor(ROOT.kMagenta)
        max = hist1.GetMaximum()
        bin_max = hist4.GetMaximumBin()
	hist4.SetBinContent(bin_max, max/4.0)
#	hist4.SetBinContent(bin_max-1, max/4.0)
	#hist4.SetBinContent(bin_max+1, max/4.0)
        hs.Add(hist4)
     
    	c = ROOT.TCanvas() 
    	hs.Draw("nostack")
    	c.SaveAs("gethiststack_%s"%todraw+"_%s"%tree.GetTitle()+"_test.png")
        return hs   
     #   c1.SaveAs("test_%d"%(n%9)+".png") 


#_______________________________________________________________________________
if __name__ == "__main__":
    
    treename = "mmctree_4204" 
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs-1M-B3-1071409.root"
    file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_100k_test_0331_B3.root"
    dir = "DiHiggsWWAna/%s"%treename
    #dir = "DiHiggsWWAna/"
    
    #test(file, dir)
    monitoringMMC(file)  
    title1 = "MMC PDF for M_{H}, Event 4204"
    h2massbins = "(30,280,370)"#for843 only
    
    pic_h2mass = "h2Mass_combined_0325_test"
    #drawh2Mass_combined(file, h2massbins, pic_h2mass)

    h2mass = "M_{H}, [GeV]"
    #cut1 = "control<2"
    tex1 = "correctly pair (muon candidates, muon lorentzVec in MMC)"
    pic1 = "evt_4204_h2Mass_0325_test"
 #   draw1D(file, dir, h2tohhmass, title1, h2massbins, h2mass, tex1, pic1)

    h2massbins_2 = "(70,270,340)"#for843 only
    cut2 = "control>1"
    tex2 = "incorrectly pair (muon candidates, muon lorentzVec in MMC)"
    pic2 = "evt843_h2Mass_incorrect_0325_test"
   # draw1D(file, dir, h2tohhmass, title1, h2massbins_2, h2mass, cut2, tex2, pic2)

    h2massbins_3 = "(70,270,340)"#for843 only
    cut3 = ""
    tex3 = "inclduing both pairs (muon candidates, muon lorentzVec in MMC)"
    pic3 = "evt843_h2Mass_both_0325_test"
  #  draw1D(file, dir, h2tohhmass, title1, h2massbins_3, h2mass, cut3, tex3, pic3)
    #draw2D(file,dir,"nu_onshellW_phi:nu_onshellW_eta","eta_nuonshellW_true","phi_nuonshellW_true","#phi Vs #eta of nu_onshell","(50,-6,6)","#eta","(30,-3.14,3.14)","#phi","tex","etaVsphi_nuonshellW_0325_test")
    #draw2D(file,dir,"nu_offshellW_phi:nu_offshellW_eta","eta_nuoffshellW_true","phi_nuoffshellW_true","#phi Vs #eta of nu_offshell","(50,-6,6)","#eta","(30,-3.14,3.14)","#phi","tex","etaVsphi_nuoffshellW_0325_test")
    """ 
    offshellWmass = "offshellW_Mass"
    titleW1 = "off-shell W mass from MMC, Event843 "

    offshellWmassbins = "(50,0,50)"
    pic_offshellW = "offshellW_combined_0325_test"
    drawoffshellWmass_combined(file, offshellWmassbins, pic_offshellW)
    
    offshellnu_eta = "nu_offshellW_eta"
    offshellnu_eta_true = "eta_nuoffshellW_true"
    nuetabins= "(100,-7,7)"
    pic_offshellnu_eta = "offshellW_nueta_0325_test"
    tex_offshellnu_eta = "ensembles of allowable solutions for #eta_{#nu}^{offshellW}"
    draw_combined(file, offshellnu_eta, offshellnu_eta_true, nuetabins, tex_offshellnu_eta, pic_offshellnu_eta)

    onshellnu_eta = "nu_onshellW_eta"
    onshellnu_eta_true = "eta_nuonshellW_true"
    pic_onshellnu_eta = "onshellW_nueta_0325_test"
    tex_onshellnu_eta = "ensembles of allowable solutions for #eta_{#nu}^{onshellW}"
    draw_combined(file, onshellnu_eta, onshellnu_eta_true, nuetabins, tex_onshellnu_eta, pic_onshellnu_eta)
 
    offshellnu_pt = "nu_offshellW_pt"
    offshellnu_pt_true = "pt_nuoffshellW_true"
    nuptbins = "(60,0,120)"
    pic_offshellnu_pt = "offshellW_nupt_0325_test"
    tex_offshellnu_pt = "ensembles of allowable solutions for p_{T}_{#nu}^{offshellW}"
    draw_combined(file, offshellnu_pt, offshellnu_pt_true, nuptbins, tex_offshellnu_pt, pic_offshellnu_pt)
    
    onshellnu_pt = "nu_onshellW_pt"
    onshellnu_pt_true = "pt_nuonshellW_true"
    nuptbins = "(60,0,120)"
    pic_onshellnu_pt = "onshellW_nupt_0325_test"
    tex_onshellnu_pt = "ensembles of allowable solutions for p_{T}_{#nu}^{onshellW}"
    draw_combined(file, onshellnu_pt, onshellnu_pt_true, nuptbins, tex_onshellnu_pt, pic_onshellnu_pt)
    

    offshellWmass_x = "off-shell W mass"
    cutW1 = "control<2"
    texW1 = "correctly pair (muon candidates, muon lorentzVec in MMC)"
    picW1 = "evt843_offshellWmass_correct_0325_test"
  #  draw1D(file, dir, offshellWmass, titleW1, offshellWmassbins, offshellWmass_x, cutW1, texW1, picW1)

    cutW2 = "control>2"
    texW2 = "incorrectly pair (muon candidates, muon lorentzVec in MMC)"
    picW2 = "evt843_offshellWmass_incorrect_0325_test"
  #  draw1D(file, dir, offshellWmass, titleW1, offshellWmassbins, offshellWmass_x, cutW2, texW2, picW2)

    cutW3 = ""
    texW3 = "incluing both pairs (muon candidates, muon lorentzVec in MMC)"
    picW3 = "evt843_offshellWmass_both_0325_test"
   # draw1D(file, dir, offshellWmass, titleW1, offshellWmassbins, offshellWmass_x, cutW1, texW3, picW3)

    """
