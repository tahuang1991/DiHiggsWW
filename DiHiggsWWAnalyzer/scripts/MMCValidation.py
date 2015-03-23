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


#used to scale hist
count = 10000.0

h2tohhmass = "h2tohh_Mass"
h2truemass = "h2_Mass_FromGen"
offshellWmass = "offshellW_Mass"
offshellWtruemass = "offshellW_Mass_FromGen"
#___________________________________________
def draw1D(file,dir,todraw,title,x_bins,x_title,cut,tex,pic_name):
    
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
    b1.SetTitle("%s"%title+" "*8 + "CMS Simulation Preliminary")
    b1.GetYaxis().SetTitle("Events")
    b1.GetXaxis().SetTitle("%s"%x_title)
    #b1.SetStats(0)
    #print "todraw ",todraw 
    t.Draw(todraw+">>b1",cut)
    b1.Draw()
    tex2 = ROOT.TLatex(0.15,0.6,"%s "%tex)
    tex2.SetTextSize(0.04)
    tex2.SetTextFont(42)
    tex2.SetNDC()
    tex2.Draw("same")
	
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
    c1.Divide(3,3)
    n = 0;
    #test = sub_key.ReadObj()
    #testhist = hist_1D(test, h2tohhmass, x_bins, cut1)
    #testhist.Draw()
    #c1.SaveAs("test.png")
    hslist = []
    while sub_key: 
	print " n", n,"n%9", n%9 
    	if n%9 == 0:
                c1.Clear
		print  "c1 divided "
        c1.cd(n%9+1)
	tree = sub_key.ReadObj()
        #tree.Print()
        hs = ROOT.THStack("hs_%d"%(n%9),"hs_%d"%(n%9));
        hist1 = hist_1D(tree, h2tohhmass, x_bins, cut1)
        hist1.SetFillColor(ROOT.kRed)
        #hist1.Draw("E3") 
	hs.Add(hist1)
        
        hist2 = hist_1D(tree, h2tohhmass, x_bins, cut2)
        hist2.SetFillColor(ROOT.kBlue)
        #hist2.Draw("sameE3") 
	hs.Add(hist2)
   
        hist3 = hist_1D(tree, h2tohhmass, x_bins, cut3)
        hist3.SetFillColor(ROOT.kGreen)
        #hist3.Draw("sameE3") 
	hs.Add(hist3)

        hist4 = hist_1D(tree, h2truemass, x_bins, cut4)
        hist4.SetFillColor(ROOT.kMagenta)
        hist4.Scale(0.1)
        hs.Add(hist4)
        
        hslist.append(hs)
     	hs.Draw("nostack")
     #   c1.SaveAs("test_%d"%(n%9)+".png") 
	if n%9 == 8:
		c1.cd()
		c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".pdf")
    		c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".png")
        n = n+1
        sub_key = sub_list.At(n+1)
    if n%9 != 0:
        c1.cd()
	c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".pdf")
    	c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".png")



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
    #test = sub_key.ReadObj()
    #testhist = hist_1D(test, h2tohhmass, x_bins, cut1)
    #testhist.Draw()
    #c1.SaveAs("test.png")
    hslist = []
    while sub_key: 
	print " n", n,"n%9", n%9 
    	if n%9 == 0:
                c1.Clear
		print  "c1 divided "
        c1.cd(n%9+1)
	tree = sub_key.ReadObj()
        #tree.Print()
        hs = ROOT.THStack("hs_%d"%(n%9),"hs_%d"%(n%9));
        hist1 = hist_1D(tree, offshellWmass, x_bins, cut1)
        hist1.SetFillColor(ROOT.kRed)
        #hist1.Draw("E3") 
	hs.Add(hist1)
        
        hist2 = hist_1D(tree, offshellWmass, x_bins, cut2)
        hist2.SetFillColor(ROOT.kBlue)
        #hist2.Draw("sameE3") 
	hs.Add(hist2)
   
        hist3 = hist_1D(tree, offshellWmass, x_bins, cut3)
        hist3.SetFillColor(ROOT.kGreen)
        #hist3.Draw("sameE3") 
	hs.Add(hist3)

        hist4 = hist_1D(tree, offshellWtruemass, x_bins, cut4)
        hist4.SetFillColor(ROOT.kMagenta)
        hist4.Scale(0.1)
       # hs.Add(hist4)
        
        hslist.append(hs)
     	hs.Draw("nostack")
     #   c1.SaveAs("test_%d"%(n%9)+".png") 
	if n%9 == 8:
		c1.cd()
		c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".pdf")
    		c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".png")
        n = n+1
        sub_key = sub_list.At(n+1)
    if n%9 != 0:
        c1.cd()
	c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".pdf")
    	c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".png")

#___________________________________________
def draw_combined(file,todraw,x_bins,pic_name):
    
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
    #test = sub_key.ReadObj()
    #testhist = hist_1D(test, h2tohhmass, x_bins, cut1)
    #testhist.Draw()
    #c1.SaveAs("test.png")
    hslist = []
    while sub_key: 
	print " n", n,"n%9", n%9 
    	if n%9 == 0:
                c1.Clear
		print  "c1 divided "
        c1.cd(n%9+1)
	tree = sub_key.ReadObj()
        #tree.Print()
        hs = ROOT.THStack("hs_%d"%(n%9),"hs_%d"%(n%9));
        hist1 = hist_1D(tree, todraw, x_bins, cut1)
        hist1.SetFillColor(ROOT.kRed)
        #hist1.Draw("E3") 
	hs.Add(hist1)
        
        hist2 = hist_1D(tree, todraw, x_bins, cut2)
        hist2.SetFillColor(ROOT.kBlue)
        #hist2.Draw("sameE3") 
	hs.Add(hist2)
   
        hist3 = hist_1D(tree, todraw, x_bins, cut3)
        hist3.SetFillColor(ROOT.kGreen)
        #hist3.Draw("sameE3") 
	hs.Add(hist3)

        hist4 = hist_1D(tree, todraw, x_bins, cut4)
        hist4.SetFillColor(ROOT.kMagenta)
        hist4.Scale(0.1)
       # hs.Add(hist4)
        
        hslist.append(hs)
     	hs.Draw("nostack")
     #   c1.SaveAs("test_%d"%(n%9)+".png") 
	if n%9 == 8:
		c1.cd()
		c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".pdf")
    		c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".png")
        n = n+1
        sub_key = sub_list.At(n+1)
    if n%9 != 0:
        c1.cd()
	c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".pdf")
    	c1.SaveAs("Dihiggs_MMC_%s"%pic_name+"_B3_%d"%(n/9)+".png")

        


        



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
	tree.Print()
        i = i+1
        sub_key = sub_list.At(i)
    print "key getclassname ", key.GetClassName()
    print "c1 ",c1
    print "obj  ", obj
    print "sub_key ", sub_key
    print "sub_list 1", sub_list.At(1)
    print "dir ", dir, "  t ", t, "files ", 


#______________________________________________________________________________
def hist_1D(tree, todraw, x_bins, cut_weight):
#constrol >1 , correct pairing 
#cut_weight should either "weight && (control<2)" or "weight" or "weight && (control>1)"
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    hist = ROOT.TH1F("hist","hist",xBins,xminBin,xmaxBin)
    tree.Draw(todraw+">>hist",cut_weight)
    hist.Scale(1.0/count)
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


#_______________________________________________________________________________
if __name__ == "__main__":
    
    treename = "mmctree_843" 
    #file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs-1M-B3-1071409.root"
    file = "/fdata/hepx/store/user/taohuang/Hhh/DiHiggs_10k_0309_B3.root"
    dir = "DiHiggsWWAna/%s"%treename
    #dir = "DiHiggsWWAna/"
    
    #test(file, dir)
      
    title1 = "heavy higgs mass from MMC, Event843 "
    h2massbins = "(30,280,370)"#for843 only
    
    pic_h2mass = "h2Mass_combined_0310"
    #drawh2Mass_combined(file, h2massbins, pic_h2mass)

    h2mass = "heavy higgs mass"
    cut1 = "control<2"
    tex1 = "correctly pair (muon candidates, muon lorentzVec in MMC)"
    pic1 = "evt843_h2Mass_correct_0310"
#    draw1D(file, dir, h2tohhmass, title1, h2massbins, h2mass, cut1, tex1, pic1)

    h2massbins_2 = "(70,270,340)"#for843 only
    cut2 = "control>1"
    tex2 = "incorrectly pair (muon candidates, muon lorentzVec in MMC)"
    pic2 = "evt843_h2Mass_incorrect_0310"
#    draw1D(file, dir, h2tohhmass, title1, h2massbins_2, h2mass, cut2, tex2, pic2)

    h2massbins_3 = "(70,270,340)"#for843 only
    cut3 = ""
    tex3 = "inclduing both pairs (muon candidates, muon lorentzVec in MMC)"
    pic3 = "evt843_h2Mass_both_0310"
#    draw1D(file, dir, h2tohhmass, title1, h2massbins_3, h2mass, cut3, tex3, pic3)

    
    offshellWmass = "offshellW_Mass"
    titleW1 = "off-shell W mass from MMC, Event843 "
    offshellWmassbins = "(50,0,50)"
    pic_offshellW = "offshellW_combined_0310"
    #drawoffshellWmass_combined(file, offshellWmassbins, pic_offshellW)
    
    offshellnu_eta = "nu_offshellW_eta"
    nuetabins= "(100,-7,7)"
    pic_offshellnu_eta = "offshellW_nueta_0310"
   # draw_combined(file, offshellnu_eta, nuetabins, pic_offshellnu_eta)
    onshellnu_eta = "nu_onshellW_eta"
    pic_onshellnu_eta = "onshellW_nueta_0310"
  #  draw_combined(file, onshellnu_eta, nuetabins, pic_onshellnu_eta)
 
    offshellnu_pt = "nu_offshellW_pt"
    nuptbins = "(60,0,120)"
    pic_offshellnu_pt = "offshellW_nupt_0319"
    draw_combined(file, offshellnu_pt, nuptbins, pic_offshellnu_pt)
    
    onshellnu_pt = "nu_onshellW_pt"
    nuptbins = "(60,0,120)"
    pic_onshellnu_pt = "onshellW_nupt_0319"
    draw_combined(file, onshellnu_pt, nuptbins, pic_onshellnu_pt)
    

    offshellWmass_x = "off-shell W mass"
    cutW1 = "control<2"
    texW1 = "correctly pair (muon candidates, muon lorentzVec in MMC)"
    picW1 = "evt843_offshellWmass_correct_0310"
  #  draw1D(file, dir, offshellWmass, titleW1, offshellWmassbins, offshellWmass_x, cutW1, texW1, picW1)

    cutW2 = "control>2"
    texW2 = "incorrectly pair (muon candidates, muon lorentzVec in MMC)"
    picW2 = "evt843_offshellWmass_incorrect_0310"
   # draw1D(file, dir, offshellWmass, titleW1, offshellWmassbins, offshellWmass_x, cutW2, texW2, picW2)

    cutW3 = ""
    texW3 = "incluing both pairs (muon candidates, muon lorentzVec in MMC)"
    picW3 = "evt843_offshellWmass_both_0310"
    #draw1D(file, dir, offshellWmass, titleW1, offshellWmassbins, offshellWmass_x, cutW1, texW3, picW3)


