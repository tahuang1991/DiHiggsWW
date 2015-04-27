import ROOT
import os


#inputprefix = "/eos/uscms/store/user/tahuang/store/user/tahuang/"
#outputprefix = "/eos/uscms/store/user/tahuang/SLHC25_GEMCSC_bendingangle_1M_2023Muon/"
#________________________________________________
def mergerootfiles(inputDir, outputfile):
        if not os.path.isdir(inputDir):
            print "ERROR: This is not a valid directory: ", inputDir
	print "inputdir ", inputDir[:]
	tmp = inputDir[:]+"tmp.root"
#	os.system("cd %s"%inputDir[:])
#	os.system("pwd")
	m = 0
	#outputfile = outp+outputfile
	if os.path.isfile(outputfile):
		os.system("rm %s"%outputfile)
        ls = os.listdir(inputDir)
	for x in ls:
		if m==0:
			x = inputDir[:]+x
			print "x ",x
			os.system("cp %s"%x+" %s"%tmp)
			print "m=0, tmp ",tmp," entries ",getEntries(x)
		if m!=0:
			x = inputDir[:]+x
			os.system("hadd -f %s"%outputfile+" %s"%tmp+" %s"%x)
			os.system("cp %s"%outputfile+" %s"%tmp)
			print "m= ",m," entries ", getEntries(x)," tmp root file enetries ",getEntries(tmp)
		m = m+1

	os.system("rm %s"%tmp)
#_________________________________________________
def getEntries(file):

     f = ROOT.TFile(file)
     tree = f.Get("htoWWAna/evtree")
     entries = tree.GetEntries()
     f.Close()
     return entries
#_______________________________________________________________________________
if __name__ == "__main__":

   inputDir = "/fdata/hepx/store/user/taohuang/Hhh/htoWWAna/"
   outputfile = "/fdata/hepx/store/user/taohuang/Hhh/htoWWAna/htoWWAna-1M-0413-mediateStates-B3-combined.root"
   mergerootfiles(inputDir, outputfile)
 #  inputDirs = ['SLHC23_patch1_2023Muon_1M_Ana_PU0_Pt7_V2/','SLHC23_patch1_2023Muon_1M_Ana_PU0_Pt10_V2/','SLHC23_patch1_2023Muon_1M_Ana_PU0_Pt15_V2/','SLHC23_patch1_2023Muon_1M_Ana_PU0_Pt30_V2/','SLHC23_patch1_2023Muon_1M_Ana_PU0_Pt40_V2/']
  # outputfiles = ['gem-csc_ana_Pt7','gem-csc_ana_Pt10','gem-csc_ana_Pt15','gem-csc_ana_Pt30','gem-csc_ana_Pt40'] 
#   for d in range(len(inputDirs)):
#	input = inputprefix+inputDirs[d]
#	output = outputprefix+outputfiles[d]+".root"
#	print "input ",input, " output ",output
#	mergerootfiles(input, output)
