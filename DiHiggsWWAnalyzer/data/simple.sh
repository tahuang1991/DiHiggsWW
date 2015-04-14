#!/bin/csh -f
#before submiting jobs, setup environment
source ~/.cshrc

set jobid=testfile
echo "jobid is $jobid"
#cd /uscms_data/d3/tahuang/CMSSW_7_2_0/src/DiHiggsWW/DiHiggsWWAnalyzer/data/
#cd ${_CONDOR_SCRATCH_DIR}
#cp /uscms_data/d3/tahuang/CMSSW_7_2_0/src/DiHiggsWW/DiHiggsWWAnalyzer/test/addGenParticles.py addGenParticles_$jobid.py
#chmod 775 addGenParticles_$jobid.py
#sed -i ""
##input file, to avoid race issue, copy input files to some dirs?? what if input file is too large?
#sed -i "s/HH-bbWW-B3\_Gen\_100k\_0215.root/$rootfile/g" runDiHiggsWWAna_$jobid_$taskid.py
#sed output filename
#sed -i "s/DiHiggs_1_618_0406_B3.root/DiHiggs-100k-0406-mediateStates-$version-$jobid.root/g" runDiHiggsWWAna_$jobid.py
#source /cvmfs/cms.cern.ch/cmsset_default.csh
#eval `scramv1 runtime -sh`
#cmsRun addGenParticles_10.py
#echo "_CONDOR_SCRATCH_DIR is ${_CONDOR_SCRATCH_DIR}"
echo "this is testscript for condor"
