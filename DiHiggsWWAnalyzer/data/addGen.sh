#!/bin/csh -f

source ~/.cshrc

set jobid=$1
echo "jobid is $jobid"
cd /uscms_data/d3/tahuang/CMSSW_7_2_0/src/DiHiggsWW/DiHiggsWWAnalyzer/data/
cp /uscms_data/d3/tahuang/CMSSW_7_2_0/src/DiHiggsWW/DiHiggsWWAnalyzer/test/addGenParticles.py addGenParticles_$jobid.py
chmod 775 addGenParticles_$jobid.py
#sed -i ""
##input file, to avoid race issue, copy input files to some dirs?? what if input file is too large?
#sed -i "s/HH-bbWW-B3\_Gen\_100k\_0215.root/$rootfile/g" runDiHiggsWWAna_$jobid_$taskid.py
#sed output filename
#sed -i "s/DiHiggs_1_618_0406_B3.root/DiHiggs-100k-0406-mediateStates-$version-$jobid.root/g" runDiHiggsWWAna_$jobid.py
source /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scramv1 runtime -sh`
cmsRun addGenParticles_$jobid.py

echo "test"
