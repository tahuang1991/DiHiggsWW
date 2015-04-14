// -*- C++ -*-
//
// Package:    ttbarAnalyzer
// Class:      ttbarAnalyzer
// 
/**\class ttbarAnalyzer ttbarAnalyzer.cc DiHiggsWW/ttbarAnalyzer/plugins/ttbarAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  tao huang
//         Created:  Wed, 26 Nov 2014 17:58:07 GMT
// $Id$
//
//


// system include files
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "math.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
//headers from root lib
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH1F.h"



#include "DiHiggsWW/DiHiggsWWAnalyzer/src/MMC.h"


class MMC;
//#define WMass 80.385   // W mass
//#define SMHMass 125.03 // SM module higgs mass
//t(6)->W^+(24) b(5)
//tbar(-6)->W^- bbar(-5)

//
// class declaration
//

using namespace reco;

typedef std::pair<float, float> EtaPhi;

class ttbarAnalyzer : public edm::EDAnalyzer {
   public:
      explicit ttbarAnalyzer(const edm::ParameterSet&);
      ~ttbarAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
    
   //runMMC
   /*private: 
      void runMMC();
      void initTree(TTree* mmctree);
        
     
      float genEtaGuass(float mean, float rms);
      float genPhiFlat();
      EtaPhi generatenu1_etaphi();
      float nu1pt_onshellW(EtaPhi nu1_etaphi, TLorentzVector* mu1lorentz, float wMass);
      bool  nulorentz_offshellW(TLorentzVector* jetlorentz, TLorentzVector* mu1lorentz, 
			       TLorentzVector* mu2lorentz, TLorentzVector* nu1lorentz, 
 			       TLorentzVector* nu2lorentz, int control, float hMass);
      bool checkSolution(TLorentzVector* jetslorentz,
                          TLorentzVector* mu1lorentz,
                          TLorentzVector* mu2lorentz,
                          TLorentzVector* nu1lorentz, int control, float hMass); 
      bool cutsCheck();
      void assignMuLorentzVec(int control);  
        */  
   private:
      //edm::ParameterSet cfg_;
      edm::ParameterSet mmcset_;
      // debuglevel constrol 
      int verbose_; 
      void print();
     // void printHtoWWChain();
      void printCandidate(const reco::Candidate* );
      void printallDecendants(const reco::Candidate* );
      void printallAncestors(const reco::Candidate* );
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
   private:
      const reco::Candidate* findmudaughter(const reco::Candidate* );
      const reco::Candidate* findnudaughter(const reco::Candidate* );
      const reco::Candidate* findmudescendants(const reco::Candidate*, int& );
      const reco::Candidate* findnudescendants(const reco::Candidate*, int& );
    private:
     bool finalStates_;
    private:
     //decendants and ancestor
     const reco::Candidate* stabledecendant(const reco::Candidate* cand, int id);
     //const reco::Candidate* stabletdecendant(Particle p, PdgId id);
     const reco::Candidate* alldecendant(const reco::Candidate* cand, int id, bool first=false);
     const reco::Candidate* allancestor(const reco::Candidate* cand, int id, bool first=false);
     bool hasMother(const reco::Candidate* cand, int id);
     bool hasDaughter(const reco::Candidate* cand, int id);
      // ---------- Candidates in signal channel ---------------------------

      // ---------- Candidates in signal channel ---------------------------
      const reco::Candidate* mu1_W1_cand;
      const reco::Candidate* nu1_W1_cand;
      const reco::Candidate* mu2_W2_cand;
      const reco::Candidate* nu2_W2_cand;
      const reco::Candidate* mu1_tbar_cand;
      const reco::Candidate* mu2_t_cand;
      const reco::Candidate* b1_t_cand;
      const reco::Candidate* b2_tbar_cand;

      const reco::Candidate* mu1cand;
      const reco::Candidate* nu1cand;
      const reco::Candidate* mu2cand;
      const reco::Candidate* nu2cand;
      const reco::Candidate* b1cand;
      const reco::Candidate* b2cand;
      const reco::Candidate* w1cand;
      const reco::Candidate* w2cand;
      const reco::Candidate* tcand;
      const reco::Candidate* tbarcand;

    private:
      
      TLorentzVector mu1_lorentz;
      TLorentzVector mu2_lorentz;
      TLorentzVector bbar_lorentz;
      TLorentzVector nu1_lorentz;
      TLorentzVector nu2_lorentz;
      TLorentzVector met_lorentz;

    private:
      TLorentzVector stableDecendantsLorentz(const reco::Candidate* cand); 
   
    private:
      TLorentzVector calculateMET(); 

    private:
      TTree *evtree;
      TFile *output;
    private:
      void fillbranches(); 
      edm::Service< TFileService > fs;
      
      //----------branches of tree ---------------------------
      int ievent;
      float mu1_energy;
      float mu1_px;
      float mu1_py;
      float mu1_pz;
      float mu1_eta;
      float mu1_phi;
      int mu1_motherid;
      float mu1_mother_energy;
      float mu1_mother_px;
      float mu1_mother_py;
      float mu1_mother_pz;
      float mu1_mother_mass;
      float nu1_energy;
      float nu1_px;
      float nu1_py;
      float nu1_pz;
      float nu1_eta;
      float nu1_phi;
      bool Wtomu1nu1;

      float mu2_energy;
      float mu2_px;
      float mu2_py;
      float mu2_pz;
      float mu2_eta;
      float mu2_phi;
      int mu2_motherid;
      float mu2_mother_energy;
      float mu2_mother_px;
      float mu2_mother_py;
      float mu2_mother_pz;
      float mu2_mother_mass;
      float nu2_energy;
      float nu2_px;
      float nu2_py;
      float nu2_pz;
      float nu2_eta;
      float nu2_phi;
      bool Wtomu2nu2;
     
    //  float w1_mass;
    //  float w2_mass;
      float b1_energy;
      float b1_px;
      float b1_py;
      float b1_pz;
      int b1_motherid;
      float b2_energy;
      float b2_px;
      float b2_py;
      float b2_pz;
      int b2_motherid;
      float jets_energy;
      float jets_px;
      float jets_py;
      float jets_pz;
      float jets_mass;
      TLorentzVector jets_lorentz;

      float tbar_energy;
      float tbar_px;
      float tbar_py;
      float tbar_pz;
      float tbar_mass;
      float t_energy;
      float t_px;
      float t_py;
      float t_pz; 
      float t_mass;
      
      float met;
      float met_phi;
      float met_px;
      float met_py;
      //cuts for higgstoWWbb
      bool mu_positive;
      bool mu_negative;
      bool nu_positive;
      bool nu_negative;
      bool bquark;
      bool bbarquark;
      bool tbartoWbbar;
      bool ttoWb;

        
    private:
      bool runMMC_;
      MMC* thismmc;
    /*
     // MMC tree branches
    private:
      float onshellWMassRandomWalk(float x0, float step, float random);
      float onshellWMassRandomWalk(float x0, float step, float random, TH1F* hist);
      float onshellWMassPDF(float wmass);
  
    private:
      TH1F* readoutonshellWMassPDF();
      TH1F* readoutoffshellWMassPDF();
      TH1F* readoutonshellnuptPDF();
 
    private:
      float weightfromhist(TH1F* pdf, float x); 
      float weightfromonshellnupt(float nupt); 
   
    private:
      bool weightfromonshellnupt_func_;
      bool weightfromonshellnupt_hist_;
      bool weightfromoffshellWmass_hist_;

    private:
      int iterations_;
      int seed_;
      std::string RefPDFfile_;
      float eta_mean;
      float eta_rms;
      float eta_gen; 
      float phi_gen;
      float wmass_gen;
      float hmass_gen;
      TLorentzVector* mu_onshellW_lorentz;
      TLorentzVector* mu_offshellW_lorentz;
      TVector2* MMCmet_vec2;
      TLorentzVector* nu_onshellW_lorentz;
      TLorentzVector* nu_offshellW_lorentz;
      TLorentzVector* offshellW_lorentz;
      TLorentzVector* onshellW_lorentz;
      TLorentzVector* t_lorentz;
      TLorentzVector* tbar_lorentz;
      TLorentzVector* h2tohh_lorentz;
      
      int control;
      float weight;
      float weight1;//extra weight
      float weight2;//extra weight
 
      float mu_onshellW_Eta;
      float mu_onshellW_Phi;
      float mu_onshellW_Pt;
      float mu_onshellW_E;
      float mu_offshellW_Eta;
      float mu_offshellW_Phi;
      float mu_offshellW_Pt;
      float mu_offshellW_E;
      float nu_onshellW_Eta;
      float nu_onshellW_Phi;
      float nu_onshellW_Pt;
      float nu_onshellW_E;
      float nu_offshellW_Eta;
      float nu_offshellW_Phi;
      float nu_offshellW_Pt;
      float nu_offshellW_E;
      
      float onshellW_Eta;
      float onshellW_Phi;
      float onshellW_Pt;
      float onshellW_E;
      float onshellW_Mass;
      float offshellW_Eta;
      float offshellW_Phi;
      float offshellW_Pt;
      float offshellW_E;
      float offshellW_Mass;
     
      float tbar_Eta;
      float tbar_Phi;
      float tbar_Pt;
      float tbar_E;
      float tbar_Mass;
      float t_Eta;
      float t_Phi;
      float t_Pt;
      float t_E;
      float t_Mass;

      float MMCmet_E;
      float MMCmet_Phi;
      float MMCmet_Px;
      float MMCmet_Py;

      float h2tohh_Eta;
      float h2tohh_Phi;
      float h2tohh_Pt;
      float h2tohh_E;
      float h2tohh_Mass;


      float eta_nuoffshellW_true;
      float phi_nuoffshellW_true;
      float pt_nuoffshellW_true;
      float eta_nuonshellW_true;
      float phi_nuonshellW_true;
      float pt_nuonshellW_true;
      float mass_offshellW_true;
      float mass_onshellW_true;
      float pt_h2tohh_true;
*/
};

    //
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ttbarAnalyzer::ttbarAnalyzer(const edm::ParameterSet& iConfig)
{
     verbose_ = iConfig.getUntrackedParameter<int>("verbose",0);
     mmcset_ = iConfig.getParameter<edm::ParameterSet>("mmcset"); 
     finalStates_ = iConfig.getParameter<bool>("finalStates");
     runMMC_ = iConfig.getParameter<bool>("runMMC");
     /*
     weightfromonshellnupt_func_ = iConfig.getParameter<bool>("weightfromonshellnupt_func");
     weightfromonshellnupt_hist_ = iConfig.getParameter<bool>("weightfromonshellnupt_hist");
     weightfromoffshellWmass_hist_ = iConfig.getParameter<bool>("weightfromoffshellWmass_hist");
     iterations_ = iConfig.getUntrackedParameter<int>("iterations",100000);
     seed_ = iConfig.getParameter<int>("seed");
     RefPDFfile_ = iConfig.getParameter<std::string>("RefPDFfile");
     std::cout <<" RefPDFfile_ " << RefPDFfile_ << std::endl;
     */
   // initilize candidates pointer



   //now do what ever initialization is needed
      ievent = 0;
      mu1_W1_cand = NULL;
      nu1_W1_cand = NULL;
      mu2_W2_cand = NULL;
      nu2_W2_cand = NULL;
      mu1_tbar_cand = NULL;
      mu2_t_cand = NULL;
      b1_t_cand = NULL;
      b2_tbar_cand = NULL;
   //  
      mu1_energy = 0.0;
      mu1_px = 0.0;
      mu1_py = 0.0;
      mu1_pz = 0.0;
      mu1_motherid = 0;
      mu1_mother_energy = 0.0;
      mu1_mother_px = 0.0;
      mu1_mother_py = 0.0;
      mu1_mother_pz = 0.0;
      mu1_mother_mass = 0.0;
      nu1_energy = 0.0;
      nu1_px = 0.0;
      nu1_py = 0.0;
      nu1_pz = 0.0;
      Wtomu1nu1 = false;

      mu2_energy = 0.0;
      mu2_px = 0.0;
      mu2_py = 0.0;
      mu2_pz = 0.0;
      mu2_motherid = 0;
      mu2_mother_energy = 0.0;
      mu2_mother_px = 0.0;
      mu2_mother_py = 0.0;
      mu2_mother_pz = 0.0;
      mu2_mother_mass = 0.0;
      nu2_energy = 0.0;
      nu2_px = 0.0;
      nu2_py = 0.0;
      nu2_pz = 0.0;
      Wtomu2nu2 = false;

      t_energy = 0.0;
      t_px = 0.0;
      t_py = 0.0;
      t_pz = 0.0;
      t_mass = 0.0;

      b1_energy = 0.0;
      b1_px = 0.0;
      b1_py = 0.0;
      b1_pz = 0.0;
      b1_motherid = 0;
      b2_energy = 0.0;
      b2_px = 0.0;
      b2_py = 0.0;
      b2_pz = 0.0;
      b2_motherid = 0;
      jets_energy = 0.0;
      jets_px = 0.0;
      jets_py = 0.0;
      jets_pz = 0.0;
      jets_mass = 0.0;
      //jets_lorentz = new TLorentzVector();

      tbar_energy = 0.0;
      tbar_px = 0.0;
      tbar_py = 0.0;
      tbar_pz = 0.0;
      tbar_mass = 0.0;
      

      mu_positive = false;
      mu_negative = false;
      bquark = false;
      bbarquark = false;
      tbartoWbbar = false;
      ttoWb = false;

     /* 
      mu1_lorentz = new TLorentzVector();
      nu1_lorentz = new TLorentzVector();
      mu2_lorentz = new TLorentzVector();
      nu2_lorentz = new TLorentzVector();
      bbar_lorentz = new TLorentzVector();
      met_lorentz = new TLorentzVector();
     // runMMC
      mu_onshellW_lorentz = new TLorentzVector();
      mu_offshellW_lorentz = new TLorentzVector();
      MMCmet_vec2 = new TVector2();
      nu_onshellW_lorentz = new TLorentzVector();
      nu_offshellW_lorentz = new TLorentzVector();
      offshellW_lorentz = new TLorentzVector();
      onshellW_lorentz = new TLorentzVector();
      t_lorentz = new TLorentzVector();
      tbar_lorentz = new TLorentzVector();
      h2tohh_lorentz = new TLorentzVector();

      mu_onshellW_Eta =0;
      mu_onshellW_Phi =0;
      mu_onshellW_Pt =0;
      mu_onshellW_E =0;
      mu_offshellW_Eta =0;
      mu_offshellW_Phi =0;
      mu_offshellW_Pt =0;
      mu_offshellW_E =0;
      nu_onshellW_Eta =0;
      nu_onshellW_Phi =0;
      nu_onshellW_Pt =0;
      nu_onshellW_E =0 ;
      nu_offshellW_Eta =0;
      nu_offshellW_Phi =0; 
      nu_offshellW_Pt =0;
      nu_offshellW_E =0;
      
      onshellW_Eta =0;
      onshellW_Phi =0;
      onshellW_Pt =0;
      onshellW_E =0;
      onshellW_Mass =0;
      offshellW_Eta =0;
      offshellW_Phi =0;
      offshellW_Pt =0;
      offshellW_E =0;
      offshellW_Mass =0;
     
      tbar_Eta =0;
      tbar_Phi =0;
      tbar_Pt =0;
      tbar_E =0;
      tbar_Mass =0;
      t_Eta =0;
      t_Phi =0;
      t_Pt =0;
      t_E =0;
      t_Mass =0;

      h2tohh_Eta =0;
      h2tohh_Phi =0;
      h2tohh_Pt =0;
      h2tohh_E =0;
      h2tohh_Mass =0;
      eta_nuoffshellW_true = 0.0;
      pt_nuoffshellW_true = 0.0;
      eta_nuonshellW_true = 0.0;
      pt_nuonshellW_true = 0.0;
      mass_offshellW_true = 0.0;
      mass_onshellW_true =0.0;
      pt_h2tohh_true = 0;
      
  */  

}


ttbarAnalyzer::~ttbarAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ttbarAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;



#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
//   if(iEvent.isRealData()) std::cout << " Not a real Data " << std::endl;
    
   std::cout << "event  " << iEvent.id().event() << std::endl;
   ievent = iEvent.id().event();
   Handle<reco::GenParticleCollection> genParticleColl;
   iEvent.getByLabel("prunedGenParticles", genParticleColl);
   Handle<std::vector<pat::Muon>> muonColl;
   iEvent.getByLabel("slimmedMuons",muonColl); 
   Handle<std::vector<pat::Jet>> jetColl;
   iEvent.getByLabel("slimmedJets",jetColl); 
   Handle<std::vector<pat::MET>> metColl;
   iEvent.getByLabel("slimmedMETs",metColl); 
      mu_positive = false;
      mu_negative = false;
      nu_positive = false;
      nu_negative = false;
      bquark = false;
      bbarquark = false;
      Wtomu1nu1 = false;
      Wtomu2nu2 = false;
      tbartoWbbar = false;
      ttoWb = false;
      
      
   for (std::vector<pat::Muon>::const_iterator it = muonColl->begin(); it != muonColl->end(); ++it){
   
       std::cout <<"muon id " << it->pdgId() <<" status " << it->status() <<" px " << it->px() << " py " << it->py() <<std::endl;
   
   }

   for (std::vector<pat::Jet>::const_iterator jet = jetColl->begin(); jet != jetColl->end(); ++jet){
   
       std::cout <<"jet mass " << jet->mass() <<" status " << jet->status() <<" px " << jet->px() << " py " << jet->py() <<std::endl;
   
   }

   for (std::vector<pat::MET>::const_iterator met = metColl->begin(); met != metColl->end(); ++met){
   
       std::cout <<"met mass " << met->mass() <<" status " << met->status() <<" px " << met->px() << " py " << met->py() <<std::endl;
   
   }

   for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {

//particle id, (muon13),(b5),(W+24),(SM higgs25)
   // particle id  it->pdgId()
   //
   //      std::cout << "Gen paticles: id " << it->pdgId() << std::endl; 
      if (it->pdgId() == 13 && it->status() == 1 && !mu_negative)
      {
	  //mu_negative = true;
          //std::cout << "find muon(-) with status 1" << std::endl;
          const reco::Candidate* tmp_mu1 = it->mother(); 
          while (tmp_mu1->pdgId() == 13 && tmp_mu1->numberOfMothers() == 1) tmp_mu1 = tmp_mu1->mother();
          if (tmp_mu1->numberOfMothers() != 1 ) std::cout << "muon has more than one mother particle" << std::endl;
          if (tmp_mu1->pdgId() == -24)  mu1_W1_cand = tmp_mu1;
           

          while (tmp_mu1->pdgId() == -24 && tmp_mu1->numberOfMothers() == 1) tmp_mu1 = tmp_mu1->mother();
          if (tmp_mu1->numberOfMothers() != 1 ) std::cout << "W- has more than one mother particle" << std::endl;
         
          if (tmp_mu1->pdgId() == -6)   
             {
                  //while (tmp_mu1->mother()->pdgId() == 25) tmp_mu1 = tmp_mu1->mother();
                  std::cout << "find muon(-) candidate" << std::endl;
	          mu_negative = true;
                  mu1_tbar_cand = tmp_mu1;
                 // std::cout << "mother of this higgs, id " << tmp_mu1->mother()->pdgId() << " energy " << tmp_mu1->mother()->energy() << std::endl;
              }
      }
      else if (it->pdgId() == -13 && it->status() == 1 && !mu_positive)
      {
         // std::cout << "find muon(+) with status 1" << std::endl;
	//  mu_positive = true;
          const reco::Candidate* tmp_mu2 = it->mother(); 
          while (tmp_mu2->pdgId() == -13 && tmp_mu2->numberOfMothers() == 1) tmp_mu2 = tmp_mu2->mother();
          if (tmp_mu2->numberOfMothers() != 1)  std::cout << "muon has more than one mother particle" << std::endl;
          if (tmp_mu2->pdgId() == 24)  mu2_W2_cand = tmp_mu2;
          while (tmp_mu2->pdgId() == 24 && tmp_mu2->numberOfMothers() == 1) tmp_mu2 = tmp_mu2->mother();
          if (tmp_mu2->numberOfMothers() != 1 ) std::cout << "W+ has more than one mother particle" << std::endl;
    //      const reco::Candidate* tmphiggs_mu2 = tmp_mu2->mother();
          if (tmp_mu2->pdgId() == 6)   
             {
                // while (tmp_mu2->mother()->pdgId() == 25)   tmp_mu2 = tmp_mu2->mother();
                 std::cout << "find muon(+) candidate" << std::endl;
	         mu_positive = true;
                 mu2_t_cand = tmp_mu2;
                // std::cout << "mother of this higgs, id " << tmp_mu2->mother()->pdgId() << " energy " << tmp_mu2->mother()->energy() << std::endl;
               }
        }
      else if (it->pdgId() == -14  && it->status() == 1 && !nu_negative )
      {
          const reco::Candidate* tmp_nu1 = it->mother();
          while (tmp_nu1->pdgId() == -14) tmp_nu1 = tmp_nu1->mother();
          if (tmp_nu1->pdgId() == -24)  nu1_W1_cand = tmp_nu1;
    //      std::cout << " the mother of nutrio" 
          while (tmp_nu1->pdgId() == -24) tmp_nu1 = tmp_nu1->mother();
          if (tmp_nu1->pdgId() == -6)
             {
            //     std::cout << "find nuetrino candidate" << std::endl;
                 nu_negative = true;
                }
         
       }
      else if (it->pdgId() == 14 && it->status() == 1 && !nu_positive )
      {
          const reco::Candidate* tmp_nu2 = it->mother(); 
          while (tmp_nu2->pdgId() == 14) tmp_nu2 = tmp_nu2->mother();
          if (tmp_nu2->pdgId() == 24)  nu2_W2_cand = tmp_nu2;
          while (tmp_nu2->pdgId() == 24) tmp_nu2 = tmp_nu2->mother();
          if (tmp_nu2->pdgId() == 6)
             {
              //   std::cout << "find antinuetrino candidate" << std::endl;
                 nu_positive = true;
                }
         
       }
      else if (it->pdgId() == 5 && it->mother()->pdgId() == 6 && !bquark)
      {
	  bquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
       //   std::cout << "find bquark candidate" << std::endl;
          b1_t_cand = it->mother();
      }
      else if (it->pdgId() == -5 && it->mother()->pdgId() == -6 && !bbarquark)
      {
	  bbarquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
         // std::cout << "find bbarquark candidate" << std::endl;
          b2_tbar_cand = it->mother();
      }

     // std::cout << "test" << std::endl;

   }// all Gen particles

    if (mu_negative and nu_negative && mu1_W1_cand == nu1_W1_cand)
       {
          Wtomu1nu1 = true;
         // std::cout << "find W1 from negative mu nu, mass " << mu1_W1_cand->mass() 
          //          << " energy " << mu1_W1_cand->energy() << std::endl;
    
          }
     else if (mu_negative && nu_negative)  std::cout << "find negative mu && nu but not find W1" << std::endl;
  
    if (mu_positive and nu_positive && mu2_W2_cand == nu2_W2_cand)
       {
          Wtomu2nu2 = true;
          //std::cout << "find W2 from positive mu nu, mass " << mu2_W2_cand->mass()
           //         << " energy " << mu2_W2_cand->energy() << std::endl;
    
          }
     else if (mu_positive && nu_positive)  std::cout << "find positive mu && nu but not find W2" << std::endl;
  
    if (Wtomu1nu1 and bbarquark and mu1_tbar_cand == b2_tbar_cand)
       {
         std::cout << "find muon and bbar, and they come from same tbar" << std::endl;
          tbarcand = mu1_tbar_cand; 
          tbartoWbbar = true;
          //if (abs(final_p4.M()-125)>5)    
         }
    else if(Wtomu1nu1 and bbarquark)   
       {
         std::cout << "find muon and bbar, and they do not come from same tbar" << std::endl;
         std::cout << "mu1_tbar energy " << mu1_tbar_cand->energy() << " px " << mu1_tbar_cand->px() << std::endl;
         std::cout << "b2_tbar energy " << b2_tbar_cand->energy() << " px " << b2_tbar_cand->px() << std::endl;

         }
    // pair b and bbar quark
    if (Wtomu2nu2 and bquark and mu2_t_cand == b1_t_cand)
       {
         std::cout << "find muon and b, and they come from same t" << std::endl;
	 tcand = mu2_t_cand;
         ttoWb = true;
          
         }
    else if(Wtomu2nu2 and bquark)   
       {
         std::cout << "find muon and b, and they do not come from same t" << std::endl;
         std::cout << "mu2_t energy " << mu2_t_cand->energy() << " px " << mu2_t_cand->px() << std::endl;
         std::cout << "b1_t energy " << b1_t_cand->energy() << " px " << b1_t_cand->px() << std::endl;

         }
//     else if()

   if (ttoWb and tbartoWbbar) {
        std::cout << "find t and tbar candidate " << std::endl;
        if (finalStates_){ 
    	mu1cand = stabledecendant(mu1_tbar_cand, 13);
        mu2cand = stabledecendant(mu2_t_cand, -13);
        nu1cand = stabledecendant(mu1_tbar_cand, -14);
        nu2cand = stabledecendant(mu2_t_cand, 14);
        b1cand = alldecendant(b1_t_cand, 5, false);
        b2cand = alldecendant(b2_tbar_cand, -5, false);
        w1cand = allancestor(mu1cand, -24, false);
	w2cand = allancestor(mu2cand, 24, false);   
        }else {
        mu1cand = findmudaughter(mu1_W1_cand);
        nu1cand = findnudaughter(mu1_W1_cand);
        mu2cand = findmudaughter(mu2_W2_cand);
        nu2cand = findnudaughter(mu2_W2_cand);

        b1cand = alldecendant(b1_t_cand, 5, true);
        b2cand = alldecendant(b2_tbar_cand, -5, true);
        w1cand = allancestor(mu1cand, -24, true);
	w2cand = allancestor(mu2cand, 24, true);   
        }
        //if ()
        std::cout <<"new mu1 "; printCandidate(mu1cand);	
        //std::cout <<"old mu1 "; printCandidate(findmudaughter(mu1_W1_cand));	
        std::cout <<"new nu1 "; printCandidate(nu1cand);
        //std::cout <<"old nu1 "; printCandidate(findnudaughter(mu1_W1_cand));	
        std::cout <<"new mu2 "; printCandidate(mu2cand);	
        //std::cout <<"old mu2 "; printCandidate(findmudaughter(mu2_W2_cand));	
        std::cout <<"new nu2 "; printCandidate(nu2cand);	
        //std::cout <<"old nu2 "; printCandidate(findnudaughter(mu2_W2_cand));	
        std::cout <<" w1 " ; printCandidate(w1cand);
        std::cout <<" w2 " ; printCandidate(w2cand);
        std::cout <<" b1(b) " ; printCandidate(b1cand);
        std::cout <<" b2(bbar) " ; printCandidate(b2cand);
        std::cout <<" t->Wb " ; printCandidate(tcand);
        std::cout <<" tbar->Wbbar " ; printCandidate(tbarcand);
        //std::cout <<"another b1 " ; printCandidate(alldecendant(tbarcand, 5, false));
        //std::cout <<"another b2 " ; printCandidate(alldecendant(tbarcand, -5, false));
       //std::cout <<" b jet"; bjetsLorentz(b1cand).Print(); 
        //std::cout <<" b jet"; bjetsLorentz(alldecendant(tbarcand, 5, false)).Print(); 
        std::cout <<" b jet"; stableDecendantsLorentz(b1cand).Print(); 
        std::cout <<" bbar jet"; stableDecendantsLorentz(b2cand).Print(); 
        //std::cout <<" t " ; printCandidate(tcand);
        //std::cout <<" tbar and its ancestors" <<std::endl; printallAncestors(tbarcand);
        //std::cout <<" t and its ancestors" <<std::endl; printallAncestors(tcand);
        //std::cout <<" tbar and its decendants" <<std::endl; printallDecendants(tbarcand);
        //std::cout <<" t and its decendants" <<std::endl; printallDecendants(tcand);
	const reco::Candidate* t_tmp = tcand->mother();
	while (t_tmp->pdgId() == 6) t_tmp = t_tmp->mother();
	const reco::Candidate* tbar_tmp = tbarcand->mother();
	while (tbar_tmp->pdgId() == -6) tbar_tmp = tbar_tmp->mother();
        std::cout <<" tbar_tmp and its decendants" <<std::endl; printallDecendants(tbar_tmp);
        std::cout <<" t_tmp and its decendants" <<std::endl; printallDecendants(t_tmp);
        
        //std::cout <<" tbar and its decendants" <<std::endl; printCandidate(tbarcand);
        jets_lorentz = stableDecendantsLorentz(b1cand) + stableDecendantsLorentz(b2cand);
        met_lorentz = calculateMET();
       
        fillbranches();
        evtree->Fill();
        //debug
     }
	//

   //if (h2tohh && runMMC_) runMMC();
   if (ttoWb and tbartoWbbar and  runMMC_){
   	mu1_lorentz.SetPtEtaPhiM(mu1cand->pt(), mu1cand->eta(), mu1cand->phi(), 0);
   	nu1_lorentz.SetPtEtaPhiM(nu1cand->pt(), nu1cand->eta(), nu1cand->phi(), 0);
        mu2_lorentz.SetPtEtaPhiM(mu2cand->pt(), mu2cand->eta(), mu2cand->phi(), 0); 
        nu2_lorentz.SetPtEtaPhiM(nu2cand->pt(), nu2cand->eta(), nu2cand->phi(), 0); 
        bbar_lorentz.SetXYZT(b1cand->px()+b2cand->px(), b1cand->py()+b2cand->py(), b1cand->pz()+b2cand->pz(), b1cand->energy()+b2cand->energy());
        int onshellMarker = -1;
        if (w1cand->mass() > w2cand->mass()) onshellMarker = 1;
        else onshellMarker = 2;
        // std::cout <<" mu1 lorenz "; mu1_lorentz.Print(); 
        //thismmc = new MMC();
        //std::cout << "onshellMarkder  " << onshellMarker << std::endl;
	thismmc = new MMC(&mu1_lorentz, &mu2_lorentz, &bbar_lorentz, &met_lorentz, &nu1_lorentz, &nu2_lorentz, onshellMarker, 
	true, ievent, mmcset_, fs, verbose_);
        //thismmc->printTrueLorentz();
        thismmc->runMMC();	
        delete thismmc;

    }    
}


// ------------ method called once each job just before starting event loop  ------------
void 
ttbarAnalyzer::beginJob()
{
   evtree = fs->make<TTree>("evtree", "evtree");
 //  output = new TFile("output.root","recreate");
  // output->cd();
   evtree->Branch("ievent",&ievent);
 //  evtree = new TTree("evtree","event tree");
   evtree->Branch("mu1_energy",&mu1_energy);
   evtree->Branch("mu1_px",&mu1_px);
   evtree->Branch("mu1_py",&mu1_py);
   evtree->Branch("mu1_pz",&mu1_pz);
   evtree->Branch("mu1_eta",&mu1_eta);
   evtree->Branch("mu1_phi",&mu1_phi);
   evtree->Branch("mu1_motherid",&mu1_motherid);
   evtree->Branch("mu1_mother_energy",&mu1_mother_energy);
   evtree->Branch("mu1_mother_px",&mu1_mother_px);
   evtree->Branch("mu1_mother_py",&mu1_mother_py);
   evtree->Branch("mu1_mother_pz",&mu1_mother_pz);
   evtree->Branch("mu1_mother_mass",&mu1_mother_mass);
   evtree->Branch("nu1_energy",&nu1_energy);
   evtree->Branch("nu1_px",&nu1_px);
   evtree->Branch("nu1_py",&nu1_py);
   evtree->Branch("nu1_pz",&nu1_pz);
   evtree->Branch("nu1_eta",&nu1_eta);
   evtree->Branch("nu1_phi",&nu1_phi);
   evtree->Branch("Wtomu1nu1",&Wtomu1nu1);

   evtree->Branch("mu2_energy",&mu2_energy);
   evtree->Branch("mu2_px",&mu2_px);
   evtree->Branch("mu2_py",&mu2_py);
   evtree->Branch("mu2_pz",&mu2_pz);
   evtree->Branch("mu2_eta",&mu2_eta);
   evtree->Branch("mu2_phi",&mu2_phi);
   evtree->Branch("mu2_mother_energy",&mu2_mother_energy);
   evtree->Branch("mu2_mother_px",&mu2_mother_px);
   evtree->Branch("mu2_mother_py",&mu2_mother_py);
   evtree->Branch("mu2_mother_pz",&mu2_mother_pz);
   evtree->Branch("mu2_mother_mass",&mu2_mother_mass);
   evtree->Branch("mu2_motherid",&mu2_motherid);
   evtree->Branch("nu2_energy",&nu2_energy);
   evtree->Branch("nu2_px",&nu2_px);
   evtree->Branch("nu2_py",&nu2_py);
   evtree->Branch("nu2_pz",&nu2_pz);
   evtree->Branch("nu2_eta",&nu2_eta);
   evtree->Branch("nu2_phi",&nu2_phi);
   evtree->Branch("Wtomu2nu2",&Wtomu2nu2);

   evtree->Branch("t_energy",&t_energy);
   evtree->Branch("t_px",&t_px);
   evtree->Branch("t_py",&t_py);
   evtree->Branch("t_pz",&t_pz);
   evtree->Branch("t_mass",&t_mass);
   
   evtree->Branch("b1_energy",&b1_energy);
   evtree->Branch("b1_px",&b1_px);
   evtree->Branch("b1_py",&b1_py);
   evtree->Branch("b1_pz",&b1_pz);
   evtree->Branch("b1_motherid",&b1_motherid);
   evtree->Branch("b2_energy",&b2_energy);
   evtree->Branch("b2_px",&b2_px);
   evtree->Branch("b2_py",&b2_py);
   evtree->Branch("b2_pz",&b2_pz);
   evtree->Branch("b2_motherid",&b2_motherid);
   evtree->Branch("jets_energy",&jets_energy);
   evtree->Branch("jets_px",&jets_px);
   evtree->Branch("jets_py",&jets_py);
   evtree->Branch("jets_pz",&jets_pz);
   evtree->Branch("jets_mass",&jets_mass);
   
   evtree->Branch("tbar_energy",&tbar_energy);
   evtree->Branch("tbar_px",&tbar_px);
   evtree->Branch("tbar_py",&tbar_py);
   evtree->Branch("tbar_pz",&tbar_pz);
   evtree->Branch("tbar_mass",&tbar_mass);
   

   evtree->Branch("met",&met);
   evtree->Branch("met_phi",&met_phi);
   evtree->Branch("met_px",&met_px);
   evtree->Branch("met_py",&met_py);
   
   evtree->Branch("tbartoWbbar",&tbartoWbbar);
   evtree->Branch("ttoWb",&ttoWb);
    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ttbarAnalyzer::endJob() 
{
//   std::cout << "endJob, ievent  " << ievent << std::endl;
    //output->Write();
   // output->Close();
/*
   // release space 
   delete mu_onshellW_lorentz;
   delete mu_offshellW_lorentz;
   delete nu_onshellW_lorentz;
   delete nu_offshellW_lorentz;
   delete onshellW_lorentz;
   delete offshellW_lorentz;
   delete t_lorentz;
   delete tbar_lorentz;
   delete h2tohh_lorentz;*/
   //delete jets_lorentz;
}


// ------------ method called when starting to processes a run  ------------
/*
void 
ttbarAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ttbarAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ttbarAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ttbarAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/


//---------- method called to fill branches --------------------------------------------
void 
ttbarAnalyzer::fillbranches(){
      mu1_energy = mu1cand->energy();
      mu1_px = mu1cand->px();
      mu1_py = mu1cand->py();
      mu1_pz = mu1cand->pz();
      mu1_mother_energy = w1cand->energy();
      mu1_mother_px = w1cand->px();
      mu1_mother_py = w1cand->py();
      mu1_mother_pz = w1cand->pz();
      mu1_mother_mass = w1cand->mass();
      nu1_energy = nu1cand->energy();
      nu1_px = nu1cand->px();
      nu1_py = nu1cand->py();
      nu1_pz = nu1cand->pz();

      mu1_eta = mu1cand->eta();
      mu1_phi = mu1cand->phi();
      nu1_eta = nu1cand->eta();
      nu1_phi = nu1cand->phi();

      mu2_energy = mu2cand->energy();
      mu2_px = mu2cand->px();
      mu2_py = mu2cand->py();
      mu2_pz = mu2cand->pz();
      mu2_mother_energy = w2cand->energy();
      mu2_mother_px = w2cand->px();
      mu2_mother_py = w2cand->py();
      mu2_mother_pz = w2cand->pz();
      mu2_mother_mass = w2cand->mass();
      nu2_energy = nu2cand->energy();
      nu2_px = nu2cand->px();
      nu2_py = nu2cand->py();
      nu2_pz = nu2cand->pz();

      mu2_eta = mu2cand->eta();
      mu2_phi = mu2cand->phi();
      nu2_eta = nu2cand->eta();
      nu2_phi = nu2cand->phi();

      t_energy = tcand->energy();
      t_px = tcand->px();
      t_py = tcand->py();
      t_pz = tcand->pz();
      t_mass = tcand->mass();

      b1_energy = b1cand->energy();
      b1_px = b1cand->px();
      b1_py = b1cand->py();
      b1_pz = b1cand->pz();
      b2_energy = b2cand->energy();
      b2_px = b2cand->px();
      b2_py = b2cand->py();
      b2_pz = b2cand->pz();

      jets_energy = jets_lorentz.Energy();
      jets_px = jets_lorentz.Px();
      jets_py = jets_lorentz.Py();
      jets_pz = jets_lorentz.Pz();
      jets_mass = jets_lorentz.M();
      
      tbar_energy = tbarcand->energy();
      tbar_px = tbarcand->px();
      tbar_py = tbarcand->py();
      tbar_pz = tbarcand->pz();
      tbar_mass = tbarcand->mass();
     
      
   
      met = met_lorentz.Energy();
      met_phi = met_lorentz.Phi();
      met_px = met_lorentz.Px();
      met_py = met_lorentz.Py();
   
}


//-------------- method called to find a muon daughter for given candidate -------------------------
const reco::Candidate* 
ttbarAnalyzer::findmudaughter(const reco::Candidate *Wcand){

    const reco::Candidate *tmpcand = NULL;
    int count = 0;
    for (unsigned int i = 0; i<Wcand->numberOfDaughters(); i++)
          if (abs(Wcand->daughter(i)->pdgId()) == 13) {
 		 count++;
                 tmpcand = Wcand->daughter(i);
    	   }
    if (count != 1) std::cout << " this candidate has more one mu daughter " << std::endl;
    return tmpcand;
}


//------------ method called to find muon descendants for given candidate -------------------------
// return one muon descendant and the number of muon descendants by "count"
const reco::Candidate*
ttbarAnalyzer::findmudescendants(const reco::Candidate *cand, int& count){

   const reco::Candidate * tmpcand = NULL;
   for (unsigned int i = 0; i<cand->numberOfDaughters(); i++){
         tmpcand = findmudescendants(cand->daughter(i), count);
         if (abs(cand->daughter(i)->pdgId()) == 13) count++;
         if (tmpcand == NULL && abs(cand->daughter(i)->pdgId()) == 13)  tmpcand = cand->daughter(i);                            
   }

   return tmpcand;
}


//-------------- method called to find a neutrino daughter for given candidate -------------------------
const reco::Candidate* 
ttbarAnalyzer::findnudaughter(const reco::Candidate *Wcand){

    const reco::Candidate *tmpcand = NULL;
    int count = 0;
    for (unsigned int i = 0; i<Wcand->numberOfDaughters(); i++)
          if (abs(Wcand->daughter(i)->pdgId()) == 14) {
 		 count++;
                 tmpcand = Wcand->daughter(i);
		
           }
    if (count != 1) std::cout << " this candidate has more one nu daughter " << std::endl;
    return tmpcand;
}

//------------ method called to find neutrino descendants for given candidate -------------------------
// return one neutrino descendant and the number of neutrno descendants by "count"
const reco::Candidate*
ttbarAnalyzer::findnudescendants(const reco::Candidate *cand, int& count){

   const reco::Candidate * tmpcand = NULL;
   for (unsigned int i = 0; i<cand->numberOfDaughters(); i++){
         tmpcand = findmudescendants(cand->daughter(i), count);
         if (abs(cand->daughter(i)->pdgId()) == 14) count++;
         if (tmpcand == NULL && abs(cand->daughter(i)->pdgId()) == 14)  tmpcand = cand->daughter(i);                            
   }

   return tmpcand;
}


//------------- method called to find stable decendant with pdgid = id
const reco::Candidate* 
ttbarAnalyzer::stabledecendant(const reco::Candidate* cand, int id){
   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++){
	if ((cand->daughter(i))->pdgId() == id && (cand->daughter(i))->status() == 1)
		return 	tmp=cand->daughter(i);
        else if (stabledecendant(cand->daughter(i),id)) 
		return tmp=stabledecendant(cand->daughter(i),id);
   }
    
    return tmp;

}


//------------ method called to find decendant with pdgid = id, 
//if first it true, then return the candidate closest to seed
//if first it false, then return the candidate farthest to seed
const reco::Candidate* 
ttbarAnalyzer::alldecendant(const reco::Candidate* cand, int id, bool first){
   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++){
        
	if ((cand->daughter(i))->pdgId() == id && first && cand->pdgId() != id)
		return 	tmp=cand->daughter(i);
	else if ((cand->daughter(i))->pdgId() == id && !first && !hasDaughter(cand->daughter(i), id)) 
		return  tmp=cand->daughter(i); // tmp does not has daughter with pdgid = id
	else if ((cand->daughter(i))->pdgId() == id && !first && (cand->daughter(i))->numberOfDaughters()>1) 
		return  tmp=cand->daughter(i);// tmp has more one daughters therefore it is final-states
        else if (alldecendant(cand->daughter(i),id, first)) 
		return tmp=alldecendant(cand->daughter(i),id);
   }
    
    return tmp;

}

//---------- method called to find a ancestor with pdgid = id, 
//if first is true, then return the candidate closest to seed
//if first is false, then return the candidate furthest to seed
const reco::Candidate*
ttbarAnalyzer::allancestor(const reco::Candidate* cand, int id, bool first){

   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfMothers(); i++){
        
	if ((cand->mother(i))->pdgId() == id && first && cand->pdgId() != id)
		return 	tmp=cand->mother(i);
	else if ((cand->mother(i))->pdgId() == id && !first && !hasMother(cand->mother(i), id)) 
		return  tmp=cand->mother(i);
        else if (allancestor(cand->mother(i),id, first)) 
		return tmp=allancestor(cand->mother(i),id, first);
   }
   return tmp;

}

//---------- method called to check whether cand has mother with pdgid = id -----------------------------------
bool
ttbarAnalyzer::hasMother(const reco::Candidate* cand, int id){

   for (unsigned int i=0; i < cand->numberOfMothers(); i++)
        if ((cand->mother(i))->pdgId() == id) return true;
   return false;

}

//-------- method called to check whether cand has daughter with pdgid = id ------------------------------------
bool
ttbarAnalyzer::hasDaughter(const reco::Candidate* cand, int id){
 
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        if ((cand->daughter(i))->pdgId() == id) return true;
   return false;

}


//----------- method called to calculate total lorentz vector from b bar jets ---------------
//
TLorentzVector
ttbarAnalyzer::stableDecendantsLorentz(const reco::Candidate *cand){
  
   TLorentzVector bjets = TLorentzVector();
   for (unsigned i = 0; i < cand->numberOfDaughters(); i++){
	if ((cand->daughter(i))->status() == 1){
		TLorentzVector tmp((cand->daughter(i))->px(), (cand->daughter(i))->py(), 
				(cand->daughter(i))->pz(),(cand->daughter(i))->energy());
                //printCandidate(cand->daughter(i));
        	bjets = bjets+tmp; 
	}else bjets = bjets+stableDecendantsLorentz(cand->daughter(i));
   }

   return bjets;
}

//-------------  method called to calculate MET in simuation ------------------ 
TLorentzVector 
ttbarAnalyzer::calculateMET(){

   TLorentzVector METlorentz = TLorentzVector();
   TVector2 met_pxpy(nu1cand->px()+nu2cand->px(), nu1cand->py()+nu2cand->py());
   METlorentz.SetPxPyPzE(nu1cand->px()+nu2cand->px(), nu1cand->py()+nu2cand->py(),0,met_pxpy.Mod());

   return METlorentz;
}




// ------------ method called for printing additional information, useful for debugging  ------------
void
ttbarAnalyzer::print() {
   //todo: combine all necessary print out for debug 
    std::cout << "print() " << std::endl;
    


}


/*
//------------- method called for printing h->WW->mumununu chain -------------------------
void 
ttbarAnalyzer::printHtoWWChain(){
          
    float h_energy = mu1_energy+mu2_energy+nu1_energy+nu2_energy;
    float h_px = mu1_px+mu2_px+nu1_px+nu2_px;
    float h_py = mu1_py+mu2_py+nu1_py+nu2_py;
    float h_pz = mu1_pz+mu2_pz+nu1_pz+nu2_pz;
          
    TLorentzVector final_p4(h_px, h_py, h_pz, h_energy);
    
    float h_energy_bb = b1_energy+b2_energy; 
    float h_px_bb = b1_px+b2_px; 
    float h_py_bb = b1_py+b2_py;
    float h_pz_bb = b1_pz+b2_pz;
    TLorentzVector h_bb_p4(h_px_bb,h_py_bb,h_pz_bb,h_energy_bb);

    TLorentzVector h2_final_tot(final_p4+h_bb_p4); 
    if (t){
             
               TLorentzVector h_p4(t_px, t_py, t_pz, t_energy);
               TLorentzVector WW_p4(mu1_W1_cand->px()+mu2_W2_cand->px(),mu1_W1_cand->py()+mu2_W2_cand->py(), 
                                    mu1_W1_cand->pz()+mu2_W2_cand->pz(),mu1_W1_cand->energy()+mu2_W2_cand->energy());
               std::cout << "invariant mass from h_p4: " << h_p4.M() 
                         << "	, from WW_p4 " << WW_p4.M() 
                         << "	, from final_p4 " << final_p4.M() << std::endl;
               if (abs(WW_p4.M()-h_p4.M())>1)  std::cout << "h->WW, invariant mass reconstruction discrepancy ? " << std::endl; 
               std::cout <<  " H -> WW " << std::endl;
               const reco::Candidate* tmp_cand1 = NULL;
               const reco::Candidate* tmp_cand2 = NULL;
               
               for (unsigned int n = 0; n<mu1_tbar_cand->numberOfDaughters(); n++){                                  
                      if ((mu1_tbar_cand->daughter(n))->pdgId()==-24) tmp_cand1 =  mu1_tbar_cand->daughter(n);
                      if ((mu1_tbar_cand->daughter(n))->pdgId()== 24) tmp_cand2 =  mu1_tbar_cand->daughter(n);
                      if (n >= 2) std::cout << "h has more 2 daughters, id " << (mu1_tbar_cand->daughter(n))->pdgId() << std::endl;
                  }
               while (tmp_cand1 != mu1_W1_cand || tmp_cand2 != mu2_W2_cand){
                TLorentzVector W1_p4(tmp_cand1->px(),tmp_cand1->py(),tmp_cand1->pz(), tmp_cand1->energy());
                TLorentzVector W2_p4(tmp_cand2->px(),tmp_cand2->py(),tmp_cand2->pz(), tmp_cand2->energy());
                TLorentzVector tmp_WW_p4(W1_p4+W2_p4);
                std::cout <<"W- num of daughters "<< tmp_cand1->numberOfDaughters() << " W-mass " << W1_p4.M() << " W1 four momentum "; W1_p4.Print();
                std::cout <<"W+ num of daughters "<< tmp_cand2->numberOfDaughters() << " W+mass " << W2_p4.M() << " W2 four momentum "; W2_p4.Print();
                std::cout << "Total invariant mass " << tmp_WW_p4.M() <<" tmp_WW four momentum "; tmp_WW_p4.Print(); 
                //if (tmp_cand1 != mu1_W1_cand) {
                     for (unsigned int i = 0; i<tmp_cand1->numberOfDaughters(); i++){
                          std::cout << " daughter of W- , id " << tmp_cand1->daughter(i)->pdgId() << "  status " << tmp_cand1->daughter(i)->status() <<std::endl; 
                          TLorentzVector dau_W1_p4(tmp_cand1->daughter(i)->px(),tmp_cand1->daughter(i)->py(),tmp_cand1->daughter(i)->pz(),tmp_cand1->daughter(i)->energy());
                          std::cout << " four momentum "; dau_W1_p4.Print();
                          if (tmp_cand1->daughter(i)->pdgId() == -24) tmp_cand1 = tmp_cand1->daughter(i);  
                          }
                // }
                //if (tmp_cand2 != mu2_W2_cand) {
                     for (unsigned int j = 0; j<tmp_cand2->numberOfDaughters(); j++){
                          std::cout << " daughter of W+ , id " << tmp_cand2->daughter(j)->pdgId() << "  status " << tmp_cand2->daughter(j)->status() <<std::endl; 
                          TLorentzVector dau_W2_p4(tmp_cand2->daughter(j)->px(),tmp_cand2->daughter(j)->py(),tmp_cand2->daughter(j)->pz(),tmp_cand2->daughter(j)->energy());
                          std::cout << " four momentum "; dau_W2_p4.Print();
                          if (tmp_cand2->daughter(j)->pdgId() == 24) tmp_cand2 = tmp_cand2->daughter(j);  
                          }
                // }
               }
               std::cout << "WW -> mumununu" << std::endl;
               while (tmp_cand1->status() != 1 || tmp_cand2->status() != 1){
                     
                     std::cout << "-------------  begin of this loop ----------------------" << std::endl;
     
                     int size1 = tmp_cand1->numberOfDaughters();
                     int size2 = tmp_cand2->numberOfDaughters();
                     int muon1 = -1;
                     int muon2 = -1;
                     float px=0;
                     float py=0;
                     float pz=0;
                     float energy=0;
                      if (tmp_cand1->pdgId() == 13)  {
                             //std::cout << "cand1 reaches final states, status of particle " << tmp_cand1->status() << std::endl;
                             px += nu1_px;
                             py += nu1_py;
                             pz += nu1_pz;
                             energy += nu1_energy;
                             }
                     std::cout << "cand1, id"<< tmp_cand1->pdgId() << " status " << tmp_cand1->status() <<" size of daughters " << size1 << std::endl; 
                     std::cout << " daughters of " << ((tmp_cand1->pdgId()==-24)?"W-  ":"muon- ") << std::endl;
                       for (int i = 0; i < size1; i++){
                             const Candidate * d1 = tmp_cand1->daughter(i); 
                             std::cout << "daughter id " << d1->pdgId() << "  status " << d1->status() << std::endl;
                             if (d1->pdgId() == 13 ) muon1 = i;
                             printCandidate(d1);
                             px += d1->px();
                             py += d1->py();
                             pz += d1->pz();
                             energy += d1->energy(); 
                   }
                      TLorentzVector cand1_lorentz(px, py, pz, energy); 
                      std::cout << " W- mass from W- Candidate " << mu1_mother_mass << " from mu-,nu " << cand1_lorentz.M() << std::endl;
                      if (muon1 != -1 && tmp_cand1->status() != 1) tmp_cand1 = tmp_cand1->daughter(muon1);
                      float px2 = 0.0;
                      float py2 = 0.0;
                      float pz2 = 0.0;
                      float energy2 = 0.0;
                      if (tmp_cand2->pdgId() == -13)  {
                           //  std::cout << "cand2 reaches final states, status of particle "<< tmp_cand2->status() << std::endl;
                             px2 += nu2_px;
                             py2 += nu2_py;
                             pz2 += nu2_pz;
                             energy2 += nu2_energy;
                             }
                     std::cout << "cand2, id" << tmp_cand2->pdgId() <<" status " << tmp_cand2->status() <<" size of daughters " << size2 << std::endl; 
                     std::cout << " daughters of " << ((tmp_cand2->pdgId()==24)?"W+  ":"muon+ ") << std::endl;
                       for (int j = 0; j < size2; j++){
                             const Candidate * d2 = tmp_cand2->daughter(j); 
                             std::cout << "daughter id " << d2->pdgId() << "  status " << d2->status() << std::endl;
                             if (d2->pdgId() == -13 ) muon2 = j;
                             printCandidate(d2);
                             px2 += d2->px();
                             py2 += d2->py();
                             pz2 += d2->pz();
                             energy2 += d2->energy(); 
                   }
                      TLorentzVector cand2_lorentz(px2, py2, pz2, energy2); 
                      std::cout << " W+ mass from W+ Candidate " << mu2_mother_mass << " from mu+,nu " << cand2_lorentz.M() << std::endl;
                      if (muon2 != -1 && tmp_cand2->status() != 1) tmp_cand2 = tmp_cand2->daughter(muon2);
                     TLorentzVector tmp = cand1_lorentz+cand2_lorentz;
                     std::cout <<"Total px " << tmp.Px() << " py " << tmp.Py() << " pz " << tmp.Pz() << " E " << tmp.Energy() << std::endl;
                     std::cout << " invariant mass from daughters of WW " << tmp.M() << std::endl;  
                     std::cout << "For Next loop status of cand1 " << tmp_cand1->status()  << "	cand2 " << tmp_cand2->status() << std::endl; 
                     std::cout << "-------------  end of this loop ----------------------" << std::endl;
              } 
            if (h2tohh){	
			std::cout <<"t invariant mass " << final_p4.M() <<" four momentum " << std::endl; 
                        final_p4.Print();
                        std::cout <<"tbar invariant mass " << h_bb_p4.M() << " four momentum " << std::endl;
                        h_bb_p4.Print();
                        std::cout <<"h2tohh invariant mass " << h2_final_tot.M() <<" total momentum " << std::endl;
                        h2_final_tot.Print();
			} 

        }//end if (t)

}

*/

//---------- method called to print candidates for debug ---------------------
void
ttbarAnalyzer::printCandidate(const reco::Candidate* cand){

   std::cout <<" Candidate id: "<< cand->pdgId() << " mass: " << cand->mass() <<" (P,E)= ("<< cand->px() <<", "<< cand->py()<<", "<< cand->pz()<<", "<< cand->energy() <<")" <<"(Pt,E) = ("<< cand->pt() <<", "<< cand->eta() <<", "<< cand->phi()<<", "<<cand->energy()<< ")" <<" status: " << cand->status() << std::endl;

}

//--------- method called to print all decendants for cand -------------------
void 
ttbarAnalyzer::printallDecendants(const reco::Candidate* cand){
   
   if (cand->status() != 1 && cand->numberOfDaughters() > 0){
        std::cout << "******************  children of id "<< cand->pdgId() <<"      *********************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        	printCandidate(cand->daughter(i));
        std::cout << "***********************************************************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
		printallDecendants(cand->daughter(i));

    }
}

//--------- method called to print all mothers for cand -------------------
void 
ttbarAnalyzer::printallAncestors(const reco::Candidate* cand){
  //status for beam particles?? 
   if (cand->status() != 0 && cand->numberOfMothers() > 0){
        std::cout << "******************  mothers of id "<< cand->pdgId() <<"      *********************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfMothers(); i++)
        	printCandidate(cand->mother(i));
        std::cout << "***********************************************************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfMothers(); i++)
		printallAncestors(cand->mother(i));

    }
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ttbarAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ttbarAnalyzer);


