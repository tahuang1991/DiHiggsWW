// -*- C++ -*-
//
// Package:    DiHiggsWWAnalyzer
// Class:      DiHiggsWWAnalyzer
// 
/**\class DiHiggsWWAnalyzer DiHiggsWWAnalyzer.cc DiHiggsWW/DiHiggsWWAnalyzer/plugins/DiHiggsWWAnalyzer.cc

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
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/GenMET.h"


//headers from root lib
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"


#include "DiHiggsWW/DiHiggsWWAnalyzer/src/MMC.h"
#include "DiHiggsWW/DiHiggsWWAnalyzer/src/deltaR.h"


class MMC;
//#define WMass 80.385   // W mass
//#define SMHMass 125.03 // SM module higgs mass


//
// class declaration
//

using namespace reco;

typedef std::pair<float, float> EtaPhi;

class DiHiggsWWAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DiHiggsWWAnalyzer(const edm::ParameterSet&);
      ~DiHiggsWWAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
    
   private:
      edm::ParameterSet cfg_;
      edm::ParameterSet mmcset_;
      // debuglevel constrol 
      int verbose_; 
      void print();
      void printHtoWWChain();
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
     std::string jetLabel_;
     std::string metLabel_;
     bool finalStates_;
     //cuts on GenJets
     double jetsPt_;
     double jetsEta_;
     double bjetsPt_;
     double bjetsEta_;
     double jetsDeltaR_;
     double jetleptonDeltaR_;
     double muonPt2_;
     double muonPt1_;
     double muonsEta_;
     double metPt_;
     double leptonIso_;
     
     bool rescalebjets_;
     bool metcorrection_;
     
    private:
     //decendants and ancestor
     const reco::Candidate* stabledecendant(const reco::Candidate* cand, int id);
     //const reco::Candidate* stablehtoWWdecendant(Particle p, PdgId id);
     const reco::Candidate* finddecendant(const reco::Candidate* cand, int id, bool first=false);
     const reco::Candidate* findancestor(const reco::Candidate* cand, int id, bool first=false);
     bool hasMother(const reco::Candidate* cand, int id);
     bool hasDaughter(const reco::Candidate* cand, int id);
      // ---------- Candidates in signal channel ---------------------------

      // ---------- Candidates in signal channel ---------------------------
      const reco::Candidate* mu1_W1_cand;
      const reco::Candidate* nu1_W1_cand;
      const reco::Candidate* mu2_W2_cand;
      const reco::Candidate* nu2_W2_cand;
      const reco::Candidate* mu1_htoWW_cand;
      const reco::Candidate* mu2_htoWW_cand;
      const reco::Candidate* b1_htobb_cand;
      const reco::Candidate* b2_htobb_cand;
      const reco::Candidate* h2tohh_cand;

      const reco::Candidate* mu1cand;
      const reco::Candidate* nu1cand;
      const reco::Candidate* mu2cand;
      const reco::Candidate* nu2cand;
      const reco::Candidate* b1cand;
      const reco::Candidate* b2cand;

      const reco::Candidate* w1cand;
      const reco::Candidate* w2cand;
      const reco::Candidate* htoWWcand;
      const reco::Candidate* htoBBcand;
      const reco::Candidate* h2tohhcand;

      const reco::GenJet* b1jet;
      const reco::GenJet* b2jet;

    private:
      
      TLorentzVector mu1_lorentz;
      TLorentzVector mu2_lorentz;
      TLorentzVector bjets_lorentz;
      TLorentzVector nu1_lorentz;
      TLorentzVector nu2_lorentz;
      TLorentzVector met_lorentz;
      TLorentzVector b1jet_lorentz;
      TLorentzVector b2jet_lorentz;
      TLorentzVector totjets_lorentz;
      TLorentzVector bjet_nu_lorentz;
      TLorentzVector bbarjet_nu_lorentz;
      TLorentzVector b_genp_lorentz;
      TLorentzVector bbar_genp_lorentz;
      TLorentzVector h2tohh_genp_lorentz;
     

    private:
      void stableDecendantsCandidates(const reco::Candidate* cand, std::vector<const reco::Candidate*> &candvec); 
      void neutrinoDecendantsCandidates(const reco::Candidate *cand, std::vector<const reco::Candidate*> &candvec);
   
    private:
      TLorentzVector calculateMET(); 

    private:
      TTree *evtree;
      TFile *output;
    private:
      void fillbranches(); 
      void init();
      edm::Service< TFileService > fs;
      
      //----------branches of tree ---------------------------
      int ievent;
      float mu1_energy;
      float mu1_px;
      float mu1_py;
      float mu1_pz;
      float mu1_eta;
      float mu1_phi;
      float w1_energy;
      float w1_px;
      float w1_py;
      float w1_pz;
      float w1_mass;
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
      float w2_energy;
      float w2_px;
      float w2_py;
      float w2_pz;
      float w2_mass;
      float nu2_energy;
      float nu2_px;
      float nu2_py;
      float nu2_pz;
      float nu2_eta;
      float nu2_phi;
      bool Wtomu2nu2;
     
      float htoWW_energy;
      float htoWW_px;
      float htoWW_py;
      float htoWW_pz; 
      float htoWW_mass;
    //  float w1_mass;
    //  float w2_mass;
      float b1_energy;
      float b1_px;
      float b1_py;
      float b1_pz;
      float b1_pt;
      float b1_eta;
      float b1_phi;
      int b1_motherid;
      float b2_energy;
      float b2_px;
      float b2_py;
      float b2_pz;
      float b2_pt;
      float b2_eta;
      float b2_phi;
      int b2_motherid;
      float bjet_energy;
      float bjet_px;
      float bjet_py;
      float bjet_pz;
      float bjet_pt;
      float bjet_eta;
      float bjet_phi;
      float bjet_mass;
      float bbarjet_energy;
      float bbarjet_px;
      float bbarjet_py;
      float bbarjet_pz;
      float bbarjet_pt;
      float bbarjet_eta;
      float bbarjet_phi;
      float bbarjet_mass;
      float dR_bjet;
      float dR_bbarjet;
      float rescalefactor;

      float dR_b1l1;
      float dR_b1l2;
      float dR_b2l1;
      float dR_b2l2;

      float bjet_energy_tot;
      float bjet_px_tot;
      float bjet_py_tot;
      float bjet_pz_tot;
      float bbarjet_energy_tot;
      float bbarjet_px_tot;
      float bbarjet_py_tot;
      float bbarjet_pz_tot;


      float bjet_decendant_energy;
      float bjet_decendant_px;
      float bjet_decendant_py;
      float bjet_decendant_pz;
      float bbarjet_decendant_energy;
      float bbarjet_decendant_px;
      float bbarjet_decendant_py;
      float bbarjet_decendant_pz;

      float bjet_nu_px;//bjet
      float bjet_nu_py;
      float bjet_nu_pz;
      float bjet_nu_energy;
      float bbarjet_nu_px;//bbarjet
      float bbarjet_nu_py;
      float bbarjet_nu_pz;
      float bbarjet_nu_energy;
      

      float totjets_energy;
      float totjets_px;
      float totjets_py;
      float totjets_pz;

      float met;
      float met_phi;
      float met_px;
      float met_py;

      float met_correction;
      float met_correction_phi;
      float met_correction_px;
      float met_correction_py;

      float htobb_energy;
      float htobb_px;
      float htobb_py;
      float htobb_pz;
      float htobb_mass;
      
      float h2tohh_energy;
      float h2tohh_px;
      float h2tohh_py;
      float h2tohh_pz;
      float h2tohh_mass;
      int h2tohh_pdgid;
      //cuts for higgstoWWbb
      bool mu_positive;
      bool mu_negative;
      bool nu_positive;
      bool nu_negative;
      bool bquark;
      bool bbarquark;
      bool htobb;
      bool htoWW;
      bool h2tohh;

      bool hasbjet;
      bool hasbbarjet;
      bool hasMET;
      bool hastwomuons;

     
    private:
      bool runMMC_;
      bool simulation_;
      MMC* thismmc;
      TTree* MMCtree;
      bool runmmc;
  
      float MMC_h2mass_prob;
      float MMC_h2massweight1_prob;
      float MMC_h2massweight4_prob;
      float MMC_h2mass_Entries;
      float MMC_h2mass_RMS;
      float MMC_h2mass_Mean;
      float MMC_h2mass_underflow;
      float MMC_h2mass_overflow;
     // MMC tree branches
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
DiHiggsWWAnalyzer::DiHiggsWWAnalyzer(const edm::ParameterSet& iConfig)
{
     verbose_ = iConfig.getUntrackedParameter<int>("verbose",0);
     mmcset_ = iConfig.getParameter<edm::ParameterSet>("mmcset"); 
     finalStates_ = iConfig.getParameter<bool>("finalStates");
     rescalebjets_ = iConfig.getParameter<bool>("rescalebjets");
     runMMC_ = iConfig.getParameter<bool>("runMMC");
     simulation_ = iConfig.getParameter<bool>("simulation");
     jetLabel_ = iConfig.getParameter<std::string>("jetLabel");
     metLabel_ = iConfig.getParameter<std::string>("metLabel");
     jetsPt_ = iConfig.getParameter<double>("jetsPt");
     jetsEta_ = iConfig.getParameter<double>("jetsEta");
     bjetsPt_ = iConfig.getParameter<double>("bjetsPt");
     bjetsEta_ = iConfig.getParameter<double>("bjetsEta");
     jetsDeltaR_ = iConfig.getParameter<double>("jetsDeltaR");
     jetleptonDeltaR_ = iConfig.getParameter<double>("jetleptonDeltaR");
     leptonIso_ = iConfig.getParameter<double>("leptonIso");
     muonPt1_ = iConfig.getParameter<double>("muonPt1");
     muonPt2_ = iConfig.getParameter<double>("muonPt2");
     muonsEta_ = iConfig.getParameter<double>("muonsEta");
     metPt_ = iConfig.getParameter<double>("metPt");
     metcorrection_ = iConfig.getParameter<bool>("metcorrection");
   // initilize candidates pointer

  

   //now do what ever initialization is needed
      ievent = 0;

}

//init tree
void
DiHiggsWWAnalyzer::init(){

      mu1_W1_cand = NULL;
      nu1_W1_cand = NULL;
      mu2_W2_cand = NULL;
      nu2_W2_cand = NULL;
      mu1_htoWW_cand = NULL;
      mu2_htoWW_cand = NULL;
      b1_htobb_cand = NULL;
      b2_htobb_cand = NULL;
      h2tohh_cand = NULL;
   //  
      mu1_energy = 0.0;
      mu1_px = 0.0;
      mu1_py = 0.0;
      mu1_pz = 0.0;
      w1_energy = 0.0;
      w1_px = 0.0;
      w1_py = 0.0;
      w1_pz = 0.0;
      w1_mass = 0.0;
      nu1_energy = 0.0;
      nu1_px = 0.0;
      nu1_py = 0.0;
      nu1_pz = 0.0;
      Wtomu1nu1 = false;

      mu2_energy = 0.0;
      mu2_px = 0.0;
      mu2_py = 0.0;
      mu2_pz = 0.0;
      w2_energy = 0.0;
      w2_px = 0.0;
      w2_py = 0.0;
      w2_pz = 0.0;
      w2_mass = 0.0;
      nu2_energy = 0.0;
      nu2_px = 0.0;
      nu2_py = 0.0;
      nu2_pz = 0.0;
      Wtomu2nu2 = false;

      htoWW_energy = 0.0;
      htoWW_px = 0.0;
      htoWW_py = 0.0;
      htoWW_pz = 0.0;
      htoWW_mass = 0.0;

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
      bjet_energy = 0.0;
      bjet_px = 0.0;
      bjet_py = 0.0;
      bjet_pz = 0.0;
      bjet_mass = 0.0;
      bbarjet_energy = 0.0;
      bbarjet_px = 0.0;
      bbarjet_py = 0.0;
      bbarjet_pz = 0.0;
      bbarjet_mass = 0.0;
      rescalefactor = 1.0;
      
      dR_b1l1 =0;
      dR_b1l2 = 0;
      dR_b2l1 =0;
      dR_b2l2 =0;

      htobb_energy = 0.0;
      htobb_px = 0.0;
      htobb_py = 0.0;
      htobb_pz = 0.0;
      htobb_mass = 0.0;
      
      h2tohh_energy = 0.0;
      h2tohh_px = 0.0;
      h2tohh_py = 0.0;
      h2tohh_pz = 0.0;
      h2tohh_mass = 0.0;

      met = 0.0;
      met_px = 0.0;
      met_py = 0.0;
      met_phi = 0.0;
     
      mu_positive = false;
      mu_negative = false;
      nu_positive = false;
      nu_negative = false;
      bquark = false;
      bbarquark = false;
      Wtomu1nu1 = false;
      Wtomu2nu2 = false;
      htobb = false;
      htoWW = false;
      h2tohh = false;
      runmmc = false;
      rescalefactor  =1.0;


   hasbjet = false;
   hasbbarjet =false;
   hasMET = false;
   hastwomuons =false;
   dR_bjet = jetsDeltaR_;
   dR_bbarjet = jetsDeltaR_;
   bjet_nu_energy = 0;
   bjet_nu_px = 0.0;
   bjet_nu_py = 0.0;
   bjet_nu_pz = 0.0;
   bbarjet_nu_energy = 0;
   bbarjet_nu_px = 0.0;
   bbarjet_nu_py = 0.0;
   bbarjet_nu_pz = 0.0;

   bjet_decendant_energy = 0;
   bjet_decendant_px = 0.0;
   bjet_decendant_py = 0.0;
   bjet_decendant_pz = 0.0;
   bbarjet_decendant_energy = 0;
   bbarjet_decendant_px = 0.0;
   bbarjet_decendant_py = 0.0;
   bbarjet_decendant_pz = 0.0;
  
   totjets_energy =0;
   totjets_px =0;
   totjets_py =0;
   totjets_pz =0;

}


DiHiggsWWAnalyzer::~DiHiggsWWAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DiHiggsWWAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
   edm::Handle<reco::GenParticleCollection> genParticleColl;
   iEvent.getByLabel("genParticles", genParticleColl);

   edm::Handle<reco::GenJetCollection> genjetColl;
   iEvent.getByLabel(jetLabel_, genjetColl);

   edm::Handle<edm::View<reco::GenMET> > genmetColl; 
   iEvent.getByLabel(metLabel_, genmetColl);

   std::vector<const reco::Candidate*> bjet_nucands;
   std::vector<const reco::Candidate*> bbarjet_nucands;
   std::vector<const reco::Candidate*> h2_nucands;

   init(); 
   //std::cout <<" genJetColl size " << genjetColl->size() <<" genmetColl size " << genmetColl->size() << std::endl;   

   for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {

//particle id, (muon13),(b5),(W+24),(SM higgs25)
   // particle id  it->pdgId()
   //
      if (it->pdgId() == 13 && it->status() == 1 && !mu_negative)
      {
	  //mu_negative = true;
          //std::cout << "find muon(-) with status 1" << std::endl;
          if (fabs(it->eta())>muonsEta_ ) continue;
          const reco::Candidate* tmp_mu1 = it->mother(); 
          while (tmp_mu1->pdgId() == 13 && tmp_mu1->numberOfMothers() == 1) tmp_mu1 = tmp_mu1->mother();
          if (tmp_mu1->numberOfMothers() != 1 ) std::cout << "muon has more than one mother particle" << std::endl;
          if (tmp_mu1->pdgId() == -24)  mu1_W1_cand = tmp_mu1;
           

          while (tmp_mu1->pdgId() == -24 && tmp_mu1->numberOfMothers() == 1) tmp_mu1 = tmp_mu1->mother();
          if (tmp_mu1->numberOfMothers() != 1 ) std::cout << "W- has more than one mother particle" << std::endl;
         
          if (tmp_mu1->pdgId() == 25)   
             {
                  //while (tmp_mu1->mother()->pdgId() == 25) tmp_mu1 = tmp_mu1->mother();
                  std::cout << "find muon(-) candidate" << std::endl;
	          mu_negative = true;
                  mu1_htoWW_cand = tmp_mu1;
		  mu1cand = it->clone();
                 // std::cout << "mother of this higgs, id " << tmp_mu1->mother()->pdgId() << " energy " << tmp_mu1->mother()->energy() << std::endl;
              }
      }
      else if (it->pdgId() == -13 && it->status() == 1 && !mu_positive)
      {
         // std::cout << "find muon(+) with status 1" << std::endl;
	//  mu_positive = true;
          if (fabs(it->eta())>muonsEta_ ) continue;
          const reco::Candidate* tmp_mu2 = it->mother(); 
          while (tmp_mu2->pdgId() == -13 && tmp_mu2->numberOfMothers() == 1) tmp_mu2 = tmp_mu2->mother();
          if (tmp_mu2->numberOfMothers() != 1)  std::cout << "muon has more than one mother particle" << std::endl;
          if (tmp_mu2->pdgId() == 24)  mu2_W2_cand = tmp_mu2;
          while (tmp_mu2->pdgId() == 24 && tmp_mu2->numberOfMothers() == 1) tmp_mu2 = tmp_mu2->mother();
          if (tmp_mu2->numberOfMothers() != 1 ) std::cout << "W+ has more than one mother particle" << std::endl;
    //      const reco::Candidate* tmphiggs_mu2 = tmp_mu2->mother();
          if (tmp_mu2->pdgId() == 25)   
             {
                // while (tmp_mu2->mother()->pdgId() == 25)   tmp_mu2 = tmp_mu2->mother();
                 std::cout << "find muon(+) candidate" << std::endl;
	         mu_positive = true;
                 mu2_htoWW_cand = tmp_mu2;
		 mu2cand = it->clone();
                // std::cout << "mother of this higgs, id " << tmp_mu2->mother()->pdgId() << " energy " << tmp_mu2->mother()->energy() << std::endl;
               }
        }
      else if (it->pdgId() == -14  && it->status() == 1 && !nu_negative )
      {
          const reco::Candidate* tmp_nu1 = it->mother();
          while (tmp_nu1->pdgId() == -14) tmp_nu1 = tmp_nu1->mother();
          if (tmp_nu1->pdgId() == -24)  nu1_W1_cand = tmp_nu1;
          while (tmp_nu1->pdgId() == -24) tmp_nu1 = tmp_nu1->mother();
          if (tmp_nu1->pdgId() == 25)
             {
            	 nu1cand = it->clone();
                 nu_negative = true;
                }
         
       }
      else if (it->pdgId() == 14 && it->status() == 1 && !nu_positive )
      {
          const reco::Candidate* tmp_nu2 = it->mother(); 
          while (tmp_nu2->pdgId() == 14) tmp_nu2 = tmp_nu2->mother();
          if (tmp_nu2->pdgId() == 24)  nu2_W2_cand = tmp_nu2;
          while (tmp_nu2->pdgId() == 24) tmp_nu2 = tmp_nu2->mother();
          if (tmp_nu2->pdgId() == 25)
             {
		 nu2cand = it->clone();
                 nu_positive = true;
                }
         
       }
      else if (it->pdgId() == 5 && it->mother()->pdgId() == 25 && !bquark)
      {
	  bquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
          b1_htobb_cand = it->mother();
      }
      else if (it->pdgId() == -5 && it->mother()->pdgId() == 25 && !bbarquark)
      {
	  bbarquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
          b2_htobb_cand = it->mother();
      }


   }// all Gen particles

    if (mu_negative and nu_negative && mu1_W1_cand == nu1_W1_cand)
       {
          Wtomu1nu1 = true;
	  w1cand = mu1_W1_cand;
         // std::cout << "find W1 from negative mu nu, mass " << mu1_W1_cand->mass() 
          //          << " energy " << mu1_W1_cand->energy() << std::endl;
    
          }
     else if (mu_negative && nu_negative)  std::cout << "find negative mu && nu but not find W1" << std::endl;
  
    if (mu_positive and nu_positive && mu2_W2_cand == nu2_W2_cand)
       {
          Wtomu2nu2 = true;
	  w2cand = mu2_W2_cand;
          //std::cout << "find W2 from positive mu nu, mass " << mu2_W2_cand->mass()
           //         << " energy " << mu2_W2_cand->energy() << std::endl;
    
          }
     else if (mu_positive && nu_positive)  std::cout << "find positive mu && nu but not find W2" << std::endl;
  
    if (Wtomu1nu1 and Wtomu2nu2 and mu1_htoWW_cand == mu2_htoWW_cand)
       {
         std::cout << "find 2 muons and they come from same higgs" << std::endl;
	 htoWWcand = mu1_htoWW_cand;
	 htoWW = true;

         }
    else if(Wtomu1nu1 and Wtomu2nu2)   
       {
         std::cout << "find 2 muons but they are not from same higgs" << std::endl;
         std::cout << "mu1_higgs energy " << mu1_htoWW_cand->energy() << " px " << mu1_htoWW_cand->px() << std::endl;
         std::cout << "mu2_higgs energy " << mu2_htoWW_cand->energy() << " px " << mu2_htoWW_cand->px() << std::endl;

         }
    // pair b and bbar quark
    if (bquark and bbarquark and b1_htobb_cand == b2_htobb_cand)
       {
        // std::cout << "find bbar and they come from same higgs" << std::endl;
         // const reco::Candidate* tmphiggs_bb = b1_htobb_cand->mother();
         // while (b1_htobb_cand->mother()->pdgId() == 25)  b1_htobb_cand = b1_htobb_cand->mother();
          
          htoBBcand = b1_htobb_cand;
          if(finalStates_){
        	b1cand = finddecendant(htoBBcand, 5, false);
        	b2cand = finddecendant(htoBBcand, -5, false);
          }else{
        	b1cand = finddecendant(htoBBcand, 5, true);
       	 	b2cand = finddecendant(htoBBcand, -5, true);
	  }
	
          htobb = true;
         }
//     else if()
          const reco::Candidate* tmp_htoWW = NULL;
          const reco::Candidate* tmp_htoBB = NULL;
       // pair htoWW and htoBB to find heavy higgs
    if (htoWW and htobb and mu1_htoWW_cand != b1_htobb_cand)
       {
          
          tmp_htoWW = mu1_htoWW_cand->mother();
          tmp_htoBB =  b1_htobb_cand->mother();
          while (tmp_htoWW->pdgId() == 25)  tmp_htoWW = tmp_htoWW->mother();
          while (tmp_htoBB->pdgId() == 25)  tmp_htoBB = tmp_htoBB->mother();
          //if (tmp_htoWW == tmp_htoBB && tmp_htoWW->pdgId() == 99927){
          if (tmp_htoWW == tmp_htoBB ){

              std::cout << "find 2 higgs and both of them come from same heavey higgs"  << std::endl;
              h2tohh = true;
              h2tohh_cand = tmp_htoWW;
              h2tohhcand = tmp_htoWW;
              if (verbose_>0) printHtoWWChain();
            }
          else {
              std::cout << "pdgId of mother htoWW " << tmp_htoWW->pdgId() << " htoBB " << tmp_htoBB->pdgId() << std::endl;  
          }
       }
    else if(htoWW and htobb and verbose_>0){
         std::cout << "find 2 higgs but they are not from same heavey higgs" << std::endl;
         std::cout << "mother of htoWW id " << tmp_htoWW->pdgId() <<" energy " << tmp_htoWW->energy() << " px " << tmp_htoWW->px() << std::endl;
         std::cout << "mother of htobb id " << tmp_htoBB->pdgId() << " energy " << tmp_htoBB->energy() << " px " << tmp_htoBB->px() << std::endl;

    }
   
   if (genmetColl->size() != 1){
	//std::cout <<" size of genmetColl " << genmetColl->size();
	} 
   reco::GenMET genMet(genmetColl->front());
   /*if (h2tohh){
   	met_lorentz = calculateMET();
   	std::cout <<"MET from two nuetrino "; met_lorentz.Print();
	}*/
   met_lorentz.SetXYZT(genMet.px(), genMet.py(), 0, genMet.pt()); 
      met = genMet.pt();
      met_phi = met_lorentz.Phi();
      met_px = genMet.px();
      met_py = genMet.py();

  // match bquark and bjet
	
   int nparts = 0;
   if ((!runMMC_ and htobb) || (runMMC_ and h2tohh)){
   for (reco::GenJetCollection::const_iterator jetit = genjetColl->begin(); jetit != genjetColl->end(); jetit++){
	//cuts on GenJets
	if (jetit->pt() >= jetsPt_ and std::fabs(jetit->eta()) <= jetsEta_){	
 		totjets_px += jetit->px();
		totjets_py += jetit->py();
		totjets_pz += jetit->pz();
		totjets_energy += jetit->energy();

	}
	if (jetit->pt() < bjetsPt_ or std::fabs(jetit->eta())> bjetsEta_) continue;
        std::vector <const reco::GenParticle*> mcparts = jetit->getGenConstituents();
  	for (unsigned i = 0; i < mcparts.size(); i++) {
    		const reco::GenParticle* mcpart = mcparts[i];
		const reco::Candidate* bcand;
		const reco::Candidate* hcand;
		if ((bcand=findancestor(mcpart, -5, false)) and bcand->mother()->pdgId() == 25){ 
			hcand = bcand->mother();
			if (hcand == htoBBcand and dR_bbarjet>deltaR(jetit->eta(), jetit->phi(), b2cand->eta(), b2cand->phi())){
				dR_bbarjet = deltaR(jetit->eta(), jetit->phi(), b2cand->eta(), b2cand->phi());
				//printCandidate(jetit->clone());
				//std::cout <<" bcand(-5) "; printCandidate(bcand);
				//std::cout <<" has h->bbar,h is the same from h->bb in genparticles flow,dR "<< dR_bbarjet <<std::endl;
				hasbbarjet = true;
				b2jet = jetit->clone();
				bbarjet_decendant_px = 0;
				bbarjet_decendant_py = 0;
				bbarjet_decendant_pz = 0;
				bbarjet_decendant_energy = 0;
				nparts = 0;
			}
			if (hasbbarjet and dR_bbarjet == deltaR(jetit->eta(), jetit->phi(), b2cand->eta(), b2cand->phi())){
				bbarjet_decendant_px += mcpart->px();
				bbarjet_decendant_py += mcpart->py();
				bbarjet_decendant_pz += mcpart->pz();
				bbarjet_decendant_energy += mcpart->energy();
				nparts++;
				continue;
			}
   			
		  }
		if ((bcand=findancestor(mcpart, 5, false)) and bcand->mother()->pdgId() == 25){
			hcand = bcand->mother();
			if (hcand == htoBBcand and dR_bjet>deltaR(jetit->eta(), jetit->phi(), b1cand->eta(), b1cand->phi())){
				dR_bjet = deltaR(jetit->eta(), jetit->phi(), b1cand->eta(), b1cand->phi());
				//printCandidate(jetit->clone());
				//std::cout <<" bcand(5) "; printCandidate(bcand);
				//std::cout <<" has h->b, h is the same from h->bb in genparticles flow,dR "<< dR_bjet <<std::endl;
				hasbjet = true;
				b1jet = jetit->clone();
				bjet_decendant_px = 0;
				bjet_decendant_py = 0;
				bjet_decendant_pz = 0;
				bjet_decendant_energy = 0;
				nparts =0;
			}
			if (hasbjet  and dR_bjet == deltaR(jetit->eta(), jetit->phi(), b1cand->eta(), b1cand->phi())){
				bjet_decendant_px += mcpart->px();
				bjet_decendant_py += mcpart->py();
				bjet_decendant_pz += mcpart->pz();
				bjet_decendant_energy += mcpart->energy();
				nparts++;
				continue;
			}
		 }	
	}
	//std::cout <<"mcparticle size " <<mcparts.size() <<"   nparts "<< nparts << std::endl;
     }
	if (!hasbbarjet or !hasbjet) std::cout <<" has h->bb, but failed to two b jets " << std::endl;
    }//htobb

   if (htobb and hasbjet and hasbbarjet and b1jet->px() == b2jet->px() and b1jet->py() == b2jet->py()) {
	std::cout <<" error: two bjets are in fact the same " << "b1jetpx "<< b1jet->px() <<" b2jetpx "<< b2jet->px()<<std::endl;
	hasbjet = false;
	hasbbarjet =false;
    }
   if (htobb and hasbjet) {
	 b1jet_lorentz.SetXYZT(b1jet->px(), b1jet->py(), b1jet->pz(), b1jet->energy());
   	 //std::cout << " neutrino from bjet " << std::endl;
	 neutrinoDecendantsCandidates(b1cand, bjet_nucands);
	 TLorentzVector tmp = TLorentzVector();
	 for (unsigned int i=0; i < bjet_nucands.size(); i++){
		printCandidate(bjet_nucands[i]);
		tmp += TLorentzVector(bjet_nucands[i]->px(), bjet_nucands[i]->py(), bjet_nucands[i]->pz(), bjet_nucands[i]->energy());
	}
	 bjet_nu_lorentz = tmp;
   }
   //else b2jet_lorentz.SetXYZT(1,0,0,1);
   else b1jet_lorentz = TLorentzVector();

   if (htobb and hasbbarjet) {
	b2jet_lorentz.SetXYZT(b2jet->px(), b2jet->py(), b2jet->pz(), b2jet->energy());
   	//std::cout << " neutrino from bbarjet " << std::endl;
	neutrinoDecendantsCandidates(b2cand, bbarjet_nucands);
	TLorentzVector tmp = TLorentzVector();
	 for (unsigned int i=0; i < bbarjet_nucands.size(); i++){
		printCandidate(bbarjet_nucands[i]);
		tmp += TLorentzVector(bbarjet_nucands[i]->px(), bbarjet_nucands[i]->py(), bbarjet_nucands[i]->pz(), bbarjet_nucands[i]->energy());
	}
	bbarjet_nu_lorentz = tmp;
   }
   //else b2jet_lorentz.SetXYZT(1,0,0,1);
   else b2jet_lorentz = TLorentzVector();

   if (htoWW){
         if (((mu1cand->pt() >= muonPt1_ and mu2cand->pt() >= muonPt2_) || (mu2cand->pt()>=muonPt1_ and mu1cand->pt()>=muonPt2_)) and
		fabs(mu1cand->eta()) <= muonsEta_ and fabs(mu2cand->eta()) <= muonsEta_) 
		hastwomuons = true;
    }
   
   if (genMet.pt() >= metPt_) hasMET =true;

   if (h2tohh) {
        std::cout << "find h2 candidate " << std::endl;
        std::cout << "h2 candidate id " << h2tohhcand->pdgId() << " mass " << h2tohhcand->mass() << std::endl;
        if (finalStates_){ 
    	//mu1cand = stabledecendant(mu1_W1_cand, 13);
        //mu2cand = stabledecendant(mu2_W2_cand, -13);
        //nu1cand = stabledecendant(mu1_W1_cand, -14);
        //nu2cand = stabledecendant(mu2_W2_cand, 14);
        w1cand = findancestor(mu1cand, -24, false);
	w2cand = findancestor(mu2cand, 24, false);   
        }else {
        mu1cand = findmudaughter(mu1_W1_cand);
        nu1cand = findnudaughter(mu1_W1_cand);
        mu2cand = findmudaughter(mu2_W2_cand);
        nu2cand = findnudaughter(mu2_W2_cand);

        w1cand = findancestor(mu1cand, -24, true);
	w2cand = findancestor(mu2cand, 24, true);   
        }
	
        
        std::cout <<" h2tohh "; printCandidate(h2tohhcand);
        std::cout <<" mu1 "; printCandidate(mu1cand);	
        std::cout <<" nu1 "; printCandidate(nu1cand);
        //std::cout <<"old mu1 "; printCandidate(findmudaughter(mu1_W1_cand));	
        //std::cout <<"old nu1 "; printCandidate(findnudaughter(mu1_W1_cand));	
        std::cout <<" mu2 "; printCandidate(mu2cand);	
        std::cout <<" nu2 "; printCandidate(nu2cand);	
	std::cout <<" nu12_px" << nu1cand->px()+nu2cand->px() <<" nu12_py " << nu1cand->py()+nu2cand->py() << std::endl;
        //std::cout <<"old mu2 "; printCandidate(findmudaughter(mu2_W2_cand));	
        //std::cout <<"old nu2 "; printCandidate(findnudaughter(mu2_W2_cand));	
        std::cout <<" w1 " ; printCandidate(w1cand);
        std::cout <<" w2 " ; printCandidate(w2cand);
        std::cout <<" b1 " ; printCandidate(b1cand);
        std::cout <<" b2 " ; printCandidate(b2cand);
 	
        //std::cout <<" htoWW " ; printCandidate(htoWWcand);
        //std::cout <<" htoBB and its decendants" <<std::endl; printCandidate(htoBBcand);
	//printallDecendants(h2tohhcand);
        //debug
         /* float h_energy = mu1_energy+mu2_energy+nu1_energy+nu2_energy;
          float h_px = mu1_px+mu2_px+nu1_px+nu2_px;
          float h_py = mu1_py+mu2_py+nu1_py+nu2_py;
          float h_pz = mu1_pz+mu2_pz+nu1_pz+nu2_pz;
          TLorentzVector final_p4(h_px, h_py, h_pz, h_energy);
        if (fabs(final_p4.M()-htoWWcand->mass())>0.15) {
		printHtoWWChain(); 
		printallDecendants(htoWWcand);
        }*/
   }

   if (hasbjet and hasbbarjet and h2tohh){
        /*
        std::cout <<" b1 " ; printCandidate(b1cand);
	std::cout<<"bjet (P,E)=( " << bjet_px << ", "<< bjet_py << ", "<<bjet_pz << ", " << bjet_energy <<")" << std::endl; 
	std::cout<<"bjet_tot (P,E)=( " << bjet_px_tot << ", "<< bjet_py_tot << ", "<<bjet_pz_tot << ", " << bjet_energy_tot <<")" << std::endl; 
	std::cout<<"bjet_bdecendant (P,E)=( " << bjet_decendant_px << ", "<< bjet_decendant_py << ", "<<bjet_decendant_pz << ", " << bjet_decendant_energy <<")" << std::endl; 
        std::cout <<" b2 " ; printCandidate(b2cand);
	std::cout<<"bbarjet (P,E)=( " << bbarjet_px << ", "<< bbarjet_py << ", "<<bbarjet_pz << ", " << bbarjet_energy <<")" << std::endl; 
	std::cout<<"bbarjet_tot (P,E)=( " << bbarjet_px_tot << ", "<< bbarjet_py_tot << ", "<<bbarjet_pz_tot << ", " << bbarjet_energy_tot <<")" << std::endl; 
	std::cout<<"bbarjet_bdecendant (P,E)=( " << bbarjet_decendant_px << ", "<< bbarjet_decendant_py << ", "<<bbarjet_decendant_pz << ", " << bbarjet_decendant_energy <<")" << std::endl; 
	*/
    	//std::cout <<" num of neutrino from bjet is " << bjet_nucands.size() <<" p4: "; bjet_nu_lorentz.Print();
    	//std::cout <<" num of neutrino from bbarjet is " << bbarjet_nucands.size() <<" p4: "; bbarjet_nu_lorentz.Print();
    	/*
    	neutrinoDecendantsCandidates(h2tohhcand, h2_nucands);
	TLorentzVector met_diff = calculateMET();
        met_diff = met_diff+bjet_nu_lorentz+bbarjet_nu_lorentz-met_lorentz;
	if (h2_nucands.size() != (bjet_nucands.size()+bbarjet_nucands.size()+2))
		std::cout <<" h2_nucands: "<< h2_nucands.size() <<"  bjets nucands+2 "<< 
			(bjet_nucands.size()+bbarjet_nucands.size()+2)<<std::endl;
	std::cout <<"nu12+bjet_nu and met diff  Px " << met_diff.Px() << " Py "<< met_diff.Py() << std::endl; 
	
        if (met_diff.Pt()>10){
	
   	for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {
		if (it->status()==1 &&(abs(it->pdgId()) == 12 || abs(it->pdgId()) == 14 || abs(it->pdgId())==16 ))
		     {unsigned int i=0; 
		      for (; i<h2_nucands.size(); i++){
			if (it->pdgId() == h2_nucands[i]->pdgId() && it->px()==h2_nucands[i]->px())
				break;
			}
		       if (i==h2_nucands.size()) {std::cout <<"extra neutrino "<< std::endl; printCandidate(it->clone());
		      			printallAncestors(it->clone());}
		      }
	} 	
	std::cout <<" neutrinos from h2 " << std::endl;
	for (unsigned int i =0; i<h2_nucands.size(); i++)
		printCandidate(h2_nucands[i]);
	}*/
	dR_b1l1 = deltaR(b1jet->eta(),b1jet->phi(),mu1cand->eta(),mu1cand->phi());
	dR_b1l2 = deltaR(b1jet->eta(),b1jet->phi(),mu2cand->eta(),mu2cand->phi());
	dR_b2l1 = deltaR(b2jet->eta(),b2jet->phi(),mu1cand->eta(),mu1cand->phi());
	dR_b2l2 = deltaR(b2jet->eta(),b2jet->phi(),mu2cand->eta(),mu2cand->phi());
        bjets_lorentz.SetXYZT(b1jet->px()+b2jet->px(), b1jet->py()+b2jet->py(), b1jet->pz()+b2jet->pz(), b1jet->energy()+b2jet->energy());
       //need to correction later
	std::cout <<"m_(jets) before rescaling " << bjets_lorentz.M(); bjets_lorentz.Print(); 
	if (rescalebjets_) rescalefactor = 125.0/bjets_lorentz.M();
        totjets_lorentz.SetXYZT(totjets_px+(rescalefactor-1)*bjets_lorentz.Px(), totjets_py+(rescalefactor-1)*bjets_lorentz.Py(), totjets_pz+(rescalefactor-1)*bjets_lorentz.Pz(), totjets_energy+(rescalefactor-1)*bjets_lorentz.Energy());
    }
    
    //met correction? remove? 
    TVector2 met_correction_vec2(met_lorentz.Px()-(rescalefactor-1)*bjets_lorentz.Px(), met_lorentz.Py()-(rescalefactor-1)*bjets_lorentz.Py());
    met_correction = met_correction_vec2.Mod();
    met_correction_px = met_correction_vec2.Px();
    met_correction_py = met_correction_vec2.Py();
    met_correction_phi = met_correction_vec2.Phi();
    if (metcorrection_){
	//met_lorentz.SetXYZT(met_correction_px,met_correction_py,0,met_correction);
	met_lorentz = met_lorentz-bjet_nu_lorentz-bbarjet_nu_lorentz;
    	std::cout <<"met due to  correction(nubjets) "; met_lorentz.Print();
	}

   //if (h2tohh && runMMC_) runMMC();
   if (h2tohh and runMMC_ and hasbjet and hasbbarjet){
   	mu1_lorentz.SetPtEtaPhiM(mu1cand->pt(), mu1cand->eta(), mu1cand->phi(), 0);
   	nu1_lorentz.SetPtEtaPhiM(nu1cand->pt(), nu1cand->eta(), nu1cand->phi(), 0);
        mu2_lorentz.SetPtEtaPhiM(mu2cand->pt(), mu2cand->eta(), mu2cand->phi(), 0); 
        nu2_lorentz.SetPtEtaPhiM(nu2cand->pt(), nu2cand->eta(), nu2cand->phi(), 0); 
	//rescale bjets
	//std::cout <<" b1jet "; b1jet_lorentz.Print();
	//std::cout <<" b2jet "; b2jet_lorentz.Print();
		
        int onshellMarker = -1;

        if (w1cand->mass() > w2cand->mass()) onshellMarker = 1;
        else onshellMarker = 2;
        b_genp_lorentz.SetXYZT(b1cand->px(), b1cand->py(), b1cand->pz(), b1cand->energy()); 
        bbar_genp_lorentz.SetXYZT(b2cand->px(), b2cand->py(), b2cand->pz(), b2cand->energy()); 
	h2tohh_genp_lorentz.SetXYZT(h2tohhcand->px(), h2tohhcand->py(), h2tohhcand->pz(), h2tohhcand->energy());
	TLorentzVector bgenp_pt1_lorentz,bgenp_pt2_lorentz,bjet_pt1_lorentz,bjet_pt2_lorentz;
	if (b1jet->pt()>b2jet->pt()) {
		bgenp_pt1_lorentz = b_genp_lorentz;//pair with large Pt
		bgenp_pt2_lorentz = bbar_genp_lorentz;
		bjet_pt1_lorentz = b1jet_lorentz;
		bjet_pt2_lorentz = b2jet_lorentz;
	}
	else{
		bgenp_pt1_lorentz = bbar_genp_lorentz;
		bgenp_pt2_lorentz = b_genp_lorentz;
		bjet_pt1_lorentz = b2jet_lorentz;
		bjet_pt2_lorentz = b1jet_lorentz;
	}
        //met_lorentz = calculateMET();
    std::cout <<"met before correction px "<< met_px  <<"  py "<< met_py << "  pt "<< met << std::endl;
    std::cout <<"met after correction px "<< met_correction_px  <<"  py "<< met_correction_py << "  pt "<< met_correction << std::endl;
    std::cout <<"NoNufromb met used in MMC  px "<< met_lorentz.Px()  <<"  py "<< met_lorentz.Py() << "  pt "<< met_lorentz.Pt() << std::endl;
        //std::cout <<" totjets lorenz "; totjets_lorentz.Print(); 
        TLorentzVector sum_lorentz = TLorentzVector();
        //std::cout <<" sum of (mu1 nu1 mu2 nu2 (corrected)bjets ) "; 
	//sum_lorentz = mu1_lorentz+mu2_lorentz+nu1_lorentz+nu2_lorentz+bjets_lorentz;
	//sum_lorentz.Print();
        //std::cout <<" sum of (mu1 nu1 mu2 nu2 (corrected)totjets ) "; 
	//sum_lorentz = mu1_lorentz+mu2_lorentz+nu1_lorentz+nu2_lorentz+totjets_lorentz;
	//sum_lorentz.Print();
        std::cout <<" sum of (nu1 nu2 (correct)MET ) "; 
 	sum_lorentz = nu1_lorentz+nu2_lorentz-met_lorentz;
	std::cout <<" px " << sum_lorentz.Px() <<" py " << sum_lorentz.Py() << std::endl;
        //thismmc = new MMC();
        //std::cout << "onshellMarkder  " << onshellMarker << std::endl;
	thismmc = new MMC(&mu1_lorentz, &mu2_lorentz, &bjet_pt1_lorentz, &bjet_pt2_lorentz, &totjets_lorentz,&met_lorentz, 
	&nu1_lorentz, &nu2_lorentz, &bgenp_pt1_lorentz, &bgenp_pt2_lorentz, &h2tohh_genp_lorentz, 
	onshellMarker, simulation_, ievent, mmcset_, fs, verbose_);
        //thismmc->printTrueLorentz();
        runmmc=thismmc->runMMC();	
        if (runmmc) {
		MMCtree =  (thismmc->getMMCTree())->CloneTree();
		std::cout <<" MMCtree entries " << MMCtree->GetEntries() << std::endl;
		TH1F* MMC_h2mass =(TH1F*)(thismmc->getMMCh2()).Clone("MMC_h2mass");
		TH1F* MMC_h2mass_weight1 =(TH1F*)(thismmc->getMMCh2weight1()).Clone("MMC_h2massweight1");
		TH1F* MMC_h2mass_weight4 =(TH1F*)(thismmc->getMMCh2weight4()).Clone("MMC_h2massweight4");
		std::cout <<" Mass_h2mass in Analyzer " << std::endl;
		MMC_h2mass_prob = (MMC_h2mass->GetXaxis())->GetBinCenter(MMC_h2mass->GetMaximumBin());
		MMC_h2massweight1_prob = (MMC_h2mass_weight1->GetXaxis())->GetBinCenter(MMC_h2mass_weight1->GetMaximumBin());
		MMC_h2massweight4_prob = (MMC_h2mass_weight4->GetXaxis())->GetBinCenter(MMC_h2mass_weight4->GetMaximumBin());
		MMC_h2mass_RMS = MMC_h2mass->GetRMS();
		MMC_h2mass_Entries = MMC_h2mass->GetEntries();
		MMC_h2mass_Mean = MMC_h2mass->GetMean();
		int nbin=(MMC_h2mass->GetXaxis())->GetNbins();
		MMC_h2mass_overflow = MMC_h2mass->GetBinContent(nbin+1);
		MMC_h2mass_underflow = MMC_h2mass->GetBinContent(-1);
                std::cout <<" most prob " << MMC_h2mass_prob <<" RMS "<< MMC_h2mass_RMS << " entries " << MMC_h2mass_Entries 
		<< " most prob weight1 "<< MMC_h2massweight1_prob <<" weight4 "<< MMC_h2massweight4_prob <<std::endl;
	
	}
        delete thismmc;

    }    
    if (htoWW or htobb){
        fillbranches();
        evtree->Fill();
    }

   /*if (h2tohh and runmmc){
   	std::cout <<" MMCtree entries " << MMCtree->GetEntries() << std::endl;
	TH1F* th1 = new TH1F("th1","th1",900,100,1000);
	//th1->SetDirectory(0);
   	//TH1::AddDirectory(kFALSE);
   	TCanvas *c1 = new TCanvas("c1","c1",600,800);
	MMCtree->Draw("h2tohh_Mass>>th1","weight4");
        c1->Print("c1test.png");
	std::cout <<"th1 entries "<< th1->GetEntries() << std::endl;
   }*/
}


// ------------ method called once each job just before starting event loop  ------------
void 
DiHiggsWWAnalyzer::beginJob()
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
   evtree->Branch("w1_energy",&w1_energy);
   evtree->Branch("w1_px",&w1_px);
   evtree->Branch("w1_py",&w1_py);
   evtree->Branch("w1_pz",&w1_pz);
   evtree->Branch("w1_mass",&w1_mass);
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
   evtree->Branch("w2_energy",&w2_energy);
   evtree->Branch("w2_px",&w2_px);
   evtree->Branch("w2_py",&w2_py);
   evtree->Branch("w2_pz",&w2_pz);
   evtree->Branch("w2_mass",&w2_mass);
   evtree->Branch("nu2_energy",&nu2_energy);
   evtree->Branch("nu2_px",&nu2_px);
   evtree->Branch("nu2_py",&nu2_py);
   evtree->Branch("nu2_pz",&nu2_pz);
   evtree->Branch("nu2_eta",&nu2_eta);
   evtree->Branch("nu2_phi",&nu2_phi);
   evtree->Branch("Wtomu2nu2",&Wtomu2nu2);

   evtree->Branch("htoWW_energy",&htoWW_energy);
   evtree->Branch("htoWW_px",&htoWW_px);
   evtree->Branch("htoWW_py",&htoWW_py);
   evtree->Branch("htoWW_pz",&htoWW_pz);
   evtree->Branch("htoWW_mass",&htoWW_mass);
   
   evtree->Branch("b1_energy",&b1_energy);
   evtree->Branch("b1_px",&b1_px);
   evtree->Branch("b1_py",&b1_py);
   evtree->Branch("b1_pz",&b1_pz);
   evtree->Branch("b1_pt",&b1_pt);
   evtree->Branch("b1_eta",&b1_eta);
   evtree->Branch("b1_phi",&b1_phi);
   evtree->Branch("b1_motherid",&b1_motherid);
   evtree->Branch("b2_energy",&b2_energy);
   evtree->Branch("b2_px",&b2_px);
   evtree->Branch("b2_py",&b2_py);
   evtree->Branch("b2_pz",&b2_pz);
   evtree->Branch("b2_pt",&b2_pt);
   evtree->Branch("b2_eta",&b2_eta);
   evtree->Branch("b2_phi",&b2_phi);
   evtree->Branch("b2_motherid",&b2_motherid);
   evtree->Branch("hasbjet",&hasbjet);
   evtree->Branch("bjet_energy",&bjet_energy);
   evtree->Branch("bjet_px",&bjet_px);
   evtree->Branch("bjet_py",&bjet_py);
   evtree->Branch("bjet_pz",&bjet_pz);
   evtree->Branch("bjet_pt",&bjet_pt);
   evtree->Branch("bjet_eta",&bjet_eta);
   evtree->Branch("bjet_phi",&bjet_phi);
   evtree->Branch("bjet_mass",&bjet_mass);
   evtree->Branch("hasbbarjet",&hasbbarjet);
   evtree->Branch("bbarjet_energy",&bbarjet_energy);
   evtree->Branch("bbarjet_px",&bbarjet_px);
   evtree->Branch("bbarjet_py",&bbarjet_py);
   evtree->Branch("bbarjet_pz",&bbarjet_pz);
   evtree->Branch("bbarjet_pt",&bbarjet_pt);
   evtree->Branch("bbarjet_eta",&bbarjet_eta);
   evtree->Branch("bbarjet_phi",&bbarjet_phi);
   evtree->Branch("bbarjet_mass",&bbarjet_mass);
   evtree->Branch("dR_bjet",&dR_bjet);
   evtree->Branch("dR_bbarjet",&dR_bbarjet);
   evtree->Branch("dR_b1l1",&dR_b1l1);
   evtree->Branch("dR_b1l2",&dR_b1l2);
   evtree->Branch("dR_b2l1",&dR_b2l1);
   evtree->Branch("dR_b2l2",&dR_b2l2);
   evtree->Branch("rescalefactor",&rescalefactor);

   evtree->Branch("bjet_nu_energy",&bjet_nu_energy);
   evtree->Branch("bjet_nu_px",&bjet_nu_px);
   evtree->Branch("bjet_nu_py",&bjet_nu_py);
   evtree->Branch("bjet_nu_pz",&bjet_nu_pz);
   evtree->Branch("bbarjet_nu_energy",&bbarjet_nu_energy);
   evtree->Branch("bbarjet_nu_px",&bbarjet_nu_px);
   evtree->Branch("bbarjet_nu_py",&bbarjet_nu_py);
   evtree->Branch("bbarjet_nu_pz",&bbarjet_nu_pz);

   evtree->Branch("totjets_energy",&totjets_energy);
   evtree->Branch("totjets_px",&totjets_px);
   evtree->Branch("totjets_py",&totjets_py);
   evtree->Branch("totjets_pz",&totjets_pz);
   
   evtree->Branch("bjet_decendant_energy",&bjet_decendant_energy);
   evtree->Branch("bjet_decendant_px",&bjet_decendant_px);
   evtree->Branch("bjet_decendant_py",&bjet_decendant_py);
   evtree->Branch("bjet_decendant_pz",&bjet_decendant_pz);
   evtree->Branch("bbarjet_decendant_energy",&bbarjet_decendant_energy);
   evtree->Branch("bbarjet_decendant_px",&bbarjet_decendant_px);
   evtree->Branch("bbarjet_decendant_py",&bbarjet_decendant_py);
   evtree->Branch("bbarjet_decendant_pz",&bbarjet_decendant_pz);
   
   evtree->Branch("htobb_energy",&htobb_energy);
   evtree->Branch("htobb_px",&htobb_px);
   evtree->Branch("htobb_py",&htobb_py);
   evtree->Branch("htobb_pz",&htobb_pz);
   evtree->Branch("htobb_mass",&htobb_mass);
   
   evtree->Branch("h2tohh_energy",&h2tohh_energy);
   evtree->Branch("h2tohh_px",&h2tohh_px);
   evtree->Branch("h2tohh_py",&h2tohh_py);
   evtree->Branch("h2tohh_pz",&h2tohh_pz);
   evtree->Branch("h2tohh_mass",&h2tohh_mass);
   evtree->Branch("h2tohh_pdgid",&h2tohh_pdgid);

   evtree->Branch("met",&met);
   evtree->Branch("met_phi",&met_phi);
   evtree->Branch("met_px",&met_px);
   evtree->Branch("met_py",&met_py);
   evtree->Branch("met_correction",&met_correction);
   evtree->Branch("met_correction_phi",&met_correction_phi);
   evtree->Branch("met_correction_px",&met_correction_px);
   evtree->Branch("met_correction_py",&met_correction_py);
   
   evtree->Branch("hasMET",&hasMET);
   evtree->Branch("hastwomuons",&hastwomuons);
   evtree->Branch("htobb",&htobb);
   evtree->Branch("htoWW",&htoWW);
   evtree->Branch("h2tohh",&h2tohh);
  
   evtree->Branch("runmmc",&runmmc); 
   evtree->Branch("MMC_h2mass_prob",&MMC_h2mass_prob);
   evtree->Branch("MMC_h2massweight1_prob",&MMC_h2massweight1_prob);
   evtree->Branch("MMC_h2massweight4_prob",&MMC_h2massweight4_prob);
   evtree->Branch("MMC_h2mass_RMS",&MMC_h2mass_RMS);
   evtree->Branch("MMC_h2mass_Mean",&MMC_h2mass_Mean);
   evtree->Branch("MMC_h2mass_Entries",&MMC_h2mass_Entries);
   evtree->Branch("MMC_h2mass_overflow",&MMC_h2mass_overflow);
   evtree->Branch("MMC_h2mass_underflow",&MMC_h2mass_underflow);
    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiHiggsWWAnalyzer::endJob() 
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
   delete htoWW_lorentz;
   delete htoBB_lorentz;
   delete h2tohh_lorentz;*/
   //delete jets_lorentz;
}


// ------------ method called when starting to processes a run  ------------
/*
void 
DiHiggsWWAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
DiHiggsWWAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
DiHiggsWWAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
DiHiggsWWAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/


//---------- method called to fill branches --------------------------------------------
void 
DiHiggsWWAnalyzer::fillbranches(){
   if(htoWW){
      mu1_energy = mu1cand->energy();
      mu1_px = mu1cand->px();
      mu1_py = mu1cand->py();
      mu1_pz = mu1cand->pz();
      w1_energy = w1cand->energy();
      w1_px = w1cand->px();
      w1_py = w1cand->py();
      w1_pz = w1cand->pz();
      w1_mass = w1cand->mass();
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
      w2_energy = w2cand->energy();
      w2_px = w2cand->px();
      w2_py = w2cand->py();
      w2_pz = w2cand->pz();
      w2_mass = w2cand->mass();
      nu2_energy = nu2cand->energy();
      nu2_px = nu2cand->px();
      nu2_py = nu2cand->py();
      nu2_pz = nu2cand->pz();

      mu2_eta = mu2cand->eta();
      mu2_phi = mu2cand->phi();
      nu2_eta = nu2cand->eta();
      nu2_phi = nu2cand->phi();

      htoWW_energy = htoWWcand->energy();
      htoWW_px = htoWWcand->px();
      htoWW_py = htoWWcand->py();
      htoWW_pz = htoWWcand->pz();
      htoWW_mass = htoWWcand->mass();
     }
    if (htobb){
      b1_energy = b1cand->energy();
      b1_px = b1cand->px();
      b1_py = b1cand->py();
      b1_pz = b1cand->pz();
      b1_pt = b1cand->pt();
      b1_eta = b1cand->eta();
      b1_phi = b1cand->phi();
      b2_energy = b2cand->energy();
      b2_px = b2cand->px();
      b2_py = b2cand->py();
      b2_pz = b2cand->pz();
      b2_pt = b2cand->pt();
      b2_eta = b2cand->eta();
      b2_phi = b2cand->phi();

      htobb_energy = htoBBcand->energy();
      htobb_px = htoBBcand->px();
      htobb_py = htoBBcand->py();
      htobb_pz = htoBBcand->pz();
      htobb_mass = htoBBcand->mass();
   
     }
     if (hasbbarjet){
      bbarjet_px = b2jet->px();
      bbarjet_py = b2jet->py();
      bbarjet_pz = b2jet->pz();
      bbarjet_pt = b2jet->pt();
      bbarjet_eta = b2jet->eta();
      bbarjet_phi = b2jet->phi();
      bbarjet_energy = b2jet->energy();
      bbarjet_mass = b2jet->mass();
      bbarjet_nu_energy = bbarjet_nu_lorentz.Energy();
      bbarjet_nu_px = bbarjet_nu_lorentz.Px();
      bbarjet_nu_py = bbarjet_nu_lorentz.Py();
      bbarjet_nu_pz = bbarjet_nu_lorentz.Pz();
      }
      if (hasbjet){
      bjet_px = b1jet->px();
      bjet_py = b1jet->py();
      bjet_pz = b1jet->pz();
      bjet_pt = b1jet->pt();
      bjet_eta = b1jet->eta();
      bjet_phi = b1jet->phi();
      bjet_energy = b1jet->energy();
      bjet_mass = b1jet->mass();
      bjet_nu_energy = bjet_nu_lorentz.Energy();
      bjet_nu_px = bjet_nu_lorentz.Px();
      bjet_nu_py = bjet_nu_lorentz.Py();
      bjet_nu_pz = bjet_nu_lorentz.Pz();
      }

    if (h2tohh){
      h2tohh_energy = h2tohhcand->energy();
      h2tohh_px = h2tohhcand->px();
      h2tohh_py = h2tohhcand->py();
      h2tohh_pz = h2tohhcand->pz();
      h2tohh_mass = h2tohhcand->mass();
      h2tohh_pdgid = h2tohhcand->pdgId();
      }
}


//-------------- method called to find a muon daughter for given candidate -------------------------
const reco::Candidate* 
DiHiggsWWAnalyzer::findmudaughter(const reco::Candidate *Wcand){

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
DiHiggsWWAnalyzer::findmudescendants(const reco::Candidate *cand, int& count){

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
DiHiggsWWAnalyzer::findnudaughter(const reco::Candidate *Wcand){

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
DiHiggsWWAnalyzer::findnudescendants(const reco::Candidate *cand, int& count){

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
DiHiggsWWAnalyzer::stabledecendant(const reco::Candidate* cand, int id){
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
DiHiggsWWAnalyzer::finddecendant(const reco::Candidate* cand, int id, bool first){
   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++){
        
	if ((cand->daughter(i))->pdgId() == id && first && cand->pdgId() != id)
		return 	tmp=cand->daughter(i);
	else if ((cand->daughter(i))->pdgId() == id && !first && !hasDaughter(cand->daughter(i), id)) 
		return  tmp=cand->daughter(i); // tmp does not has daughter with pdgid = id
	else if ((cand->daughter(i))->pdgId() == id && !first && (cand->daughter(i))->numberOfDaughters()>1) 
		return  tmp=cand->daughter(i);// tmp has more one daughters therefore it is final-states
        else if (finddecendant(cand->daughter(i),id, first)) 
		return tmp=finddecendant(cand->daughter(i),id,first);
   }
    
    return tmp;

}

//---------- method called to find a ancestor with pdgid = id, 
//if first is true, then return the candidate closest to seed
//if first is false, then return the candidate furthest to seed
const reco::Candidate*
DiHiggsWWAnalyzer::findancestor(const reco::Candidate* cand, int id, bool first){

   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfMothers(); i++){
        
	if ((cand->mother(i))->pdgId() == id && first && cand->pdgId() != id)
		return 	tmp=cand->mother(i);
	else if ((cand->mother(i))->pdgId() == id && !first && !hasMother(cand->mother(i), id)) 
		return  tmp=cand->mother(i);
        else if (findancestor(cand->mother(i),id, first)) 
		return tmp=findancestor(cand->mother(i),id, first);
   }
   return tmp;

}

//---------- method called to check whether cand has mother with pdgid = id -----------------------------------
bool
DiHiggsWWAnalyzer::hasMother(const reco::Candidate* cand, int id){

   for (unsigned int i=0; i < cand->numberOfMothers(); i++)
        if ((cand->mother(i))->pdgId() == id) return true;
   return false;

}

//-------- method called to check whether cand has daughter with pdgid = id ------------------------------------
bool
DiHiggsWWAnalyzer::hasDaughter(const reco::Candidate* cand, int id){
 
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        if ((cand->daughter(i))->pdgId() == id) return true;
   return false;

}

//----------- method called to calculate total lorentz vector for neutrinos decendants ---------------
//
void 
DiHiggsWWAnalyzer::neutrinoDecendantsCandidates(const reco::Candidate *cand, std::vector<const reco::Candidate*> &candvec){
  

   for (unsigned i = 0; i < cand->numberOfDaughters(); i++){
	if ((cand->daughter(i))->status() == 1 && (abs((cand->daughter(i))->pdgId()) == 12 
	     || abs((cand->daughter(i))->pdgId()) == 14 || abs((cand->daughter(i))->pdgId()) == 16) ){
		if (candvec.size() == 0) candvec.push_back(cand->daughter(i));
		else if (std::find(candvec.begin(), candvec.end(), cand->daughter(i)) == candvec.end()) 
			candvec.push_back(cand->daughter(i));
	}else neutrinoDecendantsCandidates(cand->daughter(i), candvec);
   }

}



//----------- method called to calculate total lorentz vector from b bar jets ---------------
//
void DiHiggsWWAnalyzer::stableDecendantsCandidates(const reco::Candidate *cand, std::vector<const reco::Candidate*> &candvec){
  
   for (unsigned i = 0; i < cand->numberOfDaughters(); i++){
	if ((cand->daughter(i))->status() == 1){
		if (candvec.size() == 0) candvec.push_back(cand->daughter(i));
		else if (std::find(candvec.begin(), candvec.end(), cand->daughter(i)) == candvec.end()) 
			candvec.push_back(cand->daughter(i));
	}else stableDecendantsCandidates(cand->daughter(i), candvec);
   }

}

//-------------  method called to calculate MET in simuation ------------------ 
TLorentzVector 
DiHiggsWWAnalyzer::calculateMET(){

   TLorentzVector METlorentz = TLorentzVector();
   TVector2 met_pxpy(nu1cand->px()+nu2cand->px(), nu1cand->py()+nu2cand->py());
   METlorentz.SetPxPyPzE(nu1cand->px()+nu2cand->px(), nu1cand->py()+nu2cand->py(),0,met_pxpy.Mod());

   return METlorentz;
}




// ------------ method called for printing additional information, useful for debugging  ------------
void
DiHiggsWWAnalyzer::print() {
   //todo: combine all necessary print out for debug 
    std::cout << "print() " << std::endl;
    if (!h2tohh)  {
       
        std::cout << "print() error " << std::endl;
        return;
    }
    
    


}

//------------- method called for printing h->WW->mumununu chain -------------------------
void 
DiHiggsWWAnalyzer::printHtoWWChain(){
          
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
    if (htoWW){
             
               TLorentzVector h_p4(htoWW_px, htoWW_py, htoWW_pz, htoWW_energy);
               TLorentzVector WW_p4(mu1_W1_cand->px()+mu2_W2_cand->px(),mu1_W1_cand->py()+mu2_W2_cand->py(), 
                                    mu1_W1_cand->pz()+mu2_W2_cand->pz(),mu1_W1_cand->energy()+mu2_W2_cand->energy());
               std::cout << "invariant mass from h_p4: " << h_p4.M() 
                         << "	, from WW_p4 " << WW_p4.M() 
                         << "	, from final_p4 " << final_p4.M() << std::endl;
               if (abs(WW_p4.M()-h_p4.M())>1)  std::cout << "h->WW, invariant mass reconstruction discrepancy ? " << std::endl; 
               std::cout <<  " H -> WW " << std::endl;
               const reco::Candidate* tmp_cand1 = NULL;
               const reco::Candidate* tmp_cand2 = NULL;
               
               for (unsigned int n = 0; n<mu1_htoWW_cand->numberOfDaughters(); n++){                                  
                      if ((mu1_htoWW_cand->daughter(n))->pdgId()==-24) tmp_cand1 =  mu1_htoWW_cand->daughter(n);
                      if ((mu1_htoWW_cand->daughter(n))->pdgId()== 24) tmp_cand2 =  mu1_htoWW_cand->daughter(n);
                      if (n >= 2) std::cout << "h has more 2 daughters, id " << (mu1_htoWW_cand->daughter(n))->pdgId() << std::endl;
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
                      std::cout << " W- mass from W- Candidate " << w1_mass << " from mu-,nu " << cand1_lorentz.M() << std::endl;
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
                      std::cout << " W+ mass from W+ Candidate " << w2_mass << " from mu+,nu " << cand2_lorentz.M() << std::endl;
                      if (muon2 != -1 && tmp_cand2->status() != 1) tmp_cand2 = tmp_cand2->daughter(muon2);
                     TLorentzVector tmp = cand1_lorentz+cand2_lorentz;
                     std::cout <<"Total px " << tmp.Px() << " py " << tmp.Py() << " pz " << tmp.Pz() << " E " << tmp.Energy() << std::endl;
                     std::cout << " invariant mass from daughters of WW " << tmp.M() << std::endl;  
                     std::cout << "For Next loop status of cand1 " << tmp_cand1->status()  << "	cand2 " << tmp_cand2->status() << std::endl; 
                     std::cout << "-------------  end of this loop ----------------------" << std::endl;
              } 
            if (h2tohh){	
			std::cout <<"htoWW invariant mass " << final_p4.M() <<" four momentum " << std::endl; 
                        final_p4.Print();
                        std::cout <<"htobb invariant mass " << h_bb_p4.M() << " four momentum " << std::endl;
                        h_bb_p4.Print();
                        std::cout <<"h2tohh invariant mass " << h2_final_tot.M() <<" total momentum " << std::endl;
                        h2_final_tot.Print();
			} 

        }//end if (htoWW)

}



//---------- method called to print candidates for debug ---------------------
void
DiHiggsWWAnalyzer::printCandidate(const reco::Candidate* cand){

   std::cout <<" Candidate id: "<< cand->pdgId() << " mass: " << cand->mass() <<" (P,E)= ("<< cand->px() <<", "<< cand->py()<<", "<< cand->pz()<<", "<< cand->energy() <<")" <<"(Pt,E) = ("<< cand->pt() <<", "<< cand->eta() <<", "<< cand->phi()<<", "<<cand->energy()<< ")" <<" status: " << cand->status() << std::endl;

}

//--------- method called to print all decendants for cand -------------------
void 
DiHiggsWWAnalyzer::printallDecendants(const reco::Candidate* cand){
   
   if (cand->status() != 0 && cand->numberOfDaughters() > 0){
        std::cout << "******************  children of id "<< cand->pdgId() <<"      *********************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        	printCandidate(cand->daughter(i));
        std::cout << "***********************************************************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
		printallDecendants(cand->daughter(i));

    }
}


//--------- method called to print all Ancestors for cand -------------------
void 
DiHiggsWWAnalyzer::printallAncestors(const reco::Candidate* cand){
   
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
DiHiggsWWAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiHiggsWWAnalyzer);


