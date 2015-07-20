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
//pat::muon...
//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/MET.h"
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



#include "DiHiggsWW/DiHiggsWWAnalyzer/src/MMC.h"
#include "DiHiggsWW/DiHiggsWWAnalyzer/src/deltaR.h"


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
     bool simulation_;
     std::string genparticleLabel_;
     std::string jetLabel_;
     std::string metLabel_;
     //cuts on GenJets
     double jetsPt_;
     double jetsEta_;
     double bjetsPt_;
     double bjetsEta_;
     double jetsDeltaR_;
     double muonPt2_;
     double muonPt1_;
     double muonsEta_;
     double metPt_;
     
    private:
     //decendants and ancestor
     const reco::Candidate* stabledecendant(const reco::Candidate* cand, int id);
     //const reco::Candidate* stabletdecendant(Particle p, PdgId id);
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

      const reco::GenJet* b1jet;
      const reco::GenJet* b2jet;
    private:
      
      TLorentzVector mu1_lorentz;
      TLorentzVector mu2_lorentz;
      TLorentzVector b_genp_lorentz;
      TLorentzVector bbar_genp_lorentz;
      TLorentzVector nu1_lorentz;
      TLorentzVector nu2_lorentz;
      TLorentzVector met_lorentz;
      TLorentzVector b1jet_lorentz;
      TLorentzVector b2jet_lorentz;
      TLorentzVector totjets_lorentz;

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
      int mu2_motherid;
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
      float bjet_energy;
      float bjet_px;
      float bjet_py;
      float bjet_pz;
      float bjet_mass;
      float bbarjet_energy;
      float bbarjet_px;
      float bbarjet_py;
      float bbarjet_pz;
      float bbarjet_mass;
      float dR_bjet;
      float dR_bbarjet;

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

      float totjets_energy;
      float totjets_px;
      float totjets_py;
      float totjets_pz;

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

      bool hasbjet;
      bool hasbbarjet;
      bool hasMET;
      bool hastwomuons;
        
    private:
      bool runMMC_;
      MMC* thismmc;
};

    //
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
     simulation_ = iConfig.getParameter<bool>("simulation");
     genparticleLabel_ = iConfig.getParameter<std::string>("genparticleLabel");
     jetLabel_ = iConfig.getParameter<std::string>("jetLabel");
     metLabel_ = iConfig.getParameter<std::string>("metLabel");
     jetsPt_ = iConfig.getParameter<double>("jetsPt");
     jetsEta_ = iConfig.getParameter<double>("jetsEta");
     bjetsPt_ = iConfig.getParameter<double>("bjetsPt");
     bjetsEta_ = iConfig.getParameter<double>("bjetsEta");
     jetsDeltaR_ = iConfig.getParameter<double>("jetsDeltaR");
     muonPt1_ = iConfig.getParameter<double>("muonPt1");
     muonPt2_ = iConfig.getParameter<double>("muonPt2");
     muonsEta_ = iConfig.getParameter<double>("muonsEta");
     metPt_ = iConfig.getParameter<double>("metPt");
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
      mu2_motherid = 0;
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
   edm::Handle<reco::GenParticleCollection> genParticleColl;
   iEvent.getByLabel(genparticleLabel_, genParticleColl);
   /*
   edm::Handle<std::vector<pat::Muon>> muonColl;
   iEvent.getByLabel("slimmedMuons",muonColl); 
   edm::Handle<std::vector<pat::Jet>> jetColl;
   iEvent.getByLabel("slimmedJets",jetColl); 
   edm::Handle<std::vector<pat::MET>> metColl;
   iEvent.getByLabel("slimmedMETs",metColl); 
*/

   edm::Handle<reco::GenJetCollection> genjetColl;
   iEvent.getByLabel(jetLabel_, genjetColl);

   edm::Handle<edm::View<reco::GenMET> > genmetColl; 
   iEvent.getByLabel(metLabel_, genmetColl);

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
      
 /*     
   for (std::vector<pat::Muon>::const_iterator it = muonColl->begin(); it != muonColl->end(); ++it){
   
   //    std::cout <<"muon id " << it->pdgId() <<" status " << it->status() <<" px " << it->px() << " py " << it->py() <<std::endl;
   
   }

   for (std::vector<pat::Jet>::const_iterator jet = jetColl->begin(); jet != jetColl->end(); ++jet){
   
     //  std::cout <<"jet mass " << jet->mass() <<" status " << jet->status() <<" px " << jet->px() << " py " << jet->py() <<std::endl;
   
   }

   for (std::vector<pat::MET>::const_iterator met = metColl->begin(); met != metColl->end(); ++met){
   
       //std::cout <<"met mass " << met->mass() <<" status " << met->status() <<" px " << met->px() << " py " << met->py() <<std::endl;
   
   }
*/
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

   }// all Gen particles

    if (mu_negative and nu_negative && mu1_W1_cand == nu1_W1_cand)
       {
          Wtomu1nu1 = true;
          }
     else if (mu_negative && nu_negative)  std::cout << "find negative mu && nu but not find W1" << std::endl;
  
    if (mu_positive and nu_positive && mu2_W2_cand == nu2_W2_cand)
       {
          Wtomu2nu2 = true;
          }
     else if (mu_positive && nu_positive)  std::cout << "find positive mu && nu but not find W2" << std::endl;
  
    if (Wtomu1nu1 and bbarquark and mu1_tbar_cand == b2_tbar_cand)
       {
         std::cout << "find muon and bbar, and they come from same tbar" << std::endl;
	 if (finalStates_)
        	b2cand = finddecendant(b2_tbar_cand, -5, false);
	 else
        	b2cand = finddecendant(b2_tbar_cand, -5, true);
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
    if (Wtomu2nu2 and bquark and mu2_t_cand == b1_t_cand)
       {
         std::cout << "find muon and b, and they come from same t" << std::endl;
	 if (finalStates_)
        	b1cand = finddecendant(b1_t_cand, 5, false);
	 else
        	b1cand = finddecendant(b1_t_cand, 5, true);
	 tcand = mu2_t_cand;
         ttoWb = true;
          
         }
    else if(Wtomu2nu2 and bquark)   
       {
         std::cout << "find muon and b, and they do not come from same t" << std::endl;
         std::cout << "mu2_t energy " << mu2_t_cand->energy() << " px " << mu2_t_cand->px() << std::endl;
         std::cout << "b1_t energy " << b1_t_cand->energy() << " px " << b1_t_cand->px() << std::endl;

         }



   if (genmetColl->size() != 1){
	std::cout <<" size of genmetColl " << genmetColl->size();
	} 
   reco::GenMET genMet(genmetColl->front());

   hasbjet = false;
   hasbbarjet = false;
   hasMET =false;
   hastwomuons = false;
  // match bquark and bjet
   dR_bjet = jetsDeltaR_;
   dR_bbarjet = jetsDeltaR_;

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

   int nparts = 0;
   if (ttoWb and tbartoWbbar){
   for (reco::GenJetCollection::const_iterator jetit = genjetColl->begin(); jetit != genjetColl->end(); jetit++){
	if (jetit->pt() >= jetsPt_ and std::fabs(jetit->eta()) <= jetsEta_){	
 		totjets_px += jetit->px();
		totjets_py += jetit->py();
		totjets_pz += jetit->pz();
		totjets_energy += jetit->energy();

	}
	//cuts on GenJets
	if (jetit->pt()<bjetsPt_ or std::fabs(jetit->eta())> bjetsEta_) continue;

	//printCandidate(jetit->clone());
        std::vector <const reco::GenParticle*> mcparts = jetit->getGenConstituents();
  	for (unsigned i = 0; i < mcparts.size(); i++) {
    		const reco::GenParticle* mcpart = mcparts[i];
		const reco::Candidate* bcand;
		const reco::Candidate* tmpcand;
		if ((bcand=findancestor(mcpart, -5, false)) and bcand->mother()->pdgId() == -6){ 
			tmpcand = bcand->mother();
			if (tmpcand == tbarcand and dR_bbarjet > deltaR(jetit->eta(), jetit->phi(), b2cand->eta(), b2cand->phi())){
				dR_bbarjet = deltaR(jetit->eta(), jetit->phi(), b2cand->eta(), b2cand->phi());
				//printCandidate(jetit->clone());
				//std::cout <<" bcand(-5) "; printCandidate(bcand);
				std::cout <<" has tbar->bbar,tbar is the same from tbar->bbar in genparticles flow,dR "<< dR_bbarjet <<std::endl;
				bbarjet_px = jetit->px();
				bbarjet_py = jetit->py();
				bbarjet_pz = jetit->pz();
				bbarjet_energy = jetit->energy();
				bbarjet_mass = jetit->mass();
   				b2jet = jetit->clone();
				hasbbarjet = true;
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
				//std::cout <<" nparts inside loop " << nparts <<" dR_bbarjet "<< dR_bbarjet << std::endl;
			}
		  }
		if ((bcand=findancestor(mcpart, 5, false)) and bcand->mother()->pdgId() == 6){
			tmpcand = bcand->mother();
			if (tmpcand == tcand and dR_bjet>deltaR(jetit->eta(), jetit->phi(), b1cand->eta(), b1cand->phi())){
				dR_bjet = deltaR(jetit->eta(), jetit->phi(), b1cand->eta(), b1cand->phi());
				//printCandidate(jetit->clone());
				//std::cout <<" bcand(5) "; printCandidate(bcand);
				std::cout <<" has t->b, t is the same as t->b in genparticles flow,dR "<< dR_bjet <<std::endl;
				bjet_px = jetit->px();
				bjet_py = jetit->py();
				bjet_pz = jetit->pz();
				bjet_energy = jetit->energy();
				bjet_mass = jetit->mass();
				b1jet = jetit->clone();
				hasbjet = true;
				bjet_decendant_px = 0;
				bjet_decendant_py = 0;
				bjet_decendant_pz = 0;
				bjet_decendant_energy = 0;
				nparts = 0;
			}
			if (hasbjet and dR_bjet == deltaR(jetit->eta(), jetit->phi(), b1cand->eta(), b1cand->phi())){
				bjet_decendant_px += mcpart->px();
				bjet_decendant_py += mcpart->py();
				bjet_decendant_pz += mcpart->pz();
				bjet_decendant_energy += mcpart->energy();
				nparts++;
				//std::cout <<" nparts inside loop " << nparts <<" dR_bjet "<< dR_bjet << std::endl;
			}
		 }	
	}
	std::cout <<"mcparticle size " <<mcparts.size() <<"   nparts "<< nparts << std::endl;
     }
	if (!hasbbarjet or !hasbjet) std::cout <<" has ttbar->WbWbbar, but failed to two b jets " << std::endl;
    }//

   if (ttoWb and hasbjet){
	b1jet_lorentz.SetXYZT(b1jet->px(), b1jet->py(), b1jet->pz(), b1jet->energy());	

     }
   if (tbartoWbbar){
	b2jet_lorentz.SetXYZT(b2jet->px(), b2jet->py(), b2jet->pz(), b2jet->energy());	

    }

   if (ttoWb and tbartoWbbar) {
        std::cout << "find t and tbar candidate " << std::endl;
        if (finalStates_){ 
    	mu1cand = stabledecendant(mu1_W1_cand, 13);
        mu2cand = stabledecendant(mu2_W2_cand, -13);
        nu1cand = stabledecendant(mu1_W1_cand, -14);
        nu2cand = stabledecendant(mu2_W2_cand, 14);
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

   	mu1_lorentz.SetPtEtaPhiM(mu1cand->pt(), mu1cand->eta(), mu1cand->phi(), 0);
   	nu1_lorentz.SetPtEtaPhiM(nu1cand->pt(), nu1cand->eta(), nu1cand->phi(), 0);
        mu2_lorentz.SetPtEtaPhiM(mu2cand->pt(), mu2cand->eta(), mu2cand->phi(), 0); 
        nu2_lorentz.SetPtEtaPhiM(nu2cand->pt(), nu2cand->eta(), nu2cand->phi(), 0); 
        totjets_lorentz.SetXYZT(totjets_px, totjets_py, totjets_pz, totjets_energy);
        met_lorentz = calculateMET();
        std::cout <<"MET from nuetrinos "; met_lorentz.Print(); 
        met_lorentz.SetXYZT(genMet.px(), genMet.py(), 0, genMet.pt()); 
        std::cout <<"MET from genMet "; met_lorentz.Print();
        std::cout <<" b1 " ; printCandidate(b1cand);
	std::cout<<"bjet (P,E)=( " << bjet_px << ", "<< bjet_py << ", "<<bjet_pz << ", " << bjet_energy <<")" << std::endl; 
	std::cout<<"bjet_tot (P,E)=( " << bjet_px_tot << ", "<< bjet_py_tot << ", "<<bjet_pz_tot << ", " << bjet_energy_tot <<")" << std::endl; 
	std::cout<<"bjet_bdecendant (P,E)=( " << bjet_decendant_px << ", "<< bjet_decendant_py << ", "<<bjet_decendant_pz << ", " << bjet_decendant_energy <<")" << std::endl; 
        std::cout <<" b2 " ; printCandidate(b2cand);
	std::cout<<"bbarjet (P,E)=( " << bbarjet_px << ", "<< bbarjet_py << ", "<<bbarjet_pz << ", " << bbarjet_energy <<")" << std::endl; 
	std::cout<<"bbarjet_tot (P,E)=( " << bbarjet_px_tot << ", "<< bbarjet_py_tot << ", "<<bbarjet_pz_tot << ", " << bbarjet_energy_tot <<")" << std::endl; 
	std::cout<<"bbarjet_bdecendant (P,E)=( " << bbarjet_decendant_px << ", "<< bbarjet_decendant_py << ", "<<bbarjet_decendant_pz << ", " << bbarjet_decendant_energy <<")" << std::endl; 

        /* 
	const reco::Candidate* t_tmp = tcand->mother();
	while (t_tmp->pdgId() == 6) t_tmp = t_tmp->mother();
	const reco::Candidate* tbar_tmp = tbarcand->mother();
	while (tbar_tmp->pdgId() == -6) tbar_tmp = tbar_tmp->mother();
        std::cout <<" tbar_tmp and its decendants" <<std::endl; printallDecendants(tbar_tmp);
        std::cout <<" t_tmp and its decendants" <<std::endl; printallDecendants(t_tmp);
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
        //std::cout <<" bbarjet, mass "<<bbarjet_lorentz.M(); bbarjet_lorentz.Print();
        std::cout <<" t->Wb " ; printCandidate(tcand);
        std::cout <<" tbar->Wbbar " ; printCandidate(tbarcand);
        //std::cout <<"another b1 " ; printCandidate(finddecendant(tbarcand, 5, false));
        //std::cout <<"another b2 " ; printCandidate(finddecendant(tbarcand, -5, false));
       //std::cout <<" b jet"; bjetsLorentz(b1cand).Print(); 
        //std::cout <<" b jet"; bjetsLorentz(finddecendant(tbarcand, 5, false)).Print(); 
        //std::cout <<" b jet"; stableDecendantsLorentz(b1cand).Print(); 
        //std::cout <<" bbar jet"; stableDecendantsLorentz(b2cand).Print(); 
        //std::cout <<" t " ; printCandidate(tcand);
        //std::cout <<" tbar and its ancestors" <<std::endl; printallAncestors(tbarcand);
        //std::cout <<" t and its ancestors" <<std::endl; printallAncestors(tcand);
        //std::cout <<" tbar and its decendants" <<std::endl; printallDecendants(tbarcand);
        //std::cout <<" t and its decendants" <<std::endl; printallDecendants(tcand);
        if (fabs(TLorentzVector(mu1_lorentz+nu1_lorentz+bbarjet_lorentz).M()-tbarcand->mass())>0.5 or  
        fabs(TLorentzVector(mu2_lorentz+nu2_lorentz+bjet_lorentz).M()-tcand->mass())>0.5  ) 
	 {
	        std::cout << "reconstructed t mass " << TLorentzVector(mu1_lorentz+nu1_lorentz+bbarjet_lorentz).M() <<std::endl;
		std::cout <<" reconstructed tbar mass"<< TLorentzVector(mu2_lorentz+nu2_lorentz+bjet_lorentz).M() << std::endl;
	        printallDecendants(t_tmp);
	 
	 }*/
        //std::cout <<" tbar and its decendants" <<std::endl; printCandidate(tbarcand);
       
        fillbranches();
        evtree->Fill();
        //debug
     }
	//

   if (genMet.pt() >= metPt_) hasMET =true;
   //if (h2tohh && runMMC_) runMMC();
   if (ttoWb and tbartoWbbar and  runMMC_){
        b_genp_lorentz.SetXYZT(b1cand->px(), b1cand->py(), b1cand->pz(), b1cand->energy()); 
        bbar_genp_lorentz.SetXYZT(b2cand->px(), b2cand->py(), b2cand->pz(), b2cand->energy()); 
        int onshellMarker = -1;
        if (w1cand->mass() > w2cand->mass()) onshellMarker = 1;
        else onshellMarker = 2;
        // std::cout <<" mu1 lorenz "; mu1_lorentz.Print(); 
        //thismmc = new MMC();
        //std::cout << "onshellMarkder  " << onshellMarker << std::endl;
        TLorentzVector tot_lorentz= mu1_lorentz+mu2_lorentz+nu1_lorentz+nu2_lorentz+b_genp_lorentz+bbar_genp_lorentz;
	thismmc = new MMC(&mu1_lorentz, &mu2_lorentz, &b1jet_lorentz, &b2jet_lorentz, &totjets_lorentz, &met_lorentz, 
	&nu1_lorentz, &nu2_lorentz,&b_genp_lorentz, &bbar_genp_lorentz, &tot_lorentz, onshellMarker, 
	simulation_, ievent, mmcset_, fs, verbose_);
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
   evtree->Branch("bjet_energy",&bjet_energy);
   evtree->Branch("bjet_px",&bjet_px);
   evtree->Branch("bjet_py",&bjet_py);
   evtree->Branch("bjet_pz",&bjet_pz);
   evtree->Branch("bjet_mass",&bjet_mass);
   evtree->Branch("bbarjet_energy",&bbarjet_energy);
   evtree->Branch("bbarjet_px",&bbarjet_px);
   evtree->Branch("bbarjet_py",&bbarjet_py);
   evtree->Branch("bbarjet_pz",&bbarjet_pz);
   evtree->Branch("bbarjet_mass",&bbarjet_mass);
   evtree->Branch("dR_bjet",&dR_bjet);
   evtree->Branch("dR_bbarjet",&dR_bbarjet);

   evtree->Branch("bjet_energy_tot",&bjet_energy_tot);
   evtree->Branch("bjet_px_tot",&bjet_px_tot);
   evtree->Branch("bjet_py_tot",&bjet_py_tot);
   evtree->Branch("bjet_pz_tot",&bjet_pz_tot);
   evtree->Branch("bbarjet_energy_tot",&bbarjet_energy_tot);
   evtree->Branch("bbarjet_px_tot",&bbarjet_px_tot);
   evtree->Branch("bbarjet_py_tot",&bbarjet_py_tot);
   evtree->Branch("bbarjet_pz_tot",&bbarjet_pz_tot);
   
   evtree->Branch("bjet_decendant_energy",&bjet_decendant_energy);
   evtree->Branch("bjet_decendant_px",&bjet_decendant_px);
   evtree->Branch("bjet_decendant_py",&bjet_decendant_py);
   evtree->Branch("bjet_decendant_pz",&bjet_decendant_pz);
   evtree->Branch("bbarjet_decendant_energy",&bbarjet_decendant_energy);
   evtree->Branch("bbarjet_decendant_px",&bbarjet_decendant_px);
   evtree->Branch("bbarjet_decendant_py",&bbarjet_decendant_py);
   evtree->Branch("bbarjet_decendant_pz",&bbarjet_decendant_pz);
   
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
ttbarAnalyzer::finddecendant(const reco::Candidate* cand, int id, bool first){
   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++){
        
	if ((cand->daughter(i))->pdgId() == id && first && cand->pdgId() != id)
		return 	tmp=cand->daughter(i);
	else if ((cand->daughter(i))->pdgId() == id && !first && !hasDaughter(cand->daughter(i), id)) 
		return  tmp=cand->daughter(i); // tmp does not has daughter with pdgid = id
	else if ((cand->daughter(i))->pdgId() == id && !first && (cand->daughter(i))->numberOfDaughters()>1) 
		return  tmp=cand->daughter(i);// tmp has more one daughters therefore it is final-states
        else if (finddecendant(cand->daughter(i),id, first)) 
		return tmp=finddecendant(cand->daughter(i),id, first);
   }
    
    return tmp;

}

//---------- method called to find a ancestor with pdgid = id, 
//if first is true, then return the candidate closest to seed
//if first is false, then return the candidate furthest to seed
const reco::Candidate*
ttbarAnalyzer::findancestor(const reco::Candidate* cand, int id, bool first){

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

//----------- method called to calculate total lorentz vector for neutrinos decendants ---------------
//
void 
ttbarAnalyzer::neutrinoDecendantsCandidates(const reco::Candidate *cand, std::vector<const reco::Candidate*> &candvec){
  

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
void 
ttbarAnalyzer::stableDecendantsCandidates(const reco::Candidate *cand, std::vector<const reco::Candidate*> &candvec){
  
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


