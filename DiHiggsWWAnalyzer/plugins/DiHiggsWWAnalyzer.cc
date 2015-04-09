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

//headers from root lib
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH1F.h"

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
    
   //runMMC
   private: 
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
          
   private:
      edm::ParameterSet cfg_;
      // debuglevel constrol 
      int verbose_; 
      void print();
      void printHtoWWChain();
      void printCandidate(const reco::Candidate* );
      void printallDecendants(const reco::Candidate* );
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
     //const reco::Candidate* stablehtoWWdecendant(Particle p, PdgId id);
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
      float virtualW_lowM;
      float virtualW_highM;

     
    private:
      bool runMMC_;
     // MMC tree branches
    private:
      TLorentzVector stableDecendantsLorentz(const reco::Candidate* cand); 
   
    private:
      TLorentzVector calculateMET(); 
    private:
      float onshellWMassRandomWalk(float x0, float step, float random);
      float onshellWMassRandomWalk(float x0, float step, float random, TH1F* hist);
      float onshellWMassPDF(float wmass);
  
    private:
      TH1F* readoutonshellWMassPDF();
      TH1F* readoutonshellnuptPDF();
 
    private:
      float weightfromonshellnupt(TH1F* pdf, float nupt); 
      float weightfromonshellnupt(float nupt); 
   
    private:
      bool weightfromonshellnupt_func_;
      bool weightfromonshellnupt_hist_;

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
      TLorentzVector* jets_lorentz;
      TLorentzVector* met_lorentz;
      TVector2* MMCmet_vec2;
      TLorentzVector* nu_onshellW_lorentz;
      TLorentzVector* nu_offshellW_lorentz;
      TLorentzVector* offshellW_lorentz;
      TLorentzVector* onshellW_lorentz;
      TLorentzVector* htoWW_lorentz;
      TLorentzVector* htoBB_lorentz;
      TLorentzVector* h2tohh_lorentz;
      
      int control;
      float weight;
      float weight1;//extra weight
 
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
     
      float htoBB_Eta;
      float htoBB_Phi;
      float htoBB_Pt;
      float htoBB_E;
      float htoBB_Mass;
      float htoWW_Eta;
      float htoWW_Phi;
      float htoWW_Pt;
      float htoWW_E;
      float htoWW_Mass;

      float MMCmet_E;
      float MMCmet_Phi;
      float MMCmet_Px;
      float MMCmet_Py;

      float h2tohh_Eta;
      float h2tohh_Phi;
      float h2tohh_Pt;
      float h2tohh_E;
      float h2tohh_Mass;

      float met;
      float met_phi;
      float met_px;
      float met_py;

      float eta_nuoffshellW_true;
      float phi_nuoffshellW_true;
      float pt_nuoffshellW_true;
      float eta_nuonshellW_true;
      float phi_nuonshellW_true;
      float pt_nuonshellW_true;
      float mass_offshellW_true;
      float mass_onshellW_true;
      float pt_h2tohh_true;

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
     runMMC_ = iConfig.getParameter<bool>("runMMC");
     finalStates_ = iConfig.getParameter<bool>("finalStates");
     weightfromonshellnupt_func_ = iConfig.getParameter<bool>("weightfromonshellnupt_func");
     weightfromonshellnupt_hist_ = iConfig.getParameter<bool>("weightfromonshellnupt_hist");
     iterations_ = iConfig.getUntrackedParameter<int>("iterations",100000);
     seed_ = iConfig.getParameter<int>("seed");
     RefPDFfile_ = iConfig.getParameter<std::string>("RefPDFfile");
     std::cout <<" RefPDFfile_ " << RefPDFfile_ << std::endl;
   // initilize candidates pointer



   //now do what ever initialization is needed
      ievent = 0;
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
      jets_energy = 0.0;
      jets_px = 0.0;
      jets_py = 0.0;
      jets_pz = 0.0;
      jets_mass = 0.0;

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

      mu_positive = false;
      mu_negative = false;
      bquark = false;
      bbarquark = false;
      htobb = false;
      htoWW = false;
      virtualW_lowM = 25;
      virtualW_highM = 45;


     // runMMC
      mu_onshellW_lorentz = new TLorentzVector();
      mu_offshellW_lorentz = new TLorentzVector();
      jets_lorentz = new TLorentzVector();
      met_lorentz = new TLorentzVector();
      MMCmet_vec2 = new TVector2();
      nu_onshellW_lorentz = new TLorentzVector();
      nu_offshellW_lorentz = new TLorentzVector();
      offshellW_lorentz = new TLorentzVector();
      onshellW_lorentz = new TLorentzVector();
      htoWW_lorentz = new TLorentzVector();
      htoBB_lorentz = new TLorentzVector();
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
     
      htoBB_Eta =0;
      htoBB_Phi =0;
      htoBB_Pt =0;
      htoBB_E =0;
      htoBB_Mass =0;
      htoWW_Eta =0;
      htoWW_Phi =0;
      htoWW_Pt =0;
      htoWW_E =0;
      htoWW_Mass =0;

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
   Handle<reco::GenParticleCollection> genParticleColl;
   iEvent.getByLabel("genParticles", genParticleColl);
    
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
         
          if (tmp_mu1->pdgId() == 25)   
             {
                  //while (tmp_mu1->mother()->pdgId() == 25) tmp_mu1 = tmp_mu1->mother();
                  std::cout << "find muon(-) candidate" << std::endl;
	          mu_negative = true;
                  mu1_htoWW_cand = tmp_mu1;
                 // std::cout << "mother of this higgs, id " << tmp_mu1->mother()->pdgId() << " energy " << tmp_mu1->mother()->energy() << std::endl;
              }
      }
    /* else if (abs(it->pdgId()) == 13) //for test
      {
          std::cout << "find a muon, id " << it->pdgId() << " energy "<< it->energy() << " status "  << it->status() << std::endl;
          const reco::Candidate* tmp=it->mother(); 
          while (abs(tmp->pdgId()) == 13) tmp = tmp->mother();
          std::cout << "above muon, id " << tmp->pdgId() << " energy " << tmp->energy()  << " status " << tmp->status() << std::endl;
          
          while (abs(tmp->pdgId()) == 24) tmp = tmp->mother();
          std::cout << "above W, id " << tmp->pdgId() << " energy " << tmp->energy() <<" status " << tmp->status() << std::endl;
         
       }*/
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
          if (tmp_mu2->pdgId() == 25)   
             {
                // while (tmp_mu2->mother()->pdgId() == 25)   tmp_mu2 = tmp_mu2->mother();
                 std::cout << "find muon(+) candidate" << std::endl;
	         mu_positive = true;
                 mu2_htoWW_cand = tmp_mu2;
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
          if (tmp_nu1->pdgId() == 25)
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
          if (tmp_nu2->pdgId() == 25)
             {
              //   std::cout << "find antinuetrino candidate" << std::endl;
                 nu_positive = true;
                }
         
       }
      else if (it->pdgId() == 5 && it->mother()->pdgId() == 25 && !bquark)
      {
	  bquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
       //   std::cout << "find bquark candidate" << std::endl;
          b1_htobb_cand = it->mother();
      }
      else if (it->pdgId() == -5 && it->mother()->pdgId() == 25 && !bbarquark)
      {
	  bbarquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
         // std::cout << "find bbarquark candidate" << std::endl;
          b2_htobb_cand = it->mother();
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
  
    if (Wtomu1nu1 and Wtomu2nu2 and mu1_htoWW_cand == mu2_htoWW_cand)
       {
         std::cout << "find 2 muons and they come from same higgs" << std::endl;
    
          htoWW = true;
          //if (abs(final_p4.M()-125)>5)    
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
          if (tmp_htoWW == tmp_htoBB && tmp_htoWW->pdgId() == 99927){

              std::cout << "find 2 higgs and both of them come from same heavey higgs"  << std::endl;
              h2tohh = true;
              h2tohh_cand = tmp_htoWW;
              //
              h2tohhcand = tmp_htoWW;
              htoWWcand = mu1_htoWW_cand;
              htoBBcand = b1_htobb_cand;
              //print();
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
   

   if (h2tohh) {
        std::cout << "find h2 candidate " << std::endl;
        std::cout << "h2 candidate id " << h2tohhcand->pdgId() << " mass " << h2tohhcand->mass() << std::endl;
        if (finalStates_){ 
    	mu1cand = stabledecendant(htoWWcand, 13);
        mu2cand = stabledecendant(htoWWcand, -13);
        nu1cand = stabledecendant(htoWWcand, -14);
        nu2cand = stabledecendant(htoWWcand, 14);
        b1cand = alldecendant(htoBBcand, 5, false);
        b2cand = alldecendant(htoBBcand, -5, false);
        w1cand = allancestor(mu1cand, -24, false);
	w2cand = allancestor(mu2cand, 24, false);   
        }else {
        mu1cand = findmudaughter(mu1_W1_cand);
        nu1cand = findnudaughter(mu1_W1_cand);
        mu2cand = findmudaughter(mu2_W2_cand);
        nu2cand = findnudaughter(mu2_W2_cand);

        b1cand = alldecendant(htoBBcand, 5, true);
        b2cand = alldecendant(htoBBcand, -5, true);
        w1cand = allancestor(mu1cand, -24, true);
	w2cand = allancestor(mu2cand, 24, true);   
        }
        //if ()
        std::cout <<" h2tohh "; printCandidate(h2tohhcand);
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
        std::cout <<" b1 " ; printCandidate(b1cand);
        std::cout <<" b2 " ; printCandidate(b2cand);
        //std::cout <<"another b1 " ; printCandidate(alldecendant(htoBBcand, 5, false));
        //std::cout <<"another b2 " ; printCandidate(alldecendant(htoBBcand, -5, false));
        //std::cout <<" b jet"; bjetsLorentz(b1cand).Print(); 
        //std::cout <<" b jet"; bjetsLorentz(alldecendant(htoBBcand, 5, false)).Print(); 
        //std::cout <<" bbar jet"; bjetsLorentz(b2cand).Print(); 
        //std::cout <<" bbar jet"; bjetsLorentz(alldecendant(htoBBcand, -5, false)).Print(); 
        //std::cout <<" htoWW " ; printCandidate(htoWWcand);
        //std::cout <<" htoBB and its decendants" <<std::endl; printCandidate(htoBBcand);
	//printallDecendants(h2tohhcand);
        *jets_lorentz = stableDecendantsLorentz(htoBBcand);
        *met_lorentz = calculateMET();
       
        fillbranches();
        evtree->Fill();
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

    
   if (h2tohh && runMMC_) runMMC();
       
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

   evtree->Branch("htoWW_energy",&htoWW_energy);
   evtree->Branch("htoWW_px",&htoWW_px);
   evtree->Branch("htoWW_py",&htoWW_py);
   evtree->Branch("htoWW_pz",&htoWW_pz);
   evtree->Branch("htoWW_mass",&htoWW_mass);
   
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

   evtree->Branch("met",&met);
   evtree->Branch("met_phi",&met_phi);
   evtree->Branch("met_px",&met_px);
   evtree->Branch("met_py",&met_py);
   
   evtree->Branch("htobb",&htobb);
   evtree->Branch("htoWW",&htoWW);
   evtree->Branch("h2tohh",&h2tohh);
    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiHiggsWWAnalyzer::endJob() 
{
//   std::cout << "endJob, ievent  " << ievent << std::endl;
    //output->Write();
   // output->Close();

   // release space 
   delete mu_onshellW_lorentz;
   delete mu_offshellW_lorentz;
   delete nu_onshellW_lorentz;
   delete nu_offshellW_lorentz;
   delete onshellW_lorentz;
   delete offshellW_lorentz;
   delete htoWW_lorentz;
   delete htoBB_lorentz;
   delete h2tohh_lorentz;
   delete jets_lorentz;
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

      htoWW_energy = htoWWcand->energy();
      htoWW_px = htoWWcand->px();
      htoWW_py = htoWWcand->py();
      htoWW_pz = htoWWcand->pz();
      htoWW_mass = htoWWcand->mass();

      b1_energy = b1cand->energy();
      b1_px = b1cand->px();
      b1_py = b1cand->py();
      b1_pz = b1cand->pz();
      b2_energy = b2cand->energy();
      b2_px = b2cand->px();
      b2_py = b2cand->py();
      b2_pz = b2cand->pz();

      jets_energy = jets_lorentz->Energy();
      jets_px = jets_lorentz->Px();
      jets_py = jets_lorentz->Py();
      jets_pz = jets_lorentz->Pz();
      jets_mass = jets_lorentz->M();
      
      htobb_energy = htoBBcand->energy();
      htobb_px = htoBBcand->px();
      htobb_py = htoBBcand->py();
      htobb_pz = htoBBcand->pz();
      htobb_mass = htoBBcand->mass();
     
      
      h2tohh_energy = h2tohhcand->energy();
      h2tohh_px = h2tohhcand->px();
      h2tohh_py = h2tohhcand->py();
      h2tohh_pz = h2tohhcand->pz();
      h2tohh_mass = h2tohhcand->mass();
   
      met = met_lorentz->Energy();
      met_phi = met_lorentz->Phi();
      met_px = met_lorentz->Px();
      met_py = met_lorentz->Py();
   
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
DiHiggsWWAnalyzer::alldecendant(const reco::Candidate* cand, int id, bool first){
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
DiHiggsWWAnalyzer::allancestor(const reco::Candidate* cand, int id, bool first){

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

//================================================================================================================
//
// runMMC algo
//
//
//================================================================================================================
//----------- method called to run MMC method for each case -------------------
// control 0 : take muon from onshellW as muon from onshell W and nu_offshellW_eta = some_eta+deltaeta
// control 1 : take muon from onshellW  as muon from onshell W and nu_offshellW_eta = some_eta-deltaeta
// control 2 : take muon from offshellW as muon from onshell W and nu_offshellW_eta = some_eta+deltaeta
// control 3 : take muon from offshellW as muon from onshell W and nu_offshellW_eta = some_eta-deltaeta
void 
DiHiggsWWAnalyzer::runMMC(){

//   TTree *mmctree = new TTree(); 
   eta_mean = 0;
   eta_rms = 0;
// genetated (eta,phi) pair
   eta_gen = 0;
   phi_gen = 0;
   std::stringstream ss;
   ss << "mmctree_" << ievent;
   const std::string name(ss.str());
   TTree *mmctree = fs->make<TTree>(name.c_str(), name.c_str());
   initTree(mmctree);
   mmctree->SetDebug(100,0,9999999);
   //int count = 100000;

   eta_mean=0;
   eta_rms=1.403;
   TRandom3 *generator = new TRandom3();
   generator->SetSeed(seed_+ievent);
   //TF1* wmasspdf = new TF1("wmasspdf","exp(x*7.87e-3+1.69)+603.47*exp(-0.5*((x-80.1)/2.0)**2)",50,90);

  // later should take into consideration both possible cases
  // int onshell_control = 0;
        
   htoBB_lorentz->SetXYZT(b1_px+b2_px,b1_py+b2_py,b1_pz+b2_pz,b1_energy+b2_energy);
  // *jets_lorentz = bjetsLorentz(htoBBcand);
  // htoBB_lorentz = jets_lorentz;
   htoBB_Eta = htoBB_lorentz->Eta();
   htoBB_Phi = htoBB_lorentz->Phi();
   htoBB_Pt = htoBB_lorentz->Pt();
   htoBB_E = htoBB_lorentz->E();
   htoBB_Mass = htoBB_lorentz->M();
   
   std::cout <<" jets_lorentz "; jets_lorentz->Print();  
   std::cout <<" htoBB_lorentz "; htoBB_lorentz->Print();  
   if (mu1_mother_mass > mu2_mother_mass){
	eta_nuoffshellW_true = nu2cand->eta();
	phi_nuoffshellW_true = nu2cand->phi();
	pt_nuoffshellW_true = nu2cand->pt();
	eta_nuonshellW_true = nu1cand->eta();
	phi_nuonshellW_true = nu1cand->phi();
	pt_nuonshellW_true = nu1cand->pt();
        mass_offshellW_true = mu2_mother_mass;
        mass_onshellW_true = mu1_mother_mass;
    }
    else{
	eta_nuoffshellW_true = nu1cand->eta();
	phi_nuoffshellW_true = nu1cand->phi();
	pt_nuoffshellW_true = nu1cand->pt();
	eta_nuonshellW_true = nu2cand->eta();
	phi_nuonshellW_true = nu2cand->phi();
	pt_nuonshellW_true = nu2cand->pt();
        mass_offshellW_true = mu1_mother_mass;
        mass_onshellW_true = mu2_mother_mass;
    }

    
   pt_h2tohh_true = h2tohhcand->pt();
   float nu_onshellW_pt =0;
   wmass_gen = 80.3;// initial value
   float step,random01;
   TH1F* wmasshist = readoutonshellWMassPDF(); 
   TH1F* onshellnupt_hist = readoutonshellnuptPDF(); 
   onshellnupt_hist->Scale(1.0/onshellnupt_hist->GetBinContent(onshellnupt_hist->GetMaximumBin()));


   for (int i = 0; i < iterations_ ; i++){

	   eta_gen = generator->Uniform(-6,6); 
	   phi_gen = generator->Uniform(-3.1415926, 3.1415926);
           //wmass_gen = generator->Gaus(80.385,0.015);
           hmass_gen = generator->Gaus(125.03,0.3);
           //generate onshell Wmass
           step = generator->Uniform(-4,4);
           //step = generator->Gaus(0,8);
           random01 = generator->Uniform(0,1);
           //wmass_gen = onshellWMassRandomWalk(wmass_gen, step, random01);
           wmass_gen = onshellWMassRandomWalk(wmass_gen, step, random01, wmasshist);
	  
           //wmass_gen = wmasspdf->GetRandom(50.0,90.0);
           //test
           /*eta_gen = eta_nuonshellW_true;
           phi_gen = phi_nuonshellW_true;
	   wmass_gen = mass_onshellW_true; 
           std::cout << "true eta phi of nuonshell ("<<eta_nuonshellW_true <<","<<phi_nuonshellW_true<<"), pt " <<pt_nuonshellW_true 
		<<" mass of onshellW " << mass_onshellW_true <<" wmass_gen "<< wmass_gen  << std::endl;
           std::cout << "true eta phi of nuoffshell ("<<eta_nuoffshellW_true <<","<<phi_nuoffshellW_true<<"), pt " <<pt_nuoffshellW_true 
		<<" mass of offshellW " << mass_offshellW_true <<  std::endl;
            */
           int solutions = 0;//count num of soluble case
           bool solution[4]={false, false, false, false}; //record whether the case is soluble or not
           for (int j = 0; j < 4; j++){
                 assignMuLorentzVec(j/2);
          	 nu_onshellW_pt = nu1pt_onshellW(std::make_pair(eta_gen, phi_gen), mu_onshellW_lorentz, wmass_gen); 
          	 nu_onshellW_lorentz->SetPtEtaPhiM(nu_onshellW_pt, eta_gen, phi_gen,0);
                  // should replace htoBB_lorentz by jets_lorentz if we have correct jets_lorentz
                 //solution[j] = nulorentz_offshellW(jets_lorentz, mu_onshellW_lorentz,
                 solution[j] = nulorentz_offshellW(htoBB_lorentz, mu_onshellW_lorentz,
					            mu_offshellW_lorentz, nu_onshellW_lorentz,
						   nu_offshellW_lorentz, j%2, hmass_gen);
                 //std::cout << " calculate nu1_pt " << nu_onshellW_pt << " eta_gen "<< eta_gen << " phi_gen " << phi_gen << std::endl; 
                 //std::cout << j << " nu_offshellW_eta " << nu_offshellW_lorentz->Eta()<<" phi " << nu_offshellW_lorentz->Phi() << std::endl; 
                 if (solution[j]) solutions++;
           }
    //       nu_offshellW_lorentz= NULL; 
           for (int j = 0; j < 4; j++){
	 //  int fill = mmctree->Fill();
	        if (!solution[j])  continue;
                // reassign muons LorentzVector
                control = j;
                assignMuLorentzVec(j/2);
                nu_onshellW_pt = nu1pt_onshellW(std::make_pair(eta_gen, phi_gen), mu_onshellW_lorentz, wmass_gen); 
          	nu_onshellW_lorentz->SetPtEtaPhiM(nu_onshellW_pt, eta_gen, phi_gen,0);
                 //nulorentz_offshellW(jets_lorentz, mu_onshellW_lorentz,
                nulorentz_offshellW(htoBB_lorentz, mu_onshellW_lorentz,
                                    mu_offshellW_lorentz, nu_onshellW_lorentz,
				    nu_offshellW_lorentz, j%2, hmass_gen);


                weight = 1.0/solutions;// change weight if we consider possibility factor  like matrix elements
                if (weightfromonshellnupt_func_) weight1 = weightfromonshellnupt(nu_onshellW_pt); 
                if (weightfromonshellnupt_hist_) weight1 = weightfromonshellnupt(onshellnupt_hist, nu_onshellW_pt); 
                weight1 = weight1*weight;
 		mu_onshellW_Eta = mu_onshellW_lorentz->Eta();
   		mu_onshellW_Phi = mu_onshellW_lorentz->Phi();
   		mu_onshellW_Pt = mu_onshellW_lorentz->Pt();
   		mu_onshellW_E = mu_onshellW_lorentz->E();

   		mu_offshellW_Eta = mu_offshellW_lorentz->Eta();
   		mu_offshellW_Phi = mu_offshellW_lorentz->Phi();
   		mu_offshellW_Pt = mu_offshellW_lorentz->Pt();
  		mu_offshellW_E = mu_offshellW_lorentz->E();

           	nu_onshellW_Eta = nu_onshellW_lorentz->Eta();
           	nu_onshellW_Phi = nu_onshellW_lorentz->Phi();
           	nu_onshellW_Pt = nu_onshellW_lorentz->Pt();
           	nu_onshellW_E = nu_onshellW_lorentz->E();

           	nu_offshellW_Eta = nu_offshellW_lorentz->Eta();
          	nu_offshellW_Phi = nu_offshellW_lorentz->Phi();
          	nu_offshellW_Pt = nu_offshellW_lorentz->Pt();
           	nu_offshellW_E = nu_offshellW_lorentz->E();

                *onshellW_lorentz = *mu_onshellW_lorentz+*nu_onshellW_lorentz;
                *offshellW_lorentz = *mu_offshellW_lorentz+*nu_offshellW_lorentz;
                *htoWW_lorentz = *onshellW_lorentz+*offshellW_lorentz;
                *h2tohh_lorentz = *htoWW_lorentz+*htoBB_lorentz;
                *MMCmet_vec2 = TVector2(nu_onshellW_lorentz->Px()+nu_offshellW_lorentz->Px(),
						nu_onshellW_lorentz->Py()+nu_offshellW_lorentz->Py());
		if (fabs(hmass_gen-htoWW_lorentz->M()) > 2) {
			std::cout << "  hmass_gen " << hmass_gen << " Higgs mass from MMC " << htoWW_lorentz->M() <<std::endl;
           		verbose_ = 2;
                 }
                if (verbose_ > 1){
                	std::cout << "mu_onshellW "; mu_onshellW_lorentz->Print();
                	std::cout << "nu_onshellW "; nu_onshellW_lorentz->Print();
                	std::cout << "mu_offshellW "; mu_offshellW_lorentz->Print();
                	std::cout << "nu_offshellW "; nu_offshellW_lorentz->Print();
                }

                if (verbose_ > 1){
  			std::cout << " onshell W mass "<< onshellW_lorentz->M();   onshellW_lorentz->Print();
  			std::cout << " offshell W mass "<< offshellW_lorentz->M(); offshellW_lorentz->Print();
  			std::cout << " htoWW mass "<< htoWW_lorentz->M(); htoWW_lorentz->Print();
  			std::cout << " htoBB mass "<< htoBB_lorentz->M(); htoBB_lorentz->Print();
                        verbose_ = 0;
                }
  		if (verbose_ > 1 && (h2tohh_lorentz->Pt()/h2tohh_lorentz->E())>0.0000001) {
 			std::cout << " h2tohh mass "<< h2tohh_lorentz->M() <<" pt " << h2tohh_lorentz->Pt();
			h2tohh_lorentz->Print();
                }
           	onshellW_Eta = onshellW_lorentz->Eta();
           	onshellW_Phi = onshellW_lorentz->Phi();
           	onshellW_Pt = onshellW_lorentz->Pt();
           	onshellW_E = onshellW_lorentz->E();
           	onshellW_Mass = onshellW_lorentz->M();
           	offshellW_Eta = offshellW_lorentz->Eta();
           	offshellW_Phi = offshellW_lorentz->Phi();
           	offshellW_Pt = offshellW_lorentz->Pt();
           	offshellW_E = offshellW_lorentz->E();
           	offshellW_Mass = offshellW_lorentz->M();
                htoWW_Eta = htoWW_lorentz->Eta();
                htoWW_Phi = htoWW_lorentz->Phi();
                htoWW_Pt = htoWW_lorentz->Pt();
                htoWW_E = htoWW_lorentz->E();
                htoWW_Mass = htoWW_lorentz->M();
                h2tohh_Pt = h2tohh_lorentz->Pt();
                h2tohh_E = h2tohh_lorentz->E();
                h2tohh_Mass = h2tohh_lorentz->M();
                MMCmet_Px = MMCmet_vec2->Px();
                MMCmet_Py = MMCmet_vec2->Py();
                MMCmet_E = MMCmet_vec2->Mod();
                MMCmet_Phi = MMCmet_vec2->Phi();
 
                if ((h2tohh_lorentz->Pt()/h2tohh_lorentz->E())>0.0000001){
                	h2tohh_Eta = h2tohh_lorentz->Eta();
               		h2tohh_Phi = h2tohh_lorentz->Phi();
                }else {//pt =0, strange case here
                        h2tohh_Eta = 1000000;
                        h2tohh_Phi = 0;
                }


             	mmctree->Fill();
           }//end controls loop,(0,1,2,3)
              	//mmctree->Fill();
   }//end of tries
   delete generator;
  // delete wmasspdf;
//   delete mmctree;
}

//------------ method called to initialize a tree for MMC for this event ------------
void
DiHiggsWWAnalyzer::initTree(TTree* mmctree){
 
   std::cout <<" init tree "<< mmctree->GetTitle() << std::endl; 
   //initial branch value if necessary
   //
   //
   weight1 = 1.0;

 
   mmctree->Branch("ievent", &ievent);
   mmctree->Branch("eta_mean", &eta_mean);
   mmctree->Branch("eta_rms", &eta_rms);
   mmctree->Branch("eta_gen",&eta_gen);
   mmctree->Branch("phi_gen",&phi_gen);
   mmctree->Branch("wmass_gen",&wmass_gen);
   mmctree->Branch("hmass_gen",&hmass_gen);

   mmctree->Branch("mu_onshellW_eta", &mu_onshellW_Eta);
   mmctree->Branch("mu_onshellW_phi", &mu_onshellW_Phi);
   mmctree->Branch("mu_onshellW_pt", &mu_onshellW_Pt);
   mmctree->Branch("mu_onshellW_E", &mu_onshellW_E);
   mmctree->Branch("mu_offshellW_eta", &mu_offshellW_Eta);
   mmctree->Branch("mu_offshellW_phi", &mu_offshellW_Phi);
   mmctree->Branch("mu_offshellW_pt", &mu_offshellW_Pt);
   mmctree->Branch("mu_offshellW_E", &mu_offshellW_E);
   mmctree->Branch("nu_onshellW_eta", &nu_onshellW_Eta);
   mmctree->Branch("nu_onshellW_phi", &nu_onshellW_Phi);
   mmctree->Branch("nu_onshellW_pt", &nu_onshellW_Pt);
   mmctree->Branch("nu_onshellW_E", &nu_onshellW_E);
   mmctree->Branch("nu_offshellW_eta", &nu_offshellW_Eta);
   mmctree->Branch("nu_offshellW_phi", &nu_offshellW_Phi);
   mmctree->Branch("nu_offshellW_pt", &nu_offshellW_Pt);
   mmctree->Branch("nu_offshellW_E", &nu_offshellW_E);
   mmctree->Branch("onshellW_eta", &onshellW_Eta);
   mmctree->Branch("onshellW_phi", &onshellW_Phi);
   mmctree->Branch("onshellW_pt", &onshellW_Pt);
   mmctree->Branch("onshellW_E", &onshellW_E);
   mmctree->Branch("onshellW_Mass", &onshellW_Mass);
   mmctree->Branch("offshellW_eta", &offshellW_Eta);
   mmctree->Branch("offshellW_phi", &offshellW_Phi);
   mmctree->Branch("offshellW_pt", &offshellW_Pt);
   mmctree->Branch("offshellW_E", &offshellW_E);
   mmctree->Branch("offshellW_Mass", &offshellW_Mass);
   mmctree->Branch("htoWW_Eta", &htoWW_Eta);
   mmctree->Branch("htoWW_Phi", &htoWW_Phi);
   mmctree->Branch("htoWW_Pt", &htoWW_Pt);
   mmctree->Branch("htoWW_E", &htoWW_E);
   mmctree->Branch("htoWW_Mass", &htoWW_Mass);
   mmctree->Branch("htoBB_Eta", &htoBB_Eta);
   mmctree->Branch("htoBB_Phi", &htoBB_Phi);
   mmctree->Branch("htoBB_Pt", &htoBB_Pt);
   mmctree->Branch("htoBB_E", &htoBB_E);
   mmctree->Branch("htoBB_Mass", &htoBB_Mass);
   mmctree->Branch("MMCmet_E",&MMCmet_E);
   mmctree->Branch("MMCmet_Phi",&MMCmet_Phi);
   mmctree->Branch("MMCmet_Px",&MMCmet_Px);
   mmctree->Branch("MMCmet_Py",&MMCmet_Py);

   mmctree->Branch("h2tohh_Eta", &h2tohh_Eta);
   mmctree->Branch("h2tohh_Phi", &h2tohh_Phi);
   mmctree->Branch("h2tohh_Pt", &h2tohh_Pt);
   mmctree->Branch("h2tohh_E", &h2tohh_E);
   mmctree->Branch("h2tohh_Mass", &h2tohh_Mass);


   mmctree->Branch("met_true",&met);
   mmctree->Branch("met_phi_true",&met_phi);
   mmctree->Branch("met_px_true",&met_px);
   mmctree->Branch("met_py_true",&met_py);
 
   mmctree->Branch("eta_nuoffshellW_true", &eta_nuoffshellW_true);
   mmctree->Branch("phi_nuoffshellW_true", &phi_nuoffshellW_true);
   mmctree->Branch("eta_nuonshellW_true", &eta_nuonshellW_true);
   mmctree->Branch("phi_nuonshellW_true", &phi_nuonshellW_true);
   mmctree->Branch("pt_nuoffshellW_true", &pt_nuoffshellW_true);
   mmctree->Branch("pt_nuonshellW_true", &pt_nuonshellW_true);
   mmctree->Branch("mass_offshellW_true", &mass_offshellW_true);
   mmctree->Branch("mass_onshellW_true", &mass_onshellW_true);
   mmctree->Branch("mass_h2_true", &h2tohh_mass);
   mmctree->Branch("pt_h2_true", &pt_h2tohh_true);
   mmctree->Branch("mass_htoWW_true", &htoWW_mass);

   mmctree->Branch("weight", &weight);
   mmctree->Branch("weight1", &weight1);
   mmctree->Branch("control", &control);
   
   
}

TLorentzVector 
DiHiggsWWAnalyzer::calculateMET(){

   TLorentzVector METlorentz = TLorentzVector();
   TVector2 met_pxpy(nu1cand->px()+nu2cand->px(), nu1cand->py()+nu2cand->py());
   METlorentz.SetPxPyPzE(nu1cand->px()+nu2cand->px(), nu1cand->py()+nu2cand->py(),0,met_pxpy.Mod());

   return METlorentz;
}


//------------ method called to assign muons lorenz vector --------------
void 
DiHiggsWWAnalyzer::assignMuLorentzVec(int control){

 //  control/2 == 0, namely control =0,1, we have correct muon lorentz Vector pair
 //
  if (mu1_mother_mass > mu2_mother_mass && control == 0){
   	mu_onshellW_lorentz->SetPtEtaPhiM(mu1cand->pt(), mu1cand->eta(), mu1cand->phi(), 0);
        mu_offshellW_lorentz->SetPtEtaPhiM(mu2cand->pt(), mu2cand->eta(), mu2cand->phi(), 0); }
   else if (mu1_mother_mass > mu2_mother_mass && control == 1){
        mu_onshellW_lorentz->SetPtEtaPhiM(mu2cand->pt(), mu2cand->eta(), mu2cand->phi(), 0);
   	mu_offshellW_lorentz->SetPtEtaPhiM(mu1cand->pt(), mu1cand->eta(), mu1cand->phi(), 0);}
   else if (mu1_mother_mass < mu2_mother_mass && control == 0){
        mu_onshellW_lorentz->SetPtEtaPhiM(mu2cand->pt(), mu2cand->eta(), mu2cand->phi(), 0);
   	mu_offshellW_lorentz->SetPtEtaPhiM(mu1cand->pt(), mu1cand->eta(), mu1cand->phi(), 0);}
   else if (mu1_mother_mass < mu2_mother_mass && control == 1){
        mu_onshellW_lorentz->SetPtEtaPhiM(mu1cand->pt(), mu1cand->eta(), mu1cand->phi(), 0);
   	mu_offshellW_lorentz->SetPtEtaPhiM(mu2cand->pt(), mu2cand->eta(), mu2cand->phi(), 0);}
}


// ------------ method called to generate a pair (eta,phi) for nuetrino1  ------------
EtaPhi 
DiHiggsWWAnalyzer::generatenu1_etaphi(){

   float eta=0.0;
   float phi=0.0;
   
   float mean=0;
   float rms=1.403;
   eta = genEtaGuass(mean, rms);
   phi = genPhiFlat();

   return std::make_pair(eta, phi);
}

// ------------ method called to generate eta from Gauss distribution  ------------
float 
DiHiggsWWAnalyzer::genEtaGuass(float mean, float rms){
	
    TRandom3 *etaGenerator = new TRandom3();
    float eta = etaGenerator->Gaus(mean, rms);
    delete etaGenerator;
     
    return eta;

}

// ------------ method called to generate phi from Flat distribution  ------------
float 
DiHiggsWWAnalyzer::genPhiFlat(){
 
   TRandom3 *phiGenerator = new TRandom3();
   float pi = 3.1415926;
   float phi = phiGenerator->Uniform(-pi, pi);
   delete phiGenerator;

   return phi;
}

//------------ method called to readout TH1F onshellWmasspdf from root file -----------------------------
//
TH1F*
DiHiggsWWAnalyzer::readoutonshellWMassPDF(){

	
   //TFile* file = new TFile("/home/taohuang/work/CMSSW_7_3_1/src/DiHiggsWW/DiHiggsWWAnalyzer/plugins/MMCRefPDF.ROOT");
   TFile* file = new TFile(RefPDFfile_.c_str());
   TH1F* onshellWmasspdf = (TH1F*)file->Get("onshellWmasspdf");
   delete file;
   return onshellWmasspdf;

}


//------------ method called to readout TH1F onshellWmasspdf from root file -----------------------------
//
TH1F*
DiHiggsWWAnalyzer::readoutonshellnuptPDF(){

	
   //TFile* file = new TFile("/home/taohuang/work/CMSSW_7_3_1/src/DiHiggsWW/DiHiggsWWAnalyzer/plugins/MMCRefPDF.ROOT");
   TFile* file = new TFile(RefPDFfile_.c_str());
   TH1F* onshellWmasspdf = (TH1F*)file->Get("onshellnuptpdf");
   delete file;
   return onshellWmasspdf;

}



//------------ method to describe onshellW mass Probability density function ------------------------------
//
float 
DiHiggsWWAnalyzer::onshellWMassPDF(float mass){

  // float sigma = 1.75;
  // float mean = 80.1;
   float p0 =7.87161e-03;
   float p1 =1.69085;
   float p2 =603.474 ;
   float p = 0;
   p = exp(mass*p0+p1)+p2*exp(-0.5*((mass-80.1)/2.00)*((mass-80.1)/2.00));
   return p;
}

//------------ use random walk to generate random onshellW mass accroding to wmass pdf --------------
//
float
DiHiggsWWAnalyzer::onshellWMassRandomWalk(float x0, float step, float random){

   float xmin = 50;
   float xmax = 90;
   float x1 = x0+step;
   while (x1 > xmax || x1 < xmin){
  	if (x1 > xmax) x1 = x1-xmax+xmin;
  	if (x1 < xmin) x1 = xmax-(xmin-x1);
   }
    //transition probability
   float w = onshellWMassPDF(x1)/onshellWMassPDF(x0);
   //std::cout <<" initial " <<x0 <<" step " << step << " x1 "<< x1 << " transition probability " << w << " random " << random << std::endl;
   if (w >= 1.00) return x1;
   if (w < 1.00 && random < w) return x1;
   else return x0;

}  


//------------ use random walk to generate random onshellW mass accroding to wmass pdf --------------
//
float
DiHiggsWWAnalyzer::onshellWMassRandomWalk(float x0, float step, float random, TH1F* hist){
   float xmin = 50;
   float xmax = 90;
   //periodic boundary codition
   while (x0 > xmax || x0 < xmin){
        if (x0 > xmax) x0 = x0-xmax+xmin;
        if (x0 < xmin) x0 = xmax-(xmin-x0);
   }

   float x1 = x0+step;
   while (x1 > xmax || x1 < xmin){
  	if (x1 > xmax) x1 = x1-xmax+xmin;
  	if (x1 < xmin) x1 = xmax-(xmin-x1);
   }
    //find
   int binx0_1,binx0_2;
   int binx1_1,binx1_2;
   double bincent0_1,bincont0_1;// center and content
   double bincent1_1,bincont1_1;
   
   binx0_1 = hist->FindBin(x0);
   binx1_1 = hist->FindBin(x1);
  
   if ((float)hist->GetBinCenter(binx0_1) < x0){
	binx0_2 = binx0_1+1;
   }
   else {
	binx0_2 = binx0_1;
	binx0_1 = binx0_1-1;
    }

   if ((float)hist->GetBinCenter(binx1_1) < x1){
	binx1_2 = binx1_1+1;
    }
   else {
	binx1_2 = binx1_1;
	binx1_1 = binx1_1-1;
    }
    bincent0_1 = hist->GetBinCenter(binx0_1);
    bincont0_1 = hist->GetBinContent(binx0_1);
    bincent1_1 = hist->GetBinCenter(binx1_1);
    bincont1_1 = hist->GetBinContent(binx1_1);
   double w0 = (x0-bincent0_1)*(bincont0_1-hist->GetBinContent(binx0_2))/(bincent0_1-hist->GetBinCenter(binx0_2))+bincont0_1;
   double w1 = (x1-bincent1_1)*(bincont1_1-hist->GetBinContent(binx1_2))/(bincent1_1-hist->GetBinCenter(binx1_2))+bincont1_1;
    //transition probability
   double w = w1/w0;
   //std::cout <<" initial " <<x0 <<" step " << step << " x1 "<< x1 << " transition probability " << w << " random " << random << std::endl;
   if (w >= 1.00) return x1;
   if (w < 1.00 && random < (float)w) return x1;
   else return x0;

}  


//---------- weight solution by nupt --------------------------------------------------------
//
float
DiHiggsWWAnalyzer::weightfromonshellnupt(TH1F* hist, float nupt){
//hist should be scaled

   float weight = 0.0;
   int bin1 = hist->FindBin(nupt);
   //first make sure that nupt is within range
   if (bin1 == 0 || bin1 == hist->GetNbinsX()+1) return weight=0;
   
   float bin1content = hist->GetBinContent(bin1);
   float bin1center = hist->GetBinCenter(bin1);
   int bin2 = 0;
   
   if ((float)hist->GetBinCenter(bin1) < nupt)  
	bin2 = bin1+1;
   else
        bin2 = bin1-1;
   //find probability of nupt and set it as weight
   weight = (nupt-bin1center)*(bin1content-hist->GetBinContent(bin2))/(bin1center-hist->GetBinCenter(bin2))+bin1content;
   return weight;
}

//---------- weight solution by nupt --------------------------------------------------------
//
float
DiHiggsWWAnalyzer::weightfromonshellnupt(float nupt){

   float weight = 0.0;
   float max = 170;
   if (nupt<0 || nupt>125) return 0.0;

   weight = -16.925+12.4066*nupt-0.2884*std::pow(nupt,2)+0.00203*std::pow(nupt,3)+7.695e-7*std::pow(nupt,4)
            -7.2191e-8*std::pow(nupt,5)+2.499e-10*std::pow(nupt,6);
   if (weight < 0 && nupt<5) return 0.0;
   if (weight < 0) std::cout << " error! nupt " << nupt << " weight " << weight << std::endl;
   weight = weight/max;
   return weight;
}

//----------- method called to calculate total lorentz vector from b bar jets ---------------
//
TLorentzVector
DiHiggsWWAnalyzer::stableDecendantsLorentz(const reco::Candidate *cand){
  
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

//------------- method called to calculate pt of nuetrinos from on-shell W decay ------------
float 
DiHiggsWWAnalyzer::nu1pt_onshellW(EtaPhi nu1_etaphi, TLorentzVector* mu1lorentz, float wMass){
  
   float nu1_pt=0.0;
//   TVector2 *numu_phi = new TVector(nu_etaphi.first(),mu1lorentz->eta());
   float deltaeta = nu1_etaphi.first - mu1lorentz->Eta();
   float deltaphi = nu1_etaphi.second - mu1lorentz->Phi();
   
   nu1_pt = wMass*wMass/(2*mu1lorentz->Pt()*(cosh(deltaeta)-cos(deltaphi)));
   return nu1_pt;

}

//------------ method called to check whether the solution in this case exist or not -------------
// not use now, may be helpful later 
bool
DiHiggsWWAnalyzer::checkSolution(TLorentzVector* jetslorentz,
                                 TLorentzVector* mu1lorentz,
                                 TLorentzVector* mu2lorentz,
                                 TLorentzVector* nu1lorentz, int control, float hMass){

    

   TLorentzVector* tmplorentz = new TLorentzVector(mu1lorentz->Px()+mu2lorentz->Px()+nu1lorentz->Px(),
                                                   mu1lorentz->Py()+mu2lorentz->Py()+nu1lorentz->Py(),
                                                   mu1lorentz->Pz()+mu2lorentz->Pz()+nu1lorentz->Pz(),
                                                   mu1lorentz->Energy()+mu2lorentz->Energy()+nu1lorentz->Energy());

   float nu_tmp_px;
   float nu_tmp_py;
   float nu_tmp_pt;
   
   nu_tmp_px = -jetslorentz->Px()-mu1lorentz->Px()-mu2lorentz->Px()-nu1lorentz->Px();
   nu_tmp_py = -jetslorentz->Py()-mu1lorentz->Py()-mu2lorentz->Py()-nu1lorentz->Py();
   TVector2 nu_pxpy(nu_tmp_px, nu_tmp_py);

   nu_tmp_pt = nu_pxpy.Mod();

   float chdeltaeta;//cosh(nu2_eta-tmp2lorenz_eta)
   TLorentzVector* tmp2lorentz = new TLorentzVector(sqrt(pow(tmplorentz->Pt(),2)+pow(tmplorentz->M(),2)),0,tmplorentz->Pz(),tmplorentz->Energy());// construct massless lorentzvector with same pz and E as tmplorentzvector
   
   chdeltaeta = (pow(hMass,2)+pow(jetslorentz->Pt(),2)-pow(tmplorentz->M(),2)-pow(tmplorentz->Pt(),2)-pow(nu_tmp_pt,2))/(2*tmp2lorentz->Pt()*nu_tmp_pt);
   
   
   delete tmplorentz;
   delete tmp2lorentz;
   
// place the cuts we may need 
//
  //at present if (|chdeltaeta|>1) return true; 
   return (fabs(chdeltaeta)>1);
}




//------------- method called to calculate lorentzvector of second nuetrinos, which is from offshell W -----------
// return true if we can get nu_offshellW_lorentz
bool 
DiHiggsWWAnalyzer::nulorentz_offshellW(TLorentzVector* jetslorentz, 
                                        TLorentzVector* mu1lorentz, 
                                        TLorentzVector* mu2lorentz, 
                                        TLorentzVector* nu1lorentz, 
                                        TLorentzVector* nu2lorentz, int control, float hMass){

   TLorentzVector* tmplorentz = new TLorentzVector(mu1lorentz->Px()+mu2lorentz->Px()+nu1lorentz->Px(),
                                                   mu1lorentz->Py()+mu2lorentz->Py()+nu1lorentz->Py(),
                                                   mu1lorentz->Pz()+mu2lorentz->Pz()+nu1lorentz->Pz(),
                                                   mu1lorentz->Energy()+mu2lorentz->Energy()+nu1lorentz->Energy());
   float nu_tmp_px;
   float nu_tmp_py;
   float nu_tmp_pt;
   
   nu_tmp_px = -jetslorentz->Px()-mu1lorentz->Px()-mu2lorentz->Px()-nu1lorentz->Px();
   nu_tmp_py = -jetslorentz->Py()-mu1lorentz->Py()-mu2lorentz->Py()-nu1lorentz->Py();
   TVector2 nu_pxpy(nu_tmp_px, nu_tmp_py);

   nu_tmp_pt = nu_pxpy.Mod();

   float chdeltaeta;//cosh(nu_offshellW_eta-tmp2lorentz_eta)
   TLorentzVector* tmp2lorentz = new TLorentzVector(sqrt(pow(tmplorentz->Pt(),2)+pow(tmplorentz->M(),2)),0,tmplorentz->Pz(),tmplorentz->Energy());//fake one massless lorentzvector with same pz and E
   
   chdeltaeta = (pow(hMass,2)+pow(jetslorentz->Pt(),2)-pow(tmplorentz->M(),2)-pow(tmplorentz->Pt(),2)-pow(nu_tmp_pt,2))/(2*tmp2lorentz->Pt()*nu_tmp_pt);
   if (verbose_>0){
        std::cout << "nu2 px: " << nu_tmp_px << " py: "<< nu_tmp_py << std::endl;
   	std::cout << "chdeltaeta " << chdeltaeta << std::endl;
        std::cout << "tmp2lorentz "; tmp2lorentz->Print();
   }
   if (chdeltaeta < 1.0) {
        delete tmplorentz;
	delete tmp2lorentz;
   	nu2lorentz->SetPtEtaPhiM(0, 0, 0, 0);
	return false;
      }
   float nu_tmp_phi = nu_pxpy.Phi_mpi_pi(nu_pxpy.Phi());
   float deltaeta = acosh(chdeltaeta);
   float nu_tmp_eta = (control == 1) ? (tmp2lorentz->Eta()-deltaeta) : (tmp2lorentz->Eta()+deltaeta);//control = j%2 
   // should check whether deltaeta > 1
  // std::cout <<"control "<< control <<" nu_tmp_px " << nu_tmp_px << "  nu_tmp_py " << nu_tmp_py << " nu_tmp_pt " << nu_tmp_pt 
    //         << " cosh(deltaeta2) " << chdeltaeta << " nu_tmp_eta " << nu_tmp_eta << " nu_tmp_phi " << nu_tmp_phi << std::endl; 
   if (fabs(nu_tmp_eta) > 7) {
        delete tmplorentz;
	delete tmp2lorentz;
   	nu2lorentz->SetPtEtaPhiM(0, 0, 0, 0);
	return false;  //from simulation, |nu_offshellW_Eta|<6
    }
   nu2lorentz->SetPtEtaPhiM(nu_tmp_pt, nu_tmp_eta, nu_tmp_phi, 0);
   TLorentzVector* htoww_tmp = new TLorentzVector(*tmplorentz+*nu2lorentz);
   if (abs(htoww_tmp->M()-hMass) >2){
   	std::cout <<" set Higgs Mass" << hMass << " MMC higgs mass" << htoww_tmp->M() << std::endl;
        htoww_tmp->Print();
        verbose_ = 1;
       }
   if (verbose_ > 0){
   	std::cout << "tmplorentz mass " << tmplorentz->M(); tmplorentz->Print();
   	std::cout << "tmp2lorentz mass " << tmp2lorentz->M(); tmp2lorentz->Print();
   	std::cout << " jets lorentz"; jetslorentz->Print(); 
   	std::cout << " mu1 lorentz "; mu1lorentz->Print();
    	std::cout << " mu2 lorentz "; mu2lorentz->Print();
   	std::cout << " nu1 lorentz "; nu1lorentz->Print();
   	std::cout << " tmp lorentz "; tmplorentz->Print();
        std::cout << " nu2 lorentz "; nu2lorentz->Print();
    }
       // std::cout << " nu_offshellW lorentz "; nu2lorentz->Print();
   delete tmplorentz;
   delete tmp2lorentz;
   delete htoww_tmp;
   
   return true; 
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
