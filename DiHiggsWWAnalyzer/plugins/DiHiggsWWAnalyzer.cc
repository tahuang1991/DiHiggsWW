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

#define WMass 80.385   // W mass
#define SMHMass 125.03 // SM module higgs mass


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
      
      EtaPhi generatenu1_etaphi();
      float nu1pt_onshellW(EtaPhi nu1_etaphi, TLorentzVector* mu1lorentz);
      TLorentzVector* nu2lorentz_offshellW(TLorentzVector* jetlorentz, TLorentzVector* mu1lorentz, TLorentzVector* mu2lorentz, TLorentzVector*       nu1lorentz);
      
      bool cutsCheck();
       
          
 
      void print();
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      const reco::Candidate* mu1_W1_cand;
      const reco::Candidate* nu1_W1_cand;
      const reco::Candidate* mu2_W2_cand;
      const reco::Candidate* nu2_W2_cand;
      const reco::Candidate* mu1_htoWW_cand;
      const reco::Candidate* mu2_htoWW_cand;
      const reco::Candidate* b1_htobb_cand;
      const reco::Candidate* b2_htobb_cand;
      const reco::Candidate* h2tohh_cand;


      TTree *evtree;
      TFile *output;
      
      edm::Service< TFileService > fs;
      
      //----------branches of tree ---------------------------
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
   //now do what ever initialization is needed
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
	  mu1_energy = it->energy();
	  mu1_px = it->px();
	  mu1_py = it->py();
	  mu1_pz = it->pz();
          mu1_eta = it->eta();
          mu1_phi = it->phi();
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
	  mu2_energy = it->energy();
	  mu2_px = it->px();
	  mu2_py = it->py();
	  mu2_pz = it->pz();
          mu2_eta = it->eta();
          mu2_phi = it->phi();
	  mu2_motherid = it->mother()->pdgId();
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
      else if (it->pdgId() == -14 && !nu_negative )
      {
          nu1_eta = it->eta();
          nu1_phi = it->phi();
          const reco::Candidate* tmp_nu1 = it->mother();
          if (tmp_nu1->pdgId() == -24)  nu1_W1_cand = tmp_nu1;
    //      std::cout << " the mother of nutrio" 
          while (tmp_nu1->pdgId() == -24) tmp_nu1 = tmp_nu1->mother();
          if (tmp_nu1->pdgId() == 25)
             {
            //     std::cout << "find nuetrino candidate" << std::endl;
	         nu1_energy = it->energy();
	         nu1_px = it->px();
	         nu1_py = it->py();
	         nu1_pz = it->pz();
                 nu_negative = true;
                }
         
       }
      else if (it->pdgId() == 14 && !nu_positive )
      {
          const reco::Candidate* tmp_nu2 = it->mother(); 
          if (tmp_nu2->pdgId() == 24)  nu2_W2_cand = tmp_nu2;
          while (tmp_nu2->pdgId() == 24) tmp_nu2 = tmp_nu2->mother();
          if (tmp_nu2->pdgId() == 25)
             {
              //   std::cout << "find antinuetrino candidate" << std::endl;
	         nu2_energy = it->energy();
	         nu2_px = it->px();
	         nu2_py = it->py();
	         nu2_pz = it->pz();
                 nu_positive = true;
                }
         
       }
      else if (it->pdgId() == 5 && it->mother()->pdgId() == 25 && !bquark)
      {
          nu2_eta = it->eta();
          nu2_phi = it->phi();
	  b1_energy = it->energy();
	  b1_px = it->px();
	  b1_py = it->py();
	  b1_pz = it->pz();
	  b1_motherid = it->mother()->pdgId();
	  bquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
       //   std::cout << "find bquark candidate" << std::endl;
          b1_htobb_cand = it->mother();
      }
      else if (it->pdgId() == -5 && it->mother()->pdgId() == 25 && !bbarquark)
      {
	  b2_energy = it->energy();
	  b2_px = it->px();
	  b2_py = it->py();
	  b2_pz = it->pz();
	  b2_motherid = it->mother()->pdgId();
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
	  mu1_motherid = mu1_W1_cand->pdgId();
          mu1_mother_px = mu1_W1_cand->px();
          mu1_mother_py = mu1_W1_cand->py();
          mu1_mother_pz = mu1_W1_cand->pz();
      //    TLorentzVector W1(mu1_mother_px,mu1_mother_py,mu1_mother_pz,mu1_mother_energy);
          mu1_mother_energy = mu1_W1_cand->energy();
          mu1_mother_mass = mu1_W1_cand->mass();
          std::cout << "find W1 from negative mu nu, mass " << mu1_W1_cand->mass() 
                    << " energy " << mu1_W1_cand->energy() << std::endl;
    
          }
     else if (mu_negative && nu_negative)  std::cout << "find negative mu && nu but not find W1" << std::endl;
  
    if (mu_positive and nu_positive && mu2_W2_cand == nu2_W2_cand)
       {
          Wtomu2nu2 = true;
	  mu2_motherid = mu2_W2_cand->pdgId();
          mu2_mother_px = mu2_W2_cand->px();
          mu2_mother_py = mu2_W2_cand->py();
          mu2_mother_pz = mu2_W2_cand->pz();
          mu2_mother_energy = mu2_W2_cand->energy();
          mu2_mother_mass = mu2_W2_cand->mass();
          std::cout << "find W2 from positive mu nu, mass " << mu2_W2_cand->mass()
                    << " energy " << mu2_W2_cand->energy() << std::endl;
    
          }
     else if (mu_positive && nu_positive)  std::cout << "find positive mu && nu but not find W2" << std::endl;
  
    if (Wtomu1nu1 and Wtomu2nu2 and mu1_htoWW_cand == mu2_htoWW_cand)
       {
         std::cout << "find 2 muons and they come from same higgs" << std::endl;
          htoWW_energy = mu1_htoWW_cand->energy();
          htoWW_px = mu1_htoWW_cand->px();
          htoWW_py = mu1_htoWW_cand->py();
          htoWW_pz = mu1_htoWW_cand->pz();
          htoWW_mass = mu1_htoWW_cand->mass();
          float h_energy = mu1_energy+mu2_energy+nu1_energy+nu2_energy;
          float h_px = mu1_px+mu2_px+nu1_px+nu2_px;
          float h_py = mu1_py+mu2_py+nu1_py+nu2_py;
          float h_pz = mu1_pz+mu2_pz+nu1_pz+nu2_pz;
          float h_mass = std::sqrt(h_energy*h_energy-h_px*h_px-h_py*h_py-h_pz*h_pz);
          if (abs(h_mass-125)>5) {
               std::cout << "    mass of higgs " << mu1_htoWW_cand->mass() << "from final state " << h_mass << std::endl;
               std::cout << "energy from higgs " << htoWW_energy << " from final state " << h_energy << std::endl;
               std::cout << "    px from higgs " << htoWW_px << " from final state " << h_px << std::endl;
               std::cout << "    py from higgs " << htoWW_py << " from final state " << h_py << std::endl;
               std::cout << "    pz from higgs " << htoWW_pz << " from final state " << h_pz << std::endl;
             }
          htoWW = true;
         }
    else if(Wtomu1nu1 and Wtomu2nu2)   
       {
         std::cout << "find 2 muons but they are not from same higgs" << std::endl;
         std::cout << "mu1_higgs energy " << mu1_htoWW_cand->energy() << " px " << mu1_htoWW_cand->px() << std::endl;
         std::cout << "mu2_higgs energy " << mu2_htoWW_cand->energy() << " px " << mu2_htoWW_cand->px() << std::endl;

         }

    if (bquark and bbarquark and b1_htobb_cand == b2_htobb_cand)
       {
         std::cout << "find bbar and they come from same higgs" << std::endl;
         // const reco::Candidate* tmphiggs_bb = b1_htobb_cand->mother();
         // while (b1_htobb_cand->mother()->pdgId() == 25)  b1_htobb_cand = b1_htobb_cand->mother();
          
          htobb_energy = b1_htobb_cand->energy();
          htobb_px = b1_htobb_cand->px();
          htobb_py = b1_htobb_cand->py();
          htobb_pz = b1_htobb_cand->pz();
          htobb_mass = b1_htobb_cand->mass();
          htobb = true;
         }
//     else if()
          const reco::Candidate* tmp_htoWW = NULL;
          const reco::Candidate* tmp_htoBB = NULL;
    if (htoWW and htobb and mu1_htoWW_cand != b1_htobb_cand)
       {
          
          tmp_htoWW = mu1_htoWW_cand->mother();
          tmp_htoBB =  b1_htobb_cand->mother();
          while (tmp_htoWW->pdgId() == 25)  tmp_htoWW = tmp_htoWW->mother();
          while (tmp_htoBB->pdgId() == 25)  tmp_htoBB = tmp_htoBB->mother();
          if (tmp_htoWW == tmp_htoBB){

              std::cout << "find 2 higgs and both of them come from same heavey higgs"  << std::endl;
              h2tohh_energy = tmp_htoWW->energy();
              h2tohh_px = tmp_htoWW->px();
              h2tohh_py = tmp_htoWW->py();
              h2tohh_pz = tmp_htoWW->pz();
              h2tohh_mass = tmp_htoWW->mass();
              h2tohh = true;
              h2tohh_cand = tmp_htoWW;
              print();
            }
          else {
              std::cout << "pdgId of mother htoWW " << tmp_htoWW->pdgId() << " htoBB " << tmp_htoBB->pdgId() << std::endl;  
                 }
         }
    else if(htoWW and htobb)   
       {
         std::cout << "find 2 higgs but they are not from same heavey higgs" << std::endl;
         std::cout << "mother of htoWW id " << tmp_htoWW->pdgId() <<" energy " << tmp_htoWW->energy() << " px " << tmp_htoWW->px() << std::endl;
         std::cout << "mother of htobb id " << tmp_htoBB->pdgId() << " energy " << tmp_htoBB->energy() << " px " << tmp_htoBB->px() << std::endl;

         }
    
   //TLorentzVector htobb;// = new TLorentzVector(0,0,0,0);
   //htobb.SetPxPyPzE(b1_px+b2_px, b1_py+b2_py, b1_pz+b2_pz, b1_energy+b2_energy);
     //  htobb_mass = htobb.M();


   //if (mu_positive and mu_negative and bquark and bbarquark and mu1_htoWW_cand != b1_htobb_cand) 
   if (htoWW or htobb)
   {
     //  std::cout << "find one event with required final state(mumubb)" << std::endl;
       evtree->Fill();
       
   }
}


// ------------ method called once each job just before starting event loop  ------------
void 
DiHiggsWWAnalyzer::beginJob()
{
   evtree = fs->make<TTree>("evtree", "evtree");
 //  output = new TFile("output.root","recreate");
  // output->cd();

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
   
   evtree->Branch("htobb",&htobb);
   evtree->Branch("htoWW",&htoWW);
   evtree->Branch("h2tohh",&h2tohh);
    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiHiggsWWAnalyzer::endJob() 
{
    //output->Write();
   // output->Close();
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

// ------------ method called for printing additional information, useful for debugging  ------------
void
DiHiggsWWAnalyzer::print() {
    
    std::cout << "print() " << std::endl;
    if (!h2tohh)  {
       
        std::cout << "print() error " << std::endl;
        return;
    }
    
    std::cout << "  particle	, id	, mass	, px	, py	, pz	, Energy	" << std::endl;
    std::cout << "  h2	"<< h2tohh_cand->pdgId() <<"	, "<< h2tohh_cand->mass() << "	, " << h2tohh_cand->px() 
              << "	, "<< h2tohh_cand->py() << "	, " << h2tohh_cand->pz() << "	, " << h2tohh_cand->energy()
              << std::endl;
    


}


// ------------ method called to generate a pair (eta,phi) for nuetrino1  ------------
EtaPhi 
DiHiggsWWAnalyzer::generatenu1_etaphi(){

   float eta=0.0;
   float phi=0.0;

   return std::make_pair(eta, phi);
}

//------------- method called to calculate pt of nuetrinos from on-shell W decay ------------
float 
DiHiggsWWAnalyzer::nu1pt_onshellW(EtaPhi nu1_etaphi, TLorentzVector* mu1lorentz){
  
   float nu1_pt=0.0;
//   TVector2 *numu_phi = new TVector(nu_etaphi.first(),mu1lorentz->eta());
   float deltaeta = nu1_etaphi.first - mu1lorentz->Eta();
   float deltaphi = nu1_etaphi.second - mu1lorentz->Phi();
   nu1_pt = WMass*WMass/(2*mu1lorentz->Pt()*(cosh(deltaeta)-cos(deltaphi)));
 
   return nu1_pt;

}

//------------- method called to calculate lorentzvector of second nuetrinos, which is from offshell W -----------
TLorentzVector* 
DiHiggsWWAnalyzer::nu2lorentz_offshellW(TLorentzVector* jetslorentz, 
                                        TLorentzVector* mu1lorentz, 
                                        TLorentzVector* mu2lorentz, 
                                        TLorentzVector* nu1lorentz){

   TLorentzVector* nu2lorentz = new TLorentzVector();
   TLorentzVector* tmplorentz = new TLorentzVector(mu1lorentz->Px()+mu2lorentz->Px()+nu1lorentz->Px(),
                                                   mu1lorentz->Py()+mu2lorentz->Py()+nu1lorentz->Py(),
                                                   mu1lorentz->Pz()+mu2lorentz->Pz()+nu1lorentz->Pz(),
                                                   mu1lorentz->Energy()+mu2lorentz->Energy()+nu1lorentz->Energy());
   float nu2_px;
   float nu2_py;
   float nu2_pt;
   
   nu2_px = -jetslorentz->Px()-mu1lorentz->Px()-mu2lorentz->Px()-nu1lorentz->Px();
   nu2_py = -jetslorentz->Py()-mu1lorentz->Py()-mu2lorentz->Py()-nu1lorentz->Py();
   nu2_pt = sqrt(nu2_px*nu2_px+nu2_py*nu2_py);

   float chdeltaeta;//cosh(nu2_eta-mu2_eta)
   
   chdeltaeta = (pow(SMHMass,2)+pow(jetslorentz->Pt(),2)-pow(tmplorentz->Pt(),2)-pow(nu2_pt,2))/(2*tmplorentz->Pt()*nu2_pt);

   // should check whether deltaeta > 1
   float nu2_eta = acosh(chdeltaeta)+tmplorentz->Eta();
   float nu2_phi = atan(nu2_py/nu2_px);
   nu2lorentz->SetPtEtaPhiM(nu2_pt, nu2_eta, nu2_phi, 0);

   return nu2lorentz;
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
