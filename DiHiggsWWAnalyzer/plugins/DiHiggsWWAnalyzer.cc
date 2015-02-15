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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//headers from root lib
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

//
// class declaration
//
using namespace reco;

class DiHiggsWWAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DiHiggsWWAnalyzer(const edm::ParameterSet&);
      ~DiHiggsWWAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      TTree *evtree;
      TFile *output;

      //branches of tree
      float mu1_energy;
      float mu1_px;
      float mu1_py;
      float mu1_pz;
      int mu1_motherid;
      float mu1_mother_energy;
      float mu1_mother_px;
      float mu1_mother_py;
      float mu1_mother_pz;
      float nu1_energy;
      float nu1_px;
      float nu1_py;
      float nu1_pz;
      bool Wtomu1nu1;

      float mu2_energy;
      float mu2_px;
      float mu2_py;
      float mu2_pz;
      int mu2_motherid;
      float mu2_mother_energy;
      float mu2_mother_px;
      float mu2_mother_py;
      float mu2_mother_pz;
      float nu2_energy;
      float nu2_px;
      float nu2_py;
      float nu2_pz;
      bool Wtomu2nu2;
     
      float htoWW_energy;
      float htoWW_px;
      float htoWW_py;
      float htoWW_pz; 
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
      
      float h2tohh_energy;
      float h2tohh_px;
      float h2tohh_py;
      float h2tohh_pz;
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
      mu1_energy = 0.0;
      mu1_px = 0.0;
      mu1_py = 0.0;
      mu1_pz = 0.0;
      mu1_motherid = 0;
      mu1_mother_energy = 0.0;
      mu1_mother_px = 0.0;
      mu1_mother_py = 0.0;
      mu1_mother_pz = 0.0;
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
      nu2_energy = 0.0;
      nu2_px = 0.0;
      nu2_py = 0.0;
      nu2_pz = 0.0;
      Wtomu2nu2 = false;

      htoWW_energy = 0.0;
      htoWW_px = 0.0;
      htoWW_py = 0.0;
      htoWW_pz = 0.0;

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
      
      h2tohh_energy = 0.0;
      h2tohh_px = 0.0;
      h2tohh_py = 0.0;
      h2tohh_pz = 0.0;

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

    const reco::Candidate* mu1_htoWW_cand=NULL;
    const reco::Candidate* mu2_htoWW_cand=NULL;
    const reco::Candidate* b1_htobb_cand=NULL;
    const reco::Candidate* b2_htobb_cand=NULL;
    std::vector<reco::GenParticle*> mu1Coll; 
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
	  //mu_negative = true;
          //std::cout << "find muon(-) with status 1" << std::endl;
          const reco::Candidate* tmp_mu1 = it->mother(); 
          while (tmp_mu1->pdgId() == 13 && tmp_mu1->numberOfMothers() == 1) tmp_mu1 = tmp_mu1->mother();
          if (tmp_mu1->numberOfMothers() != 1 ) std::cout << "muon has more than one mother particle" << std::endl;
 
	  mu1_motherid = tmp_mu1->pdgId();
          mu1_mother_px = tmp_mu1->px();
          mu1_mother_py = tmp_mu1->py();
          mu1_mother_pz = tmp_mu1->pz();
          mu1_mother_energy = tmp_mu1->energy();

          while (tmp_mu1->pdgId() == -24 && tmp_mu1->numberOfMothers() == 1) tmp_mu1 = tmp_mu1->mother();
          if (tmp_mu1->numberOfMothers() != 1 ) std::cout << "W- has more than one mother particle" << std::endl;
         
          if (tmp_mu1->pdgId() == 25)   
             {
                  while (tmp_mu1->mother()->pdgId() == 25) tmp_mu1 = tmp_mu1->mother();
                  std::cout << "find muon(-) candidate" << std::endl;
	          mu_negative = true;
                  mu1_htoWW_cand = tmp_mu1;
                  mu1Coll.push_back(it->clone());
                  std::cout << "mother of this higgs, id " << tmp_mu1->mother()->pdgId() << " energy " << tmp_mu1->mother()->energy() << std::endl;
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
	  mu2_motherid = it->mother()->pdgId();
         // std::cout << "find muon(+) with status 1" << std::endl;
	//  mu_positive = true;
          const reco::Candidate* tmp_mu2 = it->mother(); 
          while (tmp_mu2->pdgId() == -13 && tmp_mu2->numberOfMothers() == 1) tmp_mu2 = tmp_mu2->mother();
          if (tmp_mu2->numberOfMothers() != 1)  std::cout << "muon has more than one mother particle" << std::endl;
	  mu2_motherid = tmp_mu2->pdgId();
          mu2_mother_px = tmp_mu2->px();
          mu2_mother_py = tmp_mu2->py();
          mu2_mother_pz = tmp_mu2->pz();
          mu2_mother_energy = tmp_mu2->energy();
          while (tmp_mu2->pdgId() == 24 && tmp_mu2->numberOfMothers() == 1) tmp_mu2 = tmp_mu2->mother();
          if (tmp_mu2->numberOfMothers() != 1 ) std::cout << "W+ has more than one mother particle" << std::endl;
    //      const reco::Candidate* tmphiggs_mu2 = tmp_mu2->mother();
          if (tmp_mu2->pdgId() == 25)   
             {
                 while (tmp_mu2->mother()->pdgId() == 25)   tmp_mu2 = tmp_mu2->mother();
                 std::cout << "find muon(+) candidate" << std::endl;
	         mu_positive = true;
                 mu2_htoWW_cand = tmp_mu2;
                 std::cout << "mother of this higgs, id " << tmp_mu2->mother()->pdgId() << " energy " << tmp_mu2->mother()->energy() << std::endl;
               }
        }
      else if (it->pdgId() == -14 && !nu_negative )
      {
          const reco::Candidate* tmp_nu1 = it->mother();
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

    if (mu_negative and mu_positive and mu1_htoWW_cand == mu2_htoWW_cand)
       {
         std::cout << "find 2 muons and they come from same higgs" << std::endl;
          htoWW_energy = mu1_htoWW_cand->energy();
          htoWW_px = mu1_htoWW_cand->px();
          htoWW_py = mu1_htoWW_cand->py();
          htoWW_pz = mu1_htoWW_cand->pz();
          htoWW = true;
         }
    else if(mu_negative and mu_positive)   
       {
         std::cout << "find 2 muons but they are not from same higgs" << std::endl;
         std::cout << "mu1_higgs energy " << mu1_htoWW_cand->energy() << " px " << mu1_htoWW_cand->px() << std::endl;
         std::cout << "mu2_higgs energy " << mu2_htoWW_cand->energy() << " px " << mu2_htoWW_cand->px() << std::endl;

         }
    if (bquark and bbarquark and b1_htobb_cand == b2_htobb_cand)
       {
         std::cout << "find bbar and they come from same higgs" << std::endl;
         // const reco::Candidate* tmphiggs_bb = b1_htobb_cand->mother();
          while (b1_htobb_cand->mother()->pdgId() == 25)  b1_htobb_cand = b1_htobb_cand->mother();
          
          htobb_energy = b1_htobb_cand->energy();
          htobb_px = b1_htobb_cand->px();
          htobb_py = b1_htobb_cand->py();
          htobb_pz = b1_htobb_cand->pz();
          htobb = true;
         }
//     else if()
    if (htoWW and htobb and mu1_htoWW_cand != b1_htobb_cand and mu1_htoWW_cand->mother() == b1_htobb_cand->mother())
       {
         std::cout << "find 2 higgs and both of them come from same heavey higgs"  << std::endl;
          h2tohh_energy = mu1_htoWW_cand->mother()->energy();
          h2tohh_px = mu1_htoWW_cand->mother()->px();
          h2tohh_py = mu1_htoWW_cand->mother()->py();
          h2tohh_pz = mu1_htoWW_cand->mother()->pz();
          h2tohh = true;
         }
    else if(htoWW and htobb)   
       {
         std::cout << "find 2 higgs but they are not from same heavey higgs" << std::endl;
         std::cout << "mother of htoWW id " << mu1_htoWW_cand->mother()->pdgId() <<" energy " << mu1_htoWW_cand->mother()->energy() << " px " << mu1_htoWW_cand->mother()->px() << std::endl;
         std::cout << "mother of htobb id " << b1_htobb_cand->mother()->pdgId() << " energy " << b1_htobb_cand->mother()->energy() << " px " << b1_htobb_cand->mother()->px() << std::endl;

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
   output = new TFile("output.root","recreate");
   output->cd();

   evtree = new TTree("evtree","event tree");
   evtree->Branch("mu1_energy",&mu1_energy);
   evtree->Branch("mu1_px",&mu1_px);
   evtree->Branch("mu1_py",&mu1_py);
   evtree->Branch("mu1_pz",&mu1_pz);
   evtree->Branch("mu1_motherid",&mu1_motherid);
   evtree->Branch("mu1_mother_energy",&mu1_mother_energy);
   evtree->Branch("mu1_mother_px",&mu1_mother_px);
   evtree->Branch("mu1_mother_py",&mu1_mother_py);
   evtree->Branch("mu1_mother_pz",&mu1_mother_pz);
   evtree->Branch("nu1_energy",&nu1_energy);
   evtree->Branch("nu1_px",&nu1_px);
   evtree->Branch("nu1_py",&nu1_py);
   evtree->Branch("nu1_pz",&nu1_pz);
   evtree->Branch("Wtomu1nu1",&Wtomu1nu1);

   evtree->Branch("mu2_energy",&mu2_energy);
   evtree->Branch("mu2_px",&mu2_px);
   evtree->Branch("mu2_py",&mu2_py);
   evtree->Branch("mu2_pz",&mu2_pz);
   evtree->Branch("mu2_mother_energy",&mu2_mother_energy);
   evtree->Branch("mu2_mother_px",&mu2_mother_px);
   evtree->Branch("mu2_mother_py",&mu2_mother_py);
   evtree->Branch("mu2_mother_pz",&mu2_mother_pz);
   evtree->Branch("mu2_motherid",&mu2_motherid);
   evtree->Branch("nu2_energy",&nu2_energy);
   evtree->Branch("nu2_px",&nu2_px);
   evtree->Branch("nu2_py",&nu2_py);
   evtree->Branch("nu2_pz",&nu2_pz);
   evtree->Branch("Wtomu2nu2",&Wtomu2nu2);

   evtree->Branch("htoWW_energy",&htoWW_energy);
   evtree->Branch("htoWW_px",&htoWW_px);
   evtree->Branch("htoWW_py",&htoWW_py);
   evtree->Branch("htoWW_pz",&htoWW_pz);
   
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
   
   evtree->Branch("h2tohh_energy",&h2tohh_energy);
   evtree->Branch("h2tohh_px",&h2tohh_px);
   evtree->Branch("h2tohh_py",&h2tohh_py);
   evtree->Branch("h2tohh_pz",&h2tohh_pz);
   
   evtree->Branch("htobb",&htobb);
   evtree->Branch("htoWW",&htoWW);
   evtree->Branch("h2tohh",&h2tohh);
    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiHiggsWWAnalyzer::endJob() 
{
    output->Write();
    output->Close();
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
