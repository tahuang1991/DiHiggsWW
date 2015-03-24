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

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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
    private:
      std::vector<reco::GenParticle*> mu1Coll; 
      std::vector<reco::GenParticle*> mu2Coll; 
      std::vector<reco::GenParticle*> nu1Coll; 
      std::vector<reco::GenParticle*> nu2Coll; 
      std::vector<reco::GenParticle*> b1Coll; 
      std::vector<reco::GenParticle*> b2Coll;
      std::vector<const reco::Candidate*> W1Coll;
      std::vector<const reco::Candidate*> W2Coll;
      std::vector<const reco::Candidate*> htoWWColl;
      std::vector<const reco::Candidate*> htoBBColl;
      const reco::Candidate* mu1cand;//store correct mu1 in h2->hh->bbWW->bb mu1 nu1mu2nu2
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
      void clear();
    private:
     //decendants
     const reco::Candidate* stabledecendant(const reco::Candidate* cand, int id);
     //const reco::Candidate* stablehtoWWdecendant(Particle p, PdgId id);
     const reco::Candidate* alldecendant(const reco::Candidate* cand, int id, bool first=false);
     const reco::Candidate* allancestor(const reco::Candidate* cand, int id, bool first=false);
                    
    private:
      void printCandidate(const reco::Candidate* );
      

    private:
      // ----------member data ---------------------------
      TTree *evtree;
      TFile *output;
    
    private:
      void fillbranches(); 
    private: 
      edm::Service< TFileService > fs;
      //branches of tree
      float mu1_energy;
      float mu1_px;
      float mu1_py;
      float mu1_pz;
      float mu1_eta;
      float mu1_phi;
      float mu1_mother_mass;
      float mu1_mother_energy;
      float mu1_mother_px;
      float mu1_mother_py;
      float mu1_mother_pz;
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
      float mu2_mother_mass;
      float mu2_mother_energy;
      float mu2_mother_px;
      float mu2_mother_py;
      float mu2_mother_pz;
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
      mu1_energy = 0.0;
      mu1_px = 0.0;
      mu1_py = 0.0;
      mu1_pz = 0.0;
      mu1_eta = 0.0;
      mu1_phi = 0.0;
      mu1_mother_mass = 0;
      mu1_mother_energy = 0.0;
      mu1_mother_px = 0.0;
      mu1_mother_py = 0.0;
      mu1_mother_pz = 0.0;
      nu1_energy = 0.0;
      nu1_px = 0.0;
      nu1_py = 0.0;
      nu1_pz = 0.0;
      nu1_eta = 0.0;
      nu1_phi = 0.0;
      Wtomu1nu1 = false;

      mu2_energy = 0.0;
      mu2_px = 0.0;
      mu2_py = 0.0;
      mu2_pz = 0.0;
      mu2_eta = 0.0;
      mu2_phi = 0.0;
      mu2_mother_mass = 0;
      mu2_mother_energy = 0.0;
      mu2_mother_px = 0.0;
      mu2_mother_py = 0.0;
      mu2_mother_pz = 0.0;
      nu2_energy = 0.0;
      nu2_px = 0.0;
      nu2_py = 0.0;
      nu2_pz = 0.0;
      nu2_eta = 0.0;
      nu2_phi = 0.0;
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

    clear();


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
 

          while (tmp_mu1->pdgId() == -24 && tmp_mu1->numberOfMothers() == 1) tmp_mu1 = tmp_mu1->mother();
          if (tmp_mu1->numberOfMothers() != 1 ) std::cout << "W- has more than one mother particle" << std::endl;
         
          if (tmp_mu1->pdgId() == 25)   
             {
                  while (tmp_mu1->mother()->pdgId() == 25) tmp_mu1 = tmp_mu1->mother();
                  std::cout << "find muon(-) candidate" << std::endl;
	          mu_negative = true;
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
         // std::cout << "find muon(+) with status 1" << std::endl;
	//  mu_positive = true;
          const reco::Candidate* tmp_mu2 = it->mother(); 
          while (tmp_mu2->pdgId() == -13 && tmp_mu2->numberOfMothers() == 1) tmp_mu2 = tmp_mu2->mother();
          if (tmp_mu2->numberOfMothers() != 1)  std::cout << "muon has more than one mother particle" << std::endl;
          while (tmp_mu2->pdgId() == 24 && tmp_mu2->numberOfMothers() == 1) tmp_mu2 = tmp_mu2->mother();
          if (tmp_mu2->numberOfMothers() != 1 ) std::cout << "W+ has more than one mother particle" << std::endl;
    //      const reco::Candidate* tmphiggs_mu2 = tmp_mu2->mother();
          if (tmp_mu2->pdgId() == 25)   
             {
                 while (tmp_mu2->mother()->pdgId() == 25)   tmp_mu2 = tmp_mu2->mother();
                 std::cout << "find muon(+) candidate" << std::endl;
	         mu_positive = true;
                 mu2Coll.push_back(it->clone());
                 std::cout << "mother of this higgs, id " << tmp_mu2->mother()->pdgId() << " energy " << tmp_mu2->mother()->energy() << std::endl;
               }
        }
      else if (it->pdgId() == -14 && !nu_negative )
      {
          const reco::Candidate* tmp_nu1 = it->mother();
          std::cout << " the mother of nuetrino1 id " <<  tmp_nu1->pdgId() << std::endl; 
          while (tmp_nu1->pdgId() == -24) tmp_nu1 = tmp_nu1->mother();
          if (tmp_nu1->pdgId() == 25)
             {
            //     std::cout << "find nuetrino candidate" << std::endl;
                 nu_negative = true;
                 nu1Coll.push_back(it->clone());
                }
         
       }
      else if (it->pdgId() == 14 && !nu_positive )
      {
          const reco::Candidate* tmp_nu2 = it->mother(); 
          std::cout << " the mother of nuetrino2 id " <<  tmp_nu2->pdgId() << std::endl; 
          while (tmp_nu2->pdgId() == 24) tmp_nu2 = tmp_nu2->mother();
          if (tmp_nu2->pdgId() == 25)
             {
              //   std::cout << "find antinuetrino candidate" << std::endl;
                 nu_positive = true;
                 nu2Coll.push_back(it->clone());
                }
         
       }
      else if (it->pdgId() == 5 && it->mother()->pdgId() == 25 && !bquark)
      {
	  bquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
       //   std::cout << "find bquark candidate" << std::endl;
          b1Coll.push_back(it->clone());
      }
      else if (it->pdgId() == -5 && it->mother()->pdgId() == 25 && !bbarquark)
      {
	  bbarquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
         // std::cout << "find bbarquark candidate" << std::endl;
          b2Coll.push_back(it->clone());
      }

     // std::cout << "test" << std::endl;

   }// all Gen particles

    //W1
    if (mu1Coll.size() && nu1Coll.size()){
         for (auto mu1_cand : mu1Coll)
              for (auto nu1_cand : nu1Coll){
                      const reco::Candidate* mu1_mother = mu1_cand->mother();
                      const reco::Candidate* nu1_mother = nu1_cand->mother();
                      while (mu1_mother->pdgId() == 13) mu1_mother = mu1_mother->mother();
                      if (nu1_mother == mu1_mother && mu1_mother->pdgId() == -24){ 
				W1Coll.push_back(mu1_mother);
				break;
                      }        
              }
    }
    //W2
    if (mu2Coll.size() && nu2Coll.size()){
         for (auto mu2_cand : mu2Coll)
              for (auto nu2_cand : nu2Coll){
                      const reco::Candidate* mu2_mother = mu2_cand->mother();
                      const reco::Candidate* nu2_mother = nu2_cand->mother();
                      while (mu2_mother->pdgId() == -13) mu2_mother = mu2_mother->mother();
                      if (nu2_mother == mu2_mother && mu2_mother->pdgId() == 24) {
				W2Coll.push_back(mu2_mother);
				break;
		        }
                      
              }
    }

    //htoWW
    if (W1Coll.size() && W2Coll.size()){
         for (auto W1_cand : W1Coll)
              for (auto W2_cand : W2Coll){
                      const reco::Candidate* W1_mother = W1_cand->mother();
                      const reco::Candidate* W2_mother = W2_cand->mother();
                      while (W1_mother->pdgId() == -24) W1_mother = W1_mother->mother();
                      while (W2_mother->pdgId() ==  24) W2_mother = W2_mother->mother();
                      if (W1_mother == W2_mother && W1_mother->pdgId() == 25) {
				htoWWColl.push_back(W1_mother);
				break;
			}
              }
    }

     //htoBB
     if (b1Coll.size() && b2Coll.size()){
          for(auto b1_cand : b1Coll)
              for (auto b2_cand : b2Coll) {
                       const reco::Candidate* b1_mother = b1_cand->mother();
                       const reco::Candidate* b2_mother = b2_cand->mother();
                       if (b1_mother == b2_mother && b1_mother->pdgId() == 25) {
				htoBBColl.push_back(b1_mother);
				break;
			}

              }
       
     }

      //h2tohh
     if (htoWWColl.size() && htoBBColl.size()){
           for (auto htoWW_cand : htoWWColl){
               for (auto htoBB_cand : htoBBColl){
                       const reco::Candidate* htoWW_mother = htoWW_cand->mother();
                       const reco::Candidate* htoBB_mother = htoBB_cand->mother();
                       while (htoWW_mother->pdgId() == 25)  htoWW_mother = htoWW_mother->mother();
                       while (htoBB_mother->pdgId() == 25)  htoBB_mother = htoBB_mother->mother();
                       if (htoWW_mother == htoBB_mother && htoBB_mother->pdgId()==99927){ 
				h2tohhcand = htoWW_mother;
				htoWWcand = htoWW_cand;
				htoBBcand = htoBB_cand;
				h2tohh = true;
				break;
			}
               }
	   	if (h2tohh) break;
           }
     }

     if (h2tohh){
        std::cout << "find h2 candidate " << std::endl;
        std::cout << "h2 candidate id " << h2tohhcand->pdgId() << " mass " << h2tohhcand->mass() << std::endl;
    	mu1cand = stabledecendant(htoWWcand, 13);
        mu2cand = stabledecendant(htoWWcand, -13);
        nu1cand = stabledecendant(htoWWcand, -14);
        nu2cand = stabledecendant(htoWWcand, 14);
        b1cand = alldecendant(htoBBcand, 5, true);
        b2cand = alldecendant(htoBBcand, -5, true);
        w1cand = allancestor(mu1cand, -24, true);
	w2cand = allancestor(mu2cand, 24, true);   
        std::cout <<" mu1 "; printCandidate(mu1cand);	
        std::cout <<" nu1 "; printCandidate(nu1cand);	
        std::cout <<" mu2 "; printCandidate(mu2cand);	
        std::cout <<" nu2 "; printCandidate(nu2cand);	
        std::cout <<" w1 " ; printCandidate(w1cand);
        std::cout <<" w2 " ; printCandidate(w2cand);

     	fillbranches();
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
   evtree->Branch("mu1_mother_mass",&mu1_mother_mass);
   evtree->Branch("mu1_mother_energy",&mu1_mother_energy);
   evtree->Branch("mu1_mother_px",&mu1_mother_px);
   evtree->Branch("mu1_mother_py",&mu1_mother_py);
   evtree->Branch("mu1_mother_pz",&mu1_mother_pz);
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
void
DiHiggsWWAnalyzer::clear(){
   mu1Coll.clear();
   nu1Coll.clear();
   mu2Coll.clear();
   nu2Coll.clear();
   b1Coll.clear();
   b2Coll.clear();
   W1Coll.clear();
   W2Coll.clear();
   htoWWColl.clear();
   htoBBColl.clear();

   mu1cand = NULL;
   nu1cand = NULL;
   mu2cand = NULL;
   nu2cand = NULL;
   b1cand = NULL;
   b2cand = NULL;
   w1cand = NULL;
   w2cand = NULL;
   htoWWcand = NULL;
   htoBBcand = NULL;
   h2tohhcand = NULL;



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


//------------ method called to find decendant with pdgid = id, if first it true, then return the candidate closest to seed
const reco::Candidate* 
DiHiggsWWAnalyzer::alldecendant(const reco::Candidate* cand, int id, bool first){
   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++){
        
	if ((cand->daughter(i))->pdgId() == id && first && cand->pdgId() != id)
		return 	tmp=cand->daughter(i);
	else if ((cand->daughter(i))->pdgId() == id) 
		return  tmp=cand->daughter(i);
        else if (alldecendant(cand->daughter(i),id, first)) 
		return tmp=alldecendant(cand->daughter(i),id);
   }
    
    return tmp;

}

//---------- method called to find a ancestor with pdgid = id, if first is true, then return the candidate closest to seed
const reco::Candidate*
DiHiggsWWAnalyzer::allancestor(const reco::Candidate* cand, int id, bool first){

   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfMothers(); i++){
        
	if ((cand->mother(i))->pdgId() == id && first && cand->pdgId() != id)
		return 	tmp=cand->mother(i);
	else if ((cand->mother(i))->pdgId() == id) 
		return  tmp=cand->mother(i);
        else if (allancestor(cand->mother(i),id, first)) 
		return tmp=allancestor(cand->mother(i),id, first);
   }
   return tmp;

}



//---------- method called to print candidates for debug ---------------------
void
DiHiggsWWAnalyzer::printCandidate(const reco::Candidate* cand){

   std::cout <<" Candidate id: "<< cand->pdgId() << " mass: " << cand->mass() <<" (P,E)= ("<< cand->px() <<", "<< cand->py()<<", "<< cand->pz()<<", "<< cand->energy()
             <<")" << " status: " << cand->status() << std::endl;

}




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
