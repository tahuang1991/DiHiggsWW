// -*- C++ -*-
//
// Package:    ttbarto2l2BFilter
// Class:      ttbarto2l2BFilter
// 
/**\class ttbarto2l2BFilter ttbarto2l2BFilter.cc Zbb/ttbarto2l2BFilter/src/ttbarto2l2BFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Aruna Nayak
//         Created:  Thu Aug 23 11:37:45 CEST 2007
// $Id: ttbarto2l2BFilter.cc,v 1.7 2012/06/05 08:54:21 chanon Exp $
//
//




//
// constants, enums and typedefs
//

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


using namespace edm;
//
// class declaration
//

class ttbarto2l2BFilter : public edm::EDFilter {
   public:
      explicit ttbarto2l2BFilter(const edm::ParameterSet&);
      ~ttbarto2l2BFilter();

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
  // ----------member data ---------------------------
  std::string fLabel_;
  double minEtaLepton_, maxEtaLepton_;
  double minPtLepton_, maxPtLepton_;

};
//
// static data member definitions
//

//
// constructors and destructor
//
ttbarto2l2BFilter::ttbarto2l2BFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  fLabel_ = iConfig.getUntrackedParameter("moduleLabel",std::string("prunedGenParticles"));
  maxEtaLepton_ = iConfig.getUntrackedParameter<double>("MaxEtaLepton", 5);
  minEtaLepton_ = iConfig.getUntrackedParameter<double>("MinEtaLepton", 0);
  maxPtLepton_ = iConfig.getUntrackedParameter<double>("MaxPtLepton", 1000);
  minPtLepton_ = iConfig.getUntrackedParameter<double>("MinPtLepton", 2.0);
  
}


ttbarto2l2BFilter::~ttbarto2l2BFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ttbarto2l2BFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   bool accept = false;

   Handle<reco::GenParticleCollection> genParticleColl;
//iEvent.getByLabel("prunedGenParticles", genParticleColl);
   iEvent.getByLabel(fLabel_, genParticleColl);
   
   int numLeptons = 0;
   int numNeutrinos = 0;
   bool bquark = false;
   bool bbarquark = false;
   for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {

//particle id, (muon13),(b5),(W+24),(SM higgs25)
   // particle id  it->pdgId()
      if (std::abs(it->pdgId()) == 13 && it->status() == 1)
      {
          const reco::Candidate* tmp_mu1 = it->mother(); 
          while (std::abs(tmp_mu1->pdgId()) == 13 && tmp_mu1->numberOfMothers() == 1) tmp_mu1 = tmp_mu1->mother();
          if (tmp_mu1->numberOfMothers() != 1 ) std::cout << "muon has more than one mother particle" << std::endl;
          while (std::abs(tmp_mu1->pdgId()) == 24 && tmp_mu1->numberOfMothers() == 1) tmp_mu1 = tmp_mu1->mother();
          if (tmp_mu1->numberOfMothers() != 1 ) std::cout << "W- has more than one mother particle" << std::endl;
          if (std::abs(tmp_mu1->pdgId()) == 6 && it->eta() > minEtaLepton_ && it->eta() < maxEtaLepton_ && it->pt() < maxPtLepton_ && it->pt() > minPtLepton_)   
             {
		numLeptons++;
              }
      }
      else if (std::abs(it->pdgId()) == 11 && it->status() == 1)
      {
          const reco::Candidate* tmp_e = it->mother(); 
          while (std::abs(tmp_e->pdgId()) == 11 && tmp_e->numberOfMothers() == 1) tmp_e = tmp_e->mother();
          if (tmp_e->numberOfMothers() != 1 ) std::cout << "electron has more than one mother particle" << std::endl;
          while (std::abs(tmp_e->pdgId()) == 24 && tmp_e->numberOfMothers() == 1) tmp_e = tmp_e->mother();
          if (tmp_e->numberOfMothers() != 1 ) std::cout << "W- has more than one mother particle" << std::endl;
          if (std::abs(tmp_e->pdgId()) == 6 && it->eta() > minEtaLepton_ && it->eta() < maxEtaLepton_ && it->pt() < maxPtLepton_ && it->pt() > minPtLepton_)   
             {
		numLeptons++;
              }
      }
      else if (std::abs(it->pdgId()) == 14  && it->status() == 1 )
      {
          const reco::Candidate* tmp_nu1 = it->mother();
          while (std::abs(tmp_nu1->pdgId()) == 14) tmp_nu1 = tmp_nu1->mother();
          while (std::abs(tmp_nu1->pdgId()) == 24) tmp_nu1 = tmp_nu1->mother();
          if (std::abs(tmp_nu1->pdgId()) == 6)
             {
		 numNeutrinos++;
                }
       }
      else if (std::abs(it->pdgId()) == 12  && it->status() == 1 )
      {
          const reco::Candidate* tmp_nu2 = it->mother();
          while (std::abs(tmp_nu2->pdgId()) == 12) tmp_nu2 = tmp_nu2->mother();
          while (std::abs(tmp_nu2->pdgId()) == 24) tmp_nu2 = tmp_nu2->mother();
          if (std::abs(tmp_nu2->pdgId()) == 6)
             {
		 numNeutrinos++;
                }
       }

      else if (it->pdgId() == 5 && it->mother()->pdgId() == 6 && !bquark)
      {
	  bquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
      }
      else if (it->pdgId() == -5 && it->mother()->pdgId() == -6 && !bbarquark)
      {
	  bbarquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
      }

     // std::cout << "test" << std::endl;

   }// all Gen particles

  

   if(numLeptons == 2 && numNeutrinos == 2 && bquark && bbarquark){
       accept = true;
   }
   //delete evt;
   return accept;
}

// ------------ method called once each job just before starting event loop  ------------

// ------------ method called once each job just after ending the event loop  ------------
void 
ttbarto2l2BFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ttbarto2l2BFilter);
