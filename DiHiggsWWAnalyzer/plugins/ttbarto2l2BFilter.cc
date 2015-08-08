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
      
   private:
      
      void printCandidate(const reco::Candidate* );
      void printChildren(const reco::Candidate* );
      void printMothers(const reco::Candidate* );
      void printallDecendants(const reco::Candidate* );
      void printallAncestors(const reco::Candidate* );
      bool hasMother(const reco::Candidate* cand, int id);
      bool hasDaughter(const reco::Candidate* cand, int id);
  // ----------member data ---------------------------
  std::string fLabel_;
  double minEtaLepton_, maxEtaLepton_;
  double minPtLepton_, maxPtLepton_;
  bool Wtotau_;

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
  Wtotau_ = iConfig.getUntrackedParameter<bool>("Wtotau",false);
  
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
  // std::cout <<" Wtotau " << (Wtotau_ ? " is true ":" is false ")<< std::endl;
   for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {

//particle id, (muon13),(b5),(W+24),(SM higgs25)
   // particle id  it->pdgId()
      if (std::abs(it->pdgId()) == 13 && it->status() == 1)
      {
          const reco::Candidate* tmp_mu1 = it->clone(); 
          while (std::abs(tmp_mu1->pdgId()) == 13 and (hasMother(tmp_mu1, 13) or hasMother(tmp_mu1, -13))) { 
	      tmp_mu1 = tmp_mu1->mother();
	      //printMothers(tmp_mu1);
            if (tmp_mu1->numberOfMothers() != 1 )  printMothers(tmp_mu1);
	  }
          if (tmp_mu1->numberOfMothers() == 1 and std::abs(tmp_mu1->mother()->pdgId())==24) tmp_mu1 = tmp_mu1->mother();
	  else if (tmp_mu1->numberOfMothers() != 1)  printMothers(tmp_mu1);
          while (std::abs(tmp_mu1->pdgId()) == 24 and (hasMother(tmp_mu1, 24) or hasMother(tmp_mu1, -24))) 
	  { 	tmp_mu1 = tmp_mu1->mother();
	  //    printMothers(tmp_mu1);
          	if (tmp_mu1->numberOfMothers() != 1 ) printMothers(tmp_mu1);
	  }
          if (tmp_mu1->numberOfMothers() == 1 and std::abs(tmp_mu1->mother()->pdgId())==6) tmp_mu1 = tmp_mu1->mother();
	  else if (tmp_mu1->numberOfMothers() != 1)  printMothers(tmp_mu1);
          if (std::abs(tmp_mu1->pdgId()) == 6 && it->eta() > minEtaLepton_ && it->eta() < maxEtaLepton_ && it->pt() < maxPtLepton_ && it->pt() > minPtLepton_)   
             {
		numLeptons++;
              }
      }
      else if (std::abs(it->pdgId()) == 11 and std::abs(it->mother()->pdgId()) == 24)
      {
          const reco::Candidate* tmp_e = it->mother(); 
          while (std::abs(tmp_e->pdgId()) == 24 and (hasMother(tmp_e, 24) or hasMother(tmp_e, -24))) 
	  { 	tmp_e = tmp_e->mother();
          	if (tmp_e->numberOfMothers() != 1 ) printMothers(tmp_e);
	  }
          if (tmp_e->numberOfMothers() == 1 and std::abs(tmp_e->mother()->pdgId())==6 ) tmp_e = tmp_e->mother();
	  else if (tmp_e->numberOfMothers() != 1)  printMothers(tmp_e);
	  const reco::Candidate* final_e = it->clone();
	  while (final_e->numberOfDaughters()==1 and std::abs(final_e->daughter(0)->pdgId()) == 11) final_e = final_e->daughter(0);
          if (std::abs(tmp_e->pdgId()) == 6 && final_e->eta() > minEtaLepton_ && final_e->eta() < maxEtaLepton_ && final_e->pt() < maxPtLepton_ && final_e->pt() > minPtLepton_)   
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

      else if (Wtotau_ and (std::abs(it->pdgId()) == 15) and std::abs(it->mother()->pdgId()) == 24){

	const reco::Candidate* tmp_w = it->mother();
	//std::cout <<" Wtotau_ "; printCandidate(it->clone());
	//std::cout <<" its mother "; printCandidate(tmp_w);
        while (std::abs(tmp_w->pdgId()) == 24 and (hasMother(tmp_w, -24) or hasMother(tmp_w, 24))) {
	      tmp_w = tmp_w->mother();
	    if (tmp_w->numberOfMothers()>1 ) printMothers(tmp_w);
	}
          if (tmp_w->numberOfMothers() == 1 and std::abs(tmp_w->mother()->pdgId())==6) tmp_w = tmp_w->mother();
	  else if (tmp_w->numberOfMothers() != 1)  printMothers(tmp_w);
	const reco::Candidate* final_tau = it->clone();
	while (final_tau->numberOfDaughters()==1 and std::abs(final_tau->daughter(0)->pdgId()) == 15) {
	    printChildren(final_tau);
	    final_tau = final_tau->daughter(0);
	}
        if (std::abs(tmp_w->pdgId()) == 6 and final_tau->eta() > minEtaLepton_ && final_tau->eta() < maxEtaLepton_ && final_tau->pt() < maxPtLepton_ && final_tau->pt() > minPtLepton_)
		numLeptons++;
        if (std::abs(tmp_w->pdgId()) == 6) std::cout <<"t->Wb and W->tau nu" << std::endl;

       }

      else if (Wtotau_ and (std::abs(it->pdgId()) == 16) and it->status() ==1){
        //std::cout <<" Wtotau_ neutrino "; printCandidate(it->clone());
        const reco::Candidate* tmp_w = it->mother();
        while (std::abs(tmp_w->pdgId()) == 16) tmp_w = tmp_w->mother();
        while (std::abs(tmp_w->pdgId()) == 24) tmp_w = tmp_w->mother();
        if (std::abs(tmp_w->pdgId()) == 6)
                numNeutrinos++;

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
       std::cout <<" find ttbar->bbWW->bblvlv "<<(Wtotau_? " including tau":" not including tau ") << std::endl;
   }
   else std::cout <<"numLeptons " << numLeptons <<"   numNeutrinos  "<< numNeutrinos << std::endl;
   //delete evt;
   return accept;
}

// ------------ method called once each job just before starting event loop  ------------

// ------------ method called once each job just after ending the event loop  ------------
void 
ttbarto2l2BFilter::endJob() {
}


//---------- method called to print candidates for debug ---------------------
void
ttbarto2l2BFilter::printCandidate(const reco::Candidate* cand){

   std::cout <<" Candidate id: "<< cand->pdgId() << " mass: " << cand->mass() <<" (P,E)= ("<< cand->px() <<", "<< cand->py()<<", "<< cand->pz()<<", "<< cand->energy() <<")" <<"(Pt,E) = ("<< cand->pt() <<", "<< cand->eta() <<", "<< cand->phi()<<", "<<cand->energy()<< ")" <<" status: " << cand->status() << std::endl;

}

//--------- method called to print all decendants for cand -------------------
void 
ttbarto2l2BFilter::printallDecendants(const reco::Candidate* cand){
   
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
ttbarto2l2BFilter::printallAncestors(const reco::Candidate* cand){
   
   if (cand->status() != 0 && cand->numberOfMothers() > 0){
        std::cout << "******************  mothers of id "<< cand->pdgId() <<"      *********************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfMothers(); i++)
        	printCandidate(cand->mother(i));
        std::cout << "***********************************************************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfMothers(); i++)
		printallAncestors(cand->mother(i));

    }
}


//--------- method called to print children for cand -------------------
void 
ttbarto2l2BFilter::printChildren(const reco::Candidate* cand){
   
   if (cand->status() != 0 && cand->numberOfDaughters() > 0){
        std::cout << "******************  children of id "<< cand->pdgId() <<"      *********************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        	printCandidate(cand->daughter(i));
        std::cout << "***********************************************************" << std::endl;


    }
}


//--------- method called to print all Ancestors for cand -------------------
void 
ttbarto2l2BFilter::printMothers(const reco::Candidate* cand){
   
   if (cand->status() != 0 && cand->numberOfMothers() > 0){
        std::cout << "******************  mothers of id "<< cand->pdgId() <<"      *********************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfMothers(); i++)
        	printCandidate(cand->mother(i));
        std::cout << "***********************************************************" << std::endl;

    }
}

//---------- method called to check whether cand has mother with pdgid = id -----------------------------------
bool
ttbarto2l2BFilter::hasMother(const reco::Candidate* cand, int id){

   for (unsigned int i=0; i < cand->numberOfMothers(); i++)
        if ((cand->mother(i))->pdgId() == id) return true;
   return false;

}

//-------- method called to check whether cand has daughter with pdgid = id ------------------------------------
bool
ttbarto2l2BFilter::hasDaughter(const reco::Candidate* cand, int id){
 
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        if ((cand->daughter(i))->pdgId() == id) return true;
   return false;

}



//define this as a plug-in
DEFINE_FWK_MODULE(ttbarto2l2BFilter);
