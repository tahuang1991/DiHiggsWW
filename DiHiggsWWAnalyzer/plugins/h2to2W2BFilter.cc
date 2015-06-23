// -*- C++ -*-
//
// Package:    h2to2W2BFilter
// Class:      h2to2W2BFilter
// 
/**\class h2to2W2BFilter h2to2W2BFilter.cc Zbb/h2to2W2BFilter/src/h2to2W2BFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Aruna Nayak
//         Created:  Thu Aug 23 11:37:45 CEST 2007
// $Id: h2to2W2BFilter.cc,v 1.7 2012/06/05 08:54:21 chanon Exp $
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

class h2to2W2BFilter : public edm::EDFilter {
   public:
      explicit h2to2W2BFilter(const edm::ParameterSet&);
      ~h2to2W2BFilter();

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
   private:
      std::vector<reco::GenParticle*> b1Coll; 
      std::vector<reco::GenParticle*> b2Coll;
      std::vector<reco::GenParticle*> W1Coll;
      std::vector<reco::GenParticle*> W2Coll;
      std::vector<const reco::Candidate*> htoWWColl;
      std::vector<const reco::Candidate*> htoBBColl;

      const reco::Candidate* htoWWcand;
      const reco::Candidate* htoBBcand;

    //private function:
    private:
      void clear();
      void printCandidate(const reco::Candidate* );
      void printallDecendants(const reco::Candidate* );
      bool leptonicaldecay(const reco::Candidate* );
      const reco::Candidate* finddecendant(const reco::Candidate* cand, int id, bool first=false);
      bool hasDaughter(const reco::Candidate* cand,int id );
  // ----------member data ---------------------------
  std::string fLabel_;

};
//
// static data member definitions
//

//
// constructors and destructor
//
h2to2W2BFilter::h2to2W2BFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  fLabel_ = iConfig.getUntrackedParameter("moduleLabel",std::string("prunedGenParticles"));
  //minPtLepton_ = iConfig.getUntrackedParameter<double>("MinPtLepton", 2.0);
  
}


h2to2W2BFilter::~h2to2W2BFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
h2to2W2BFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   bool accept = false;

   Handle<reco::GenParticleCollection> genParticleColl;
//iEvent.getByLabel("prunedGenParticles", genParticleColl);
   iEvent.getByLabel(fLabel_, genParticleColl);
   
   bool h2tohh = false;
   clear();
   for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {

//particle id, (muon13),(b5),(W+24),(SM higgs25)
   // particle id  it->pdgId()
      if (it->pdgId() == 24 ){
        const reco::Candidate* tmpw1 = it->mother();
        while (tmpw1->pdgId() == 24 && tmpw1->numberOfMothers() == 1) tmpw1 = tmpw1->mother();
        if (tmpw1->pdgId() == 25)  W1Coll.push_back(it->clone());
	}
      else if (it->pdgId() == -24 )
      {
        const reco::Candidate* tmpw2 = it->mother();
        while (tmpw2->pdgId() == -24 && tmpw2->numberOfMothers() == 1) tmpw2 = tmpw2->mother();
        if (tmpw2->pdgId() == 25)  W2Coll.push_back(it->clone());
        }

      else if (it->pdgId() == 5 && it->mother()->pdgId() == 25 )
      {
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
          b1Coll.push_back(it->clone());
      }
      else if (it->pdgId() == -5 && it->mother()->pdgId() == 25 )
      {
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
          b2Coll.push_back(it->clone());
      }


   }// all Gen particles


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
                       if (htoWW_mother == htoBB_mother){ 
				htoWWcand = htoWW_cand;
				htoBBcand = htoBB_cand;
				h2tohh = true;
				break;
			}
               }
	   	if (h2tohh) break;
           }
     }
   
    if(h2tohh){
    	const reco::Candidate* W1cand = finddecendant(htoWWcand, 24, false);
    	const reco::Candidate* W2cand = finddecendant(htoWWcand, -24, false);
        //if (leptonicaldecay(W1cand)) std::cout <<" W1 is leptonical decay " << std::endl;     
        //if (leptonicaldecay(W2cand)) std::cout <<" W2 is leptonical decay " << std::endl;     
        accept = (leptonicaldecay(W1cand) and leptonicaldecay(W2cand));
        if (accept) std::cout <<" both W are leptonical decay " << std::endl;

	} 
   //delete evt;
   return accept;
}


void
h2to2W2BFilter::clear(){
   b1Coll.clear();
   b2Coll.clear();
   W1Coll.clear();
   W2Coll.clear();
   htoWWColl.clear();
   htoBBColl.clear();

   htoWWcand = NULL;
   htoBBcand = NULL;

}


//------------ method called to find decendant with pdgid = id, 
//if first it true, then return the candidate closest to seed
//if first it false, then return the candidate farthest to seed
const reco::Candidate* 
h2to2W2BFilter::finddecendant(const reco::Candidate* cand, int id, bool first){
   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++){
        
	if ((cand->daughter(i))->pdgId() == id && first && cand->pdgId() != id)
		return 	tmp=cand->daughter(i);
	else if ((cand->daughter(i))->pdgId() == id && !first && !hasDaughter(cand->daughter(i), id)) 
		return  tmp=cand->daughter(i);
	else if ((cand->daughter(i))->pdgId() == id && !first && (cand->daughter(i))->numberOfDaughters()>1) 
		return  tmp=cand->daughter(i);// tmp has more than one daughters therefore it is final-states
        else if (finddecendant(cand->daughter(i),id, first)) 
		return tmp=finddecendant(cand->daughter(i),id, first);
   }
    
    return tmp;
}

bool
h2to2W2BFilter::leptonicaldecay(const reco::Candidate* cand){
   
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++){
       if (abs((cand->daughter(i))->pdgId())==13 || abs((cand->daughter(i))->pdgId()==11)) return true;
       //printCandidate(cand->daughter(i));
   }
   return false;
}



//---------- method called to print candidates for debug ---------------------
void
h2to2W2BFilter::printCandidate(const reco::Candidate* cand){

   if (!cand) return;
   std::cout <<" Candidate id: "<< cand->pdgId() << " mass: " << cand->mass() <<" (P,E)= ("<< cand->px() <<", "<< cand->py()<<", "<< cand->pz()<<", "<< cand->energy()
             <<")" << " status: " << cand->status() << std::endl;

}


//--------- method called to print all decendants for cand -------------------
void 
h2to2W2BFilter::printallDecendants(const reco::Candidate* cand){
   if (!cand) return;

   if (cand->status() != 0 && cand->numberOfDaughters() > 0){
        std::cout << "******************  children of id "<< cand->pdgId() <<"      *********************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        	printCandidate(cand->daughter(i));
        std::cout << "***********************************************************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
		printallDecendants(cand->daughter(i));

    }
}


//-------- method called to check whether cand has daughter with pdgid = id ------------------------------------
bool
h2to2W2BFilter::hasDaughter(const reco::Candidate* cand, int id){
   
   if (cand->status() ==1 ) return false; 
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        if ((cand->daughter(i))->pdgId() == id) return true;
   return false;

}





// ------------ method called once each job just before starting event loop  ------------

// ------------ method called once each job just after ending the event loop  ------------
void 
h2to2W2BFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(h2to2W2BFilter);
