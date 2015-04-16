// -*- C++ -*-
//
// Package:    htoWWAnalyzer
// Class:      htoWWAnalyzer
// 
/**\class htoWWAnalyzer htoWWAnalyzer.cc htoWW/htoWWAnalyzer/plugins/htoWWAnalyzer.cc

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

class htoWWAnalyzer : public edm::EDAnalyzer {
   public:
      explicit htoWWAnalyzer(const edm::ParameterSet&);
      ~htoWWAnalyzer();

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
      std::vector<reco::GenParticle*> b1Coll; 
      std::vector<reco::GenParticle*> b2Coll;
      std::vector<reco::GenParticle*> W1Coll;
      std::vector<reco::GenParticle*> W2Coll;
      std::vector<const reco::Candidate*> htoWWColl;
      std::vector<const reco::Candidate*> htoBBColl;
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
     const reco::Candidate* finddecendant(const reco::Candidate* cand, int id, bool first=false);
     const reco::Candidate* findancestor(const reco::Candidate* cand, int id, bool first=false);
     bool hasMother(const reco::Candidate* cand, int id);
     bool hasDaughter(const reco::Candidate* cand, int id);
                    
    private:
      void printCandidate(const reco::Candidate* );
      void printallDecendants(const reco::Candidate* );
      



    private:
      // ----------member data ---------------------------
      TTree *evtree;
      TFile *output;

    private:
      bool finalStates_;    
    private:
      void fillbranches(); 
    private: 
      edm::Service< TFileService > fs;
      //branches of tree
      float w1_energy;
      float w1_px;
      float w1_py;
      float w1_pz;
      float w1_mass;
      float w2_energy;
      float w2_px;
      float w2_py;
      float w2_pz;
      float w2_mass;
      
     
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
      float b2_energy;
      float b2_px;
      float b2_py;
      float b2_pz;

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
      bool bquark;
      bool bbarquark;
      bool htobb;
      bool htoWW;
      bool h2tohh;




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
htoWWAnalyzer::htoWWAnalyzer(const edm::ParameterSet& iConfig)

{


     finalStates_ = iConfig.getParameter<bool>("finalStates");
   //now do what ever initialization is needed
      w1_energy = 0.0;
      w1_px = 0.0;
      w1_py = 0.0;
      w1_pz = 0.0;
      w1_mass = 0.0;

      w2_energy = 0.0;
      w2_px = 0.0;
      w2_py = 0.0;
      w2_pz = 0.0;
      w2_mass = 0.0;

      htoWW_energy = 0.0;
      htoWW_px = 0.0;
      htoWW_py = 0.0;
      htoWW_pz = 0.0;
      htoWW_mass = 0.0;

      b1_energy = 0.0;
      b1_px = 0.0;
      b1_py = 0.0;
      b1_pz = 0.0;
      b2_energy = 0.0;
      b2_px = 0.0;
      b2_py = 0.0;
      b2_pz = 0.0;
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

      bquark = false;
      bbarquark = false;
      htobb = false;
      htoWW = false;

}


htoWWAnalyzer::~htoWWAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
htoWWAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    
      bquark = false;
      bbarquark = false;
      htobb = false;
      htoWW = false;
      h2tohh = false;

    clear();


   for (reco::GenParticleCollection::const_iterator it = genParticleColl->begin(); it != genParticleColl->end(); ++it) {

//particle id, (muon13),(b5),(W+24),(SM higgs25)
   // particle id  it->pdgId()
   //
//      std::cout << "Gen paticles: id " << it->pdgId() << std::endl; 
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
	  bquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
       //   std::cout << "find bquark candidate" << std::endl;
          b1Coll.push_back(it->clone());
      }
      else if (it->pdgId() == -5 && it->mother()->pdgId() == 25 )
      {
	  bbarquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
         // std::cout << "find bbarquark candidate" << std::endl;
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
	if (finalStates_){
       		b1cand = finddecendant(htoBBcand, 5, false);
        	b2cand = finddecendant(htoBBcand, -5, false);
        	w1cand = finddecendant(htoWWcand, 24, false);
		w2cand = finddecendant(htoWWcand, -24, false);   
	 }else{
       		b1cand = finddecendant(htoBBcand, 5, true);
        	b2cand = finddecendant(htoBBcand, -5, true);
        	w1cand = finddecendant(htoWWcand, 24, true);
		w2cand = finddecendant(htoWWcand, -24, true);   

	}
        std::cout <<" w1 " ; printCandidate(w1cand);
        std::cout <<" w2 " ; printCandidate(w2cand);
        std::cout <<" b1 " ; printCandidate(b1cand);
        std::cout <<" b2 " ; printCandidate(b2cand);

     	fillbranches();
        evtree->Fill();
     }

}


// ------------ method called once each job just before starting event loop  ------------
void 
htoWWAnalyzer::beginJob()
{
   evtree = fs->make<TTree>("evtree", "evtree");
 //  output = new TFile("output.root","recreate");
  // output->cd();

 //  evtree = new TTree("evtree","event tree");
   evtree->Branch("w1_energy",&w1_energy);
   evtree->Branch("w1_px",&w1_px);
   evtree->Branch("w1_py",&w1_py);
   evtree->Branch("w1_pz",&w1_pz);
   evtree->Branch("w1_mass",&w1_mass);
   evtree->Branch("w2_energy",&w2_energy);
   evtree->Branch("w2_px",&w2_px);
   evtree->Branch("w2_py",&w2_py);
   evtree->Branch("w2_pz",&w2_pz);
   evtree->Branch("w2_mass",&w2_mass);

   evtree->Branch("htoWW_energy",&htoWW_energy);
   evtree->Branch("htoWW_px",&htoWW_px);
   evtree->Branch("htoWW_py",&htoWW_py);
   evtree->Branch("htoWW_pz",&htoWW_pz);
   evtree->Branch("htoWW_mass",&htoWW_mass);
   
   evtree->Branch("b1_energy",&b1_energy);
   evtree->Branch("b1_px",&b1_px);
   evtree->Branch("b1_py",&b1_py);
   evtree->Branch("b1_pz",&b1_pz);
   evtree->Branch("b2_energy",&b2_energy);
   evtree->Branch("b2_px",&b2_px);
   evtree->Branch("b2_py",&b2_py);
   evtree->Branch("b2_pz",&b2_pz);
   
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
htoWWAnalyzer::endJob() 
{
    //output->Write();
   // output->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
htoWWAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
htoWWAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
htoWWAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
htoWWAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
void
htoWWAnalyzer::clear(){
   b1Coll.clear();
   b2Coll.clear();
   W1Coll.clear();
   W2Coll.clear();
   htoWWColl.clear();
   htoBBColl.clear();

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
htoWWAnalyzer::stabledecendant(const reco::Candidate* cand, int id){
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
htoWWAnalyzer::finddecendant(const reco::Candidate* cand, int id, bool first){
   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++){
        
	if ((cand->daughter(i))->pdgId() == id && first && cand->pdgId() != id)
		return 	tmp=cand->daughter(i);
	else if ((cand->daughter(i))->pdgId() == id && !first && !hasDaughter(cand->daughter(i), id)) 
		return  tmp=cand->daughter(i);
	else if ((cand->daughter(i))->pdgId() == id && !first && (cand->daughter(i))->numberOfDaughters()>1) 
		return  tmp=cand->daughter(i);// tmp has more one daughters therefore it is final-states
        else if (finddecendant(cand->daughter(i),id, first)) 
		return tmp=finddecendant(cand->daughter(i),id);
   }
    
    return tmp;

}

//---------- method called to find a ancestor with pdgid = id, 
//if first is true, then return the candidate closest to seed
//if first is false, then return the candidate furthest to seed
const reco::Candidate*
htoWWAnalyzer::findancestor(const reco::Candidate* cand, int id, bool first){

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
htoWWAnalyzer::hasMother(const reco::Candidate* cand, int id){

   for (unsigned int i=0; i < cand->numberOfMothers(); i++)
        if ((cand->mother(i))->pdgId() == id) return true;
   return false;

}

//-------- method called to check whether cand has daughter with pdgid = id ------------------------------------
bool
htoWWAnalyzer::hasDaughter(const reco::Candidate* cand, int id){
 
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        if ((cand->daughter(i))->pdgId() == id) return true;
   return false;

}



//---------- method called to print candidates for debug ---------------------
void
htoWWAnalyzer::printCandidate(const reco::Candidate* cand){

   std::cout <<" Candidate id: "<< cand->pdgId() << " mass: " << cand->mass() <<" (P,E)= ("<< cand->px() <<", "<< cand->py()<<", "<< cand->pz()<<", "<< cand->energy()
             <<")" << " status: " << cand->status() << std::endl;

}


//--------- method called to print all decendants for cand -------------------
void 
htoWWAnalyzer::printallDecendants(const reco::Candidate* cand){
   
   if (cand->status() != 0 && cand->numberOfDaughters() > 0){
        std::cout << "******************  children of id "<< cand->pdgId() <<"      *********************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        	printCandidate(cand->daughter(i));
        std::cout << "***********************************************************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
		printallDecendants(cand->daughter(i));

    }
}



//---------- method called to fill branches --------------------------------------------
void 
htoWWAnalyzer::fillbranches(){
      w1_energy = w1cand->energy();
      w1_px = w1cand->px();
      w1_py = w1cand->py();
      w1_pz = w1cand->pz();
      w1_mass = w1cand->mass();
      w2_energy = w2cand->energy();
      w2_px = w2cand->px();
      w2_py = w2cand->py();
      w2_pz = w2cand->pz();
      w2_mass = w2cand->mass();

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
htoWWAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(htoWWAnalyzer);
