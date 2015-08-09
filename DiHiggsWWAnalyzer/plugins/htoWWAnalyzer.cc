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
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/GenMET.h"

//headers from root lib
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "DiHiggsWW/DiHiggsWWAnalyzer/src/deltaR.h"
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
      const reco::Candidate* l1cand;
      const reco::Candidate* l2cand;
      const reco::Candidate* nu1cand;
      const reco::Candidate* nu2cand;
      const reco::Candidate* htoWWcand;
      const reco::Candidate* htoBBcand;
      const reco::Candidate* h2tohhcand;

      const reco::GenJet* b1jet;
      const reco::GenJet* b2jet;
      const reco::GenMET* genMet;
    private:
      TLorentzVector bjets_lorentz;
      TLorentzVector ll_lorentz;
    private:
      void clear();
      void init();
    private:
     //decendants
     const reco::Candidate* stabledecendant(const reco::Candidate* cand, int id);
     //const reco::Candidate* stablehtoWWdecendant(Particle p, PdgId id);
     const reco::Candidate* finddecendant(const reco::Candidate* cand, int id, bool first=false);
     const reco::Candidate* findancestor(const reco::Candidate* cand, int id, bool first=false);
     bool hasMother(const reco::Candidate* cand, int id);
     bool hasDaughter(const reco::Candidate* cand, int id);
     bool Wtolepton(const reco::Candidate* Wcand);
                    
    private:
      void printCandidate(const reco::Candidate* );
      void printallDecendants(const reco::Candidate* );
      void printChildren(const reco::Candidate* );
      void printMothers(const reco::Candidate* );
      



    private:
      // ----------member data ---------------------------
      TTree *evtree;
      TFile *output;

    private:
      bool finalStates_;    
     std::string jetLabel_;
     std::string metLabel_;
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
      bool Wtotau_;    
    private:
      void fillbranches(); 
    private: 
      edm::Service< TFileService > fs;
      //branches of tree
      float l1_energy;//e, mu or tau, not including neutrino
      float l1_px;
      float l1_py;
      float l1_pz;
      float l1_eta;
      float l1_phi;
      float l1_mass;
      float l1_id;
      float nu1_energy;
      float nu1_px;
      float nu1_py;
      float nu1_pz;
      float nu1_eta;
      float nu1_phi;
      float nu1_mass;
      float nu1_id;
      float l2_energy;
      float l2_px;
      float l2_py;
      float l2_pz;
      float l2_eta;
      float l2_phi;
      float l2_mass;
      float l2_id;
      float nu2_energy;
      float nu2_px;
      float nu2_py;
      float nu2_pz;
      float nu2_eta;
      float nu2_phi;
      float nu2_mass;
      float nu2_id;
      
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

      float b1jet_energy;
      float b1jet_px;
      float b1jet_py;
      float b1jet_pz;
      float b1jet_pt;
      float b1jet_eta;
      float b1jet_phi;
      float b1jet_mass;
      float b2jet_energy;
      float b2jet_px;
      float b2jet_py;
      float b2jet_pz;
      float b2jet_pt;
      float b2jet_eta;
      float b2jet_phi;
      float b2jet_mass;
      float dR_bjet;
      float dR_bbarjet;

      float dR_b1l1;
      float dR_b1l2;
      float dR_b2l1;
      float dR_b2l2;
      float dR_l1l2;
      float dR_b1b2;
      float mass_l1l2;
      float mass_b1b2;

      float h2tohh_energy;
      float h2tohh_px;
      float h2tohh_py;
      float h2tohh_pz;
      float h2tohh_mass;

      float met;
      float met_phi;
      float met_px;
      float met_py;

      //cuts for higgstoWWbb
      bool bquark;
      bool bbarquark;
      bool htobb;
      bool htoWW;
      bool h2tohh;
      bool W1tolepton;
      bool W2tolepton;

      bool hasbjet;
      bool hasbbarjet;
      bool hasMET;
      bool hastwomuons;
      bool hasdRljet;
      bool hasdRjets;



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
     Wtotau_ = iConfig.getParameter<bool>("Wtotau");
   //now do what ever initialization is needed
   // initilize candidates pointer
 }


//init tree
void 
htoWWAnalyzer::init(){
      l1_energy = 0;
      l1_px =0;
      l1_py =0;
      l1_pz =0;
      l1_eta =0;
      l1_phi =0;
      l1_mass =0;
      l1_id =0;
      l2_energy = 0;
      l2_px =0;
      l2_py =0;
      l2_pz =0;
      l2_eta =0;
      l2_phi =0;
      l2_mass =0;
      l2_id =0;

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

      b1jet_energy = 0.0;
      b1jet_px = 0.0;
      b1jet_py = 0.0;
      b1jet_pz = 0.0;
      b1jet_mass = 0.0;
      b2jet_energy = 0.0;
      b2jet_px = 0.0;
      b2jet_py = 0.0;
      b2jet_pz = 0.0;
      b2jet_mass = 0.0;
       
      dR_bjet = jetsDeltaR_; 
      dR_bbarjet = jetsDeltaR_; 
      dR_b1l1 =0;
      dR_b1l2 = 0;
      dR_b2l1 =0;
      dR_b2l2 =0;
      dR_l1l2 =0;
      dR_b1b2 =0;
      mass_l1l2 =0;
      mass_b1b2 =0;


      h2tohh_energy = 0.0;
      h2tohh_px = 0.0;
      h2tohh_py = 0.0;
      h2tohh_pz = 0.0;
      h2tohh_mass = 0.0;

      met = 0.0;
      met_px = 0.0;
      met_py = 0.0;
      met_phi = 0.0;

      bquark = false;
      bbarquark = false;
      htobb = false;
      htoWW = false;
      h2tohh = false;
      W1tolepton =false;
      W2tolepton =false;
      b1cand = 0;
      b2cand = 0;
      w1cand = 0;
      w2cand = 0;
      l1cand = 0;
      l2cand = 0;
      nu1cand = 0;
      nu2cand = 0;
      htoBBcand = 0;
      htoWWcand = 0;
      h2tohhcand = 0;
      genMet = 0;
      b1jet = 0;
      b2jet =0;

      hasbjet =false;
      hasbbarjet =false;
      hasMET =false;
      hastwomuons = false;
      hasdRljet = false;
      hasdRjets =false;

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
    
   edm::Handle<reco::GenJetCollection> genjetColl;
   iEvent.getByLabel(jetLabel_, genjetColl);

   edm::Handle<edm::View<reco::GenMET> > genmetColl; 
   iEvent.getByLabel(metLabel_, genmetColl);


   clear();
   init();
   
   genMet = &(genmetColl->front());

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
          b1Coll.push_back(it->clone());
      }
      else if (it->pdgId() == -5 && it->mother()->pdgId() == 25 )
      {
	  bbarquark = true;
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
				htoWW = true;
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
				htobb = true;
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
                       if (htoWW_mother == htoBB_mother && (htoBB_mother->pdgId()==99927 or htoBB_mother->pdgId()==99926)){ 
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
//	printallDecendants(h2tohhcand);
	if (Wtolepton(w1cand)) {//w1 pdgid = 24
		if (finddecendant(w1cand, -11)) {
			l1cand = finddecendant(w1cand, -11, false); nu1cand = finddecendant(w1cand, 12, false);}
		if (finddecendant(w1cand, -13)) {
			l1cand = finddecendant(w1cand, -13, false);nu1cand = finddecendant(w1cand, 14, false);}
		if (Wtotau_ and finddecendant(w1cand, -15)) {
			l1cand = finddecendant(w1cand, -15, false);nu1cand = finddecendant(w1cand, 16, false);
			std::cout <<" find tau from W+ ";}
		if (l1cand and nu1cand) W1tolepton=true;
		else W1tolepton =false;
		if (not W1tolepton) std::cout <<" w1 has lepton decendant but program failed to find it " << std::endl;
		else printCandidate(l1cand);
	}
	if (Wtolepton(w2cand)) {//w2 pdgid = -24
		if (finddecendant(w2cand, 11)) {
			l2cand = finddecendant(w2cand, 11, false); nu2cand = finddecendant(w2cand, -12, false);}
		if (finddecendant(w2cand, 13)) {
			l2cand = finddecendant(w2cand, 13, false);nu2cand = finddecendant(w2cand, -14, false);}
		if (Wtotau_ and finddecendant(w2cand, 15)) {
			l2cand = finddecendant(w2cand, 15, false); nu2cand = finddecendant(w2cand, -16, false);
			std::cout <<" find tau from W- ";}
		if (l2cand and nu2cand) W2tolepton=true;
		else W2tolepton =false;
		if (not W2tolepton) std::cout <<" w2 has lepton decendant but program failed to find it " << std::endl;
		else printCandidate(l2cand);
	}
	if (W1tolepton and W2tolepton) std::cout <<"both Ws decay leptonically " << std::endl;

      }

  // match bquark and bjet
	
   if (h2tohh){
   
   for (reco::GenJetCollection::const_iterator jetit = genjetColl->begin(); jetit != genjetColl->end(); jetit++){
	//cuts on GenJets
	/*if (jetit->pt() >= jetsPt_ and std::fabs(jetit->eta()) <= jetsEta_){	
 		totjets_px += jetit->px();
		totjets_py += jetit->py();
		totjets_pz += jetit->pz();
		totjets_energy += jetit->energy();

	}*/
	if (jetit->pt() < bjetsPt_ or std::fabs(jetit->eta())> bjetsEta_) continue;
        std::vector <const reco::GenParticle*> mcparts = jetit->getGenConstituents();
  	for (unsigned i = 0; i < mcparts.size(); i++) {
    		const reco::GenParticle* mcpart = mcparts[i];
		const reco::Candidate* bcand;
		const reco::Candidate* hcand;
		//std::cout <<"GenP id "<< mcpart->pdgId() <<" mass "<< mcpart->mass()  <<" energy "<< mcpart->energy() <<std::endl;
		if  ( dR_bbarjet>deltaR(jetit->eta(), jetit->phi(), b2cand->eta(), b2cand->phi())
		 	and (bcand=findancestor(mcpart, -5, false)) and bcand->mother()->pdgId() == 25){ 
			
			hcand = bcand->mother();
			if (hcand == htoBBcand){
				dR_bbarjet = deltaR(jetit->eta(), jetit->phi(), b2cand->eta(), b2cand->phi());
				//printCandidate(jetit->clone());
				//std::cout <<" bcand(-5) "; printCandidate(bcand);
				//std::cout <<" has h->bbar,h is the same from h->bb in genparticles flow,dR "<< dR_bbarjet <<std::endl;
				hasbbarjet = true;
				b2jet = jetit->clone();
			}
   			
		  }
		if ( dR_bjet>deltaR(jetit->eta(), jetit->phi(), b1cand->eta(), b1cand->phi()) 
			and (bcand=findancestor(mcpart, 5, false)) and bcand->mother()->pdgId() == 25){
			hcand = bcand->mother();
			if (hcand == htoBBcand){
				dR_bjet = deltaR(jetit->eta(), jetit->phi(), b1cand->eta(), b1cand->phi());
				//printCandidate(jetit->clone());
				//std::cout <<" bcand(5) "; printCandidate(bcand);
				//std::cout <<" has h->b, h is the same from h->bb in genparticles flow,dR "<< dR_bjet <<std::endl;
				hasbjet = true;
				b1jet = jetit->clone();
			}
		 }	
	}//mcparts
	//std::cout <<"mcparticle size " <<mcparts.size() <<"   nparts "<< nparts << std::endl;
     }// genjetColl
	if (!hasbbarjet or !hasbjet) std::cout <<" has h->bb, but failed to two b jets " << std::endl;
    }//htobb

   if (htobb and hasbjet and hasbbarjet and b1jet->px() == b2jet->px() and b1jet->py() == b2jet->py()) {
	std::cout <<" error: two bjets are in fact the same " << "b1jetpx "<< b1jet->px() <<" b2jetpx "<< b2jet->px()<<std::endl;
	hasbjet = false;
	hasbbarjet =false;
    }

   if (genMet->pt()>= metPt_) hasMET = true;
   else hasMET = false;
   
   if (h2tohh and W1tolepton and W2tolepton){
         if (((l1cand->pt() >= muonPt1_ and l2cand->pt() >= muonPt2_) || (l2cand->pt()>=muonPt1_ and l1cand->pt()>=muonPt2_)) and
		fabs(l1cand->eta()) <= muonsEta_ and fabs(l2cand->eta()) <= muonsEta_) 
		hastwomuons = true;
    }
    
   if (h2tohh and W1tolepton and W2tolepton and hasbjet and hasbbarjet){
	dR_b1l1 = deltaR(b1jet->eta(),b1jet->phi(),l1cand->eta(),l1cand->phi());
	dR_b1l2 = deltaR(b1jet->eta(),b1jet->phi(),l2cand->eta(),l2cand->phi());
	dR_b2l1 = deltaR(b2jet->eta(),b2jet->phi(),l1cand->eta(),l1cand->phi());
	dR_b2l2 = deltaR(b2jet->eta(),b2jet->phi(),l2cand->eta(),l2cand->phi());
	dR_l1l2 = deltaR(l1cand->eta(),l1cand->phi(),l2cand->eta(),l2cand->phi());
        dR_b1b2 = deltaR(b1jet->eta(),b1jet->phi(),b2jet->eta(),b2jet->phi());
        bjets_lorentz = TLorentzVector(b1jet->px()+b2jet->px(),b1jet->py()+b2jet->py(),b1jet->pz()+b2jet->pz(),b1jet->energy()+b2jet->energy());
	mass_b1b2 = bjets_lorentz.M();
        ll_lorentz = TLorentzVector(l1cand->px()+l2cand->px(),l1cand->py()+l2cand->py(),l1cand->pz()+l2cand->pz(),l1cand->energy()+l2cand->energy());
	mass_l1l2 = ll_lorentz.M();
        if (dR_b1l1>jetleptonDeltaR_ and dR_b1l2>jetleptonDeltaR_ and dR_b2l1>jetleptonDeltaR_ and dR_b2l2>jetleptonDeltaR_) hasdRljet =true;
        if (dR_b1b2>jetsDeltaR_) hasdRjets = true;

   }
        //std::cout <<" w1 " ; printCandidate(w1cand);
        //std::cout <<" w2 " ; printCandidate(w2cand);
        //std::cout <<" b1 " ; printCandidate(b1cand);
        //std::cout <<" b2 " ; printCandidate(b2cand);
   if (h2tohh){
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
   evtree->Branch("l1_energy",&l1_energy);
   evtree->Branch("l1_px",&l1_px);
   evtree->Branch("l1_py",&l1_py);
   evtree->Branch("l1_pz",&l1_pz);
   evtree->Branch("l1_eta",&l1_eta);
   evtree->Branch("l1_phi",&l1_phi);
   evtree->Branch("l1_mass",&l1_mass);
   evtree->Branch("l1_id",&l1_id);
   evtree->Branch("nu1_energy",&nu1_energy);
   evtree->Branch("nu1_px",&nu1_px);
   evtree->Branch("nu1_py",&nu1_py);
   evtree->Branch("nu1_pz",&nu1_pz);
   evtree->Branch("nu1_eta",&nu1_eta);
   evtree->Branch("nu1_phi",&nu1_phi);
   evtree->Branch("nu1_mass",&nu1_mass);
   evtree->Branch("nu1_id",&nu1_id);

   evtree->Branch("l2_energy",&l2_energy);
   evtree->Branch("l2_px",&l2_px);
   evtree->Branch("l2_py",&l2_py);
   evtree->Branch("l2_pz",&l2_pz);
   evtree->Branch("l2_eta",&l2_eta);
   evtree->Branch("l2_phi",&l2_phi);
   evtree->Branch("l2_mass",&l2_mass);
   evtree->Branch("l2_id",&l2_id);
   evtree->Branch("nu2_energy",&nu2_energy);
   evtree->Branch("nu2_px",&nu2_px);
   evtree->Branch("nu2_py",&nu2_py);
   evtree->Branch("nu2_pz",&nu2_pz);
   evtree->Branch("nu2_eta",&nu2_eta);
   evtree->Branch("nu2_phi",&nu2_phi);
   evtree->Branch("nu2_mass",&nu2_mass);
   evtree->Branch("nu2_id",&nu2_id);


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
   
   evtree->Branch("hasbjet",&hasbjet);
   evtree->Branch("b1jet_energy",&b1jet_energy);
   evtree->Branch("b1jet_px",&b1jet_px);
   evtree->Branch("b1jet_py",&b1jet_py);
   evtree->Branch("b1jet_pz",&b1jet_pz);
   evtree->Branch("b1jet_pt",&b1jet_pt);
   evtree->Branch("b1jet_eta",&b1jet_eta);
   evtree->Branch("b1jet_phi",&b1jet_phi);
   evtree->Branch("b1jet_mass",&b1jet_mass);
   evtree->Branch("hasbbarjet",&hasbbarjet);
   evtree->Branch("b2jet_energy",&b2jet_energy);
   evtree->Branch("b2jet_px",&b2jet_px);
   evtree->Branch("b2jet_py",&b2jet_py);
   evtree->Branch("b2jet_pz",&b2jet_pz);
   evtree->Branch("b2jet_pt",&b2jet_pt);
   evtree->Branch("b2jet_eta",&b2jet_eta);
   evtree->Branch("b2jet_phi",&b2jet_phi);
   evtree->Branch("b2jet_mass",&b2jet_mass);
   evtree->Branch("dR_bjet",&dR_bjet);
   evtree->Branch("dR_bbarjet",&dR_bbarjet);
   evtree->Branch("dR_b1l1",&dR_b1l1);
   evtree->Branch("dR_b1l2",&dR_b1l2);
   evtree->Branch("dR_b2l1",&dR_b2l1);
   evtree->Branch("dR_b2l2",&dR_b2l2);
   evtree->Branch("dR_l1l2",&dR_l1l2);
   evtree->Branch("dR_b1b2",&dR_b1b2);
   evtree->Branch("mass_l1l2",&mass_l1l2);
   evtree->Branch("mass_b1b2",&mass_b1b2);

   evtree->Branch("h2tohh_energy",&h2tohh_energy);
   evtree->Branch("h2tohh_px",&h2tohh_px);
   evtree->Branch("h2tohh_py",&h2tohh_py);
   evtree->Branch("h2tohh_pz",&h2tohh_pz);
   evtree->Branch("h2tohh_mass",&h2tohh_mass);
   
   evtree->Branch("met",&met);
   evtree->Branch("met_phi",&met_phi);
   evtree->Branch("met_px",&met_px);
   evtree->Branch("met_py",&met_py);
   
   evtree->Branch("hasMET",&hasMET);
   evtree->Branch("hastwomuons",&hastwomuons);
   evtree->Branch("hasdRljet",&hasdRljet);
   evtree->Branch("hasdRjets",&hasdRjets);
  
   evtree->Branch("htobb",&htobb);
   evtree->Branch("htoWW",&htoWW);
   evtree->Branch("W1tolepton",&W1tolepton);
   evtree->Branch("W2tolepton",&W2tolepton);
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
        //std::cout <<"to find decendant id "<< id <<" daughter "<< i <<" child id " << (cand->daughter(i))->pdgId()
	//	<< " numberOfdaghter for this child " << (cand->daughter(i))->numberOfDaughters() << std::endl; 
	if ((cand->daughter(i))->pdgId() == id && first && cand->pdgId() != id)
		return 	tmp=cand->daughter(i);
	else if ((cand->daughter(i))->pdgId() == id && !first && !hasDaughter(cand->daughter(i), id)) 
		return  tmp=cand->daughter(i);
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
htoWWAnalyzer::findancestor(const reco::Candidate* cand, int id, bool first){

   const reco::Candidate* tmp = NULL;
   //std::cout <<"find ancestor id "<<id ; printCandidate(cand);
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


//-------- method called to check whether W decays into lepton and if so, return the 
bool
htoWWAnalyzer::Wtolepton(const reco::Candidate* Wcand){
   
   if (not hasDaughter(Wcand, Wcand->pdgId())){
     for (unsigned int i=0; i < Wcand->numberOfDaughters(); i++)
	if(abs((Wcand->daughter(i))->pdgId()) == 11 or abs((Wcand->daughter(i))->pdgId()) == 13 or 
	   (Wtotau_ and abs((Wcand->daughter(i))->pdgId()) == 15))
		return true;
    }
   else {
  //no lepton found, is this a medium W
   std::cout <<"W has a W child " << std::endl;
    for (unsigned int i=0; i < Wcand->numberOfDaughters(); i++)
	if ((Wcand->daughter(i))->pdgId() == Wcand->pdgId()) return Wtolepton(Wcand->daughter(i));
  }
   
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

//--------- method called to print children for cand -------------------
void 
htoWWAnalyzer::printChildren(const reco::Candidate* cand){
   
   if (cand->status() != 0 && cand->numberOfDaughters() > 0){
        std::cout << "******************  children of id "<< cand->pdgId() <<"      *********************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        	printCandidate(cand->daughter(i));
        std::cout << "***********************************************************" << std::endl;


    }
}


//--------- method called to print all Ancestors for cand -------------------
void 
htoWWAnalyzer::printMothers(const reco::Candidate* cand){
   
   if (cand->status() != 0 && cand->numberOfMothers() > 0){
        std::cout << "******************  mothers of id "<< cand->pdgId() <<"      *********************" << std::endl;
   	for (unsigned int i=0; i < cand->numberOfMothers(); i++)
        	printCandidate(cand->mother(i));
        std::cout << "***********************************************************" << std::endl;

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
	
      if (W1tolepton){
      l1_energy = l1cand->energy();
      l1_px = l1cand->px();
      l1_py = l1cand->py();
      l1_pz = l1cand->pz();
      l1_eta = l1cand->eta();
      l1_phi = l1cand->phi();
      l1_mass = l1cand->mass();
      l1_id = l1cand->pdgId();
      nu1_energy = nu1cand->energy();
      nu1_px = nu1cand->px();
      nu1_py = nu1cand->py();
      nu1_pz = nu1cand->pz();
      nu1_eta = nu1cand->eta();
      nu1_phi = nu1cand->phi();
      nu1_mass = nu1cand->mass();
      nu1_id = nu1cand->pdgId();
	}
      if (W2tolepton){
      l2_energy = l2cand->energy();
      l2_px = l2cand->px();
      l2_py = l2cand->py();
      l2_pz = l2cand->pz();
      l2_eta = l2cand->eta();
      l2_phi = l2cand->phi();
      l2_mass = l2cand->mass();
      l2_id = l2cand->pdgId();
      nu2_energy = nu2cand->energy();
      nu2_px = nu2cand->px();
      nu2_py = nu2cand->py();
      nu2_pz = nu2cand->pz();
      nu2_eta = nu2cand->eta();
      nu2_phi = nu2cand->phi();
      nu2_mass = nu2cand->mass();
      nu2_id = nu2cand->pdgId();
     }
    
     if (hasbbarjet){
      b2jet_px = b2jet->px();
      b2jet_py = b2jet->py();
      b2jet_pz = b2jet->pz();
      b2jet_pt = b2jet->pt();
      b2jet_eta = b2jet->eta();
      b2jet_phi = b2jet->phi();
      b2jet_energy = b2jet->energy();
      b2jet_mass = b2jet->mass();
      }
      if (hasbjet){
      b1jet_px = b1jet->px();
      b1jet_py = b1jet->py();
      b1jet_pz = b1jet->pz();
      b1jet_pt = b1jet->pt();
      b1jet_eta = b1jet->eta();
      b1jet_phi = b1jet->phi();
      b1jet_energy = b1jet->energy();
      b1jet_mass = b1jet->mass();
      }

      if (genMet){
      met = genMet->pt();
      met_px = genMet->px();
      met_py = genMet->py();
      met_phi = genMet->phi();
      }

   
   
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
