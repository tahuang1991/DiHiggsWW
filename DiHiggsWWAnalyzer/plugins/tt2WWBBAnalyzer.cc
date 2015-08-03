// -*- C++ -*-
//
// Package:    tt2WWBBAnalyzer
// Class:      tt2WWBBAnalyzer
// 
/**\class tt2WWBBAnalyzer tt2WWBBAnalyzer.cc htoWW/tt2WWBBAnalyzer/plugins/tt2WWBBAnalyzer.cc

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

class tt2WWBBAnalyzer : public edm::EDAnalyzer {
   public:
      explicit tt2WWBBAnalyzer(const edm::ParameterSet&);
      ~tt2WWBBAnalyzer();

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
	//t->bW+ , 6->5,24
      std::vector<reco::GenParticle*> b1Coll; //5
      std::vector<reco::GenParticle*> b2Coll; //-5
      std::vector<reco::GenParticle*> W1Coll;
      std::vector<reco::GenParticle*> W2Coll;
      std::vector<const reco::Candidate*> tColl;
      std::vector<const reco::Candidate*> tbarColl;
      const reco::Candidate* b1cand;
      const reco::Candidate* b2cand;
      const reco::Candidate* w1cand;
      const reco::Candidate* w2cand;
      const reco::Candidate* l1cand;
      const reco::Candidate* l2cand;
      const reco::Candidate* nu1cand;
      const reco::Candidate* nu2cand;
      const reco::Candidate* t1cand;
      const reco::Candidate* t2cand;

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
      
     
    //  float w2_mass;
      float b1_energy;
      float b1_px;
      float b1_py;
      float b1_pz;
      float b2_energy;
      float b2_px;
      float b2_py;
      float b2_pz;


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


      float t1_energy;
      float t1_px;
      float t1_py;
      float t1_pz; 
      float t1_mass;
      float t2_energy;
      float t2_px;
      float t2_py;
      float t2_pz; 
      float t2_mass;
      
      

      float met;
      float met_phi;
      float met_px;
      float met_py;

      //cuts for higgstoWWbb
      bool bquark;
      bool bbarquark;
      bool ttbar;
      bool ttoWb;
      bool tbartoWbbar;
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
tt2WWBBAnalyzer::tt2WWBBAnalyzer(const edm::ParameterSet& iConfig)

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
tt2WWBBAnalyzer::init(){
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


      b1_energy = 0.0;
      b1_px = 0.0;
      b1_py = 0.0;
      b1_pz = 0.0;
      b2_energy = 0.0;
      b2_px = 0.0;
      b2_py = 0.0;
      b2_pz = 0.0;

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

      t1_energy = 0.0;
      t1_px = 0.0;
      t1_py = 0.0;
      t1_pz = 0.0;
      t1_mass = 0.0;
      t2_energy = 0.0;
      t2_px = 0.0;
      t2_py = 0.0;
      t2_pz = 0.0;
      t2_mass = 0.0;


      met = 0.0;
      met_px = 0.0;
      met_py = 0.0;
      met_phi = 0.0;

      bquark = false;
      bbarquark = false;
      ttbar =false;
      ttoWb = false;
      tbartoWbbar = false;
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
      t1cand = 0;
      t2cand = 0;
      genMet = 0;
      b1jet = 0;
      b2jet =0;

      hasbjet =false;
      hasbbarjet =false;
      hasMET =false;
      hastwomuons = false;
      hasdRjets=false;
      hasdRljet =false;
}


tt2WWBBAnalyzer::~tt2WWBBAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
tt2WWBBAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
        //while (tmpw1->pdgId() == 24 && tmpw1->numberOfMothers() == 1) tmpw1 = tmpw1->mother();
        if (tmpw1->pdgId() == 6)  W1Coll.push_back(it->clone());
	}
      else if (it->pdgId() == -24 )
      {
        const reco::Candidate* tmpw2 = it->mother();
        //while (tmpw2->pdgId() == -24 && tmpw2->numberOfMothers() == 1) tmpw2 = tmpw2->mother();
        if (tmpw2->pdgId() == -6)  W2Coll.push_back(it->clone());
        }

      else if (it->pdgId() == 5 && it->mother()->pdgId() == 6 )
      {
	  bquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
          b1Coll.push_back(it->clone());
      }
      else if (it->pdgId() == -5 && it->mother()->pdgId() == -6 )
      {
	  bbarquark = true;
          if (it->numberOfMothers() != 1) std::cout << "bquark has more than one mother particle" << std::endl;
          b2Coll.push_back(it->clone());
      }


   }// all Gen particles


    //t->bW+ 6->5,24
    if (W1Coll.size() && b1Coll.size()){
         for (auto W1_cand : W1Coll)
              for (auto b1_cand : b1Coll){
                      const reco::Candidate* b1_mother = b1_cand->mother();
                      const reco::Candidate* W1_mother = W1_cand->mother();
              //        while (W1_mother->pdgId() == 24) W1_mother = W1_mother->mother();
                      if (W1_mother == b1_mother && W1_mother->pdgId() == 6) {
				tColl.push_back(W1_mother);
				ttoWb = true;
				std::cout <<" find t->bW+" << std::endl;
				break;
			}
              }
    }

     //htoBB
     if (W2Coll.size() && b2Coll.size()){
          for(auto W2_cand : W2Coll)
              for (auto b2_cand : b2Coll) {
                      const reco::Candidate* W2_mother = W2_cand->mother();
                      const reco::Candidate* b2_mother = b2_cand->mother();
                   //   while (W2_mother->pdgId() == -24) W2_mother = W2_mother->mother();
                       if (W2_mother == b2_mother && W2_mother->pdgId() == -6) {
				tbarColl.push_back(W2_mother);
				tbartoWbbar = true;
				std::cout <<" find tbar->bbarW- " << std::endl;
				break;
			}

              }
       
     }
     /*std::cout <<"tColl size " << tColl.size() <<" tbarColl " << tbarColl.size() << std::endl;
     for (unsigned int i = 0; i<tColl.size(); i++){
	std::cout <<" t1cand " ; printCandidate(tColl.at(i));
	}
     for (unsigned int j = 0; j<tbarColl.size(); j++){
	std::cout <<" t2cand " ; printCandidate(tbarColl.at(j));
	}*/

     if (tColl.size()==1 and tbarColl.size() ==1){
	t1cand = tColl.at(0);
	t2cand = tbarColl.at(0);
	while (t1cand->numberOfMothers()==1 and (t1cand->mother())->pdgId()==6) t1cand = t1cand->mother();
	while (t2cand->numberOfMothers()==1 and (t2cand->mother())->pdgId()==-6) t2cand = t2cand->mother();
	std::cout <<"t1 mothers " << t1cand->numberOfMothers() <<" t2 mothers " << t2cand->numberOfMothers() << std::endl; 
	if (t1cand->numberOfMothers() == t2cand->numberOfMothers()) ttbar =true; 
	if (ttbar){
	//	std::cout <<" t1 t2 has same number of mother " << std::endl;
		
		for (unsigned int i=0; i<t1cand->numberOfMothers(); i++){
			bool matchttbar=false;
			for (unsigned int j=0; j<t2cand->numberOfMothers(); j++){
				if ((t1cand->mother(i))->mass()==(t2cand->mother(j))->mass()) {	
					matchttbar = true; 
			//		std::cout <<"t1 mother can match to t2 mother" << std::endl;
					break;
				}
			}
			if (not matchttbar) {	
			 std::cout <<"fail to match t1 mother to t2 mother " << std::endl;
			 ttbar =false;
			 break;
                        }
		}
	}
      }
    
     if (ttbar){
        std::cout << "find t and tbar candidate " << std::endl;
	t1cand = tColl.at(0);
	t2cand = tbarColl.at(0);
	
	if (finalStates_){
       		b1cand = finddecendant(t1cand, 5, false);
        	w1cand = finddecendant(t1cand, 24, false);
        	b2cand = finddecendant(t2cand, -5, false);
		w2cand = finddecendant(t2cand, -24, false);   
	 }else{
       		b1cand = finddecendant(t1cand, 5, true);
        	w1cand = finddecendant(t1cand, 24, true);
        	b2cand = finddecendant(t2cand, -5, true);
		w2cand = finddecendant(t2cand, -24, true);   

	}

	if (Wtolepton(w1cand)) {//w1 pdgid = 24
		if (finddecendant(w1cand, -11)) {
			l1cand = finddecendant(w1cand, -11, false); nu1cand = finddecendant(w1cand, 12, false);}
		if (finddecendant(w1cand, -13)) {
			l1cand = finddecendant(w1cand, -13, false);nu1cand = finddecendant(w1cand, 14, false);}
		if (Wtotau_ and finddecendant(w1cand, -15)) {
			l1cand = finddecendant(w1cand, -15, false);nu1cand = finddecendant(w1cand, 16, false);}
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
			l2cand = finddecendant(w1cand, 15, false); nu2cand = finddecendant(w2cand, -16, false);}
		if (l2cand and nu2cand) W2tolepton=true;
		else W2tolepton =false;
		if (not W2tolepton) std::cout <<" w2 has lepton decendant but program failed to find it " << std::endl;
		else printCandidate(l2cand);
	}
	if (W1tolepton and W2tolepton) std::cout <<"both Ws decay leptonically " << std::endl;

      }

  // match bquark and bjet
	
   if (ttbar){
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
		const reco::Candidate* tcand;
		if  ( dR_bbarjet>deltaR(jetit->eta(), jetit->phi(), b2cand->eta(), b2cand->phi())
		 	and (bcand=findancestor(mcpart, -5, false)) and bcand->mother()->pdgId() == -6){ 
			tcand = bcand->mother();
			if (tcand == t2cand and dR_bbarjet>deltaR(jetit->eta(), jetit->phi(), b2cand->eta(), b2cand->phi())){
				dR_bbarjet = deltaR(jetit->eta(), jetit->phi(), b2cand->eta(), b2cand->phi());
				//printCandidate(jetit->clone());
				//std::cout <<" bcand(-5) "; printCandidate(bcand);
			     std::cout <<" has tbar->Wbbar,h is the same from tbar->Wbbar in genparticles flow,dR "<< dR_bbarjet <<std::endl;
				hasbbarjet = true;
				b2jet = jetit->clone();
			}
   			
		  }
		if ( dR_bjet>deltaR(jetit->eta(), jetit->phi(), b1cand->eta(), b1cand->phi()) 
			and (bcand=findancestor(mcpart, 5, false)) and bcand->mother()->pdgId() == 6){
			tcand = bcand->mother();
			if (tcand == t1cand and dR_bjet>deltaR(jetit->eta(), jetit->phi(), b1cand->eta(), b1cand->phi())){
				dR_bjet = deltaR(jetit->eta(), jetit->phi(), b1cand->eta(), b1cand->phi());
				//printCandidate(jetit->clone());
				//std::cout <<" bcand(5) "; printCandidate(bcand);
				std::cout <<" has t->Wb, t is the same from t->Wb in genparticles flow,dR "<< dR_bjet <<std::endl;
				hasbjet = true;
				b1jet = jetit->clone();
			}
		 }	
	}//mcparts
	//std::cout <<"mcparticle size " <<mcparts.size() <<"   nparts "<< nparts << std::endl;
     }// genjetColl
	if (!hasbbarjet or !hasbjet) std::cout <<" has ttbar,  failed to two b jets " << std::endl;
    }//ttbar

   if (ttbar and hasbjet and hasbbarjet and b1jet->px() == b2jet->px() and b1jet->py() == b2jet->py()) {
	std::cout <<" error: two bjets are in fact the same " << "b1jetpx "<< b1jet->px() <<" b2jetpx "<< b2jet->px()<<std::endl;
	hasbjet = false;
	hasbbarjet =false;
    }

   if (genMet->pt()>= metPt_) hasMET = true;
   else hasMET = false;
   
   if (ttbar and W1tolepton and W2tolepton){
         if (((l1cand->pt() >= muonPt1_ and l2cand->pt() >= muonPt2_) || (l2cand->pt()>=muonPt1_ and l1cand->pt()>=muonPt2_)) and
		fabs(l1cand->eta()) <= muonsEta_ and fabs(l2cand->eta()) <= muonsEta_) 
		hastwomuons = true;
    }
    
   if (ttbar and W1tolepton and W2tolepton and hasbjet and hasbbarjet){
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
   if (ttbar){
     	fillbranches();
        evtree->Fill();
    }

}


// ------------ method called once each job just before starting event loop  ------------
void 
tt2WWBBAnalyzer::beginJob()
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

   
   evtree->Branch("b1_energy",&b1_energy);
   evtree->Branch("b1_px",&b1_px);
   evtree->Branch("b1_py",&b1_py);
   evtree->Branch("b1_pz",&b1_pz);
   evtree->Branch("b2_energy",&b2_energy);
   evtree->Branch("b2_px",&b2_px);
   evtree->Branch("b2_py",&b2_py);
   evtree->Branch("b2_pz",&b2_pz);
   
   evtree->Branch("t1_energy",&t1_energy);
   evtree->Branch("t1_px",&t1_px);
   evtree->Branch("t1_py",&t1_py);
   evtree->Branch("t1_pz",&t1_pz);
   evtree->Branch("t1_mass",&t1_mass);
   evtree->Branch("t2_energy",&t2_energy);
   evtree->Branch("t2_px",&t2_px);
   evtree->Branch("t2_py",&t2_py);
   evtree->Branch("t2_pz",&t2_pz);
   evtree->Branch("t2_mass",&t2_mass);
   
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

   evtree->Branch("met",&met);
   evtree->Branch("met_phi",&met_phi);
   evtree->Branch("met_px",&met_px);
   evtree->Branch("met_py",&met_py);
   
   evtree->Branch("hasMET",&hasMET);
   evtree->Branch("hastwomuons",&hastwomuons);
   evtree->Branch("hasdRljet",&hasdRljet);
   evtree->Branch("hasdRjets",&hasdRjets);
    
   evtree->Branch("ttbar",&ttbar);
   evtree->Branch("ttoWb",&ttoWb);
   evtree->Branch("tbartoWbbar",&tbartoWbbar);
   evtree->Branch("W1tolepton",&W1tolepton);
   evtree->Branch("W2tolepton",&W2tolepton);
    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
tt2WWBBAnalyzer::endJob() 
{
    //output->Write();
   // output->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
tt2WWBBAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
tt2WWBBAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
tt2WWBBAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
tt2WWBBAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
void
tt2WWBBAnalyzer::clear(){
   b1Coll.clear();
   b2Coll.clear();
   W1Coll.clear();
   W2Coll.clear();
   tColl.clear();
   tbarColl.clear();

   b1cand = NULL;
   b2cand = NULL;
   w1cand = NULL;
   w2cand = NULL;
   t1cand = NULL;
   t2cand = NULL;



}



//------------- method called to find stable decendant with pdgid = id
const reco::Candidate* 
tt2WWBBAnalyzer::stabledecendant(const reco::Candidate* cand, int id){
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
tt2WWBBAnalyzer::finddecendant(const reco::Candidate* cand, int id, bool first){
   const reco::Candidate* tmp = NULL;
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++){
        
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
tt2WWBBAnalyzer::findancestor(const reco::Candidate* cand, int id, bool first){

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
tt2WWBBAnalyzer::hasMother(const reco::Candidate* cand, int id){

   for (unsigned int i=0; i < cand->numberOfMothers(); i++)
        if ((cand->mother(i))->pdgId() == id) return true;
   return false;

}

//-------- method called to check whether cand has daughter with pdgid = id ------------------------------------
bool
tt2WWBBAnalyzer::hasDaughter(const reco::Candidate* cand, int id){
 
   for (unsigned int i=0; i < cand->numberOfDaughters(); i++)
        if ((cand->daughter(i))->pdgId() == id) return true;
   return false;

}


//-------- method called to check whether W decays into lepton and if so, return the 
bool
tt2WWBBAnalyzer::Wtolepton(const reco::Candidate* Wcand){
   
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
tt2WWBBAnalyzer::printCandidate(const reco::Candidate* cand){

   std::cout <<" Candidate id: "<< cand->pdgId() << " mass: " << cand->mass() <<" (P,E)= ("<< cand->px() <<", "<< cand->py()<<", "<< cand->pz()<<", "<< cand->energy()
             <<")" << " status: " << cand->status() << std::endl;

}


//--------- method called to print all decendants for cand -------------------
void 
tt2WWBBAnalyzer::printallDecendants(const reco::Candidate* cand){
   
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
tt2WWBBAnalyzer::fillbranches(){
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


      b1_energy = b1cand->energy();
      b1_px = b1cand->px();
      b1_py = b1cand->py();
      b1_pz = b1cand->pz();
      b2_energy = b2cand->energy();
      b2_px = b2cand->px();
      b2_py = b2cand->py();
      b2_pz = b2cand->pz();
	
      t1_energy = t1cand->energy();
      t1_px = t1cand->px();
      t1_py = t1cand->py();
      t1_pz = t1cand->pz();
      t1_mass = t1cand->mass();
      t2_energy = t2cand->energy();
      t2_px = t2cand->px();
      t2_py = t2cand->py();
      t2_pz = t2cand->pz();
      t2_mass = t2cand->mass();

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
tt2WWBBAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(tt2WWBBAnalyzer);
