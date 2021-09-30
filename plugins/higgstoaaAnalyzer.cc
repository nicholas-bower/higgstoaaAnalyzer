// -*- C++ -*-
//
// Package:    higgstoaaAnalyzer/higgstoaaAnalyzer
// Class:      higgstoaaAnalyzer
//
/**\class higgstoaaAnalyzer higgstoaaAnalyzer.cc higgstoaaAnalyzer/higgstoaaAnalyzer/plugins/higgstoaaAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nicholas Bower
//         Created:  Tue, 21 Sep 2021 19:14:44 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "TLorentzVector.h"

#include "higgstoaaAnalyzer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <string>
#include <vector>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



class higgstoaaAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit higgstoaaAnalyzer(const edm::ParameterSet&);
      ~higgstoaaAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      const edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticles_;
      const edm::EDGetTokenT< std::vector<reco::GsfElectron> > electrons_;
      const edm::EDGetTokenT<std::vector<reco::Muon>> muons_;
      const edm::EDGetTokenT<std::vector<reco::PFJet>> jets_;
      const edm::EDGetTokenT<std::vector<reco::PFMET>> mets_;
      const edm::EDGetTokenT<std::vector<reco::PFTau>> taus_;
      const edm::EDGetTokenT<double> eventrho_;
      const edm::EDGetTokenT<GenEventInfoProduct> genInfoProduct_;
      const edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      
      std::vector<std::string> DoubleElectron {"HLT_DoubleEle33_CaloIdL_MW_v",
         "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"
         };
      std::vector<std::string> DoubleMuon {"HLT_DoubleMu20_7_Mass0to30_Photon23_v","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",
         "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v"
      };
      std::vector<std::string> MuEG {"HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v",
      "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v","HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v",
      "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v"
      };
      std::vector<std::string> singleEle {"HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v",
      "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1_v","HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_CrossL1_v",
      "HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1_v","HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v",
      "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1_v",
      "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v","HLT_Ele32_WPTight_Gsf_L1DoubleEG_v",
      "HLT_Ele35_WPTight_Gsf_v","HLT_Ele38_WPTight_Gsf_v",
      "HLT_Ele40_WPTight_Gsf_v"
      };

      std::vector<std::string> singleMu {"HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v",
         "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v","HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v",
         "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v","HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v",
         "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v","HLT_IsoMu27_v",
         "HLT_Mu50_v"
      };
      std::vector<std::string> SingleTau {"HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v","HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v",
      "HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v","HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v",
      "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v","HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v",
      "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v","HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v"
      };


      TH1F *h_GenParticleID;
      TH1F *h_ElectronPt;
      TH1F *h_GenElectronPt;
      TH1F *h_GenBJetPt;
      TH1F *h_JetPt;
      TH1F *h_TauPt;
      TH1F *h_GenTauPt;
      TH1F *h_GenMuonPt;
      TH1F *h_MuonPt;
      TH1F *h_MET; 
      TH1F *h_nElectrons;
      TH1F *h_nJets;
      TH1F *h_nMuons;
      TH1F *h_nTaus;
      TH1F *h_nbbee;
      TH1F *h_nbbmue;
      TH1F *h_nbbmumu;

      TH1F *h_ehad_tau;
      TH1F *h_ehad_singleE;
      TH2F *h_ehad_singleETau;

      TH1F *h_muhad_tau;
      TH1F *h_muhad_singleMu;
      TH2F *h_muhad_singleMuTau;

      TH1F *h_ee_DiEle;
      TH1F *h_ee_singleE;

      TH1F *h_eMu_singleE;
      TH1F *h_eMu_singleMu;
      TH2F *h_eMu_singleESingleMu;
      TH1F *h_eMu_MuEG;

      TH1F *h_MuMu_singleMu;
      TH1F *h_MuMu_diMu;

      TH1F *h_nbbehad;
      TH1F *h_nbbmuhad;
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
higgstoaaAnalyzer::higgstoaaAnalyzer(const edm::ParameterSet& iConfig) :
  genParticles_{consumes<std::vector<reco::GenParticle> >(edm::InputTag(std::string("genParticles")))},
  electrons_{consumes<std::vector<reco::GsfElectron>> (edm::InputTag(std::string("gedGsfElectrons")))},
  muons_{consumes<std::vector<reco::Muon>> (edm::InputTag(string("muons")))},
  jets_{consumes<std::vector<reco::PFJet>  >(edm::InputTag(std::string("ak4PFJets")))},
  mets_{consumes<std::vector<reco::PFMET>> (edm::InputTag(std::string("pfMet")))},
  taus_{consumes<std::vector<reco::PFTau>> (edm::InputTag(std::string("hpsPFTauProducer")))},
  eventrho_{consumes<double>(edm::InputTag(std::string("fixedGridRhoFastjetAll")))},
  genInfoProduct_{consumes<GenEventInfoProduct>(edm::InputTag(std::string("generator")))},
  triggerBits_{consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"))}
  

{
   //now do what ever initialization is needed
      edm::Service<TFileService> fs;

      h_GenParticleID=fs->make<TH1F>("h_GenParticleID", "Gen PDG ID", 37,0,37);
      h_GenParticleID->Sumw2();
      h_ElectronPt=fs->make<TH1F>("h_ElectronPt", "Reco Electron PT (GeV)", 40,0,160);
      h_ElectronPt->Sumw2();
      h_GenElectronPt=fs->make<TH1F>("h_GenElectronPt", "Gen Electron PT (GeV)", 40,0,160);
      h_GenElectronPt->Sumw2();
      h_GenBJetPt=fs->make<TH1F>("h_GenBJetPt", "Gen BJet PT (GeV)", 40,0,200);
      h_GenBJetPt->Sumw2();
      h_JetPt=fs->make<TH1F>("h_JetPt", "Reco Jet PT (GeV)", 40,0,200);
      h_JetPt->Sumw2();
      h_TauPt=fs->make<TH1F>("h_TauPt", "Reco TauHad PT (GeV)", 50,0,250);
      h_TauPt->Sumw2();
      h_GenTauPt=fs->make<TH1F>("h_GenTauPt", "Gen Tau PT (GeV)", 50,0,250);
      h_GenTauPt->Sumw2();
      h_GenMuonPt=fs->make<TH1F>("h_GenMuonPt", "Gen Muon PT (GeV)", 25,0,125);
      h_GenMuonPt->Sumw2();
      h_MuonPt=fs->make<TH1F>("h_MuonPt", "Reco Muon PT (GeV)", 25,0,125);
      h_MuonPt->Sumw2();
      h_MET=fs->make<TH1F>("h_MET", "MET (GeV)", 25,0,250);
      h_MET->Sumw2();
      h_nElectrons=fs->make<TH1F>("h_nElectrons", "Electron Cuts", 4,0,4);
      h_nElectrons->Sumw2();
      h_nMuons=fs->make<TH1F>("h_nMuons", "Muon Cuts", 3,0,3);
      h_nMuons->Sumw2();
      h_nTaus=fs->make<TH1F>("h_nTaus", "Tau Cuts", 2,0,2);
      h_nTaus->Sumw2();
      h_nJets=fs->make<TH1F>("h_nJets", "Jet Cuts", 9,0,9);
      h_nJets->Sumw2();

///////BB EHAD PLOTS
      h_ehad_tau=fs->make<TH1F>("h_ehad_tau", "EHAD Tau Triggers",3,0,3);
      h_ehad_tau->SetCanExtend(TH1::kAllAxes);
      h_ehad_tau->Sumw2();

      h_ehad_singleE=fs->make<TH1F>("h_ehad_singleE", "EHAD Ele Triggers",3,0,3);
      h_ehad_singleE->SetCanExtend(TH1::kAllAxes);
      h_ehad_singleE->Sumw2();

      h_ehad_singleETau=fs->make<TH2F>("h_ehad_singleETau", "EHad Single Ele single+Tau Triggers",3,0,3,3,0,3);
      h_ehad_singleETau->SetCanExtend(TH1::kAllAxes);
      h_ehad_singleETau->Sumw2();

/////// BBMuHAD PLOTS

      h_muhad_tau=fs->make<TH1F>("h_ehad_tau", "EHAD Tau Triggers",3,0,3);
      h_muhad_tau->SetCanExtend(TH1::kAllAxes);
      h_muhad_tau->Sumw2();
      
      h_muhad_singleMu=fs->make<TH1F>("h_muhad_singleMu", "MuHAD Single Mu Triggers",3,0,3);
      h_muhad_singleMu->SetCanExtend(TH1::kAllAxes);
      h_muhad_singleMu->Sumw2();

      h_muhad_singleMuTau=fs->make<TH2F>("h_muhad_singleMuTau", "EMu Single Mu single+Tau Triggers",3,0,3,3,0,3);
      h_muhad_singleMuTau->SetCanExtend(TH1::kAllAxes);
      h_muhad_singleMuTau->Sumw2();
/////// BBEMu PLOTS

      h_eMu_singleMu=fs->make<TH1F>("h_eMu_singleMu", "EMu Single Mu Triggers",3,0,3);
      h_eMu_singleMu->SetCanExtend(TH1::kAllAxes);
      h_eMu_singleMu->Sumw2();
      
      h_eMu_singleE=fs->make<TH1F>("h_eMu_singleE", "EMu Single Ele Triggers",3,0,3);
      h_eMu_singleE->SetCanExtend(TH1::kAllAxes);
      h_eMu_singleE->Sumw2();
      
      h_eMu_singleESingleMu=fs->make<TH2F>("h_eMu_singleESingleMu", "EMu Single Ele+single Mu Triggers",3,0,3,3,0,3);
      h_eMu_singleESingleMu->SetCanExtend(TH1::kAllAxes);
      h_eMu_singleESingleMu->Sumw2();
      
      h_eMu_MuEG=fs->make<TH1F>("h_eMu_MuEG", "EMu MuonEG Triggers",3,0,3);
      h_eMu_MuEG->SetCanExtend(TH1::kAllAxes);
      h_eMu_MuEG->Sumw2();
/////// BBee PLOTS

      h_ee_DiEle=fs->make<TH1F>("h_ee_DiEle", "ee DiEle Triggers",3,0,3);
      h_ee_DiEle->SetCanExtend(TH1::kAllAxes);
      h_ee_DiEle->Sumw2();

      h_ee_singleE=fs->make<TH1F>("h_ee_singleE", "ee Single Ele Triggers",3,0,3);
      h_ee_singleE->SetCanExtend(TH1::kAllAxes);
      h_ee_singleE->Sumw2();
/////// BBMumu PLOTS

      h_MuMu_diMu=fs->make<TH1F>("h_MuMu_diMu", "MuMu Double Mu Triggers",3,0,3);
      h_MuMu_diMu->SetCanExtend(TH1::kAllAxes);
      h_MuMu_diMu->Sumw2();

      h_MuMu_singleMu=fs->make<TH1F>("h_MuMu_singleMu", "MuMu Single Mu Triggers",3,0,3);
      h_MuMu_singleMu->SetCanExtend(TH1::kAllAxes);
      h_MuMu_singleMu->Sumw2();

///////Selection Plots   

      h_nbbee=fs->make<TH1F>("h_nbbee", "NEvent ee",3,0,3);
      h_nbbee->SetCanExtend(TH1::kAllAxes);
      h_nbbee->Sumw2();

      h_nbbmue=fs->make<TH1F>("h_nbbmue", "NEvent eMu",3,0,3);
      h_nbbmue->SetCanExtend(TH1::kAllAxes);
      h_nbbmue->Sumw2();

      h_nbbmumu=fs->make<TH1F>("h_nbbmumu", "N Mu Mu Events",3,0,3);
      h_nbbmumu->SetCanExtend(TH1::kAllAxes);
      h_nbbmumu->Sumw2();

      h_nbbehad=fs->make<TH1F>("h_nbbehad", "n eHad Events",3,0,3);
      h_nbbehad->SetCanExtend(TH1::kAllAxes);
      h_nbbehad->Sumw2();

      h_nbbmuhad=fs->make<TH1F>("h_nbbmuhad", "N Mu Had Events",3,0,3);
      h_nbbmuhad->SetCanExtend(TH1::kAllAxes);
      h_nbbmuhad->Sumw2();





}


higgstoaaAnalyzer::~higgstoaaAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
higgstoaaAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;


   std::vector<reco::PFTau> taus = h2AA::getObject<std::vector<reco::PFTau>>(taus_, iEvent);
   std::vector<reco::PFMET> mets = h2AA::getObject<std::vector<reco::PFMET>>(mets_, iEvent);
   std::vector<reco::GenParticle> genParticles = h2AA::getObject<std::vector<reco::GenParticle>>(genParticles_, iEvent);
   std::vector<reco::Muon> muons = h2AA::getObject<std::vector<reco::Muon>>(muons_, iEvent);
   std::vector<reco::PFJet> jets = h2AA::getObject<std::vector<reco::PFJet>>(jets_, iEvent);
   std::vector<reco::GsfElectron>electrons = h2AA::getObject<std::vector<reco::GsfElectron>>(electrons_, iEvent);
///Get Trigger objects for out handles
   edm::InputTag trigResultsTag("TriggerResults","","HLT"); //make sure have correct process on MC
   edm::Handle<edm::TriggerResults>  triggerBits;
   iEvent.getByLabel(trigResultsTag,triggerBits);
   
   const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerBits);   

   double rho = h2AA::getObject<double>(eventrho_,iEvent);
   edm::Handle<GenEventInfoProduct>  genInfoHandle;
   try{iEvent.getByToken(genInfoProduct_, genInfoHandle);}
   catch(...){;}
   double genWeight=genInfoHandle->weight();

   h_MET->Fill(mets[0].pt(),genWeight);
   //////Gen Particle Sorting
   std::vector<reco::GenParticle*> genElectrons;
   std::vector<reco::GenParticle*> genMuons;
   std::vector<reco::GenParticle*> genTaus;
   reco::GenParticle *genNuTau=NULL; //////These are just the tau neutrinos We know there is a maximum of 2 per event either one anti, one normal, both or none.
   reco::GenParticle *genNuTauBar=NULL; 
   for(auto& gen : genParticles){
      h_GenParticleID->Fill(abs(gen.pdgId())+.5,genWeight);
      if(abs(gen.pdgId())==11 &&
         gen.isDirectHardProcessTauDecayProductFinalState()) {
         h_GenElectronPt->Fill(gen.pt(),genWeight);
         genElectrons.push_back(&gen);

      }
      if (abs(gen.pdgId())==13 && gen.isDirectHardProcessTauDecayProductFinalState()){
         h_GenMuonPt->Fill(gen.pt(),genWeight);
         genMuons.push_back(&gen);
      }
      if (abs(gen.pdgId())==15 && gen.isHardProcess()){
         h_GenTauPt->Fill(gen.pt(),genWeight);
         genTaus.push_back(&gen);
      }
      if (gen.pdgId()==5){

         h_GenBJetPt->Fill(gen.pt(),genWeight);
      }
      
      if (gen.pdgId()==16 && gen.isDirectHardProcessTauDecayProductFinalState()){
         genNuTau = &gen;
      }
      if (gen.pdgId()==-16 && gen.isDirectHardProcessTauDecayProductFinalState()){
         genNuTauBar = &gen;
      } 
   }
   //Determin which is the neutrino associated with 
   reco::GenParticle *genTauHad=NULL;
   reco::GenParticle *genNuTauHad=NULL;
   bool DiTauHad=false;


   for (auto& gen : genTaus){
      if ((genElectrons.size()==1&& genMuons.size()==0)){
         if ((*genElectrons[0]).charge()<0&&(*gen).charge()>0){
               //std::cout<<"e-\n";

         genTauHad = gen;
         genNuTauHad = genNuTauBar;
         }
         if ((*genElectrons[0]).charge()>0&&(*gen).charge()<0){
                     //std::cout<<"e+\n";
         genTauHad = gen;
         genNuTauHad = genNuTau;
         }
      }
      else if((genElectrons.size()==0&& genMuons.size()==1)){
         if ((*genMuons[0]).charge()<0&&(*gen).charge()>0){
         //std::cout<<"m-\n";
         genTauHad = gen;
         genNuTauHad = genNuTauBar;
         }
         if ((*genMuons[0]).charge()>0&&(*gen).charge()<0){
         //std::cout<<"m+\n";
         genTauHad= gen;
         genNuTauHad = genNuTau;
         }

      }
      else{
         DiTauHad=true;
      }
   }
   cout<<"Gen Sorting";
   std::vector<reco::GsfElectron> selected_electrons;
   for(auto& ele : electrons){
      h_nElectrons->Fill(0.5,genWeight);
      if (ele.pt()>7 && abs(ele.eta())<2.5){
         h_nElectrons->Fill(1.5,genWeight);
         float E_c=ele.superCluster()->energy();
         if(h2AA::checkID(ele,1,rho, E_c)){
            h_nElectrons->Fill(2.5,genWeight);
            h_ElectronPt->Fill(ele.pt(),genWeight);
            selected_electrons.push_back(ele);
            bool matched= false;
            for(unsigned int iGenE=0; iGenE<genElectrons.size(); iGenE++){
               reco::GenParticle gen = *genElectrons[iGenE];
               TLorentzVector e;
               TLorentzVector genE;
               e.SetPtEtaPhiM(ele.pt(), ele.eta(), ele.phi(), ele.mass());
               genE.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), gen.mass());
               float dr = genE.DeltaR(e);
               if(dr <.1&& ele.pt()>1){
                  matched = true;
               }
            if (matched){
               h_ElectronPt->Fill(ele.pt(),genWeight);
               h_nElectrons->Fill(3.5,genWeight);
               selected_electrons.push_back(ele);}
            
            }
         }
      }
   }
   cout<<"Electron selectrion";
   std::vector<reco::PFJet> selected_jets;
   for (unsigned int iJet=0;iJet<jets.size(); iJet++){
      reco::PFJet jet = jets[iJet];
      h_nJets->Fill(.5,genWeight);
      if (jet.pt()>10 && abs(jet.eta())<2.5 ){
         h_nJets->Fill(1.5,genWeight);

         float NHF  = jet.neutralHadronEnergyFraction();
         float NEMF = jet.neutralEmEnergyFraction();
         float CHF  = jet.chargedHadronEnergyFraction();
         float MUF  = jet.muonEnergyFraction();
         float CEMF = jet.chargedEmEnergyFraction();
         float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
         float CHM = jet.chargedMultiplicity();
         if (CEMF>0.8){continue;}
         h_nJets->Fill(2.5,genWeight);
         if (CHM<0){continue;}
         h_nJets->Fill(3.5,genWeight);
         if (CHF <0){continue;}
         h_nJets->Fill(4.5,genWeight);

         if (NumConst<1){continue;}
         h_nJets->Fill(5.5,genWeight);

         if (NEMF>0.9){continue;}
         h_nJets->Fill(6.5,genWeight);
         if (MUF>0.8){continue;}
         h_nJets->Fill(7.5,genWeight);
         if (NHF>.9){continue;}
         h_nJets->Fill(8.5,genWeight);
         h_JetPt->Fill(jet.pt(),genWeight);
         selected_jets.push_back(jet);
      }
   }
   cout<<"jet selection";
   std::vector<reco::Muon> selected_muons;
   for (unsigned int iMu=0; iMu<muons.size(); ++iMu){
      reco::Muon muon = muons[iMu];
      h_nMuons->Fill(.5,genWeight);
      if (muon::isLooseMuon(muon)==false){continue;}////////////Loose Muon ID
      h_nMuons->Fill(1.5,genWeight);
      if((muon.pt())<3|| muon.eta()>2.4){continue;} 
      h_nMuons->Fill(2.5,genWeight);
      bool matched=false;
      for (unsigned int iGen; iGen<genMuons.size();iGen++){
         reco::GenParticle gen = *genMuons[iGen];
         TLorentzVector mu;
         TLorentzVector genMu;
         mu.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), muon.mass());
         genMu.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), gen.mass());
         float dr = genMu.DeltaR(mu);
         if(dr <.1){
            matched = true;
         }
      }
      if (!matched){continue;}
      h_MuonPt->Fill(muon.pt(),genWeight);
      selected_muons.push_back(muon);
   }
   cout<<"Mu selection";
   std::vector<reco::PFTau> selected_taus;
   for (unsigned int iTau = 0; iTau<taus.size();iTau++){
      reco::PFTau tau = taus[iTau];
      h_nTaus->Fill(.5, genWeight);
      if (tau.pt()<10 || abs(tau.eta())>2.4){continue;}
      //if (tau.tauID("decayModeFinding")==false){continue;}
      h_nTaus->Fill(1.5, genWeight);
      if(genTauHad==nullptr||genNuTauHad==nullptr){continue;}
      if(DiTauHad==true){continue;}
      reco::GenParticle thad =*genTauHad;
      reco::GenParticle nt = *genNuTauHad;
      TLorentzVector t;
      TLorentzVector genT;
      TLorentzVector genNu;
      t.SetPtEtaPhiM(tau.pt(), tau.eta(), tau.phi(), tau.mass());
      genT.SetPtEtaPhiM(thad.pt(), thad.eta(), thad.phi(), thad.mass());
      genNu.SetPtEtaPhiM(nt.pt(), nt.eta(), nt.phi(), nt.mass());
      if(t.DeltaR(genT-genNu)>.1){continue;}
      h_nTaus->Fill(2.5, genWeight);
      h_TauPt->Fill(tau.pt(),genWeight);
      selected_taus.push_back(tau);
   }
   
   cout<<"tau selection";
   ///////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////
   //////////////////////EVENT SELECTION///////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////
   sort(selected_electrons.begin(),selected_electrons.end(), h2AA::sortByPt<reco::GsfElectron>);
   sort(selected_muons.begin(),selected_muons.end(), h2AA::sortByPt<reco::Muon>);
   sort(selected_taus.begin(),selected_taus.end(), h2AA::sortByPt<reco::PFTau>);
   sort(selected_jets.begin(),selected_jets.end(), h2AA::sortByPt<reco::PFJet>);
   cout<<"Sort";
   h_nbbmue->Fill(0.5);
   if (selected_jets.size()>1){
      h_nbbmue->Fill(1.5);
      if(selected_muons.size()>0){
         h_nbbmue->Fill(2.5);
         if(selected_electrons.size()>0){
            h_nbbmue->Fill(3.5);
            for (unsigned int iTrigger = 0; iTrigger<MuEG.size();iTrigger++){
               std::string Trigger = MuEG[iTrigger];
               bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
               if (passTrig){h_eMu_MuEG->Fill(Trigger.c_str(), genWeight);}
            }
            for (unsigned int iTrigger = 0; iTrigger<singleEle.size();iTrigger++){
               std::string Trigger = singleEle[iTrigger];
               bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
               if (passTrig){h_eMu_singleE->Fill(Trigger.c_str(), genWeight);}
            }
            for (unsigned int iTrigger = 0; iTrigger<singleMu.size();iTrigger++){
               std::string Trigger = singleMu[iTrigger];
               bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
               if (passTrig){h_eMu_singleE->Fill(Trigger.c_str(), genWeight);}
            }
            for (unsigned int iTrigger1 = 0; iTrigger1<singleMu.size();iTrigger1++){
               std::string Trigger1 = singleMu[iTrigger1];
               for (unsigned int iTrigger2 = 0; iTrigger2<singleEle.size();iTrigger2++){
                  std::string Trigger2 = singleEle[iTrigger2];
                  bool passTrig1=triggerBits->accept(trigNames.triggerIndex(Trigger1));
                  bool passTrig2=triggerBits->accept(trigNames.triggerIndex(Trigger2));
                  if(passTrig1&&passTrig2){h_eMu_singleESingleMu->Fill(Trigger1.c_str(),Trigger2.c_str(), genWeight);}
               }
            }
         }
      }
   }
   cout<<"bbmue";
   h_nbbee->Fill(0.5);
   if (selected_jets.size()>1){
      h_nbbee->Fill(1.5);
      if(selected_electrons.size()>1){
         h_nbbee->Fill(2.5);
         for (unsigned int iTrigger = 0; iTrigger<DoubleElectron.size();iTrigger++){
            std::string Trigger = DoubleElectron[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_ee_DiEle->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<singleEle.size();iTrigger++){
            std::string Trigger = singleEle[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_ee_singleE->Fill(Trigger.c_str(), genWeight);}
         } 
      }
   }
   cout<<"bbee";
   h_nbbmumu->Fill(0.5);
   if (selected_jets.size()>1){
      h_nbbmumu->Fill(1.5);
      if(selected_muons.size()>1){
         h_nbbmumu->Fill(2.5);
         for (unsigned int iTrigger = 0; iTrigger<DoubleMuon.size();iTrigger++){
            std::string Trigger = DoubleMuon[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_MuMu_diMu->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<singleMu.size();iTrigger++){
            std::string Trigger = singleMu[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_MuMu_singleMu->Fill(Trigger.c_str(), genWeight);}
         }
      }
   }
   cout<<"bbmumu";
   h_nbbehad->Fill(0.5);
   if (selected_jets.size()>1){
      h_nbbehad->Fill(1.5);
      cout<<"|Event Test1|";
      if(selected_taus.size()>0){
         h_nbbehad->Fill(2.5);
         cout<<"|Event Test2|";
         if(selected_electrons.size()>0){
            h_nbbehad->Fill(3.5);
            cout<<"|Event Test3|";
            for (unsigned int iTrigger = 0; iTrigger<singleEle.size();iTrigger++){
               cout<<"|ETrigger|";
               std::string Trigger = singleEle[iTrigger];
               cout<<Trigger;
               bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger.c_str()));
               cout << "PASTRIGGERTEST";
               if (passTrig){h_ehad_singleE->Fill(Trigger.c_str(), genWeight);}
            }
            for (unsigned int iTrigger = 0; iTrigger<SingleTau.size();iTrigger++){
               std::string Trigger = SingleTau[iTrigger];
               bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
               if (passTrig){cout<<"|PassTauTrigger|";h_ehad_tau->Fill(Trigger.c_str(), genWeight);}
            }
            for (unsigned int iTrigger1 = 0; iTrigger1<SingleTau.size();iTrigger1++){
               std::string Trigger1 = SingleTau[iTrigger1];
               for (unsigned int iTrigger2 = 0; iTrigger2<singleEle.size();iTrigger2++){
                  std::string Trigger2 = singleEle[iTrigger2];
                  bool passTrig1=triggerBits->accept(trigNames.triggerIndex(Trigger1));
                  bool passTrig2=triggerBits->accept(trigNames.triggerIndex(Trigger2));
                  if(passTrig1&&passTrig2){
                     cout<<"|DualTrigger|";
                     const char *T1 = Trigger1.c_str();
                     const char *T2 = Trigger2.c_str();
                     h_ehad_singleETau->Fill(T1,T2, genWeight);
                  }
               }
            }
         }
      }
   }
   cout<<"bbehad";
   h_nbbmuhad->Fill(0.5);
   if (selected_jets.size()>1){
      h_nbbmuhad->Fill(1.5);
      if(selected_taus.size()>0){
         h_nbbmuhad->Fill(2.5);
         if(selected_electrons.size()>0){
            h_nbbmuhad->Fill(3.5);

            for (unsigned int iTrigger = 0; iTrigger<singleMu.size();iTrigger++){
               std::string Trigger = singleMu[iTrigger];
               bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
               if (passTrig){
                  const char *T = Trigger.c_str();
                  h_muhad_singleMu->Fill(T, genWeight);
               }
            }
            for (unsigned int iTrigger = 0; iTrigger<SingleTau.size();iTrigger++){
               std::string Trigger = SingleTau[iTrigger];
               bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
               if (passTrig){
                  const char *T = Trigger.c_str();
                  h_muhad_tau->Fill(T, genWeight);
               }
            }
            for (unsigned int iTrigger1 = 0; iTrigger1<SingleTau.size();iTrigger1++){
               std::string Trigger1 = SingleTau[iTrigger1];
               for (unsigned int iTrigger2 = 0; iTrigger2<singleMu.size();iTrigger2++){
                  std::string Trigger2 = singleMu[iTrigger2];
                  bool passTrig1=triggerBits->accept(trigNames.triggerIndex(Trigger1));
                  bool passTrig2=triggerBits->accept(trigNames.triggerIndex(Trigger2));
                  if(passTrig1&&passTrig2){
                     const char *T1 = Trigger1.c_str();
                     const char *T2 = Trigger2.c_str();
                     h_muhad_singleMuTau->Fill(T1,T2, genWeight);
                  }
               }
            }
         }
      }
   }

}
// ------------ method called once each job just before starting event loop  ------------
void
higgstoaaAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
higgstoaaAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
higgstoaaAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

}

//define this as a plug-in
DEFINE_FWK_MODULE(higgstoaaAnalyzer);
