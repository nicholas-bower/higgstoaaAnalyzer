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
#include "TLorentzVector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
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
      const edm::EDGetTokenT<GenEventInfoProduct> genInfoProduct_;
      const edm::EDGetTokenT< std::vector<pat::Electron> > electrons_;
      const edm::EDGetTokenT< std::vector<pat::Muon> > muons_;
      const edm::EDGetTokenT<double> eventrho_;
      const edm::EDGetTokenT<std::vector<reco::Vertex>> vertex_;
      const edm::EDGetTokenT<std::vector<pat::Jet>> jets_;
      const edm::EDGetTokenT<std::vector<reco::GenJet>> genJets_;
      const edm::EDGetTokenT<std::vector<pat::MET>> mets_;
      const edm::EDGetTokenT<std::vector<reco::GenMET>>  genMets_;
      const edm::EDGetTokenT<std::vector<pat::Tau>> taus_;
      const edm::EDGetTokenT<reco::ConversionCollection> convs_;
      const edm::EDGetTokenT<reco::BeamSpot> thebs_;
      const edm::EDGetTokenT< std::vector<pat::Electron> > std_electrons_;
      const edm::EDGetTokenT< std::vector<pat::Electron> > lowPt_electrons_;
      const edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      const edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavorMatching_;

      
      std::vector<std::string> DoubleElectron {"HLT_DoubleEle33_CaloIdL_MW_v15","HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v17",
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v17","HLT_ECALHT800_v9"
         };
      std::vector<std::string> DoubleMuon {"HLT_DoubleMu20_7_Mass0to30_Photon23_v5","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v3",
      "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v12"
      };
      std::vector<std::string> MuEG {"HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v14","HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v12",
      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12"
      };
      std::vector<std::string> singleEle {"HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v9",
      "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_TightID_CrossL1_v9","HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTau30_eta2p1_TightID_CrossL1_v9",
      "HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_CrossL1_v9","HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTau30_eta2p1_TightID_CrossL1_v9",
      "HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v9","HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v9",
      "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v7","HLT_Ele35_WPTight_Gsf_v7","HLT_Ele38_WPTight_Gsf_v7","HLT_Ele40_WPTight_Gsf_v7"

      };

      std::vector<std::string> singleMu {"HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v8","HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_TightID_CrossL1_v8",
         "HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_CrossL1_v8","HLT_IsoMu20_eta2p1_MediumChargedIsoPFTau27_eta2p1_TightID_CrossL1_v8",
         "HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_CrossL1_v8","HLT_IsoMu20_eta2p1_TightChargedIsoPFTau27_eta2p1_TightID_CrossL1_v8",
         "HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v8","HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8",
         "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v8","HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8",
         "HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1_v8","HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1_v8",
         "HLT_IsoMu27_v13","HLT_Mu50_v11"
         };
      std::vector<std::string> SingleTau {"HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v8","HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v8",
      "HLT_DoubleTightChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v8","HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v8",
      "HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v8","HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v8",
      "HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v8","HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v8"
      };


      TH1F *h_GenParticleID;
      TH1F *h_std_TauElectronPt;
      TH1F *h_GenTauElectronPt;
      TH1F *h_GenbPt;
      TH1F *h_JetPt;
      TH1F *h_JetPt_m;
      TH1F *h_JetPt_f;
      TH1F *h_TauPt;
      TH1F *h_GenTauPt;
      TH1F *h_GenTauMuonPt;
      TH1F *h_TauMuonPt;
      
      TH1F *h_TauMuonIsolation;
      TH1F *h_MET; 
      TH1F *h_nStdElectrons;
      TH1F *h_nLowPtElectrons;

      TH1F *h_nJets;
      TH1F *h_bDiscriminator;
      TH1F *h_bDiscriminator_m;
      TH1F *h_bDiscriminator_f;
      TH1F *h_nMuons;
      TH1F *h_nTaus;
      TH1F *h_DoubleElectron;


      TH1F *h_singleE;
      TH1F *h_singleMu;
      
      TH1F *h_DoubleMuon;
      TH1F *h_MuonEGamma;

      TH1F *h_TauElectronID;
      TH1F *h_TauElectronIDandIso;
      TH1F *h_TauElectronID_FULL;
//B Daughter studies
      TH1F *h_nEStatus1_NotfromTau_Gen;
      TH1F *h_nMuStatus1_NotfromTau_Gen;
      TH1F *h_nEMotherB;
      TH1F *h_nMuMotherB;
      TH1F *h_genS1Muon_NonPrimary_Pt;
      TH1F *h_genbMu_Pt;
      TH1F *h_genbE_Pt;

//BDaughter Reco
      TH1F *h_LowPtElectron_bmatched_Pt;
      TH1F *h_LowPtElectron_bmatched_ID;

      TH1F *h_LowPtElectron_Taumatched_Pt;
      TH1F *h_LowPtElectron_Taumatched_ID;


      TH1F *h_StandardElectron_bmatched_Pt;

      TH1F *h_Muon_bmatched_Pt;
      TH1F *h_Muon_bmatched_Iso;
////Event Selection
      TH1F *h_nEvent_tauE_tauE;
      TH1F *h_nEvent_tauMu_tauMu;
      TH1F *h_nEvent_tauE_tauMu;

      TH1F *h_tauE_tauE_leadtauEpt;
      TH1F *h_tauE_tauE_subtauEpt;

      TH1F *h_tauMu_tauE_tauEpt;
      TH1F *h_tauMu_tauE_tauMupt;


      TH1F *h_tauMu_tauMu_leadtauMupt;
      TH1F *h_tauMu_tauMu_subtauMupt;


////Trigger plots
      TH1F *h_Trigger_tauE_tauE_DoubleElectron;
      TH1F *h_Trigger_tauE_tauE_singleEle;
      TH1F *h_Trigger_tauE_tauE_MuEG;

      TH1F *h_tauE_TauE_Mvis_EG;
      TH1F *h_tauE_TauE_Mvis_EGorDoubleEG;
      TH1F *h_tauE_TauE_Mvis_DoubleEG;
      TH1F *h_tauE_TauE_Mvis_noTrig;

      TH1F *h_tauE_TauE_bMass_EG;
      TH1F *h_tauE_TauE_bMass_EGorDoubleEG;
      TH1F *h_tauE_TauE_bMass_DoubleEG;
      TH1F *h_tauE_TauE_bMass_noTrig;

      TH1F *h_tauE_TauE_bDRe_EG;
      TH1F *h_tauE_TauE_bDRe_EGorDoubleEG;
      TH1F *h_tauE_TauE_bDRe_DoubleEG;
      TH1F *h_tauE_TauE_bDRe_noTrig;

      TH1F *h_tauE_TauE_DiTauPt_EG;
      TH1F *h_tauE_TauE_DiTauPt_EGorDoubleEG;
      TH1F *h_tauE_TauE_DiTauPt_DoubleEG;
      TH1F *h_tauE_TauE_DiTauPt_noTrig;

      TH1F *h_tauE_TauE_eDRe_EG;
      TH1F *h_tauE_TauE_eDRe_EGorDoubleEG;
      TH1F *h_tauE_TauE_eDRe_DoubleEG;
      TH1F *h_tauE_TauE_eDRe_noTrig;

      TH1F *h_tauE_TauE_MET_EG;
      TH1F *h_tauE_TauE_MET_EGorDoubleEG;
      TH1F *h_tauE_TauE_MET_DoubleEG;
      TH1F *h_tauE_TauE_MET_noTrig;

//Tau E Tau Mu
      TH1F *h_Trigger_tauE_tauMu_singleMu;
      TH1F *h_Trigger_tauE_tauMu_singleEle;
      TH1F *h_Trigger_tauE_tauMu_MuEG;
      TH1F *h_Trigger_tauE_tauMu_DoubleElectron;
      TH1F *h_Trigger_tauE_tauMu_DoubleMuon;

      TH1F *h_tauMu_TauE_Mvis_SingleMu;
      TH1F *h_tauMu_TauE_Mvis_EGorSingleMu;
      TH1F *h_tauMu_TauE_Mvis_MuEGorSingleMu;
      TH1F *h_tauMu_TauE_Mvis_MuEGorSingleMuorSingleE;
      TH1F *h_tauMu_TauE_Mvis_EG;
      TH1F *h_tauMu_TauE_Mvis_MuEG;
      TH1F *h_tauMu_TauE_Mvis_noTrig;

      TH1F *h_tauMu_TauE_bMass_SingleMu;
      TH1F *h_tauMu_TauE_bMass_EGorSingleMu;
      TH1F *h_tauMu_TauE_bMass_MuEGorSingleMu;
      TH1F *h_tauMu_TauE_bMass_MuEGorSingleMuorSingleE;
      TH1F *h_tauMu_TauE_bMass_EG;
      TH1F *h_tauMu_TauE_bMass_MuEG;
      TH1F *h_tauMu_TauE_bMass_noTrig;

      TH1F *h_tauMu_TauE_bDRmu_SingleMu;
      TH1F *h_tauMu_TauE_bDRmu_EGorSingleMu;
      TH1F *h_tauMu_TauE_bDRmu_MuEGorSingleMu;
      TH1F *h_tauMu_TauE_bDRmu_MuEGorSingleMuorSingleE;
      TH1F *h_tauMu_TauE_bDRmu_EG;
      TH1F *h_tauMu_TauE_bDRmu_MuEG;
      TH1F *h_tauMu_TauE_bDRmu_noTrig;

      TH1F *h_tauMu_TauE_DiTauPt_SingleMu;
      TH1F *h_tauMu_TauE_DiTauPt_EGorSingleMu;
      TH1F *h_tauMu_TauE_DiTauPt_MuEGorSingleMu;
      TH1F *h_tauMu_TauE_DiTauPt_MuEGorSingleMuorSingleE;
      TH1F *h_tauMu_TauE_DiTauPt_EG;
      TH1F *h_tauMu_TauE_DiTauPt_MuEG;
      TH1F *h_tauMu_TauE_DiTauPt_noTrig;

      TH1F *h_tauMu_TauE_muDRe_SingleMu;
      TH1F *h_tauMu_TauE_muDRe_EGorSingleMu;
      TH1F *h_tauMu_TauE_muDRe_MuEGorSingleMu;
      TH1F *h_tauMu_TauE_muDRe_MuEGorSingleMuorSingleE;
      TH1F *h_tauMu_TauE_muDRe_EG;
      TH1F *h_tauMu_TauE_muDRe_MuEG;
      TH1F *h_tauMu_TauE_muDRe_noTrig;

      TH1F *h_tauMu_TauE_muIso_SingleMu;
      TH1F *h_tauMu_TauE_muIso_EGorSingleMu;
      TH1F *h_tauMu_TauE_muIso_MuEGorSingleMu;
      TH1F *h_tauMu_TauE_muIso_MuEGorSingleMuorSingleE;
      TH1F *h_tauMu_TauE_muIso_EG;
      TH1F *h_tauMu_TauE_muIso_MuEG;
      TH1F *h_tauMu_TauE_muIso_noTrig;

      TH1F *h_tauMu_TauE_MET_SingleMu;
      TH1F *h_tauMu_TauE_MET_EGorSingleMu;
      TH1F *h_tauMu_TauE_MET_MuEGorSingleMu;
      TH1F *h_tauMu_TauE_MET_MuEGorSingleMuorSingleE;
      TH1F *h_tauMu_TauE_MET_EG;
      TH1F *h_tauMu_TauE_MET_MuEG;
      TH1F *h_tauMu_TauE_MET_noTrig;
      ///TauMu Tau Mu

      TH1F *h_Trigger_tauMu_tauMu_DoubleMuon;
      TH1F *h_Trigger_tauMu_tauMu_singleMu;

      TH1F *h_tauMu_TauMu_Mvis_SingleMu;
      TH1F *h_tauMu_TauMu_Mvis_SingleMuorDoubleMu;
      TH1F *h_tauMu_TauMu_Mvis_DoubleMu;
      TH1F *h_tauMu_TauMu_Mvis_noTrig;

      TH1F *h_tauMu_TauMu_bMass_SingleMu;
      TH1F *h_tauMu_TauMu_bMass_SingleMuorDoubleMu;
      TH1F *h_tauMu_TauMu_bMass_DoubleMu;
      TH1F *h_tauMu_TauMu_bMass_noTrig;

      TH1F *h_tauMu_TauMu_bDRmu_SingleMu;
      TH1F *h_tauMu_TauMu_bDRmu_SingleMuorDoubleMu;
      TH1F *h_tauMu_TauMu_bDRmu_DoubleMu;
      TH1F *h_tauMu_TauMu_bDRmu_noTrig;

      TH1F *h_tauMu_TauMu_DiTauPt_SingleMu;
      TH1F *h_tauMu_TauMu_DiTauPt_SingleMuorDoubleMu;
      TH1F *h_tauMu_TauMu_DiTauPt_DoubleMu;
      TH1F *h_tauMu_TauMu_DiTauPt_noTrig;

      TH1F *h_tauMu_TauMu_muDRmu_SingleMu;
      TH1F *h_tauMu_TauMu_muDRmu_SingleMuorDoubleMu;
      TH1F *h_tauMu_TauMu_muDRmu_DoubleMu;
      TH1F *h_tauMu_TauMu_muDRmu_noTrig;

      TH1F *h_tauMu_TauMu_muIso_SingleMu;
      TH1F *h_tauMu_TauMu_muIso_SingleMuorDoubleMu;
      TH1F *h_tauMu_TauMu_muIso_DoubleMu;
      TH1F *h_tauMu_TauMu_muIso_noTrig;

      TH1F *h_tauMu_TauMu_MET_SingleMu;
      TH1F *h_tauMu_TauMu_MET_SingleMuorDoubleMu;
      TH1F *h_tauMu_TauMu_MET_DoubleMu;
      TH1F *h_tauMu_TauMu_MET_noTrig;

//Trigger improvement checks
      TH1F *h_Trigger_tauE_TauE_DiEandEG;
      TH1F *h_Trigger_tauE_TauE_EGandMuEG;

      TH1F *h_Trigger_tauMu_TauE_MuandEG;
      TH1F *h_Trigger_tauMu_TauE_EGandMuEG;

      TH1F *h_Trigger_tauMu_TauE_MuandMuEG;

      TH1F *h_Trigger_tauMu_TauMu_DiMuandMu;
      TH1F *h_Trigger_tauMu_TauMu_MuandMuEG;

      TH1F *h_Trigger_general_MuandDiMu;
      TH1F *h_Trigger_general_MuandMuEG;
      TH1F *h_Trigger_general_EGandMuEG;
      TH1F *h_Trigger_general_EGandDiEG;
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
   genParticles_{consumes<std::vector<reco::GenParticle> >(edm::InputTag(std::string("prunedGenParticles")))},
   genInfoProduct_{consumes<GenEventInfoProduct> (edm::InputTag(std::string("generator")))},
   electrons_{consumes< std::vector<pat::Electron> >(edm::InputTag(std::string("slimmedElectrons")))},  
   muons_{consumes< std::vector<pat::Muon> >(edm::InputTag(std::string("slimmedMuons")))},  
   eventrho_{consumes<double>(edm::InputTag(std::string("fixedGridRhoFastjetAll")))},
   vertex_{consumes<std::vector<reco::Vertex>  >(edm::InputTag(std::string("offlineSlimmedPrimaryVertices")))},
   jets_{consumes<std::vector<pat::Jet>  >(edm::InputTag(std::string("slimmedJets")))},
   genJets_{consumes<std::vector<reco::GenJet>  >(edm::InputTag(std::string("slimmedGenJets")))},
   mets_{consumes<std::vector<pat::MET>> (edm::InputTag(std::string("slimmedMETs")))},
   //taus_{consumes<std::vector<pat::Tau>> (edm::InputTag(std::string("slimmedTausUnCleaned")))},
   taus_{consumes<std::vector<pat::Tau>> (edm::InputTag(std::string("slimmedTaus")))},
   std_electrons_{consumes< std::vector<pat::Electron> >(edm::InputTag(std::string("slimmedElectrons")))},
   lowPt_electrons_{consumes< std::vector<pat::Electron> >(edm::InputTag(std::string("slimmedLowPtElectrons")))},
   triggerBits_{consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"))},
   jetFlavorMatching_{consumes<reco::JetFlavourInfoMatchingCollection>(edm::InputTag(std::string("slimmedGenJetsFlavourInfos")))}
{
   //now do what ever initialization is needed
      edm::Service<TFileService> fs;

      h_GenParticleID=fs->make<TH1F>("h_GenParticleID", "Gen PDG ID(Status=1)", 37,0,37);
      h_GenParticleID->Sumw2();
      h_std_TauElectronPt=fs->make<TH1F>("h_std_TauElectronPt", "Std Reco Tau Electron PT (GeV)", 45,0,90);
      h_std_TauElectronPt->Sumw2();
      h_GenTauElectronPt=fs->make<TH1F>("h_GenTauElectronPt", "Gen Tau Electron PT (GeV)", 45,0,90);
      h_GenTauElectronPt->Sumw2();
      h_GenbPt=fs->make<TH1F>("h_GenbPt", "Gen b PT (GeV)", 40,0,200);
      h_GenbPt->Sumw2();
      h_JetPt=fs->make<TH1F>("h_JetPt", "Reco Jet PT (GeV)", 40,0,200);
      h_JetPt->Sumw2();
      h_JetPt_m=fs->make<TH1F>("h_JetPt_m", "RecoMatched Jet PT (GeV)", 40,0,200);
      h_JetPt_m->Sumw2();      
      h_JetPt_f=fs->make<TH1F>("h_JetPt_f", "Reco Fake Jet PT (GeV)", 40,0,200);
      h_JetPt_f->Sumw2();
      h_TauPt=fs->make<TH1F>("h_TauPt", "Reco TauHad PT (GeV)", 50,0,250);
      h_TauPt->Sumw2();
      h_GenTauPt=fs->make<TH1F>("h_GenTauPt", "Gen Tau PT (GeV)", 50,0,250);
      h_GenTauPt->Sumw2();
      h_GenTauMuonPt=fs->make<TH1F>("h_GenTauMuonPt", "Gen Tau Muon PT (GeV)", 45,0,90);
      h_GenTauMuonPt->Sumw2();
      h_TauMuonPt=fs->make<TH1F>("h_TauMuonPt", "Reco Tau Muon PT (GeV)", 45,0,90);
      h_TauMuonPt->Sumw2();
      h_TauMuonIsolation=fs->make<TH1F>("h_TauMuonIsolation", "Reco Muon Isolation", 25,0,1);
      h_TauMuonIsolation->Sumw2();
      h_MET=fs->make<TH1F>("h_MET", "MET (GeV)", 25,0,250);
      h_MET->Sumw2();
      h_nStdElectrons=fs->make<TH1F>("h_nStdElectrons", "StdElectron Cuts", 5,0,5);
      h_nStdElectrons->Sumw2();
      h_bDiscriminator=fs->make<TH1F>("h_bDiscriminator", "pfDeepFlavourJetTags:probb + pfDeepFlavourJetTags:probbb + pfDeepFlavourJetTags:problepb", 25,0,1);
      h_bDiscriminator->Sumw2();
      h_bDiscriminator_m=fs->make<TH1F>("h_bDiscriminator_m", "Matched BDiscriminator", 25,0,1);
      h_bDiscriminator_m->Sumw2();
      h_bDiscriminator_f=fs->make<TH1F>("h_bDiscriminator_f", "Fake BDiscriminator", 25,0,1);
      h_bDiscriminator_f->Sumw2();
      h_nLowPtElectrons=fs->make<TH1F>("h_nLowPtElectrons", "StdElectron Cuts", 5,0,5);
      h_nLowPtElectrons->Sumw2();
      h_nMuons=fs->make<TH1F>("h_nMuons", "Muon Cuts", 5,0,5);
      h_nMuons->Sumw2();
      h_nTaus=fs->make<TH1F>("h_nTaus", "Tau Cuts", 3,0,3);
      h_nTaus->Sumw2();
      h_nJets=fs->make<TH1F>("h_nJets", "Jet Cuts", 11,0,11);
      h_nJets->Sumw2();
      h_TauElectronID=fs->make<TH1F>("h_TauElectronID", "Tau Matched Electron ID(noIso)", 3,0,3);
      h_TauElectronID->SetCanExtend(TH1::kAllAxes);
      h_TauElectronID->Sumw2();

      h_TauElectronIDandIso=fs->make<TH1F>("h_TauElectronIDandIso", "Tau Matched Electron ID And Iso", 3,0,3);
      h_TauElectronIDandIso->SetCanExtend(TH1::kAllAxes);
      h_TauElectronIDandIso->Sumw2();


      h_TauElectronID_FULL=fs->make<TH1F>("h_TauElectronID_FULL", "Tau Matched Electron FullReq", 3,0,3);
      h_TauElectronID_FULL->SetCanExtend(TH1::kAllAxes);
      h_TauElectronID_FULL->Sumw2();
///Event Selection Plots
      h_nEvent_tauE_tauE=fs->make<TH1F>("h_nEvent_tauE_tauE", "NEvent TauE TauE",3,0,3);
      h_nEvent_tauE_tauE->SetCanExtend(TH1::kAllAxes);
      h_nEvent_tauE_tauE->Sumw2();

      h_nEvent_tauE_tauMu=fs->make<TH1F>("h_nEvent_tauE_tauMu", "NEvent TauE TauMu",4,0,4);
      h_nEvent_tauE_tauMu->SetCanExtend(TH1::kAllAxes);
      h_nEvent_tauE_tauMu->Sumw2();

      h_nEvent_tauMu_tauMu=fs->make<TH1F>("h_nEvent_tauMu_tauMu", "NEvent TauMu TauMu",3,0,3);
      h_nEvent_tauMu_tauMu->SetCanExtend(TH1::kAllAxes);

      h_nEvent_tauMu_tauMu->Sumw2();
///Event PT Plots

      h_tauE_tauE_leadtauEpt=fs->make<TH1F>("h_tauE_leadtauE_tauEpt", "tauE|TauE Leading Tau Ele PT (GeV)", 45,0,90);
      h_tauE_tauE_leadtauEpt->Sumw2();

      h_tauE_tauE_subtauEpt=fs->make<TH1F>("h_tauE_tauE_subtauEpt", "tauE|TauE Subleading tau Ele PT (GeV)", 45,0,90);
      h_tauE_tauE_subtauEpt->Sumw2();


      h_tauMu_tauE_tauEpt=fs->make<TH1F>("h_tauMu_tauE_tauEpt", "tauMu|TauE Tau Muon PT (GeV)", 45,0,90);
      h_tauMu_tauE_tauEpt->Sumw2();

      h_tauMu_tauE_tauMupt=fs->make<TH1F>("h_tauMu_tauE_tauMupt", "tauMu|TauE Tau Ele PT (GeV)", 45,0,90);
      h_tauMu_tauE_tauMupt->Sumw2();


      h_tauMu_tauMu_leadtauMupt=fs->make<TH1F>("h_tauMu_tauMu_leadtauMupt", "tauMu|TauMu Leading Tau Muon PT (GeV)", 45,0,90);
      h_tauMu_tauMu_leadtauMupt->Sumw2();

      h_tauMu_tauMu_subtauMupt=fs->make<TH1F>("h_tauMu_tauMu_subtauMupt", "tauMu|TauMu  SubLeading Tau Muon PT (GeV)", 45,0,90);
      h_tauMu_tauMu_subtauMupt->Sumw2();





///////B Parantage Study


      h_nEStatus1_NotfromTau_Gen=fs->make<TH1F>("h_nEStatus1_NotfromTau_Gen", "N Electrons Status 1 (!fromtau) (Gen)", 8,0,8);
      h_nEStatus1_NotfromTau_Gen->Sumw2(); 

      h_nMuStatus1_NotfromTau_Gen=fs->make<TH1F>("h_nMuStatus1_NotfromTau_Gen", "N Muons Status (!fromtau) 1 (Gen)", 8,0,8);
      h_nMuStatus1_NotfromTau_Gen->Sumw2(); 




      h_nEMotherB=fs->make<TH1F>("h_nEMotherB", "N E with b Mother Per Event", 8,0,8);
      h_nEMotherB->Sumw2(); 

      h_nMuMotherB=fs->make<TH1F>("h_nMuMotherB", "N Mu with b Mother Per Event", 5,0,5);
      h_nMuMotherB->Sumw2(); 

      h_genbMu_Pt=fs->make<TH1F>("h_genbMu_Pt", "Gen S1 bMuons PT from A (GeV)", 45,0,90);
      h_genbMu_Pt->Sumw2();

      h_genbE_Pt=fs->make<TH1F>("h_genbE_Pt", "Gen S1 bElectrons from A PT (GeV)", 45,0,90);
      h_genbE_Pt->Sumw2();

      h_LowPtElectron_bmatched_Pt=fs->make<TH1F>("h_LowPtElectron_bmatched_Pt", "b-matched LowPtElectrons PT (GeV)", 45,0,90);
      h_LowPtElectron_bmatched_Pt->Sumw2();

      h_LowPtElectron_bmatched_ID=fs->make<TH1F>("h_LowPtElectron_bmatched_ID", "b-matched LowPtElectrons ID", 45,0,9);
      h_LowPtElectron_bmatched_ID->Sumw2();
      h_LowPtElectron_Taumatched_Pt=fs->make<TH1F>("h_LowPtElectron_Taumatched_Pt", "Tau-matched LowPtElectrons PT (GeV)", 45,0,90);
      h_LowPtElectron_Taumatched_Pt->Sumw2();

      h_LowPtElectron_Taumatched_ID=fs->make<TH1F>("h_LowPtElectron_Taumatched_ID", "Tau-matched LowPtElectrons ID", 45,0,9);
      h_LowPtElectron_Taumatched_ID->Sumw2();
      h_StandardElectron_bmatched_Pt=fs->make<TH1F>("h_StandardElectron_bmatched_Pt", "b-matched Standard Electrons PT (GeV)", 45,0,90);
      h_StandardElectron_bmatched_Pt->Sumw2();

      h_Muon_bmatched_Pt=fs->make<TH1F>("h_Muon_bmatched_Pt", "b-matched Muon PT (GeV)", 45,0,90);
      h_Muon_bmatched_Pt->Sumw2();

      h_Muon_bmatched_Iso=fs->make<TH1F>("h_Muon_bmatched_Iso", "b-matched Muon Isolation", 25,0,1);
      h_Muon_bmatched_Iso->Sumw2();

///Trigger EE

      h_Trigger_tauE_tauE_DoubleElectron=fs->make<TH1F>("h_Trigger_tauE_tauE_DoubleElectron", "tauE|TauE Double Electron Triggers",3,0,3);
      h_Trigger_tauE_tauE_DoubleElectron->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauE_tauE_DoubleElectron->Sumw2();

      h_Trigger_tauE_tauE_singleEle=fs->make<TH1F>("h_Trigger_tauE_tauE_singleEle", "tauE|TauE Single Ele Triggers",3,0,3);
      h_Trigger_tauE_tauE_singleEle->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauE_tauE_singleEle->Sumw2();

      h_Trigger_tauE_tauE_MuEG=fs->make<TH1F>("h_Trigger_tauE_tauE_MuEG", "tauE|TauE MuEgamma Triggers",3,0,3);
      h_Trigger_tauE_tauE_MuEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauE_tauE_MuEG->Sumw2();



//Triger tauETauMu
      h_Trigger_tauE_tauMu_DoubleElectron=fs->make<TH1F>("h_Trigger_tauE_tauMu_DoubleElectron", "tauE|TauMuDouble Electron Triggers",3,0,3);
      h_Trigger_tauE_tauMu_DoubleElectron->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauE_tauMu_DoubleElectron->Sumw2();

      h_Trigger_tauE_tauMu_singleEle=fs->make<TH1F>("h_Trigger_tauE_tauMu_singleEle", "tauE|TauMu Single Ele Triggers",3,0,3);
      h_Trigger_tauE_tauMu_singleEle->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauE_tauMu_singleEle->Sumw2();

      h_Trigger_tauE_tauMu_singleMu=fs->make<TH1F>("h_Trigger_tauE_tauMu_singleMu", "tauE|TauMu Single Mu Triggers",3,0,3);
      h_Trigger_tauE_tauMu_singleMu->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauE_tauMu_singleMu->Sumw2();

      h_Trigger_tauE_tauMu_DoubleMuon=fs->make<TH1F>("h_Trigger_tauE_tauMu_DoubleMuon", "tauE|TauMu Double Mu Triggers",3,0,3);
      h_Trigger_tauE_tauMu_DoubleMuon->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauE_tauMu_DoubleMuon->Sumw2();

      h_Trigger_tauE_tauMu_MuEG=fs->make<TH1F>("h_Trigger_tauE_tauMu_MuEG", "tauE|TauMu MuEgamma Triggers",3,0,3);
      h_Trigger_tauE_tauMu_MuEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauE_tauMu_MuEG->Sumw2();



//Trigger bmuTauMu




//Trigger improvment
//EE
      h_Trigger_tauE_TauE_DiEandEG=fs->make<TH1F>("h_Trigger_tauE_TauE_DiEandEG", "tauE|TauE DiEandEG Triggers",3,0,3);
      h_Trigger_tauE_TauE_DiEandEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauE_TauE_DiEandEG->Sumw2();


      h_Trigger_tauE_TauE_EGandMuEG=fs->make<TH1F>("h_Trigger_tauE_TauE_EGandMuEG", "tauE|TauE EGandMuEG Triggers",3,0,3);
      h_Trigger_tauE_TauE_EGandMuEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauE_TauE_EGandMuEG->Sumw2();
//Invariant Mass
      h_tauE_TauE_Mvis_EG=fs->make<TH1F>("h_tauE_TauE_Mvis_EG", "tauE|TauE Mvis(EG Trigger)",45,0,20);
      h_tauE_TauE_Mvis_EG->Sumw2();

      h_tauE_TauE_Mvis_EGorDoubleEG=fs->make<TH1F>("h_tauE_TauE_Mvis_EGorDoubleEG", "tauE|TauE  Mvis(EG or Double Trigger)",45,0,20);
      h_tauE_TauE_Mvis_EGorDoubleEG->Sumw2();

      h_tauE_TauE_Mvis_noTrig=fs->make<TH1F>("h_tauE_TauE_Mvis_noTrig", "tauE|TauE  Mvis (No Trigger)",45,0,20);
      h_tauE_TauE_Mvis_noTrig->Sumw2();

      h_tauE_TauE_Mvis_DoubleEG=fs->make<TH1F>("h_tauE_TauE_Mvis_DoubleEG", "tauE|TauE  Mvis(Double Trigger)",45,0,20);
      h_tauE_TauE_Mvis_DoubleEG->Sumw2();


//bMass
      h_tauE_TauE_bMass_EG=fs->make<TH1F>("h_tauE_TauE_bMass_EG", "tauE|TauE bMass(EG Trigger)",40,0,20);
      h_tauE_TauE_bMass_EG->Sumw2();

      h_tauE_TauE_bMass_EGorDoubleEG=fs->make<TH1F>("h_tauE_TauE_bMass_EGorDoubleEG", "tauE|TauE  bMass(EG or Double Trigger)",40,0,20);
      h_tauE_TauE_bMass_EGorDoubleEG->Sumw2();

      h_tauE_TauE_bMass_noTrig=fs->make<TH1F>("h_tauE_TauE_bMass_noTrig", "tauE|TauE  bMass (No Trigger)",40,0,20);
      h_tauE_TauE_bMass_noTrig->Sumw2();

      h_tauE_TauE_bMass_DoubleEG=fs->make<TH1F>("h_tauE_TauE_bMass_DoubleEG", "tauE|TauE  bMass(Double Trigger)",40,0,20);
      h_tauE_TauE_bMass_DoubleEG->Sumw2();

//bDRe
      h_tauE_TauE_bDRe_EG=fs->make<TH1F>("h_tauE_TauE_bDRe_EG", "tauE|TauE bDRe(EG Trigger)",40,0,5);
      h_tauE_TauE_bDRe_EG->Sumw2();

      h_tauE_TauE_bDRe_EGorDoubleEG=fs->make<TH1F>("h_tauE_TauE_bDRe_EGorDoubleEG", "tauE|TauE  bDRe(EG or Double Trigger)",40,0,5);
      h_tauE_TauE_bDRe_EGorDoubleEG->Sumw2();

      h_tauE_TauE_bDRe_noTrig=fs->make<TH1F>("h_tauE_TauE_bDRe_noTrig", "tauE|TauE  bDRe (No Trigger)",40,0,5);
      h_tauE_TauE_bDRe_noTrig->Sumw2();

      h_tauE_TauE_bDRe_DoubleEG=fs->make<TH1F>("h_tauE_TauE_bDRe_DoubleEG", "tauE|TauE  bDRe(Double Trigger)",40,0,5);
      h_tauE_TauE_bDRe_DoubleEG->Sumw2();
//eDRe
      h_tauE_TauE_eDRe_EG=fs->make<TH1F>("h_tauE_TauE_eDRe_EG", "tauE|TauE eDRe(EG Trigger)",40,0,4);
      h_tauE_TauE_eDRe_EG->Sumw2();

      h_tauE_TauE_eDRe_EGorDoubleEG=fs->make<TH1F>("h_tauE_TauE_eDRe_EGorDoubleEG", "tauE|TauE  eDRe(EG or Double Trigger)",40,0,4);
      h_tauE_TauE_eDRe_EGorDoubleEG->Sumw2();

      h_tauE_TauE_eDRe_noTrig=fs->make<TH1F>("h_tauE_TauE_eDRe_noTrig", "tauE|TauE  eDRe (No Trigger)",40,0,4);
      h_tauE_TauE_eDRe_noTrig->Sumw2();

      h_tauE_TauE_eDRe_DoubleEG=fs->make<TH1F>("h_tauE_TauE_eDRe_DoubleEG", "tauE|TauE  eDRe(Double Trigger)",40,0,4);
      h_tauE_TauE_eDRe_DoubleEG->Sumw2();
//DiTauPt
      h_tauE_TauE_DiTauPt_EG=fs->make<TH1F>("h_tauE_TauE_DiTauPt_EG", "tauE|TauE DiTauPt(EG Trigger)",50,0,200);
      h_tauE_TauE_DiTauPt_EG->Sumw2();

      h_tauE_TauE_DiTauPt_EGorDoubleEG=fs->make<TH1F>("h_tauE_TauE_DiTauPt_EGorDoubleEG", "tauE|TauE  DiTauPt(EG or Double Trigger)",50,0,200);
      h_tauE_TauE_DiTauPt_EGorDoubleEG->Sumw2();

      h_tauE_TauE_DiTauPt_noTrig=fs->make<TH1F>("h_tauE_TauE_DiTauPt_noTrig", "tauE|TauE  DiTauPt (No Trigger)",50,0,200);
      h_tauE_TauE_DiTauPt_noTrig->Sumw2();

      h_tauE_TauE_DiTauPt_DoubleEG=fs->make<TH1F>("h_tauE_TauE_DiTauPt_DoubleEG", "tauE|TauE  DiTauPt(Double Trigger)",50,0,200);
      h_tauE_TauE_DiTauPt_DoubleEG->Sumw2();
//MET
      h_tauE_TauE_MET_EG=fs->make<TH1F>("h_tauE_TauE_MET_EG", "tauE|TauE MET(EG Trigger)",40,0,250);
      h_tauE_TauE_MET_EG->Sumw2();

      h_tauE_TauE_MET_EGorDoubleEG=fs->make<TH1F>("h_tauE_TauE_MET_EGorDoubleEG", "tauE|TauE  MET(EG or Double Trigger)",40,0,250);
      h_tauE_TauE_MET_EGorDoubleEG->Sumw2();

      h_tauE_TauE_MET_noTrig=fs->make<TH1F>("h_tauE_TauE_MET_noTrig", "tauE|TauE  MET (No Trigger)",40,0,250);
      h_tauE_TauE_MET_noTrig->Sumw2();

      h_tauE_TauE_MET_DoubleEG=fs->make<TH1F>("h_tauE_TauE_MET_DoubleEG", "tauE|TauE  MET(Double Trigger)",40,0,250);
      h_tauE_TauE_MET_DoubleEG->Sumw2();
////// Mu E
      h_Trigger_tauMu_TauE_MuandEG=fs->make<TH1F>("h_Trigger_tauMu_TauE_MuandEG", "tauMu|TauE MuandEG Triggers",3,0,3);
      h_Trigger_tauMu_TauE_MuandEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauMu_TauE_MuandEG->Sumw2();

      h_Trigger_tauMu_TauE_EGandMuEG=fs->make<TH1F>("h_Trigger_tauMu_TauE_EGandMuEG", "tauMu|TauE EGandMuEG Triggers",3,0,3);
      h_Trigger_tauMu_TauE_EGandMuEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauMu_TauE_EGandMuEG->Sumw2();



      h_Trigger_tauMu_TauE_MuandMuEG=fs->make<TH1F>("h_Trigger_tauMu_TauE_MuandMuEG", "tauMu|TauE MuandMuEG Triggers",3,0,3);
      h_Trigger_tauMu_TauE_MuandMuEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauMu_TauE_MuandMuEG->Sumw2();

      h_Trigger_tauMu_tauMu_singleMu=fs->make<TH1F>("h_Trigger_tauMu_tauMu_singleMu", "tauMu|TauMu Single Mu Triggers",3,0,3);
      h_Trigger_tauMu_tauMu_singleMu->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauMu_tauMu_singleMu->Sumw2();

      h_Trigger_tauMu_tauMu_DoubleMuon=fs->make<TH1F>("h_Trigger_tauMu_tauMu_DoubleMuon", "tauMu|TauMu Double Mu Triggers",3,0,3);
      h_Trigger_tauMu_tauMu_DoubleMuon->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauMu_tauMu_DoubleMuon->Sumw2();
//Mvis
      h_tauMu_TauE_Mvis_SingleMu=fs->make<TH1F>("h_tauMu_TauE_Mvis_SingleMu", "tauMu|TauE Mvis(SingleMu Trigger)",45,0,20);
      h_tauMu_TauE_Mvis_SingleMu->Sumw2();

      h_tauMu_TauE_Mvis_EGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_Mvis_EGorSingleMu", "tauMu|TauE Mvis(Single Mu or Single E Trigger)",45,0,20);
      h_tauMu_TauE_Mvis_EGorSingleMu->Sumw2();

      h_tauMu_TauE_Mvis_MuEGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_Mvis_MuEGorSingleMu", "tauMu|TauE Mvis(Single Mu or MuEG Trigger)",45,0,20);
      h_tauMu_TauE_Mvis_MuEGorSingleMu->Sumw2();

      h_tauMu_TauE_Mvis_MuEGorSingleMuorSingleE=fs->make<TH1F>("h_tauMu_TauE_Mvis_MuEGorSingleMuorSingleE", "tauMu|TauE Mvis(Single Mu or MuEG or E Trigger)",45,0,20);
      h_tauMu_TauE_Mvis_MuEGorSingleMuorSingleE->Sumw2();

      h_tauMu_TauE_Mvis_MuEG=fs->make<TH1F>("h_tauMu_TauE_Mvis_MuEG", "tauMu|TauE Mvis(MuEG Trigger)",45,0,20);
      h_tauMu_TauE_Mvis_MuEG->Sumw2();

      h_tauMu_TauE_Mvis_EG=fs->make<TH1F>("h_tauMu_TauE_Mvis_EG", "tauMu|TauE Mvis(Single E Trigger)",45,0,20);
      h_tauMu_TauE_Mvis_EG->Sumw2();

      h_tauMu_TauE_Mvis_noTrig=fs->make<TH1F>("h_tauMu_TauE_Mvis_noTrig", "tauMu|TauE  Mvis (No Trigger)",45,0,20);
      h_tauMu_TauE_Mvis_noTrig->Sumw2();

//bMass
      h_tauMu_TauE_bMass_SingleMu=fs->make<TH1F>("h_tauMu_TauE_bMass_SingleMu", "tauMu|TauE bMass(SingleMu Trigger)",45,0,20);
      h_tauMu_TauE_bMass_SingleMu->Sumw2();

      h_tauMu_TauE_bMass_EGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_bMass_EGorSingleMu", "tauMu|TauE bMass(Single Mu or Single E Trigger)",45,0,20);
      h_tauMu_TauE_bMass_EGorSingleMu->Sumw2();

      h_tauMu_TauE_bMass_MuEGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_bMass_MuEGorSingleMu", "tauMu|TauE bMass(Single Mu or MuEG Trigger)",45,0,20);
      h_tauMu_TauE_bMass_MuEGorSingleMu->Sumw2();

      h_tauMu_TauE_bMass_MuEGorSingleMuorSingleE=fs->make<TH1F>("h_tauMu_TauE_bMass_MuEGorSingleMuorSingleE", "tauMu|TauE bMass(Single Mu or MuEG or E Trigger)",45,0,20);
      h_tauMu_TauE_bMass_MuEGorSingleMuorSingleE->Sumw2();

      h_tauMu_TauE_bMass_MuEG=fs->make<TH1F>("h_tauMu_TauE_bMass_MuEG", "tauMu|TauE bMass(MuEG Trigger)",45,0,20);
      h_tauMu_TauE_bMass_MuEG->Sumw2();

      h_tauMu_TauE_bMass_EG=fs->make<TH1F>("h_tauMu_TauE_bMass_EG", "tauMu|TauE bMass(Single E Trigger)",45,0,20);
      h_tauMu_TauE_bMass_EG->Sumw2();

      h_tauMu_TauE_bMass_noTrig=fs->make<TH1F>("h_tauMu_TauE_bMass_noTrig", "tauMu|TauE  bMass (No Trigger)",45,0,20);
      h_tauMu_TauE_bMass_noTrig->Sumw2();

//bDRmu
      h_tauMu_TauE_bDRmu_SingleMu=fs->make<TH1F>("h_tauMu_TauE_bDRmu_SingleMu", "tauMu|TauE bDRmu(SingleMu Trigger)",45,0,5);
      h_tauMu_TauE_bDRmu_SingleMu->Sumw2();

      h_tauMu_TauE_bDRmu_EGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_bDRmu_EGorSingleMu", "tauMu|TauE bDRmu(Single Mu or Single E Trigger)",45,0,5);
      h_tauMu_TauE_bDRmu_EGorSingleMu->Sumw2();

      h_tauMu_TauE_bDRmu_MuEGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_bDRmu_MuEGorSingleMu", "tauMu|TauE bDRmu(Single Mu or MuEG Trigger)",45,0,5);
      h_tauMu_TauE_bDRmu_MuEGorSingleMu->Sumw2();

      h_tauMu_TauE_bDRmu_MuEGorSingleMuorSingleE=fs->make<TH1F>("h_tauMu_TauE_bDRmu_MuEGorSingleMuorSingleE", "tauMu|TauE bDRmu(Single Mu or MuEG or E Trigger)",45,0,5);
      h_tauMu_TauE_bDRmu_MuEGorSingleMuorSingleE->Sumw2();

      h_tauMu_TauE_bDRmu_MuEG=fs->make<TH1F>("h_tauMu_TauE_bDRmu_MuEG", "tauMu|TauE bDRmu(MuEG Trigger)",45,0,5);
      h_tauMu_TauE_bDRmu_MuEG->Sumw2();

      h_tauMu_TauE_bDRmu_EG=fs->make<TH1F>("h_tauMu_TauE_bDRmu_EG", "tauMu|TauE bDRmu(Single E Trigger)",45,0,5);
      h_tauMu_TauE_bDRmu_EG->Sumw2();

      h_tauMu_TauE_bDRmu_noTrig=fs->make<TH1F>("h_tauMu_TauE_bDRmu_noTrig", "tauMu|TauE  bDRmu (No Trigger)",45,0,5);
      h_tauMu_TauE_bDRmu_noTrig->Sumw2();

//DiTauPt
      h_tauMu_TauE_DiTauPt_SingleMu=fs->make<TH1F>("h_tauMu_TauE_DiTauPt_SingleMu", "tauMu|TauE DiTauPt(SingleMu Trigger)",50,0,200);
      h_tauMu_TauE_DiTauPt_SingleMu->Sumw2();

      h_tauMu_TauE_DiTauPt_EGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_DiTauPt_EGorSingleMu", "tauMu|TauE DiTauPt(Single Mu or Single E Trigger)",50,0,200);
      h_tauMu_TauE_DiTauPt_EGorSingleMu->Sumw2();

      h_tauMu_TauE_DiTauPt_MuEGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_DiTauPt_MuEGorSingleMu", "tauMu|TauE DiTauPt(Single Mu or MuEG Trigger)",50,0,200);
      h_tauMu_TauE_DiTauPt_MuEGorSingleMu->Sumw2();

      h_tauMu_TauE_DiTauPt_MuEGorSingleMuorSingleE=fs->make<TH1F>("h_tauMu_TauE_DiTauPt_MuEGorSingleMuorSingleE", "tauMu|TauE DiTauPt(Single Mu or MuEG or E Trigger)",50,0,200);
      h_tauMu_TauE_DiTauPt_MuEGorSingleMuorSingleE->Sumw2();

      h_tauMu_TauE_DiTauPt_MuEG=fs->make<TH1F>("h_tauMu_TauE_DiTauPt_MuEG", "tauMu|TauE DiTauPt(MuEG Trigger)",50,0,200);
      h_tauMu_TauE_DiTauPt_MuEG->Sumw2();

      h_tauMu_TauE_DiTauPt_EG=fs->make<TH1F>("h_tauMu_TauE_DiTauPt_EG", "tauMu|TauE DiTauPt(Single E Trigger)",50,0,200);
      h_tauMu_TauE_DiTauPt_EG->Sumw2();

      h_tauMu_TauE_DiTauPt_noTrig=fs->make<TH1F>("h_tauMu_TauE_DiTauPt_noTrig", "tauMu|TauE  DiTauPt (No Trigger)",50,0,200);
      h_tauMu_TauE_DiTauPt_noTrig->Sumw2();

//muDRe
      h_tauMu_TauE_muDRe_SingleMu=fs->make<TH1F>("h_tauMu_TauE_muDRe_SingleMu", "tauMu|TauE muDRe(SingleMu Trigger)",40,0,4);
      h_tauMu_TauE_muDRe_SingleMu->Sumw2();

      h_tauMu_TauE_muDRe_EGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_muDRe_EGorSingleMu", "tauMu|TauE muDRe(Single Mu or Single E Trigger)",40,0,4);
      h_tauMu_TauE_muDRe_EGorSingleMu->Sumw2();

      h_tauMu_TauE_muDRe_MuEGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_muDRe_MuEGorSingleMu", "tauMu|TauE muDRe(Single Mu or MuEG Trigger)",40,0,4);
      h_tauMu_TauE_muDRe_MuEGorSingleMu->Sumw2();

      h_tauMu_TauE_muDRe_MuEGorSingleMuorSingleE=fs->make<TH1F>("h_tauMu_TauE_muDRe_MuEGorSingleMuorSingleE", "tauMu|TauE muDRe(Single Mu or MuEG or E Trigger)",40,0,4);
      h_tauMu_TauE_muDRe_MuEGorSingleMuorSingleE->Sumw2();

      h_tauMu_TauE_muDRe_MuEG=fs->make<TH1F>("h_tauMu_TauE_muDRe_MuEG", "tauMu|TauE muDRe(MuEG Trigger)",40,0,4);
      h_tauMu_TauE_muDRe_MuEG->Sumw2();

      h_tauMu_TauE_muDRe_EG=fs->make<TH1F>("h_tauMu_TauE_muDRe_EG", "tauMu|TauE muDRe(Single E Trigger)",40,0,4);
      h_tauMu_TauE_muDRe_EG->Sumw2();

      h_tauMu_TauE_muDRe_noTrig=fs->make<TH1F>("h_tauMu_TauE_muDRe_noTrig", "tauMu|TauE  muDRe (No Trigger)",40,0,4);
      h_tauMu_TauE_muDRe_noTrig->Sumw2();


//muIso
      h_tauMu_TauE_muIso_SingleMu=fs->make<TH1F>("h_tauMu_TauE_muIso_SingleMu", "tauMu|TauE muIso(SingleMu Trigger)",20,0,1);
      h_tauMu_TauE_muIso_SingleMu->Sumw2();

      h_tauMu_TauE_muIso_EGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_muIso_EGorSingleMu", "tauMu|TauE muIso(Single Mu or Single E Trigger)",20,0,1);
      h_tauMu_TauE_muIso_EGorSingleMu->Sumw2();

      h_tauMu_TauE_muIso_MuEGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_muIso_MuEGorSingleMu", "tauMu|TauE muIso(Single Mu or MuEG Trigger)",20,0,1);
      h_tauMu_TauE_muIso_MuEGorSingleMu->Sumw2();

      h_tauMu_TauE_muIso_MuEGorSingleMuorSingleE=fs->make<TH1F>("h_tauMu_TauE_muIso_MuEGorSingleMuorSingleE", "tauMu|TauE muIso(Single Mu or MuEG or E Trigger)",20,0,1);
      h_tauMu_TauE_muIso_MuEGorSingleMuorSingleE->Sumw2();

      h_tauMu_TauE_muIso_MuEG=fs->make<TH1F>("h_tauMu_TauE_muIso_MuEG", "tauMu|TauE muIso(MuEG Trigger)",20,0,1);
      h_tauMu_TauE_muIso_MuEG->Sumw2();

      h_tauMu_TauE_muIso_EG=fs->make<TH1F>("h_tauMu_TauE_muIso_EG", "tauMu|TauE muIso(Single E Trigger)",20,0,1);
      h_tauMu_TauE_muIso_EG->Sumw2();

      h_tauMu_TauE_muIso_noTrig=fs->make<TH1F>("h_tauMu_TauE_muIso_noTrig", "tauMu|TauE  muIso (No Trigger)",20,0,1);
      h_tauMu_TauE_muIso_noTrig->Sumw2();
//MET
      h_tauMu_TauE_MET_SingleMu=fs->make<TH1F>("h_tauMu_TauE_MET_SingleMu", "tauMu|TauE MET(SingleMu Trigger)",40,0,250);
      h_tauMu_TauE_MET_SingleMu->Sumw2();

      h_tauMu_TauE_MET_EGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_MET_EGorSingleMu", "tauMu|TauE MET(Single Mu or Single E Trigger)",40,0,250);
      h_tauMu_TauE_MET_EGorSingleMu->Sumw2();

      h_tauMu_TauE_MET_MuEGorSingleMu=fs->make<TH1F>("h_tauMu_TauE_MET_MuEGorSingleMu", "tauMu|TauE MET(Single Mu or MuEG Trigger)",40,0,250);
      h_tauMu_TauE_MET_MuEGorSingleMu->Sumw2();

      h_tauMu_TauE_MET_MuEGorSingleMuorSingleE=fs->make<TH1F>("h_tauMu_TauE_MET_MuEGorSingleMuorSingleE", "tauMu|TauE MET(Single Mu or MuEG or E Trigger)",40,0,250);
      h_tauMu_TauE_MET_MuEGorSingleMuorSingleE->Sumw2();

      h_tauMu_TauE_MET_MuEG=fs->make<TH1F>("h_tauMu_TauE_MET_MuEG", "tauMu|TauE MET(MuEG Trigger)",40,0,250);
      h_tauMu_TauE_MET_MuEG->Sumw2();

      h_tauMu_TauE_MET_EG=fs->make<TH1F>("h_tauMu_TauE_MET_EG", "tauMu|TauE MET(Single E Trigger)",40,0,250);
      h_tauMu_TauE_MET_EG->Sumw2();

      h_tauMu_TauE_MET_noTrig=fs->make<TH1F>("h_tauMu_TauE_MET_noTrig", "tauMu|TauE  MET (No Trigger)",40,0,250);
      h_tauMu_TauE_MET_noTrig->Sumw2();
// Mu Mu
      h_Trigger_tauMu_TauMu_DiMuandMu=fs->make<TH1F>("h_Trigger_tauMu_TauMu_DiMuandMu", "tauMu|TauMu DiMuandMu Triggers",3,0,3);
      h_Trigger_tauMu_TauMu_DiMuandMu->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauMu_TauMu_DiMuandMu->Sumw2();

      h_Trigger_tauMu_TauMu_MuandMuEG=fs->make<TH1F>("h_Trigger_tauMu_TauMu_MuandMuEG", "tauMu|TauMu MuandMuEG Triggers",3,0,3);
      h_Trigger_tauMu_TauMu_MuandMuEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_tauMu_TauMu_MuandMuEG->Sumw2();

//Mvis
      h_tauMu_TauMu_Mvis_SingleMu=fs->make<TH1F>("h_tauMu_TauMu_Mvis_SingleMu", "tauMu|tauMu Mvis(SingleMu Trigger)",45,0,20);
      h_tauMu_TauMu_Mvis_SingleMu->Sumw2();

      h_tauMu_TauMu_Mvis_SingleMuorDoubleMu=fs->make<TH1F>("h_tauMu_TauMu_Mvis_SingleMuorDoubleMu", "tauMu|tauMu Mvis(DoubleMu or Single Trigger)",45,0,20);
      h_tauMu_TauMu_Mvis_SingleMuorDoubleMu->Sumw2();

      h_tauMu_TauMu_Mvis_DoubleMu=fs->make<TH1F>("h_tauMu_TauMu_Mvis_DoubleMu", "tauMu|tauMu Mvis(DoubleMu)",45,0,20);
      h_tauMu_TauMu_Mvis_DoubleMu->Sumw2();

      h_tauMu_TauMu_Mvis_noTrig=fs->make<TH1F>("h_tauMu_TauMu_Mvis_noTrig", "tauMu|tauMu Mvis(No Trigger)",45,0,20);
      h_tauMu_TauMu_Mvis_noTrig->Sumw2();
//bMass
      h_tauMu_TauMu_bMass_SingleMu=fs->make<TH1F>("h_tauMu_TauMu_bMass_SingleMu", "tauMu|tauMu bMass(SingleMu Trigger)",45,0,20);
      h_tauMu_TauMu_bMass_SingleMu->Sumw2();

      h_tauMu_TauMu_bMass_SingleMuorDoubleMu=fs->make<TH1F>("h_tauMu_TauMu_bMass_SingleMuorDoubleMu", "tauMu|tauMu bMass(DoubleMu or Single Trigger)",45,0,20);
      h_tauMu_TauMu_bMass_SingleMuorDoubleMu->Sumw2();

      h_tauMu_TauMu_bMass_DoubleMu=fs->make<TH1F>("h_tauMu_TauMu_bMass_DoubleMu", "tauMu|tauMu bMass(DoubleMu)",45,0,20);
      h_tauMu_TauMu_bMass_DoubleMu->Sumw2();

      h_tauMu_TauMu_bMass_noTrig=fs->make<TH1F>("h_tauMu_TauMu_bMass_noTrig", "tauMu|tauMu bMass(No Trigger)",45,0,20);
      h_tauMu_TauMu_bMass_noTrig->Sumw2();
//bDRmu
      h_tauMu_TauMu_bDRmu_SingleMu=fs->make<TH1F>("h_tauMu_TauMu_bDRmu_SingleMu", "tauMu|tauMu bDRmu(SingleMu Trigger)",45,0,5);
      h_tauMu_TauMu_bDRmu_SingleMu->Sumw2();

      h_tauMu_TauMu_bDRmu_SingleMuorDoubleMu=fs->make<TH1F>("h_tauMu_TauMu_bDRmu_SingleMuorDoubleMu", "tauMu|tauMu bDRmu(DoubleMu or Single Trigger)",45,0,5);
      h_tauMu_TauMu_bDRmu_SingleMuorDoubleMu->Sumw2();

      h_tauMu_TauMu_bDRmu_DoubleMu=fs->make<TH1F>("h_tauMu_TauMu_bDRmu_DoubleMu", "tauMu|tauMu bDRmu(DoubleMu)",45,0,5);
      h_tauMu_TauMu_bDRmu_DoubleMu->Sumw2();

      h_tauMu_TauMu_bDRmu_noTrig=fs->make<TH1F>("h_tauMu_TauMu_bDRmu_noTrig", "tauMu|tauMu bDRmu(No Trigger)",45,0,5);
      h_tauMu_TauMu_bDRmu_noTrig->Sumw2();

//DiTauPt
      h_tauMu_TauMu_DiTauPt_SingleMu=fs->make<TH1F>("h_tauMu_TauMu_DiTauPt_SingleMu", "tauMu|tauMu DiTauPt(SingleMu Trigger)",50,0,200);
      h_tauMu_TauMu_DiTauPt_SingleMu->Sumw2();

      h_tauMu_TauMu_DiTauPt_SingleMuorDoubleMu=fs->make<TH1F>("h_tauMu_TauMu_DiTauPt_SingleMuorDoubleMu", "tauMu|tauMu DiTauPt(DoubleMu or Single Trigger)",50,0,200);
      h_tauMu_TauMu_DiTauPt_SingleMuorDoubleMu->Sumw2();

      h_tauMu_TauMu_DiTauPt_DoubleMu=fs->make<TH1F>("h_tauMu_TauMu_DiTauPt_DoubleMu", "tauMu|tauMu DiTauPt(DoubleMu)",50,0,200);
      h_tauMu_TauMu_DiTauPt_DoubleMu->Sumw2();

      h_tauMu_TauMu_DiTauPt_noTrig=fs->make<TH1F>("h_tauMu_TauMu_DiTauPt_noTrig", "tauMu|tauMu DiTauPt(No Trigger)",50,0,200);
      h_tauMu_TauMu_DiTauPt_noTrig->Sumw2();
//muDRmu
      h_tauMu_TauMu_muDRmu_SingleMu=fs->make<TH1F>("h_tauMu_TauMu_muDRmu_SingleMu", "tauMu|tauMu muDRmu(SingleMu Trigger)",40,0,4);
      h_tauMu_TauMu_muDRmu_SingleMu->Sumw2();

      h_tauMu_TauMu_muDRmu_SingleMuorDoubleMu=fs->make<TH1F>("h_tauMu_TauMu_muDRmu_SingleMuorDoubleMu", "tauMu|tauMu muDRmu(DoubleMu or Single Trigger)",40,0,4);
      h_tauMu_TauMu_muDRmu_SingleMuorDoubleMu->Sumw2();

      h_tauMu_TauMu_muDRmu_DoubleMu=fs->make<TH1F>("h_tauMu_TauMu_muDRmu_DoubleMu", "tauMu|tauMu muDRmu(DoubleMu)",40,0,4);
      h_tauMu_TauMu_muDRmu_DoubleMu->Sumw2();

      h_tauMu_TauMu_muDRmu_noTrig=fs->make<TH1F>("h_tauMu_TauMu_muDRmu_noTrig", "tauMu|tauMu muDRmu(No Trigger)",40,0,4);
      h_tauMu_TauMu_muDRmu_noTrig->Sumw2();

//muIso
      h_tauMu_TauMu_muIso_SingleMu=fs->make<TH1F>("h_tauMu_TauMu_muIso_SingleMu", "tauMu|tauMu muIso(SingleMu Trigger)",20,0,1);
      h_tauMu_TauMu_muIso_SingleMu->Sumw2();

      h_tauMu_TauMu_muIso_SingleMuorDoubleMu=fs->make<TH1F>("h_tauMu_TauMu_muIso_SingleMuorDoubleMu", "tauMu|tauMu muIso(DoubleMu or Single Trigger)",20,0,1);
      h_tauMu_TauMu_muIso_SingleMuorDoubleMu->Sumw2();

      h_tauMu_TauMu_muIso_DoubleMu=fs->make<TH1F>("h_tauMu_TauMu_muIso_DoubleMu", "tauMu|tauMu muIso(DoubleMu)",20,0,1);
      h_tauMu_TauMu_muIso_DoubleMu->Sumw2();

      h_tauMu_TauMu_muIso_noTrig=fs->make<TH1F>("h_tauMu_TauMu_muIso_noTrig", "tauMu|tauMu muIso(No Trigger)",20,0,1);
      h_tauMu_TauMu_muIso_noTrig->Sumw2();

//MET
      h_tauMu_TauMu_MET_SingleMu=fs->make<TH1F>("h_tauMu_TauMu_MET_SingleMu", "tauMu|tauMu MET(SingleMu Trigger)",40,0,250);
      h_tauMu_TauMu_MET_SingleMu->Sumw2();

      h_tauMu_TauMu_MET_SingleMuorDoubleMu=fs->make<TH1F>("h_tauMu_TauMu_MET_SingleMuorDoubleMu", "tauMu|tauMu MET(DoubleMu or Single Trigger)",40,0,250);
      h_tauMu_TauMu_MET_SingleMuorDoubleMu->Sumw2();

      h_tauMu_TauMu_MET_DoubleMu=fs->make<TH1F>("h_tauMu_TauMu_MET_DoubleMu", "tauMu|tauMu MET(DoubleMu)",40,0,250);
      h_tauMu_TauMu_MET_DoubleMu->Sumw2();

      h_tauMu_TauMu_MET_noTrig=fs->make<TH1F>("h_tauMu_TauMu_MET_noTrig", "tauMu|tauMu MET(No Trigger)",40,0,250);
      h_tauMu_TauMu_MET_noTrig->Sumw2();
//General Mu

      h_Trigger_general_MuandDiMu=fs->make<TH1F>("h_Trigger_general_MuandDiMu", "General Mu&DiMu Triggers",3,0,3);
      h_Trigger_general_MuandDiMu->SetCanExtend(TH1::kAllAxes);
      h_Trigger_general_MuandDiMu->Sumw2();

      h_Trigger_general_MuandMuEG=fs->make<TH1F>("h_Trigger_general_MuandMuEG", "General Mu&MuEG Triggers",3,0,3);
      h_Trigger_general_MuandMuEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_general_MuandMuEG->Sumw2();

      h_Trigger_general_EGandMuEG=fs->make<TH1F>("h_Trigger_general_EGandMuEG", "General EG&MuEG Triggers",3,0,3);
      h_Trigger_general_EGandMuEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_general_EGandMuEG->Sumw2();

      h_Trigger_general_EGandDiEG=fs->make<TH1F>("h_Trigger_general_EGandDiEG", "General Eg&DiEG Triggers",3,0,3);
      h_Trigger_general_EGandDiEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_general_EGandDiEG->Sumw2();

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

   edm::Handle<std::vector<pat::Tau>> taus;
	iEvent.getByToken(taus_, taus);
   edm::Handle<std::vector<pat::MET>> mets;
	iEvent.getByToken(mets_, mets);
   edm::Handle<std::vector<pat::Muon>> muons;
	iEvent.getByToken(muons_, muons);
   edm::Handle<std::vector<pat::Jet>> jets;
	iEvent.getByToken(jets_, jets);
   edm::Handle<std::vector<pat::Electron>> std_electrons;
	iEvent.getByToken(std_electrons_, std_electrons);
   edm::Handle<std::vector<pat::Electron>> lowPt_electrons;
	iEvent.getByToken(lowPt_electrons_, lowPt_electrons);
   edm::Handle<reco::JetFlavourInfoMatchingCollection> jetFlavorMatching;
   iEvent.getByToken(jetFlavorMatching_,jetFlavorMatching);
   
   //std::vector<reco::GenParticle> genParticles = h2AA::getObject<std::vector<reco::GenParticle>>(genParticles_, iEvent);
   edm::Handle<std::vector<reco::GenParticle>> genParticles;
	iEvent.getByToken(genParticles_, genParticles);

///Get Trigger objects for out handles
   edm::InputTag trigResultsTag("TriggerResults","","HLT"); //make sure have correct process on MC
   edm::Handle<edm::TriggerResults>  triggerBits;
   iEvent.getByLabel(trigResultsTag,triggerBits);

   edm::TriggerResults triggerResults = *triggerBits.product();
   cout<<to_string(triggerResults.size());
   const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerBits);  

   //for (unsigned int i = 0; i<trigNames.size();i++){cout<< to_string(i)<<" "<< trigNames.triggerName(i)<<" "<< triggerResults.accept(i)<<"\n";}
   double rho = h2AA::getObject<double>(eventrho_,iEvent);
   edm::Handle<GenEventInfoProduct>  genInfoHandle;
   try{iEvent.getByToken(genInfoProduct_, genInfoHandle);}
   catch(...){;}
   double genWeight=genInfoHandle->weight();
   pat::MET met=(*mets.product())[0];
   h_MET->Fill(met.pt(),genWeight);
   float MET = met.pt();



   
   //////Gen Particle Sorting
   std::vector<const reco::GenParticle*> genTauElectrons;
   std::vector<const reco::GenParticle*> genTauMuons;
   std::vector<const reco::GenParticle*> genTaus;
   std::vector<const reco::GenParticle*> genBs;
   std::vector<const reco::GenParticle*> genS1Electrons_NonTau;
   std::vector<const reco::GenParticle*> genS1Muons_NonTau;
   std::vector<const reco::GenParticle*> genNuTau; //////These are just the tau neutrinos We know there is a maximum of 2 per event either one anti, one normal, both or none.
   std::vector<const reco::GenParticle*> genNuTauBar; 
   for(const auto& gen : *genParticles){
      if (gen.status()==1){
         h_GenParticleID->Fill(abs(gen.pdgId())+.5,genWeight);
      }
      if(abs(gen.pdgId())==11 &&
         gen.isDirectHardProcessTauDecayProductFinalState()) {
         h_GenTauElectronPt->Fill(gen.pt(),genWeight);
         genTauElectrons.push_back(&gen);

      }
      if(abs(gen.pdgId())==11&&gen.status()==1) {
         if (!gen.isDirectHardProcessTauDecayProductFinalState()){
            genS1Electrons_NonTau.push_back(&gen);

         }
      }
      if(abs(gen.pdgId())==13 &&gen.status()==1) {
         if (!gen.isDirectHardProcessTauDecayProductFinalState()){genS1Muons_NonTau.push_back(&gen);}
      }
      if (abs(gen.pdgId())==13 && gen.isDirectHardProcessTauDecayProductFinalState()){
         h_GenTauMuonPt->Fill(gen.pt(),genWeight);
         genTauMuons.push_back(&gen);
      }
      if (abs(gen.pdgId())==15 && gen.isHardProcess()){
         h_GenTauPt->Fill(gen.pt(),genWeight);
         genTaus.push_back(&gen);
      }

      if (gen.pdgId()==16 && gen.isDirectHardProcessTauDecayProductFinalState()){
         genNuTau.push_back(&gen);

      }
      if (gen.pdgId()==-16 && gen.isDirectHardProcessTauDecayProductFinalState()){
         genNuTauBar.push_back(&gen);
      } 
      if (abs(gen.pdgId())==5&&gen.isHardProcess()){
         h_GenbPt->Fill(gen.pt(),genWeight);
         genBs.push_back(&gen);

      
      }

   }

   ////////////////////BJet ELECTRONS///////////////////////////////////////////
   h_nEStatus1_NotfromTau_Gen->Fill(genS1Electrons_NonTau.size(),genWeight);
   h_nMuStatus1_NotfromTau_Gen->Fill(genS1Muons_NonTau.size(),genWeight);

   std::vector<const reco::GenParticle*> genEfromB_S1;

   for(auto& ele : genS1Electrons_NonTau){
      int n = 0;
      int particleID = ele->pdgId();
      auto *mom = ele->mother();
      bool isfromb=false;
      bool isfroma=false;
      int steps_afterb= 0;
      TLorentzVector e;
      int nbinChian=0;
      e.SetPtEtaPhiM(ele->pt(),ele->eta(),ele->phi(),ele->mass());
      TLorentzVector b;
      while(mom!=NULL){
         particleID = mom->pdgId();//get the particle id of the nth mother
         mom=mom->mother();
         if (abs(particleID)==5){
            if (!isfromb){
               b.SetPtEtaPhiM(mom->pt(),mom->eta(),mom->phi(),mom->mass());
               
            }

            isfromb = true;
            nbinChian++;
         }
         if (abs(particleID)==36){isfroma=true;}
         if(!isfromb){n++;}
         if (isfromb){steps_afterb++;}

      }
      if (isfromb&&isfroma){
         genEfromB_S1.push_back(ele);
         h_genbE_Pt->Fill(ele->pt(),genWeight);
      }

   }
   h_nEMotherB->Fill(genEfromB_S1.size(),genWeight);
   /////////////////////////////////B Muon Sorting/////////////////////////

   std::vector<const reco::GenParticle*> genMufromB_S1;

   for(auto& muon : genS1Muons_NonTau){
      int n = 0;
      int particleID = muon->pdgId();
      auto *mom = muon->mother();
      bool isfromb=false;
      bool isfroma=false;
      int steps_afterb= 0;
      TLorentzVector mu;
      int nbinChian=0;
      mu.SetPtEtaPhiM(muon->pt(),muon->eta(),muon->phi(),muon->mass());
      TLorentzVector b;
      while(mom!=NULL){
         particleID = mom->pdgId();//get the particle id of the nth mother
         mom=mom->mother();
         if (abs(particleID)==5){
            if (!isfromb){
               b.SetPtEtaPhiM(mom->pt(),mom->eta(),mom->phi(),mom->mass());
               
            }
            isfromb = true;
            nbinChian++;
         }
         if (abs(particleID)==36){isfroma=true;}
         if(!isfromb){n++;}
         if (isfromb){steps_afterb++;}

      }
      if (isfromb&&isfroma){
         genMufromB_S1.push_back(muon);
         h_genbMu_Pt->Fill(muon->pt(),genWeight);
      }

   }
   h_nMuMotherB->Fill(genMufromB_S1.size(),genWeight);
   //std::cout<<"GenParticle Sorting\n";
   //Determin which is the neutrino associated with 


   //std::cout<<"Tau DR plotting\n";
   bool DiTauHad=false;
   std::vector<const reco::GenParticle*> genTauHad;
   std::vector<const reco::GenParticle*> genNuTauHad;
   //std::cout<<to_string(genTaus.size())<<"=size\n";
   for (auto& gen : genTaus){
      if ((genTauElectrons.size()==1&& genTauMuons.size()==0)){
         if ((*genTauElectrons[0]).charge()<0&&(*gen).charge()>0){
            //std::cout<<"e-\n";

            genTauHad.push_back(gen);
            //std::cout<<"e-test1\n";
            genNuTauHad.push_back(genNuTauBar[0]);
            //std::cout<<"e-test2\n";
         }
         //std::cout<<"test3\n";
         if ((*genTauElectrons[0]).charge()>0&&(*gen).charge()<0){
            //std::cout<<"e+\n";
            genTauHad.push_back(gen);
            genNuTauHad.push_back(genNuTau[0]);
         }
      }
      else if((genTauElectrons.size()==0&& genTauMuons.size()==1)){
            //std::cout<<"test4\n";
            if ((*genTauMuons[0]).charge()<0&&(*gen).charge()>0){
            //std::cout<<"m-\n";
            genTauHad.push_back(gen);
            genNuTauHad.push_back(genNuTauBar[0]);
         }
         if ((*genTauMuons[0]).charge()>0&&(*gen).charge()<0){
            //std::cout<<"m+\n";
            genTauHad.push_back(gen);
            genNuTauHad.push_back(genNuTauBar[0]);
         }

      }
      
      else{
         //std::cout<<"test5\n";
         DiTauHad=true;
      }
      //std::cout<<"test6\n";
   }
   

   //std::cout<<"GenTaus\n";
   std::vector<pat::Electron> selected_std_tau_electrons;
   std::vector<pat::Electron> selected_std_b_electrons;
   std::vector<pat::Electron> selected_std_electrons;
   std::vector<int> matchedTauEleGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level
   std::vector<int> matchedbEleGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level

   for(auto& ele : *std_electrons){
      h_nStdElectrons->Fill(0.5,genWeight);
      if (ele.pt()>7 && abs(ele.eta())<2.5){
         h_nStdElectrons->Fill(1.5,genWeight);
         float E_c=ele.superCluster()->energy();
         if(ele.electronID("cutBasedElectronID-Fall17-94X-V2-loose")==1){
            h_nStdElectrons->Fill(2.5,genWeight);
            selected_std_electrons.push_back(ele);
            bool TauEmatched= false;
            float tauDr_min=.1;
            float tauMatchedGen_index=-1;
            for(unsigned int iGenE=0; iGenE<genTauElectrons.size(); iGenE++){
               if (std::count(matchedTauEleGenIndexs.begin(), matchedTauEleGenIndexs.end(), iGenE)==0){
                  reco::GenParticle gen = *genTauElectrons[iGenE];
                  TLorentzVector e;
                  TLorentzVector genE;
                  e.SetPtEtaPhiM(ele.pt(), ele.eta(), ele.phi(), ele.mass());
                  genE.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), gen.mass());
                  float dr = genE.DeltaR(e);
                  if(dr <tauDr_min){//By setting this we make sure we always match to the closest Dr Particle
                     TauEmatched = true;
                     tauDr_min=dr;
                     tauMatchedGen_index=iGenE;///Keep track of the gen particle that is closest to the reco candidate
                  }
               }
            }
            if (TauEmatched){
               h_std_TauElectronPt->Fill(ele.pt(),genWeight);
               h_nStdElectrons->Fill(3.5,genWeight);
               selected_std_tau_electrons.push_back(ele);
               matchedTauEleGenIndexs.push_back(tauMatchedGen_index);
            }

            h_nStdElectrons->Fill(2.5,genWeight);
            bool bEmatched= false;
            float bDr_min=.1;
            float bMatchedGen_index=-1;
            for(unsigned int iGenE=0; iGenE<genEfromB_S1.size(); iGenE++){
               if (std::count(matchedbEleGenIndexs.begin(), matchedbEleGenIndexs.end(), iGenE)==0){
                  reco::GenParticle gen = *genEfromB_S1[iGenE];
                  TLorentzVector e;
                  TLorentzVector genE;
                  e.SetPtEtaPhiM(ele.pt(), ele.eta(), ele.phi(), ele.mass());
                  genE.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), gen.mass());
                  float dr = genE.DeltaR(e);
                  if(dr <bDr_min){//By setting this we make sure we always match to the closest Dr Particle
                     bEmatched = true;
                     bDr_min=dr;
                     bMatchedGen_index=iGenE;///Keep track of the gen particle that is closest to the reco candidate
                  }
               }
            }
            if (bEmatched){
               h_StandardElectron_bmatched_Pt->Fill(ele.pt(),genWeight);
               h_nStdElectrons->Fill(4.5,genWeight);
               selected_std_b_electrons.push_back(ele);
               matchedbEleGenIndexs.push_back(bMatchedGen_index);
           
            }
         }
      }
   }
   //std::cout<<"stdEle\n";
   std::vector<pat::ElectronRef> selected_LowPt_tau_electrons;
   std::vector<pat::ElectronRef> selected_LowPt_b_electrons;
   std::vector<int> matchedLowPtTauEleGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level
   std::vector<int> matchedLowPtbEleGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level

   for(unsigned int iElec = 0; iElec<lowPt_electrons->size(); iElec++){
      pat::ElectronRef ele (lowPt_electrons, iElec);
      h_nLowPtElectrons->Fill(0.5,genWeight);
      if (ele->pt()>1 && abs(ele->eta())<2.5){
         h_nLowPtElectrons->Fill(1.5,genWeight);
         float idVal = ele->electronID("ID");
         bool tauEmatched= false;
         float tauDr_min=.1;
         float tauMatchedGen_index=-1;
         for(unsigned int iGenE=0; iGenE<genTauElectrons.size(); iGenE++){
            if (std::count(matchedLowPtTauEleGenIndexs.begin(), matchedLowPtTauEleGenIndexs.end(), iGenE)==0){
               reco::GenParticle gen = *genTauElectrons[iGenE];
               TLorentzVector e;
               TLorentzVector genE;
               e.SetPtEtaPhiM(ele->pt(), ele->eta(), ele->phi(), ele->mass());
               genE.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), gen.mass());
               float dr = genE.DeltaR(e);
               if(dr <tauDr_min){//By setting this we make sure we always match to the closest Dr Particle
                  tauEmatched = true;
                  tauDr_min=dr;
                  tauMatchedGen_index=iGenE;///Keep track of the gen particle that is closest to the reco candidate
               }
            }
         }
         if (tauEmatched){
            h_LowPtElectron_Taumatched_Pt->Fill(ele->pt(),genWeight);
            h_LowPtElectron_Taumatched_ID->Fill(idVal,genWeight);
            h_nLowPtElectrons->Fill(2.5,genWeight);
            selected_LowPt_tau_electrons.push_back(ele);
            matchedLowPtTauEleGenIndexs.push_back(tauMatchedGen_index);
         
         }
         bool bEmatched= false;
         float bDr_min=.1;
         float bMatchedGen_index=-1;
         for(unsigned int iGenE=0; iGenE<genEfromB_S1.size(); iGenE++){
            if (std::count(matchedLowPtbEleGenIndexs.begin(), matchedLowPtbEleGenIndexs.end(), iGenE)==0){
               reco::GenParticle gen = *genEfromB_S1[iGenE];
               TLorentzVector e;
               TLorentzVector genE;
               e.SetPtEtaPhiM(ele->pt(), ele->eta(), ele->phi(), ele->mass());
               genE.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), gen.mass());
               float dr = genE.DeltaR(e);
               if(dr <bDr_min){//By setting this we make sure we always match to the closest Dr Particle
                  bEmatched = true;
                  bDr_min=dr;
                  bMatchedGen_index=iGenE;///Keep track of the gen particle that is closest to the reco candidate
               }
            }
         }
         if (bEmatched){
            h_LowPtElectron_bmatched_Pt->Fill(ele->pt(),genWeight);
            h_LowPtElectron_bmatched_ID->Fill(idVal,genWeight);
            h_nLowPtElectrons->Fill(3.5,genWeight);
            selected_LowPt_b_electrons.push_back(ele);
            matchedLowPtbEleGenIndexs.push_back(bMatchedGen_index);
         }
      }
   }
   std::vector<pat::Muon> selected_tau_Muons;
   std::vector<pat::Muon> selected_b_Muons;
   std::vector<pat::Muon> selected_Muons;
   std::vector<int> matchedTauMuGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level
   std::vector<int> matchedbMuGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level

   for (auto& muon: *muons){
      h_nMuons->Fill(.5,genWeight);
      if (muon::isLooseMuon(muon)==false){continue;}////////////Loose Muon ID
      h_nMuons->Fill(1.5,genWeight);
      if((muon.pt())<3|| muon.eta()>2.4){continue;} 
      if(h2AA::muonIsolation(muon)>.25){continue;}
      h_nMuons->Fill(2.5,genWeight);
      bool tauMumatched= false;
      float tauDr_min=.1;
      float tauMatchedGen_index=-1;
      selected_Muons.push_back(muon);
      cout<<"Muon phase1\n";
      for (unsigned int iGen=0; iGen<genTauMuons.size();iGen++){
         if (std::count(matchedTauMuGenIndexs.begin(), matchedTauMuGenIndexs.end(), iGen)==0){

            reco::GenParticle gen = *genTauMuons[iGen];
            TLorentzVector mu;
            TLorentzVector genMu;
            mu.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), muon.mass());
            genMu.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), gen.mass());
            float dr = genMu.DeltaR(mu);

            if(dr <tauDr_min){
               tauMumatched = true;
               tauDr_min=dr;
               tauMatchedGen_index=iGen;///Keep track of the gen particle that is closest to the reco candidate
            }
         }
      }
      cout<<"TauMuon phase1\n";
      if (tauMumatched==true){
         h_nMuons->Fill(3.5,genWeight);
         matchedTauMuGenIndexs.push_back(tauMatchedGen_index);

         h_TauMuonPt->Fill(muon.pt(),genWeight);
         selected_tau_Muons.push_back(muon);
         h_TauMuonIsolation->Fill(h2AA::muonIsolation(muon),genWeight);
         cout<<"TauMuon matched phase2\n";
      }
      bool bMumatched= false;
      float bDr_min=.1;
      float bMatchedGen_index=-1;
      cout<<"btest1\n";
      for (unsigned int iGen=0; iGen<genMufromB_S1.size();iGen++){
         cout<<"btest2\n";
         if (std::count(matchedbMuGenIndexs.begin(), matchedbMuGenIndexs.end(), iGen)==0){
            cout<<"btest3\n";
            reco::GenParticle gen = *genMufromB_S1[iGen];
            TLorentzVector mu;
            TLorentzVector genMu;
            mu.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), muon.mass());
            genMu.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), gen.mass());
            float dr = genMu.DeltaR(mu);
            //std::cout<<"DR = "<<to_string(dr)<<"\n";
            if(dr <bDr_min){
               //std::cout<<"BMuon matched phase1\n";
               bMumatched = true;
               bDr_min=dr;
               bMatchedGen_index=iGen;///Keep track of the gen particle that is closest to the reco candidate
            }
         }
      }
      
      if (bMumatched==true){
         cout<<"BMuon matched phase2\n";
         h_nMuons->Fill(4.5,genWeight);
         matchedbMuGenIndexs.push_back(bMatchedGen_index);

         h_Muon_bmatched_Pt->Fill(muon.pt(),genWeight);
         selected_b_Muons.push_back(muon);
         h_Muon_bmatched_Iso->Fill(h2AA::muonIsolation(muon),genWeight);
      }
   }
   std::sort(selected_std_electrons.begin(), selected_std_electrons.end(), h2AA::sortByPt<pat::Electron>);
   std::sort(selected_std_tau_electrons.begin(), selected_std_tau_electrons.end(), h2AA::sortByPt<pat::Electron>);
   std::sort(selected_std_b_electrons.begin(), selected_std_b_electrons.end(), h2AA::sortByPt<pat::Electron>);
   std::sort(selected_Muons.begin(), selected_Muons.end(), h2AA::sortByPt<pat::Muon>);
   std::sort(selected_b_Muons.begin(), selected_b_Muons.end(), h2AA::sortByPt<pat::Muon>);
   std::sort(selected_tau_Muons.begin(), selected_tau_Muons.end(), h2AA::sortByPt<pat::Muon>);
   std::sort(genEfromB_S1.begin(), genEfromB_S1.end(), h2AA::sortGenByPt);
   std::sort(genTauElectrons.begin(), genTauElectrons.end(), h2AA::sortGenByPt);
   std::sort(genTauMuons.begin(), genTauMuons.end(), h2AA::sortGenByPt);
   cout<<"Mu selection\n";
   std::vector<pat::Tau> selected_taus;
   for (auto& tau:*taus){
      h_nTaus->Fill(.5, genWeight);
      if (tau.pt()<10 || abs(tau.eta())>2.4){continue;}
      //if (tau.tauID("decayModeFinding")==false){continue;}
      h_nTaus->Fill(1.5, genWeight);
      if(genTauHad.size()==0||genNuTauHad.size()==0){continue;}
      if(DiTauHad==true){continue;}
      reco::GenParticle thad =*genTauHad[0];
      reco::GenParticle nt = *genNuTauHad[0];
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
   //std::cout<<"LowPt Electron selection\n";
   std::vector<pat::Jet> selected_jets;
   std::vector<int> matchedBIDs;
   std::vector<bool> jetLocks(jetFlavorMatching->size(),false);
   for (auto& jet:*jets){
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
         bool match=false;
         float dRMin=9999.;
         int matchedIdx=-1;
         h_bDiscriminator->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb")  + jet.bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
         for (reco::JetFlavourInfoMatchingCollection::const_iterator j  = jetFlavorMatching->begin(); j != jetFlavorMatching->end(); ++j){
            

            if( jetLocks.at(j - jetFlavorMatching->begin()) ) continue;
            if ((*j).second.getbHadrons().size()<1)continue;
            double dR = reco::deltaR( j->first->eta(), j->first->phi(), jet.eta(), jet.phi() );
            if (dR<.4&&dR<dRMin){
               match=true;
               dRMin=dR;
               matchedIdx=j - jetFlavorMatching->begin();
            }
         }
         if (match==true){
            h_bDiscriminator_m->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb")  + jet.bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
            h_JetPt_m->Fill(jet.pt(),genWeight);
            //selected_jets.push_back(jet);
            jetLocks.at(matchedIdx) = true;
         }
         else{
            h_bDiscriminator_f->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb")  + jet.bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
            h_JetPt_f->Fill(jet.pt(),genWeight);
         }
         if (jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb")  + jet.bDiscriminator("pfDeepFlavourJetTags:problepb")<.7476){continue;}

         h_nJets->Fill(9.5,genWeight);
         bool ematched=false;
         //if (selected_std_electrons.size()>0){
         //   double dR = reco::deltaR( selected_std_electrons[0].eta(), selected_std_electrons[0].phi(), jet.eta(), jet.phi());
         //   if (dR<.4){ematched=true;}
         //} 
         //bool mumatched=false;

         //if (selected_Muons.size()>0){
         //   double dR = reco::deltaR( selected_Muons[0].eta(), selected_Muons[0].phi(), jet.eta(), jet.phi());
         //   if (dR<.4){mumatched=true;}
         //} 
         //if (mumatched||ematched){continue;}
         h_nJets->Fill(10.5,genWeight);
         selected_jets.push_back(jet);
         h_JetPt->Fill(jet.pt(),genWeight);


      }
   }
   cout<<"jet selection\n";
   std::sort(selected_jets.begin(), selected_jets.end(), h2AA::sortByPt<pat::Jet>);

   
   cout<<"tau selection\n";
   for (auto& ele:selected_std_tau_electrons){
      h_TauElectronID->Fill("No ID",genWeight);
      if (ele.electronID("cutBasedElectronID-Fall17-94X-V2-veto")==1){h_TauElectronID->Fill("Veto",genWeight);}
      if (ele.electronID("cutBasedElectronID-Fall17-94X-V2-loose")==1){h_TauElectronID->Fill("Loose",genWeight);}
      if (ele.electronID("cutBasedElectronID-Fall17-94X-V2-medium")==1){h_TauElectronID->Fill("Medium",genWeight);}
      if (ele.electronID("cutBasedElectronID-Fall17-94X-V2-tight")==1){h_TauElectronID->Fill("Tight",genWeight);}

      h_TauElectronIDandIso->Fill("No ID",genWeight);
      if (ele.electronID("cutBasedElectronID-Fall17-94X-V2-veto")==2){h_TauElectronIDandIso->Fill("Veto",genWeight);}
      if (ele.electronID("cutBasedElectronID-Fall17-94X-V2-loose")==2){h_TauElectronIDandIso->Fill("Loose",genWeight);}
      if (ele.electronID("cutBasedElectronID-Fall17-94X-V2-medium")==2){h_TauElectronIDandIso->Fill("Medium",genWeight);}
      if (ele.electronID("cutBasedElectronID-Fall17-94X-V2-tight")==2){h_TauElectronIDandIso->Fill("Tight",genWeight);}

      h_TauElectronID_FULL->Fill("No ID",genWeight);
      if (ele.electronID("cutBasedElectronID-Fall17-94X-V2-veto")==7){h_TauElectronID_FULL->Fill("Veto",genWeight);}
      if (ele.electronID("cutBasedElectronID-Fall17-94X-V2-loose")==7){h_TauElectronID_FULL->Fill("Loose",genWeight);}
      if (ele.electronID("cutBasedElectronID-Fall17-94X-V2-medium")==7){h_TauElectronID_FULL->Fill("Medium",genWeight);}
      if (ele.electronID("cutBasedElectronID-Fall17-94X-V2-tight")==7){h_TauElectronID_FULL->Fill("Tight",genWeight);}
   }




//////////////////////Event Selection////////////////////////
///////////////////
   std::sort(genMufromB_S1.begin(), genMufromB_S1.end(), h2AA::sortGenByPt);
   std::sort(genEfromB_S1.begin(), genEfromB_S1.end(), h2AA::sortGenByPt);
   std::sort(genTauElectrons.begin(), genTauElectrons.end(), h2AA::sortGenByPt);
   std::sort(genTauMuons.begin(), genTauMuons.end(), h2AA::sortGenByPt);
   std::string EGLoose = "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v7";
   std::string singleMuLoose="HLT_IsoMu27_v13";
   std::string doubleMuonLoose="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v3";
   std::string MuEGLoose="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12";
   std::string DoubleElectronLoose="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v17";
   bool passTrig_EGLoose=triggerBits->accept(trigNames.triggerIndex(EGLoose));
   bool passTrig_singleMuLoose=triggerBits->accept(trigNames.triggerIndex(singleMuLoose));
   bool passTrig_doubleMuonLoose=triggerBits->accept(trigNames.triggerIndex(doubleMuonLoose));
   bool passTrig_DoubleElectronLoose=triggerBits->accept(trigNames.triggerIndex(DoubleElectronLoose));
   bool passTrig_MuEGLoose=triggerBits->accept(trigNames.triggerIndex(MuEGLoose));


//tauE TauE
   h_nEvent_tauE_tauE->Fill("Total Events",genWeight);
   if (selected_std_electrons.size()>1){
      h_nEvent_tauE_tauE->Fill("nStdTauE>1",genWeight);
      if (selected_jets.size()>0){
         h_nEvent_tauE_tauE->Fill("nbJets>0",genWeight);
         if(selected_std_electrons[0].charge()!=selected_std_electrons[1].charge()){
            h_nEvent_tauE_tauMu->Fill("LeptonCharge",genWeight);
            h_tauE_tauE_leadtauEpt->Fill(selected_std_electrons[0].pt(),genWeight);
            h_tauE_tauE_subtauEpt->Fill(selected_std_electrons[1].pt(),genWeight);
            TLorentzVector e1;
            e1.SetPtEtaPhiM(selected_std_electrons[0].pt(),selected_std_electrons[0].eta(),selected_std_electrons[0].phi(),selected_std_electrons[0].mass());
            TLorentzVector e2;
            e2.SetPtEtaPhiM(selected_std_electrons[1].pt(),selected_std_electrons[1].eta(),selected_std_electrons[1].phi(),selected_std_electrons[1].mass());
            TLorentzVector b;
            b.SetPtEtaPhiM(selected_jets[0].pt(),selected_jets[0].eta(),selected_jets[0].phi(),selected_jets[0].mass());
            float mvis = (e1+e2).M();
            float DiTauPt=(e2+e1).Pt();
            float bMass=b.M();
            float bDRe=b.DeltaR(e1);
            float eDRe=e1.DeltaR(e2);
            if (bMass<20){
               h_nEvent_tauE_tauMu->Fill("bMass<20",genWeight);
               if (mvis<17&&mvis>2){
                  h_nEvent_tauE_tauMu->Fill("2<mvis<17",genWeight);
                  for (unsigned int iTrigger = 0; iTrigger<DoubleElectron.size();iTrigger++){
                     std::string Trigger = DoubleElectron[iTrigger];
                     bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
                     if (passTrig){h_Trigger_tauE_tauE_DoubleElectron->Fill(Trigger.c_str(), genWeight);}
                  }
                  for (unsigned int iTrigger = 0; iTrigger<singleEle.size();iTrigger++){
                     std::string Trigger = singleEle[iTrigger];
                     bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
                     if (passTrig){h_Trigger_tauE_tauE_singleEle->Fill(Trigger.c_str(), genWeight);}
                  }
                  for (unsigned int iTrigger = 0; iTrigger<MuEG.size();iTrigger++){
                     std::string Trigger = MuEG[iTrigger];
                     bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
                     if (passTrig){h_Trigger_tauE_tauE_MuEG->Fill(Trigger.c_str(), genWeight);}
                  }
                  if (passTrig_DoubleElectronLoose){
                     h_tauE_TauE_Mvis_DoubleEG->Fill(mvis, genWeight);
                     h_Trigger_tauE_TauE_DiEandEG->Fill("Double Electron", genWeight);
                     h_tauE_TauE_Mvis_DoubleEG->Fill(mvis, genWeight);
                     h_tauE_TauE_bMass_DoubleEG->Fill(bMass, genWeight);
                     h_tauE_TauE_DiTauPt_DoubleEG->Fill(DiTauPt, genWeight);
                     h_tauE_TauE_bDRe_DoubleEG->Fill(bDRe, genWeight);
                     h_tauE_TauE_eDRe_DoubleEG->Fill(eDRe, genWeight);
                     h_tauE_TauE_MET_DoubleEG->Fill(MET, genWeight);
                  }
                  if (passTrig_EGLoose){
                     h_Trigger_tauE_TauE_EGandMuEG->Fill("Single Electron", genWeight);
                     h_tauE_TauE_Mvis_EG->Fill(mvis, genWeight);
                     h_tauE_TauE_bMass_EG->Fill(bMass, genWeight);
                     h_tauE_TauE_DiTauPt_EG->Fill(DiTauPt, genWeight);
                     h_tauE_TauE_bDRe_EG->Fill(bDRe, genWeight);
                     h_tauE_TauE_eDRe_EG->Fill(eDRe, genWeight);
                     h_tauE_TauE_MET_EG->Fill(MET, genWeight);
                     h_Trigger_tauE_TauE_DiEandEG->Fill("Single Electron", genWeight);
                  }
                  if (passTrig_DoubleElectronLoose||passTrig_EGLoose){
                     h_tauE_TauE_Mvis_EGorDoubleEG->Fill(mvis, genWeight);
                     h_tauE_TauE_bMass_EGorDoubleEG->Fill(bMass, genWeight);
                     h_tauE_TauE_DiTauPt_EGorDoubleEG->Fill(DiTauPt, genWeight);
                     h_tauE_TauE_bDRe_EGorDoubleEG->Fill(bDRe, genWeight);
                     h_tauE_TauE_eDRe_EGorDoubleEG->Fill(eDRe, genWeight);
                     h_tauE_TauE_MET_EGorDoubleEG->Fill(MET, genWeight);
                     h_Trigger_tauE_TauE_DiEandEG->Fill("Double or Single", genWeight);
                  }
                  h_tauE_TauE_Mvis_noTrig->Fill(mvis, genWeight);
                  h_tauE_TauE_bMass_noTrig->Fill(bMass, genWeight);
                  h_tauE_TauE_DiTauPt_noTrig->Fill(DiTauPt, genWeight);
                  h_tauE_TauE_bDRe_noTrig->Fill(bDRe, genWeight);
                  h_tauE_TauE_eDRe_noTrig->Fill(eDRe, genWeight);
                  h_tauE_TauE_MET_noTrig->Fill(MET, genWeight);

                  if (passTrig_MuEGLoose){
                     h_Trigger_tauE_TauE_EGandMuEG->Fill("Muon+EG", genWeight);
                  }
                  if (passTrig_EGLoose||passTrig_MuEGLoose){h_Trigger_tauE_TauE_EGandMuEG->Fill("Single or MuEG", genWeight);
                  }
               }
            }
         }
      }
   }


   h_nEvent_tauE_tauMu->Fill("Total Events",genWeight);
   if (selected_std_electrons.size()>0){
      h_nEvent_tauE_tauMu->Fill("Std Tau Ele>0",genWeight);
      if (selected_Muons.size()>0){
         h_nEvent_tauE_tauMu->Fill("Tau Muons>0",genWeight);
         if (selected_jets.size()>0){
            h_nEvent_tauE_tauMu->Fill("bjets>0",genWeight);
            if(selected_std_electrons[0].charge()!=selected_Muons[0].charge()){
               h_nEvent_tauE_tauMu->Fill("LeptonCharge",genWeight);
               h_tauMu_tauE_tauEpt->Fill(selected_std_electrons[0].pt(),genWeight);
               h_tauMu_tauE_tauMupt->Fill(selected_Muons[0].pt(),genWeight);
               TLorentzVector b;
               b.SetPtEtaPhiM(selected_jets[0].pt(),selected_jets[0].eta(),selected_jets[0].phi(),selected_jets[0].mass());

               TLorentzVector e;
               e.SetPtEtaPhiM(selected_std_electrons[0].pt(),selected_std_electrons[0].eta(),selected_std_electrons[0].phi(),selected_std_electrons[0].mass());
               TLorentzVector mu;
               mu.SetPtEtaPhiM(selected_Muons[0].pt(),selected_Muons[0].eta(),selected_Muons[0].phi(),selected_Muons[0].mass());
               float mvis = (e+mu).M();
               float DiTauPt=(e+mu).Pt();
               float bMass=b.M();
               float bDRmu=b.DeltaR(mu);
               float muDRe=mu.DeltaR(e);
               float muIso=h2AA::muonIsolation(selected_Muons[0]);
               if (bMass<20){
                  h_nEvent_tauE_tauMu->Fill("bMass<20",genWeight);
                  if (mvis<17&&mvis>2){
                     h_nEvent_tauE_tauMu->Fill("2<mvis<17",genWeight);
                     for (unsigned int iTrigger = 0; iTrigger<DoubleElectron.size();iTrigger++){
                        std::string Trigger = DoubleElectron[iTrigger];
                        bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
                        if (passTrig){h_Trigger_tauE_tauMu_DoubleElectron->Fill(Trigger.c_str(), genWeight);}
                     }
                     for (unsigned int iTrigger = 0; iTrigger<singleEle.size();iTrigger++){
                        std::string Trigger = singleEle[iTrigger];
                        bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
                        if (passTrig){h_Trigger_tauE_tauMu_singleEle->Fill(Trigger.c_str(), genWeight);}
                     }
                     for (unsigned int iTrigger = 0; iTrigger<MuEG.size();iTrigger++){
                        std::string Trigger = MuEG[iTrigger];
                        bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
                        if (passTrig){h_Trigger_tauE_tauMu_MuEG->Fill(Trigger.c_str(), genWeight);}
                     }
                     for (unsigned int iTrigger = 0; iTrigger<DoubleMuon.size();iTrigger++){
                        std::string Trigger = DoubleMuon[iTrigger];
                        bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
                        if (passTrig){h_Trigger_tauE_tauMu_DoubleMuon->Fill(Trigger.c_str(), genWeight);}
                     }
                     for (unsigned int iTrigger = 0; iTrigger<singleMu.size();iTrigger++){
                        std::string Trigger = singleMu[iTrigger];
                        bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
                        if (passTrig){h_Trigger_tauE_tauMu_singleMu->Fill(Trigger.c_str(), genWeight);}
                     }

                     if (passTrig_EGLoose){
                        h_tauMu_TauE_Mvis_EG->Fill(mvis,genWeight);
                        h_tauMu_TauE_bDRmu_EG->Fill(bDRmu,genWeight);
                        h_tauMu_TauE_bMass_EG->Fill(bMass,genWeight);
                        h_tauMu_TauE_DiTauPt_EG->Fill(DiTauPt,genWeight);
                        h_tauMu_TauE_muDRe_EG->Fill(muDRe,genWeight);
                        h_tauMu_TauE_muIso_EG->Fill(muIso,genWeight);
                        h_tauMu_TauE_MET_EG->Fill(MET,genWeight);
                     }
                     if (passTrig_singleMuLoose||passTrig_EGLoose){
                        h_Trigger_tauMu_TauE_MuandEG->Fill("Single Mu or Single E", genWeight);
                        h_tauMu_TauE_Mvis_EGorSingleMu->Fill(mvis,genWeight);
                        h_tauMu_TauE_bDRmu_EGorSingleMu->Fill(bDRmu,genWeight);
                        h_tauMu_TauE_bMass_EGorSingleMu->Fill(bMass,genWeight);
                        h_tauMu_TauE_DiTauPt_EGorSingleMu ->Fill(DiTauPt,genWeight);
                        h_tauMu_TauE_muDRe_EGorSingleMu->Fill(muDRe,genWeight);
                        h_tauMu_TauE_muIso_EGorSingleMu->Fill(muIso,genWeight);
                        h_tauMu_TauE_MET_EGorSingleMu->Fill(MET,genWeight);

                     }
                     if (passTrig_EGLoose||passTrig_MuEGLoose||passTrig_singleMuLoose){
                        h_Trigger_tauMu_TauE_EGandMuEG->Fill("MuEG or Single E or Single Mu", genWeight);
                        h_tauMu_TauE_Mvis_MuEGorSingleMuorSingleE->Fill(mvis,genWeight);
                        h_tauMu_TauE_bDRmu_MuEGorSingleMuorSingleE->Fill(bDRmu,genWeight);
                        h_tauMu_TauE_bMass_MuEGorSingleMuorSingleE->Fill(bMass,genWeight);
                        h_tauMu_TauE_DiTauPt_MuEGorSingleMuorSingleE->Fill(DiTauPt,genWeight);
                        h_tauMu_TauE_muDRe_MuEGorSingleMuorSingleE->Fill(muDRe,genWeight);
                        h_tauMu_TauE_muIso_MuEGorSingleMuorSingleE->Fill(muIso,genWeight);
                        h_tauMu_TauE_MET_MuEGorSingleMuorSingleE->Fill(MET,genWeight);

                     }

                     if (passTrig_singleMuLoose){
                        h_Trigger_tauMu_TauE_MuandMuEG->Fill("Single Mu", genWeight);
                        h_Trigger_tauMu_TauE_MuandEG->Fill("Single E", genWeight);
                        h_tauMu_TauE_bDRmu_SingleMu->Fill(bDRmu,genWeight);
                        h_tauMu_TauE_Mvis_SingleMu->Fill(mvis,genWeight);
                        h_tauMu_TauE_bMass_SingleMu->Fill(bMass,genWeight);
                        h_tauMu_TauE_DiTauPt_SingleMu->Fill(DiTauPt,genWeight);
                        h_tauMu_TauE_muDRe_SingleMu->Fill(muDRe,genWeight);
                        h_tauMu_TauE_muIso_SingleMu->Fill(muIso,genWeight);
                        h_tauMu_TauE_MET_SingleMu->Fill(MET,genWeight);
                     }
                     if (passTrig_MuEGLoose){
                        h_Trigger_tauMu_TauE_MuandMuEG->Fill("MuEG", genWeight);
                        h_tauMu_TauE_Mvis_MuEG->Fill(mvis,genWeight);
                        h_tauMu_TauE_DiTauPt_MuEG->Fill(DiTauPt,genWeight);
                        h_tauMu_TauE_bMass_MuEG->Fill(bMass,genWeight);
                        h_tauMu_TauE_bDRmu_MuEG->Fill(bDRmu,genWeight);
                        h_tauMu_TauE_muDRe_MuEG->Fill(muDRe,genWeight);
                        h_tauMu_TauE_muIso_MuEG->Fill(muIso,genWeight);
                        h_tauMu_TauE_MET_MuEG->Fill(MET,genWeight);
                     }
                     if (passTrig_singleMuLoose||passTrig_MuEGLoose){
                        h_Trigger_tauMu_TauE_MuandMuEG->Fill("Mu EG or single Mu", genWeight);
                        h_tauMu_TauE_Mvis_MuEGorSingleMu->Fill(mvis,genWeight);
                        h_tauMu_TauE_DiTauPt_MuEGorSingleMu->Fill(DiTauPt,genWeight);
                        h_tauMu_TauE_bMass_MuEGorSingleMu->Fill(bMass,genWeight);
                        h_tauMu_TauE_bDRmu_MuEGorSingleMu->Fill(bDRmu,genWeight);
                        h_tauMu_TauE_muDRe_MuEGorSingleMu->Fill(muDRe,genWeight);
                        h_tauMu_TauE_muIso_MuEGorSingleMu->Fill(muIso,genWeight);
                        h_tauMu_TauE_MET_MuEGorSingleMu->Fill(MET,genWeight);
                     }
                     h_tauMu_TauE_Mvis_noTrig->Fill(mvis,genWeight);
                     h_tauMu_TauE_bMass_noTrig->Fill(bMass,genWeight);
                     h_tauMu_TauE_DiTauPt_noTrig->Fill(DiTauPt,genWeight);
                     h_tauMu_TauE_bDRmu_noTrig->Fill(bDRmu,genWeight);
                     h_tauMu_TauE_muDRe_noTrig->Fill(muDRe,genWeight);
                     h_tauMu_TauE_muIso_noTrig->Fill(muIso,genWeight);
                     h_tauMu_TauE_MET_noTrig->Fill(MET,genWeight);
                  }
               }
            }
                  
         }
      }
   }

   
   if (passTrig_DoubleElectronLoose){h_Trigger_general_EGandDiEG->Fill("Double Electron", genWeight);}
   if (passTrig_EGLoose){h_Trigger_general_EGandDiEG->Fill("Single Electron", genWeight);}
   if (passTrig_EGLoose||passTrig_DoubleElectronLoose){h_Trigger_general_EGandDiEG->Fill("Single or Double Electron", genWeight);}

   if (passTrig_EGLoose){h_Trigger_general_EGandMuEG->Fill("Single Electron", genWeight);}
   if (passTrig_MuEGLoose){h_Trigger_general_EGandMuEG->Fill("Muon EG", genWeight);}
   if (passTrig_EGLoose||passTrig_MuEGLoose){h_Trigger_general_EGandMuEG->Fill("Single Ele or MuonEG", genWeight);}

   if (passTrig_singleMuLoose){h_Trigger_general_MuandMuEG->Fill("Single Muon", genWeight);}
   if (passTrig_MuEGLoose){h_Trigger_general_MuandMuEG->Fill("Muon EG", genWeight);}
   if (passTrig_singleMuLoose||passTrig_MuEGLoose){h_Trigger_general_MuandMuEG->Fill("Single Mu or Muon EG", genWeight);}

   if (passTrig_singleMuLoose){h_Trigger_general_MuandDiMu->Fill("Single Muon", genWeight);}
   if (passTrig_doubleMuonLoose){h_Trigger_general_MuandDiMu->Fill("DiMuon", genWeight);}
   if (passTrig_singleMuLoose||passTrig_doubleMuonLoose){h_Trigger_general_MuandDiMu->Fill("Single Muon or DiMuon", genWeight);}

   h_nEvent_tauMu_tauMu->Fill("Total Events",genWeight);
   if (selected_Muons.size()>1){
      h_nEvent_tauMu_tauMu->Fill("Std Tau Mu>1",genWeight);
      if (selected_jets.size()>0){
         TLorentzVector b;
         b.SetPtEtaPhiM(selected_jets[0].pt(),selected_jets[0].eta(),selected_jets[0].phi(),selected_jets[0].mass());

         h_nEvent_tauMu_tauMu->Fill("Njets>0",genWeight);
            if(selected_Muons[1].charge()!=selected_Muons[0].charge()){
            h_nEvent_tauE_tauMu->Fill("LeptonCharge",genWeight);
            h_tauMu_tauMu_leadtauMupt->Fill(selected_Muons[0].pt(),genWeight);
            h_tauMu_tauMu_subtauMupt->Fill(selected_Muons[1].pt(),genWeight);
            TLorentzVector mu1;
            mu1.SetPtEtaPhiM(selected_Muons[0].pt(),selected_Muons[0].eta(),selected_Muons[0].phi(),selected_Muons[0].mass());
            TLorentzVector mu2;
            mu2.SetPtEtaPhiM(selected_Muons[1].pt(),selected_Muons[1].eta(),selected_Muons[1].phi(),selected_Muons[1].mass());
            float mvis = (mu1+mu2).M();
            float DiTauPt=(mu2+mu1).Pt();
            float bMass=b.M();
            float bDRmu=b.DeltaR(mu1);
            float muDRmu=mu1.DeltaR(mu2);
            float muIso=h2AA::muonIsolation(selected_Muons[0]);

            if (bMass<20){
               h_nEvent_tauMu_tauMu->Fill("bMass<20",genWeight);
               if (mvis<17&&mvis>2){
                  h_nEvent_tauMu_tauMu->Fill("2<mvis<17",genWeight);
                  for (unsigned int iTrigger = 0; iTrigger<DoubleMuon.size();iTrigger++){
                     std::string Trigger = DoubleMuon[iTrigger];
                     bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
                     if (passTrig){h_Trigger_tauMu_tauMu_DoubleMuon->Fill(Trigger.c_str(), genWeight);}
                  }
                  for (unsigned int iTrigger = 0; iTrigger<singleMu.size();iTrigger++){
                     std::string Trigger = singleMu[iTrigger];
                     bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
                     if (passTrig){h_Trigger_tauMu_tauMu_singleMu->Fill(Trigger.c_str(), genWeight);}
                  }
                  if (passTrig_MuEGLoose){
                     h_Trigger_tauMu_TauMu_MuandMuEG->Fill("Muon EGamma", genWeight);
                  }
                  if (passTrig_singleMuLoose||passTrig_MuEGLoose){
                     h_Trigger_tauMu_TauMu_MuandMuEG->Fill("Single Mu OR MuEG", genWeight);
                  }


                  if (passTrig_doubleMuonLoose){
                     h_Trigger_tauMu_TauMu_DiMuandMu->Fill("Double Muon", genWeight);
                     h_tauMu_TauMu_Mvis_DoubleMu->Fill(mvis,genWeight);
                     h_tauMu_TauMu_bMass_DoubleMu->Fill(bMass,genWeight);
                     h_tauMu_TauMu_bDRmu_DoubleMu->Fill(bDRmu,genWeight);
                     h_tauMu_TauMu_DiTauPt_DoubleMu->Fill(DiTauPt,genWeight);
                     h_tauMu_TauMu_muDRmu_DoubleMu->Fill(muDRmu,genWeight);
                     h_tauMu_TauMu_MET_DoubleMu->Fill(MET,genWeight);
                     h_tauMu_TauMu_muIso_DoubleMu->Fill(muIso,genWeight);

                  }
                  if (passTrig_singleMuLoose){
                     h_Trigger_tauMu_TauMu_DiMuandMu->Fill("Single Muon", genWeight);
                     h_Trigger_tauMu_TauMu_MuandMuEG->Fill("Single Muon", genWeight);
                     h_tauMu_TauMu_Mvis_SingleMu->Fill(mvis,genWeight);
                     h_tauMu_TauMu_bMass_SingleMu->Fill(bMass,genWeight);
                     h_tauMu_TauMu_bDRmu_SingleMu->Fill(bDRmu,genWeight);
                     h_tauMu_TauMu_DiTauPt_SingleMu->Fill(DiTauPt,genWeight);
                     h_tauMu_TauMu_muDRmu_SingleMu->Fill(muDRmu,genWeight);
                     h_tauMu_TauMu_MET_SingleMu->Fill(MET,genWeight);
                     h_tauMu_TauMu_muIso_SingleMu->Fill(muIso,genWeight);


                  }
                  if (passTrig_doubleMuonLoose||passTrig_singleMuLoose){
                     h_Trigger_tauMu_TauMu_DiMuandMu->Fill("Single OR Double Muon", genWeight);
                     h_tauMu_TauMu_Mvis_SingleMuorDoubleMu->Fill(mvis,genWeight);
                     h_tauMu_TauMu_bMass_SingleMuorDoubleMu->Fill(bMass,genWeight);
                     h_tauMu_TauMu_bDRmu_SingleMuorDoubleMu->Fill(bDRmu,genWeight);
                     h_tauMu_TauMu_DiTauPt_SingleMuorDoubleMu->Fill(DiTauPt,genWeight);
                     h_tauMu_TauMu_muDRmu_SingleMuorDoubleMu->Fill(muDRmu,genWeight);
                     h_tauMu_TauMu_MET_SingleMuorDoubleMu->Fill(MET,genWeight);
                     h_tauMu_TauMu_muIso_SingleMuorDoubleMu->Fill(muIso,genWeight);

                  }
                  h_tauMu_TauMu_Mvis_noTrig->Fill(mvis,genWeight);
                  h_tauMu_TauMu_bMass_noTrig->Fill(bMass,genWeight);
                  h_tauMu_TauMu_bDRmu_noTrig->Fill(bDRmu,genWeight);
                  h_tauMu_TauMu_DiTauPt_noTrig->Fill(DiTauPt,genWeight);
                  h_tauMu_TauMu_muDRmu_noTrig->Fill(muDRmu,genWeight);
                  h_tauMu_TauMu_MET_noTrig->Fill(MET,genWeight);
                  h_tauMu_TauMu_muIso_noTrig->Fill(muIso,genWeight);
               }
            }
         }
      }
   }

///////////End
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
