// -*- C++ -*-
//
// Package:    higgstoaaAnalyzer_RECO/higgstoaaAnalyzer_RECO
// Class:      higgstoaaAnalyzer_RECO
//
/**\class higgstoaaAnalyzer_RECO higgstoaaAnalyzer_RECO.cc higgstoaaAnalyzer_RECO/higgstoaaAnalyzer_RECO/plugins/higgstoaaAnalyzer_RECO.cc

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

#include "higgstoaaAnalyzer_RECO.h"
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

#include "DataFormats/Candidate/interface/Candidate.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.



class higgstoaaAnalyzer_RECO : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit higgstoaaAnalyzer_RECO(const edm::ParameterSet&);
      ~higgstoaaAnalyzer_RECO();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      const edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticles_;
      const edm::EDGetTokenT< std::vector<reco::GsfElectron> > std_electrons_;
      const edm::EDGetTokenT< std::vector<reco::GsfElectron> > lowPt_electrons_;
      const edm::EDGetTokenT< edm::ValueMap<float> >  mvaIds_;
      const edm::EDGetTokenT<std::vector<reco::Muon>> muons_;
      const edm::EDGetTokenT<std::vector<reco::PFJet>> jets_;
      const edm::EDGetTokenT<std::vector<reco::PFMET>> mets_;
      const edm::EDGetTokenT<std::vector<reco::PFTau>> taus_;
      const edm::EDGetTokenT<double> eventrho_;
      const edm::EDGetTokenT<GenEventInfoProduct> genInfoProduct_;
      const edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      
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
      TH1F *h_TauPt;
      TH1F *h_GenTauPt;
      TH1F *h_GenTauMuonPt;
      TH1F *h_TauMuonPt;
      
      TH1F *h_TauMuonIsolation;
      TH1F *h_MET; 
      TH1F *h_nStdElectrons;
      TH1F *h_nLowPtElectrons;

      TH1F *h_nJets;
      TH1F *h_nMuons;
      TH1F *h_nTaus;




      TH1F *h_DoubleElectron;


      TH1F *h_singleE;
      TH1F *h_singleMu;
      
      TH1F *h_DoubleMuon;
      TH1F *h_MuonEGamma;

//B Daughter studies
      TH1F *h_nEStatus1_NotfromTau_Gen;
      TH1F *h_nMuStatus1_NotfromTau_Gen;
      TH1F *h_nEMotherB;
      TH1F *h_nMuMotherB;
      TH1F *h_MunStepsToFinal;
      TH1F *h_EnStepsToFinal;
      TH1F *h_EnBinChain;
      TH1F *h_MunBinChain;
      TH1F *h_ELastMother_b;
      TH1F *h_ELastMother_Nonb;
      TH1F *h_MuLastMother_Nonb;
      TH1F *h_bE_notfromA_Pt;
      TH1F *h_bMu_notfromA_Pt;

      TH1F *h_MuLastMother_b;

      TH1F *h_genS1Muon_NonPrimary_Pt;
      TH1F *h_genS1Electrons_NonPrimary_Pt;

      TH1F *h_nESteps_afterB;
      TH1F *h_nMuSteps_afterB;
      TH1F *h_bDrMu;
      TH1F *h_bDrE;

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
      TH1F *h_nEvent_bE_tauE;
      TH1F *h_nEvent_bMu_tauE;
      TH1F *h_nEvent_bMu_tauMu;
      TH1F *h_nEvent_bE_tauMu;

      TH1F *h_bE_tauE_gentauEpt;
      TH1F *h_bE_tauE_genbEpt;

      TH1F *h_bMu_tauE_gentauEpt;
      TH1F *h_bMu_tauE_genbMupt;

      TH1F *h_bE_tauMu_gentauMupt;
      TH1F *h_bE_tauMu_genbEpt;

      TH1F *h_bMu_tauMu_gentauMupt;
      TH1F *h_bMu_tauMu_genbMupt;


////Trigger plots
      TH1F *h_Trigger_bE_tauE_DoubleElectron;
      TH1F *h_Trigger_bE_tauE_singleEle;
      TH1F *h_Trigger_bE_tauE_MuEG;


      TH1F *h_Trigger_bE_tauMu_singleMu;
      TH1F *h_Trigger_bE_tauMu_singleEle;
      TH1F *h_Trigger_bE_tauMu_MuEG;
      TH1F *h_Trigger_bE_tauMu_DoubleElectron;
      TH1F *h_Trigger_bE_tauMu_DoubleMuon;


      TH1F *h_Trigger_bMu_tauE_singleMu;
      TH1F *h_Trigger_bMu_tauE_singleEle;
      TH1F *h_Trigger_bMu_tauE_MuEG;
      TH1F *h_Trigger_bMu_tauE_DoubleElectron;
      TH1F *h_Trigger_bMu_tauE_DoubleMuon;

      TH1F *h_Trigger_bMu_tauMu_DoubleMuon;
      TH1F *h_Trigger_bMu_tauMu_singleMu;
      TH1F *h_Trigger_bMu_tauMu_MuEG;
//Trigger improvement checks
      TH1F *h_Trigger_bE_TauE_DiEandEG;
      TH1F *h_Trigger_bE_TauE_EGandMuEG;

      TH1F *h_Trigger_bMu_TauE_MuandEG;
      TH1F *h_Trigger_bMu_TauE_EGandMuEG;
      TH1F *h_Trigger_bMu_TauE_MuandMuEG;

      TH1F *h_Trigger_bE_TauMu_MuandEG;
      TH1F *h_Trigger_bE_TauMu_EGandMuEG;
      TH1F *h_Trigger_bE_TauMu_MuandMuEG;

      TH1F *h_Trigger_bMu_TauMu_DiMuandMu;
      TH1F *h_Trigger_bMu_TauMu_MuandMuEG;

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
higgstoaaAnalyzer_RECO::higgstoaaAnalyzer_RECO(const edm::ParameterSet& iConfig) :
   genParticles_{consumes<std::vector<reco::GenParticle> >(edm::InputTag(std::string("genParticles")))},
   std_electrons_{consumes<std::vector<reco::GsfElectron>> (edm::InputTag(std::string("gedGsfElectrons")))},
   lowPt_electrons_{consumes<std::vector<reco::GsfElectron>> (edm::InputTag(std::string("lowPtGsfElectrons")))},
   mvaIds_{ consumes< edm::ValueMap<float>>(edm::InputTag(std::string("lowPtGsfElectronID")))},
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
      
      h_nLowPtElectrons=fs->make<TH1F>("h_nLowPtElectrons", "StdElectron Cuts", 5,0,5);
      h_nLowPtElectrons->Sumw2();
      h_nMuons=fs->make<TH1F>("h_nMuons", "Muon Cuts", 5,0,5);
      h_nMuons->Sumw2();
      h_nTaus=fs->make<TH1F>("h_nTaus", "Tau Cuts", 3,0,3);
      h_nTaus->Sumw2();
      h_nJets=fs->make<TH1F>("h_nJets", "Jet Cuts", 9,0,9);
      h_nJets->Sumw2();
///Event Selection Plots
      h_nEvent_bE_tauE=fs->make<TH1F>("h_nEvent_bE_tauE", "NEvent BE TauE",3,0,3);
      h_nEvent_bE_tauE->Sumw2();

      h_nEvent_bMu_tauE=fs->make<TH1F>("h_nEvent_bMu_tauE", "NEvent BMu TauE",3,0,3);
      h_nEvent_bMu_tauE->Sumw2();

      h_nEvent_bE_tauMu=fs->make<TH1F>("h_nEvent_bE_tauMu", "NEvent BE TauMu",3,0,3);
      h_nEvent_bE_tauMu->Sumw2();

      h_nEvent_bMu_tauMu=fs->make<TH1F>("h_nEvent_bMu_tauMu", "NEvent BMu TauMu",3,0,3);
      h_nEvent_bMu_tauMu->Sumw2();
///Event PT Plots

      h_bE_tauE_gentauEpt=fs->make<TH1F>("h_bE_tauE_gentauEpt", "bE|TauE Gen Tau Ele PT (GeV)", 45,0,90);
      h_bE_tauE_gentauEpt->Sumw2();

      h_bE_tauE_genbEpt=fs->make<TH1F>("h_bE_tauE_genbEpt", "bE|TauE Gen b Ele PT (GeV)", 45,0,90);
      h_bE_tauE_genbEpt->Sumw2();


      h_bMu_tauE_gentauEpt=fs->make<TH1F>("h_bMu_tauE_gentauEpt", "bMu|TauE Gen b Muon PT (GeV)", 45,0,90);
      h_bMu_tauE_gentauEpt->Sumw2();

      h_bMu_tauE_genbMupt=fs->make<TH1F>("h_bMu_tauE_genbMupt", "bMu|TauE Gen Ele PT (GeV)", 45,0,90);
      h_bMu_tauE_genbMupt->Sumw2();

      h_bE_tauMu_gentauMupt=fs->make<TH1F>("h_bE_tauMu_gentauMupt", "bE|TauMu Gen Tau Muon PT (GeV)", 45,0,90);
      h_bE_tauMu_gentauMupt->Sumw2();

      h_bE_tauMu_genbEpt=fs->make<TH1F>("h_bE_tauMu_genbEpt", "bE|TauMu Gen b Ele PT (GeV)", 45,0,90);
      h_bE_tauMu_genbEpt->Sumw2();

      h_bMu_tauMu_gentauMupt=fs->make<TH1F>("h_bMu_tauMu_gentauMupt", "bMu|TauMu Gen Tau Muon PT (GeV)", 45,0,90);
      h_bMu_tauMu_gentauMupt->Sumw2();

      h_bMu_tauMu_genbMupt=fs->make<TH1F>("h_bMu_tauMu_genbMupt", "bMu|TauMu Gen b Muon PT (GeV)", 45,0,90);
      h_bMu_tauMu_genbMupt->Sumw2();

//////// General Trigger plots
      h_DoubleElectron=fs->make<TH1F>("h_DoubleElectron", "Double Electron Triggers",3,0,3);
      h_DoubleElectron->SetCanExtend(TH1::kAllAxes);
      h_DoubleElectron->Sumw2();

      h_singleE=fs->make<TH1F>("h_singleE", "General Single Ele Triggers",3,0,3);
      h_singleE->SetCanExtend(TH1::kAllAxes);
      h_singleE->Sumw2();

      h_singleMu=fs->make<TH1F>("h_singleMu", "General Single Mu Triggers",3,0,3);
      h_singleMu->SetCanExtend(TH1::kAllAxes);
      h_singleMu->Sumw2();

      h_DoubleMuon=fs->make<TH1F>("h_DoubleMuon", "General Double Mu Triggers",3,0,3);
      h_DoubleMuon->SetCanExtend(TH1::kAllAxes);
      h_DoubleMuon->Sumw2();

      h_MuonEGamma=fs->make<TH1F>("h_MuonEGamma", "General MuEgamma Triggers",3,0,3);
      h_MuonEGamma->SetCanExtend(TH1::kAllAxes);
      h_MuonEGamma->Sumw2();




///////B Parantage Study


      h_nEStatus1_NotfromTau_Gen=fs->make<TH1F>("h_nEStatus1_NotfromTau_Gen", "N Electrons Status 1 (!fromtau) (Gen)", 8,0,8);
      h_nEStatus1_NotfromTau_Gen->Sumw2(); 

      h_nMuStatus1_NotfromTau_Gen=fs->make<TH1F>("h_nMuStatus1_NotfromTau_Gen", "N Muons Status (!fromtau) 1 (Gen)", 8,0,8);
      h_nMuStatus1_NotfromTau_Gen->Sumw2(); 

      
      h_EnStepsToFinal=fs->make<TH1F>("h_EnStepsToFinal", "N E steps to Final", 12,0,12);
      h_EnStepsToFinal->Sumw2(); 


      h_MunStepsToFinal=fs->make<TH1F>("h_MunStepsToFinal", "N Mu steps to Final", 12,0,12);
      h_MunStepsToFinal->Sumw2(); 


      h_ELastMother_b=fs->make<TH1F>("h_ELastMother_b", "E Last Mother PDG ID", 37,0,37);
      h_ELastMother_b->Sumw2();

      h_EnBinChain=fs->make<TH1F>("h_EnBinChain", "Electrons: n b in chain", 12,0,12);
      h_EnBinChain->Sumw2(); 

      h_MunBinChain=fs->make<TH1F>("h_MunBinChain", "Muons: n b in chain", 12,0,12);
      h_MunBinChain->Sumw2(); 

      h_MuLastMother_b=fs->make<TH1F>("h_MuLastMother_b", "Mu Last Mother PDG ID", 37,0,37);
      h_MuLastMother_b->Sumw2();

      h_ELastMother_Nonb=fs->make<TH1F>("h_ELastMother_Nonb", "Ele Last Mother PDG ID(b)", 37,0,37);
      h_ELastMother_Nonb->Sumw2();



      h_MuLastMother_Nonb=fs->make<TH1F>("h_MuLastMother_Nonb", "Mu Last Mother PDG ID(b)", 37,0,37);
      h_MuLastMother_Nonb->Sumw2();


      h_nEMotherB=fs->make<TH1F>("h_nEMotherB", "N E with b Mother Per Event", 8,0,8);
      h_nEMotherB->Sumw2(); 

      h_nMuMotherB=fs->make<TH1F>("h_nMuMotherB", "N Mu with b Mother Per Event", 5,0,5);
      h_nMuMotherB->Sumw2(); 


      h_nESteps_afterB=fs->make<TH1F>("h_nESteps_afterB", "Ele Number of steps after encountering first b to end", 12,0,12);
      h_nESteps_afterB->Sumw2(); 

      h_nMuSteps_afterB=fs->make<TH1F>("h_nMuSteps_afterB", "Mu Number of steps after encountering first b to end", 12,0,12);
      h_nMuSteps_afterB->Sumw2(); 

      h_genS1Electrons_NonPrimary_Pt=fs->make<TH1F>("h_genS1Electrons_NonPrimary_Pt", "Gen S1 Electron(!Tau!b) Pt (GeV)", 45,0,90);
      h_genS1Electrons_NonPrimary_Pt->Sumw2();

      h_genS1Muon_NonPrimary_Pt=fs->make<TH1F>("h_genS1Muon_NonPrimary_Pt", "Gen S1 Muons(!Tau!b) PT (GeV)", 45,0,90);
      h_genS1Muon_NonPrimary_Pt->Sumw2();

      h_bDrMu=fs->make<TH1F>("h_bDrMu", "Gen DeltaR (last b,Muon) [daughter]", 100,0,1);
      h_bDrMu->Sumw2();

      h_bDrE=fs->make<TH1F>("h_bDrE", "Gen DeltaR (last b,Ele) [daughter]", 100,0,1);
      h_bDrE->Sumw2();

      h_genbMu_Pt=fs->make<TH1F>("h_genbMu_Pt", "Gen S1 bMuons PT from A (GeV)", 45,0,90);
      h_genbMu_Pt->Sumw2();

      h_genbE_Pt=fs->make<TH1F>("h_genbE_Pt", "Gen S1 bElectrons from A PT (GeV)", 45,0,90);
      h_genbE_Pt->Sumw2();

      h_bMu_notfromA_Pt=fs->make<TH1F>("h_bMu_notfromA_Pt", "Gen S1 bMuons PT !from A (GeV)", 45,0,90);
      h_bMu_notfromA_Pt->Sumw2();

      h_bE_notfromA_Pt=fs->make<TH1F>("h_bE_notfromA_Pt", "Gen S1 bElectrons !from A PT (GeV)", 45,0,90);
      h_bE_notfromA_Pt->Sumw2();

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

      h_Trigger_bE_tauE_DoubleElectron=fs->make<TH1F>("h_Trigger_bE_tauE_DoubleElectron", "bE|TauE Double Electron Triggers",3,0,3);
      h_Trigger_bE_tauE_DoubleElectron->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bE_tauE_DoubleElectron->Sumw2();

      h_Trigger_bE_tauE_singleEle=fs->make<TH1F>("h_Trigger_bE_tauE_singleEle", "bE|TauE Single Ele Triggers",3,0,3);
      h_Trigger_bE_tauE_singleEle->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bE_tauE_singleEle->Sumw2();

      h_Trigger_bE_tauE_MuEG=fs->make<TH1F>("h_Trigger_bE_tauE_MuEG", "bE|TauE MuEgamma Triggers",3,0,3);
      h_Trigger_bE_tauE_MuEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bE_tauE_MuEG->Sumw2();



//Triger bETauMu
      h_Trigger_bE_tauMu_DoubleElectron=fs->make<TH1F>("h_Trigger_bE_tauMu_DoubleElectron", "bE|TauMuDouble Electron Triggers",3,0,3);
      h_Trigger_bE_tauMu_DoubleElectron->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bE_tauMu_DoubleElectron->Sumw2();

      h_Trigger_bE_tauMu_singleEle=fs->make<TH1F>("h_Trigger_bE_tauMu_singleEle", "bE|TauMu Single Ele Triggers",3,0,3);
      h_Trigger_bE_tauMu_singleEle->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bE_tauMu_singleEle->Sumw2();

      h_Trigger_bE_tauMu_singleMu=fs->make<TH1F>("h_Trigger_bE_tauMu_singleMu", "bE|TauMu Single Mu Triggers",3,0,3);
      h_Trigger_bE_tauMu_singleMu->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bE_tauMu_singleMu->Sumw2();

      h_Trigger_bE_tauMu_DoubleMuon=fs->make<TH1F>("h_Trigger_bE_tauMu_DoubleMuon", "bE|TauMu Double Mu Triggers",3,0,3);
      h_Trigger_bE_tauMu_DoubleMuon->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bE_tauMu_DoubleMuon->Sumw2();

      h_Trigger_bE_tauMu_MuEG=fs->make<TH1F>("h_Trigger_bE_tauMu_MuEG", "bE|TauMu MuEgamma Triggers",3,0,3);
      h_Trigger_bE_tauMu_MuEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bE_tauMu_MuEG->Sumw2();

//Triger bMuTauE
      h_Trigger_bMu_tauE_DoubleElectron=fs->make<TH1F>("h_Trigger_bMu_tauE_DoubleElectron", "bMu|TauEDouble Electron Triggers",3,0,3);
      h_Trigger_bMu_tauE_DoubleElectron->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bMu_tauE_DoubleElectron->Sumw2();

      h_Trigger_bMu_tauE_singleEle=fs->make<TH1F>("h_Trigger_bMu_tauE_singleEle", "bMu|TauE Single Ele Triggers",3,0,3);
      h_Trigger_bMu_tauE_singleEle->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bMu_tauE_singleEle->Sumw2();

      h_Trigger_bMu_tauE_singleMu=fs->make<TH1F>("h_Trigger_bMu_tauE_singleMu", "bMu|TauE Single Mu Triggers",3,0,3);
      h_Trigger_bMu_tauE_singleMu->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bMu_tauE_singleMu->Sumw2();

      h_Trigger_bMu_tauE_DoubleMuon=fs->make<TH1F>("h_Trigger_bMu_tauE_DoubleMuon", "bMu|TauE Double Mu Triggers",3,0,3);
      h_Trigger_bMu_tauE_DoubleMuon->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bMu_tauE_DoubleMuon->Sumw2();

      h_Trigger_bMu_tauE_MuEG=fs->make<TH1F>("h_Trigger_bMu_tauE_MuEG", "bMu|TauE MuEgamma Triggers",3,0,3);
      h_Trigger_bMu_tauE_MuEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bMu_tauE_MuEG->Sumw2();

//Trigger bmuTauMu

      h_Trigger_bMu_tauMu_singleMu=fs->make<TH1F>("h_Trigger_bMu_tauMu_singleMu", "bMu|TauMu Single Mu Triggers",3,0,3);
      h_Trigger_bMu_tauMu_singleMu->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bMu_tauMu_singleMu->Sumw2();

      h_Trigger_bMu_tauMu_DoubleMuon=fs->make<TH1F>("h_Trigger_bMu_tauMu_DoubleMuon", "bMu|TauMu Double Mu Triggers",3,0,3);
      h_Trigger_bMu_tauMu_DoubleMuon->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bMu_tauMu_DoubleMuon->Sumw2();

      h_Trigger_bMu_tauMu_MuEG=fs->make<TH1F>("h_Trigger_bMu_tauMu_MuEG", "bMu|TauMu MuEgamma Triggers",3,0,3);
      h_Trigger_bMu_tauMu_MuEG->SetCanExtend(TH1::kAllAxes);
      h_Trigger_bMu_tauMu_MuEG->Sumw2();
//Trigger improvment
      h_Trigger_bE_TauE_DiEandEG=fs->make<TH1F>("h_Trigger_bE_TauE_DiEandEG", "bE|TauE DiEandEG Triggers",3,0,3);
      h_Trigger_bE_TauE_DiEandEG->Sumw2();


      h_Trigger_bE_TauE_EGandMuEG=fs->make<TH1F>("h_Trigger_bE_TauE_EGandMuEG", "bE|TauE EGandMuEG Triggers",3,0,3);
      h_Trigger_bE_TauE_EGandMuEG->Sumw2();

      h_Trigger_bMu_TauE_MuandEG=fs->make<TH1F>("h_Trigger_bMu_TauE_MuandEG", "bMu|TauE MuandEG Triggers",3,0,3);
      h_Trigger_bMu_TauE_MuandEG->Sumw2();

      h_Trigger_bMu_TauE_EGandMuEG=fs->make<TH1F>("h_Trigger_bMu_TauE_EGandMuEG", "bMu|TauE EGandMuEG Triggers",3,0,3);
      h_Trigger_bMu_TauE_EGandMuEG->Sumw2();

      h_Trigger_bMu_TauE_MuandMuEG=fs->make<TH1F>("h_Trigger_bMu_TauE_MuandMuEG", "bMu|TauE MuandMuEG Triggers",3,0,3);
      h_Trigger_bMu_TauE_MuandMuEG->Sumw2();


      h_Trigger_bE_TauMu_MuandEG=fs->make<TH1F>("h_Trigger_bE_TauMu_MuandEG", "bE|TauMu MuandEG Triggers",3,0,3);
      h_Trigger_bE_TauMu_MuandEG->Sumw2();

      h_Trigger_bE_TauMu_EGandMuEG=fs->make<TH1F>("h_Trigger_bE_TauMu_EGandMuEG", "bE|TauMu EGandMuEG Triggers",3,0,3);
      h_Trigger_bE_TauMu_EGandMuEG->Sumw2();

      h_Trigger_bE_TauMu_MuandMuEG=fs->make<TH1F>("h_Trigger_bE_TauMu_MuandMuEG", "bE|TauMu MuandMuEG Triggers",3,0,3);
      h_Trigger_bE_TauMu_MuandMuEG->Sumw2();



      h_Trigger_bMu_TauMu_DiMuandMu=fs->make<TH1F>("h_Trigger_bMu_TauMu_DiMuandMu", "bMu|TauMu DiMuandMu Triggers",3,0,3);
      h_Trigger_bMu_TauMu_DiMuandMu->Sumw2();

      h_Trigger_bMu_TauMu_MuandMuEG=fs->make<TH1F>("h_Trigger_bMu_TauMu_MuandMuEG", "bMu|TauMu MuandMuEG Triggers",3,0,3);
      h_Trigger_bMu_TauMu_MuandMuEG->Sumw2();

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


higgstoaaAnalyzer_RECO::~higgstoaaAnalyzer_RECO()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
higgstoaaAnalyzer_RECO::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;

   edm::Handle<std::vector<reco::PFTau>> taus;
	iEvent.getByToken(taus_, taus);
   edm::Handle<std::vector<reco::PFMET>> mets;
	iEvent.getByToken(mets_, mets);
   edm::Handle<std::vector<reco::Muon>> muons;
	iEvent.getByToken(muons_, muons);
   edm::Handle<std::vector<reco::PFJet>> jets;
	iEvent.getByToken(jets_, jets);
   edm::Handle<std::vector<reco::GsfElectron>> std_electrons;
	iEvent.getByToken(std_electrons_, std_electrons);
   edm::Handle<std::vector<reco::GsfElectron>> lowPt_electrons;
	iEvent.getByToken(lowPt_electrons_, lowPt_electrons);
   //std::vector<reco::GenParticle> genParticles = ::getObject<std::vector<reco::GenParticle>>(genParticles_, iEvent);
   edm::Handle<std::vector<reco::GenParticle>> genParticles;
	iEvent.getByToken(genParticles_, genParticles);
   edm::Handle< edm::ValueMap<float> > mvaIds;
   try { iEvent.getByToken(mvaIds_, mvaIds); }
   catch (...) {;}
   edm::ValueMap<float> EleMVAID = *mvaIds.product();
///Get Trigger objects for out handles
   edm::InputTag trigResultsTag("TriggerResults","","HLT"); //make sure have correct process on MC
   edm::Handle<edm::TriggerResults>  triggerBits;
   iEvent.getByLabel(trigResultsTag,triggerBits);

   edm::TriggerResults triggerResults = *triggerBits.product();
   cout<<to_string(triggerResults.size());
   const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerBits);  

   //for (unsigned int i = 0; i<trigNames.size();i++){cout<< to_string(i)<<" "<< trigNames.triggerName(i)<<" "<< triggerResults.accept(i)<<"\n";}
   double rho = h2AA_RECO::getObject_RECO<double>(eventrho_,iEvent);
   edm::Handle<GenEventInfoProduct>  genInfoHandle;
   try{iEvent.getByToken(genInfoProduct_, genInfoHandle);}
   catch(...){;}
   double genWeight=genInfoHandle->weight();
   reco::PFMET met=(*mets.product())[0];
   h_MET->Fill(met.pt(),genWeight);
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
      if (isfromb&&!isfroma){h_bE_notfromA_Pt->Fill(ele->pt(),genWeight);}
      if (isfromb&&isfroma){
         genEfromB_S1.push_back(ele);
         h_EnStepsToFinal->Fill(n,genWeight);
         h_ELastMother_b->Fill(particleID,genWeight);
         h_nESteps_afterB->Fill(steps_afterb,genWeight);
         h_EnBinChain->Fill(nbinChian,genWeight);
         h_genbE_Pt->Fill(ele->pt(),genWeight);
         h_bDrE->Fill(b.DeltaR(e),genWeight);
      }
      else{
         h_genS1Electrons_NonPrimary_Pt->Fill(ele->pt(),genWeight);
         h_ELastMother_Nonb->Fill(particleID,genWeight);
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
      if (isfromb&&!isfroma){h_bMu_notfromA_Pt->Fill(muon->pt(),genWeight);}
      if (isfromb&&isfroma){
         genMufromB_S1.push_back(muon);
         h_MunStepsToFinal->Fill(n,genWeight);
         h_MuLastMother_b->Fill(particleID,genWeight);
         h_nMuSteps_afterB->Fill(steps_afterb,genWeight);
         h_MunBinChain->Fill(nbinChian,genWeight);
         h_genbMu_Pt->Fill(muon->pt(),genWeight);
         h_bDrMu->Fill(b.DeltaR(mu),genWeight);
      }
      else{
         h_genS1Muon_NonPrimary_Pt->Fill(muon->pt(),genWeight);
         h_MuLastMother_Nonb->Fill(particleID,genWeight);
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
   std::vector<reco::GsfElectron> selected_std_tau_electrons;
   std::vector<reco::GsfElectron> selected_std_b_electrons;
   std::vector<int> matchedTauEleGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level
   std::vector<int> matchedbEleGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level

   for(auto& ele : *std_electrons){
      h_nStdElectrons->Fill(0.5,genWeight);
      if (ele.pt()>7 && abs(ele.eta())<2.5){
         h_nStdElectrons->Fill(1.5,genWeight);
         float E_c=ele.superCluster()->energy();
         if(h2AA_RECO::checkID_RECO(ele,1,rho, E_c)){
            h_nStdElectrons->Fill(2.5,genWeight);
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
   std::vector<reco::GsfElectronRef> selected_LowPt_tau_electrons;
   std::vector<reco::GsfElectronRef> selected_LowPt_b_electrons;
   std::vector<int> matchedLowPtTauEleGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level
   std::vector<int> matchedLowPtbEleGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level

   for(unsigned int iElec = 0; iElec<lowPt_electrons->size(); iElec++){
      reco::GsfElectronRef ele (lowPt_electrons, iElec);
      h_nLowPtElectrons->Fill(0.5,genWeight);
      if (ele->pt()>1 && abs(ele->eta())<2.5){
         h_nLowPtElectrons->Fill(1.5,genWeight);
         float idVal = (EleMVAID)[ele];
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
   //std::cout<<"LowPt Electron selection\n";
   std::vector<reco::PFJet> selected_jets;
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
         h_JetPt->Fill(jet.pt(),genWeight);
         selected_jets.push_back(jet);
      }
   }
   cout<<"jet selection\n";
   std::vector<reco::Muon> selected_tau_Muons;
   std::vector<reco::Muon> selected_b_Muons;
   std::vector<int> matchedTauMuGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level
   std::vector<int> matchedbMuGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level

   for (auto& muon: *muons){
      h_nMuons->Fill(.5,genWeight);
      if (muon::isLooseMuon(muon)==false){continue;}////////////Loose Muon ID
      h_nMuons->Fill(1.5,genWeight);
      if((muon.pt())<3|| muon.eta()>2.4){continue;} 
      h_nMuons->Fill(2.5,genWeight);
      bool tauMumatched= false;
      float tauDr_min=.1;
      float tauMatchedGen_index=-1;
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
         h_TauMuonIsolation->Fill(h2AA_RECO::muonIsolation_RECO(muon),genWeight);
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
         h_Muon_bmatched_Iso->Fill(h2AA_RECO::muonIsolation_RECO(muon),genWeight);
      }
   }
   cout<<"Mu selection\n";
   std::vector<reco::PFTau> selected_taus;
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
   
   cout<<"tau selection\n";



   ////////////////TRIGGER///////////////

   for (unsigned int iTrigger = 0; iTrigger<MuEG.size();iTrigger++){
      std::string Trigger = MuEG[iTrigger];
      bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
      if (passTrig){h_MuonEGamma->Fill(Trigger.c_str(), genWeight);}
   }
   for (unsigned int iTrigger = 0; iTrigger<singleMu.size();iTrigger++){
      std::string Trigger = singleMu[iTrigger];
      bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
      if (passTrig){h_singleMu->Fill(Trigger.c_str(), genWeight);}
   }
   for (unsigned int iTrigger = 0; iTrigger<singleEle.size();iTrigger++){
      std::string Trigger = singleEle[iTrigger];
      bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
      if (passTrig){h_singleE->Fill(Trigger.c_str(), genWeight);}
   }
   for (unsigned int iTrigger = 0; iTrigger<DoubleMuon.size();iTrigger++){
      std::string Trigger = DoubleMuon[iTrigger];
      bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
      if (passTrig){h_DoubleMuon->Fill(Trigger.c_str(), genWeight);}
   }
   for (unsigned int iTrigger = 0; iTrigger<DoubleElectron.size();iTrigger++){
      std::string Trigger = DoubleElectron[iTrigger];
      bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
      if (passTrig){h_DoubleElectron->Fill(Trigger.c_str(), genWeight);}
   }

//////////////////////Event Selection////////////////////////
///////////////////
   std::sort(genMufromB_S1.begin(), genMufromB_S1.end(), h2AA_RECO::sortGenByPt_RECO);
   std::sort(genEfromB_S1.begin(), genEfromB_S1.end(), h2AA_RECO::sortGenByPt_RECO);
   std::sort(genTauElectrons.begin(), genTauElectrons.end(), h2AA_RECO::sortGenByPt_RECO);
   std::sort(genTauMuons.begin(), genTauMuons.end(), h2AA_RECO::sortGenByPt_RECO);
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


//bE TauE
   h_nEvent_bE_tauE->Fill(.5,genWeight);
   if (genEfromB_S1.size()>0){
      h_nEvent_bE_tauE->Fill(1.5,genWeight);
      if (genTauElectrons.size()>0){
         h_nEvent_bE_tauE->Fill(2.5,genWeight);
         h_bE_tauE_genbEpt->Fill(genEfromB_S1[0]->pt(),genWeight);
         h_bE_tauE_gentauEpt->Fill(genTauElectrons[0]->pt(),genWeight);
         for (unsigned int iTrigger = 0; iTrigger<DoubleElectron.size();iTrigger++){
            std::string Trigger = DoubleElectron[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bE_tauE_DoubleElectron->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<singleEle.size();iTrigger++){
            std::string Trigger = singleEle[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bE_tauE_singleEle->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<MuEG.size();iTrigger++){
            std::string Trigger = MuEG[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bE_tauE_MuEG->Fill(Trigger.c_str(), genWeight);}
         }
         if (passTrig_DoubleElectronLoose&&!passTrig_EGLoose){h_Trigger_bE_TauE_DiEandEG->Fill(0.5, genWeight);}
         if (!passTrig_DoubleElectronLoose&&passTrig_EGLoose){h_Trigger_bE_TauE_DiEandEG->Fill(1.5, genWeight);}
         if (passTrig_DoubleElectronLoose&&passTrig_EGLoose){h_Trigger_bE_TauE_DiEandEG->Fill(2.5, genWeight);}

         if (passTrig_EGLoose&&!passTrig_MuEGLoose){h_Trigger_bE_TauE_EGandMuEG->Fill(0.5, genWeight);}
         if (!passTrig_EGLoose&&passTrig_MuEGLoose){h_Trigger_bE_TauE_EGandMuEG->Fill(1.5, genWeight);}
         if (passTrig_EGLoose&&passTrig_MuEGLoose){h_Trigger_bE_TauE_EGandMuEG->Fill(2.5, genWeight);}
      }
   }


   h_nEvent_bMu_tauE->Fill(.5,genWeight);
   if (genMufromB_S1.size()>0){
      h_nEvent_bMu_tauE->Fill(1.5,genWeight);
      if (genTauElectrons.size()>0){
         h_nEvent_bMu_tauE->Fill(2.5,genWeight);

         h_bMu_tauE_genbMupt->Fill(genMufromB_S1[0]->pt(),genWeight);
         h_bMu_tauE_gentauEpt->Fill(genTauElectrons[0]->pt(),genWeight);
         for (unsigned int iTrigger = 0; iTrigger<DoubleElectron.size();iTrigger++){
            std::string Trigger = DoubleElectron[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bMu_tauE_DoubleElectron->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<singleEle.size();iTrigger++){
            std::string Trigger = singleEle[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bMu_tauE_singleEle->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<MuEG.size();iTrigger++){
            std::string Trigger = MuEG[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bMu_tauE_MuEG->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<DoubleMuon.size();iTrigger++){
            std::string Trigger = DoubleMuon[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bMu_tauE_DoubleMuon->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<singleMu.size();iTrigger++){
            std::string Trigger = singleMu[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bMu_tauE_singleMu->Fill(Trigger.c_str(), genWeight);}
         }
         if (passTrig_singleMuLoose&&!passTrig_EGLoose){h_Trigger_bMu_TauE_MuandEG->Fill(0.5, genWeight);}
         if (!passTrig_singleMuLoose&&passTrig_EGLoose){h_Trigger_bMu_TauE_MuandEG->Fill(1.5, genWeight);}
         if (passTrig_singleMuLoose&&passTrig_EGLoose){h_Trigger_bMu_TauE_MuandEG->Fill(2.5, genWeight);}

         if (passTrig_EGLoose&&!passTrig_MuEGLoose){h_Trigger_bMu_TauE_EGandMuEG->Fill(0.5, genWeight);}
         if (!passTrig_EGLoose&&passTrig_MuEGLoose){h_Trigger_bMu_TauE_EGandMuEG->Fill(1.5, genWeight);}
         if (passTrig_EGLoose&&passTrig_MuEGLoose){h_Trigger_bMu_TauE_EGandMuEG->Fill(2.5, genWeight);}

         if (passTrig_singleMuLoose&&!passTrig_MuEGLoose){h_Trigger_bMu_TauE_MuandMuEG->Fill(0.5, genWeight);}
         if (!passTrig_singleMuLoose&&passTrig_MuEGLoose){h_Trigger_bMu_TauE_MuandMuEG->Fill(1.5, genWeight);}
         if (passTrig_singleMuLoose&&passTrig_MuEGLoose){h_Trigger_bMu_TauE_MuandMuEG->Fill(2.5, genWeight);}

         
      }
   }

   
   if (!passTrig_EGLoose&&passTrig_DoubleElectronLoose){h_Trigger_general_EGandDiEG->Fill(1.5, genWeight);}
   if (passTrig_EGLoose&&!passTrig_DoubleElectronLoose){h_Trigger_general_EGandDiEG->Fill(0.5, genWeight);}
   if (passTrig_EGLoose&&passTrig_DoubleElectronLoose){h_Trigger_general_EGandDiEG->Fill(2.5, genWeight);}

   if (passTrig_EGLoose&&!passTrig_MuEGLoose){h_Trigger_general_EGandMuEG->Fill(0.5, genWeight);}
   if (!passTrig_EGLoose&&passTrig_MuEGLoose){h_Trigger_general_EGandMuEG->Fill(1.5, genWeight);}
   if (passTrig_EGLoose&&passTrig_MuEGLoose){h_Trigger_general_EGandMuEG->Fill(2.5, genWeight);}

   if (passTrig_singleMuLoose&&!passTrig_MuEGLoose){h_Trigger_general_MuandMuEG->Fill(0.5, genWeight);}
   if (!passTrig_singleMuLoose&&passTrig_MuEGLoose){h_Trigger_general_MuandMuEG->Fill(1.5, genWeight);}
   if (passTrig_singleMuLoose&&passTrig_MuEGLoose){h_Trigger_general_MuandMuEG->Fill(2.5, genWeight);}

   if (passTrig_singleMuLoose&&!passTrig_doubleMuonLoose){h_Trigger_general_MuandDiMu->Fill(0.5, genWeight);}
   if (!passTrig_singleMuLoose&&passTrig_doubleMuonLoose){h_Trigger_general_MuandDiMu->Fill(1.5, genWeight);}
   if (passTrig_singleMuLoose&&passTrig_doubleMuonLoose){h_Trigger_general_MuandDiMu->Fill(2.5, genWeight);}

   h_nEvent_bE_tauMu->Fill(.5,genWeight);
   if (genEfromB_S1.size()>0){
      h_nEvent_bE_tauMu->Fill(1.5,genWeight);
      if (genTauMuons.size()>0){
         h_nEvent_bE_tauMu->Fill(2.5,genWeight);

         h_bE_tauMu_gentauMupt->Fill(genTauMuons[0]->pt(),genWeight);
         h_bE_tauMu_genbEpt->Fill(genEfromB_S1[0]->pt(),genWeight);
         for (unsigned int iTrigger = 0; iTrigger<DoubleElectron.size();iTrigger++){
            std::string Trigger = DoubleElectron[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bE_tauMu_DoubleElectron->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<singleEle.size();iTrigger++){
            std::string Trigger = singleEle[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bE_tauMu_singleEle->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<MuEG.size();iTrigger++){
            std::string Trigger = MuEG[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bE_tauMu_MuEG->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<DoubleMuon.size();iTrigger++){
            std::string Trigger = DoubleMuon[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bE_tauMu_DoubleMuon->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<singleMu.size();iTrigger++){
            std::string Trigger = singleMu[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bE_tauMu_singleMu->Fill(Trigger.c_str(), genWeight);}
         }
         if (passTrig_singleMuLoose&&!passTrig_EGLoose){h_Trigger_bE_TauMu_MuandEG->Fill(0.5, genWeight);}
         if (!passTrig_singleMuLoose&&passTrig_EGLoose){h_Trigger_bE_TauMu_MuandEG->Fill(1.5, genWeight);}
         if (passTrig_singleMuLoose&&passTrig_EGLoose){h_Trigger_bE_TauMu_MuandEG->Fill(2.5, genWeight);}

         if (passTrig_EGLoose&&!passTrig_MuEGLoose){h_Trigger_bE_TauMu_EGandMuEG->Fill(0.5, genWeight);}
         if (!passTrig_EGLoose&&passTrig_MuEGLoose){h_Trigger_bE_TauMu_EGandMuEG->Fill(1.5, genWeight);}
         if (passTrig_EGLoose&&passTrig_MuEGLoose){h_Trigger_bE_TauMu_EGandMuEG->Fill(2.5, genWeight);}

         if (passTrig_singleMuLoose&&!passTrig_MuEGLoose){h_Trigger_bE_TauMu_MuandMuEG->Fill(0.5, genWeight);}
         if (!passTrig_singleMuLoose&&passTrig_MuEGLoose){h_Trigger_bE_TauMu_MuandMuEG->Fill(1.5, genWeight);}
         if (passTrig_singleMuLoose&&passTrig_MuEGLoose){h_Trigger_bE_TauMu_MuandMuEG->Fill(2.5, genWeight);}
      }
   }
   h_nEvent_bMu_tauMu->Fill(.5,genWeight);
   if (genMufromB_S1.size()>0){
      h_nEvent_bMu_tauMu->Fill(1.5,genWeight);
      if (genTauMuons.size()>0){
         h_nEvent_bMu_tauMu->Fill(2.5,genWeight);
         h_bMu_tauMu_genbMupt->Fill(genMufromB_S1[0]->pt(),genWeight);
         h_bMu_tauMu_gentauMupt->Fill(genTauMuons[0]->pt(),genWeight);
         for (unsigned int iTrigger = 0; iTrigger<DoubleMuon.size();iTrigger++){
            std::string Trigger = DoubleMuon[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bMu_tauMu_DoubleMuon->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<singleMu.size();iTrigger++){
            std::string Trigger = singleMu[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bMu_tauMu_singleMu->Fill(Trigger.c_str(), genWeight);}
         }
         for (unsigned int iTrigger = 0; iTrigger<MuEG.size();iTrigger++){
            std::string Trigger = MuEG[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_Trigger_bMu_tauMu_MuEG->Fill(Trigger.c_str(), genWeight);}
         }
         if (passTrig_doubleMuonLoose&&!passTrig_singleMuLoose){h_Trigger_bMu_TauMu_DiMuandMu->Fill(0.5, genWeight);}
         if (!passTrig_doubleMuonLoose&&passTrig_singleMuLoose){h_Trigger_bMu_TauMu_DiMuandMu->Fill(1.5, genWeight);}
         if (passTrig_doubleMuonLoose&&passTrig_singleMuLoose){h_Trigger_bMu_TauMu_DiMuandMu->Fill(2.5, genWeight);}

         if (passTrig_singleMuLoose&&!passTrig_MuEGLoose){h_Trigger_bMu_TauMu_MuandMuEG->Fill(0.5, genWeight);}
         if (!passTrig_singleMuLoose&&passTrig_MuEGLoose){h_Trigger_bMu_TauMu_MuandMuEG->Fill(1.5, genWeight);}
         if (passTrig_singleMuLoose&&passTrig_MuEGLoose){h_Trigger_bMu_TauMu_MuandMuEG->Fill(2.5, genWeight);}
      }
   }

///////////End
}
// ------------ method called once each job just before starting event loop  ------------
void
higgstoaaAnalyzer_RECO::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
higgstoaaAnalyzer_RECO::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
higgstoaaAnalyzer_RECO::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

}

//define this as a plug-in
DEFINE_FWK_MODULE(higgstoaaAnalyzer_RECO);
