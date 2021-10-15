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

#include "DataFormats/Candidate/interface/Candidate.h"
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
      TH1F *h_ElectronPt;
      TH1F *h_GenTauElectronPt;
      TH1F *h_GenBElectronPt;
      TH1F *h_GenBJetPt;
      TH1F *h_JetPt;
      TH1F *h_TauPt;
      TH1F *h_GenTauPt;
      TH1F *h_GenTauMuonPt;
      TH1F *h_GenBMuonPt;
      TH1F *h_MuonPt;
      
      TH1F *h_MuonIsolation;
      TH1F *h_MET; 
      TH1F *h_nElectrons;
      TH1F *h_nJets;
      TH1F *h_nMuons;
      TH1F *h_nTaus;
      TH1F *h_neebb;
      TH1F *h_nemubb;
      TH1F *h_nmumubb;

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

      TH1F *h_nehadbb;
      TH1F *h_nmuhadbb;
      TH1F *h_nGenMu;

      TH1F *h_eebb_e1pt;
      TH1F *h_eebb_e2pt;

      TH1F *h_emubb_ept;
      TH1F *h_emubb_mupt;
      TH1F *h_ehadbb_ept;
      TH1F *h_ehadbb_hadpt;
      TH1F *h_mumubb_mu1pt;
      TH1F *h_mumubb_mu2pt;
      TH1F *h_muhadbb_mupt;
      TH1F *h_muhadbb_hadpt;

      TH1F *h_muhadbb_DeltaRMuHad;
      TH1F *h_muhadbb_DeltaRMuj1;

      TH1F *h_mumubb_DeltaRMuMu;
      TH1F *h_mumubb_DeltaRMuj1;

      TH1F *h_emubb_DeltaRMuE;
      TH1F *h_emubb_DeltaRMuj1;

      TH1F *h_eebb_DeltaREE;
      TH1F *h_eebb_DeltaREj1;

      TH1F *h_ehadbb_DeltaREHad;
      TH1F *h_ehadbb_DeltaREJ1;
      TH1F *h_DeltaR_GenTauGenTau;

      TH1F *h_nEStatus1_Gen;
      TH1F *h_nMuStatus1_Gen;

      TH1F *h_nEStatus1_NotfromTau_Gen;
      TH1F *h_nMuStatus1_NotfromTau_Gen;

      TH1F *h_nEHPFinal_Gen;
      TH1F *h_nMuHPFinal_Gen;

      TH1F *h_nEHPFinal_NotfromTau_Gen;
      TH1F *h_nMuHPFinal_NotfromTau_Gen;

      TH1F *h_DeltaR_GenTauGenB;
      TH1F *h_DeltaR_GenBGenB;

      TH1F *h_DeltaR_GenBGenE_OS;
      TH1F *h_DeltaR_GenBGenMu_OS;


      TH1F *h_DeltaR_GenBGenE_SS;
      TH1F *h_DeltaR_GenBGenMu_SS;
      
      TH1F *h_BDaughter_GenParticleID;
      TH1F *h_BDaughter_Status1_GenParticleID;

      TH1F *h_DoubleElectron;


      TH1F *h_singleE;
      TH1F *h_singleMu;
      
      TH1F *h_DoubleMuon;
      TH1F *h_MuonEGamma;

      TH1F *h_nBDaughters;
      TH1F *h_ENullMother;
      TH1F *h_EMotherB;
      TH1F *h_EnStepsToFinal;
      TH1F *h_EfirstMother;
      TH1F *h_ELastMother;

      TH1F *h_nTauDaughters;
      TH1F *h_TauDaughtersIDs;
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
      h_GenTauElectronPt=fs->make<TH1F>("h_GenTauElectronPt", "Gen Tau Electron PT (GeV)", 40,0,160);
      h_GenTauElectronPt->Sumw2();
      h_GenBJetPt=fs->make<TH1F>("h_GenBJetPt", "Gen BJet PT (GeV)", 40,0,200);
      h_GenBJetPt->Sumw2();
      h_JetPt=fs->make<TH1F>("h_JetPt", "Reco Jet PT (GeV)", 40,0,200);
      h_JetPt->Sumw2();
      h_TauPt=fs->make<TH1F>("h_TauPt", "Reco TauHad PT (GeV)", 50,0,250);
      h_TauPt->Sumw2();
      h_GenTauPt=fs->make<TH1F>("h_GenTauPt", "Gen Tau PT (GeV)", 50,0,250);
      h_GenTauPt->Sumw2();
      h_GenTauMuonPt=fs->make<TH1F>("h_GenTauMuonPt", "Gen Tau Muon PT (GeV)", 25,0,125);
      h_GenTauMuonPt->Sumw2();
      h_GenBMuonPt=fs->make<TH1F>("h_GenBMuonPt", "Gen B Muon PT (GeV)", 25,0,125);
      h_GenBMuonPt->Sumw2();
      h_GenBElectronPt=fs->make<TH1F>("h_GenBElectronPt", "Gen B Muon PT (GeV)", 25,0,125);
      h_GenBElectronPt->Sumw2();
      h_MuonPt=fs->make<TH1F>("h_MuonPt", "Reco Muon PT (GeV)", 25,0,125);
      h_MuonPt->Sumw2();
      h_MuonIsolation=fs->make<TH1F>("h_MuonIsolation", "Reco Muon Isolation", 25,0,1);
      h_MuonIsolation->Sumw2();
      h_MET=fs->make<TH1F>("h_MET", "MET (GeV)", 25,0,250);
      h_MET->Sumw2();
      h_nElectrons=fs->make<TH1F>("h_nElectrons", "Electron Cuts", 4,0,4);
      h_nElectrons->Sumw2();
      h_nMuons=fs->make<TH1F>("h_nMuons", "Muon Cuts", 4,0,4);
      h_nMuons->Sumw2();
      h_nTaus=fs->make<TH1F>("h_nTaus", "Tau Cuts", 3,0,3);
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

      h_muhad_tau=fs->make<TH1F>("h_muhad_tau", "EHAD Tau Triggers",3,0,3);
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
///////Selection Plots   

      h_neebb=fs->make<TH1F>("h_neebb", "NEvent ee",3,0,3);
      h_neebb->Sumw2();

      h_nemubb=fs->make<TH1F>("h_nemubb", "NEvent eMu",4,0,4);
      h_nemubb->Sumw2();

      h_nmumubb=fs->make<TH1F>("h_nmumubb", "N Mu Mu Events",3,0,3);
      h_nmumubb->Sumw2();

      h_nehadbb=fs->make<TH1F>("h_nehadbb", "n eHad Events",4,0,4);
      h_nehadbb->Sumw2();

      h_nmuhadbb=fs->make<TH1F>("h_nmuhadbb", "N Mu Had Events",4,0,4);
      h_nmuhadbb->Sumw2();

      h_nGenMu=fs->make<TH1F>("h_nGenMu", "genMuon",3,0,3);
      h_nGenMu->SetCanExtend(TH1::kAllAxes);
      h_nGenMu->Sumw2();

      h_eebb_e1pt=fs->make<TH1F>("h_eebb_e1pt", "eebb Electron1 PT (GeV)", 40,0,160);
      h_eebb_e1pt->Sumw2();

      h_eebb_e2pt=fs->make<TH1F>("h_eebb_e2pt", "eebb Electron2 PT (GeV)", 40,0,160);
      h_eebb_e2pt->Sumw2();

      h_emubb_mupt=fs->make<TH1F>("h_emubb_mupt", "emubb Muon PT (GeV)", 25,0,125);
      h_emubb_mupt->Sumw2();

      h_emubb_ept=fs->make<TH1F>("h_emubb_ept", "eebb Electron PT (GeV)", 40,0,160);
      h_emubb_ept->Sumw2();

      h_ehadbb_hadpt=fs->make<TH1F>("h_ehadbb_hadpt", "ehadbb TauHad PT (GeV)", 50,0,250);
      h_ehadbb_hadpt->Sumw2();

      h_ehadbb_ept=fs->make<TH1F>("h_ehadbb_ept", "eebb Electron PT (GeV)", 40,0,160);
      h_ehadbb_ept->Sumw2();

      h_muhadbb_hadpt=fs->make<TH1F>("h_muhadbb_hadpt", "muhadbb TauHad PT (GeV)", 50,0,250);
      h_muhadbb_hadpt->Sumw2();

      h_muhadbb_mupt=fs->make<TH1F>("h_muhadbb_mupt", "hadmubb Muon PT (GeV)", 25,0,125);
      h_muhadbb_mupt->Sumw2();
      h_mumubb_mu1pt=fs->make<TH1F>("h_mumubb_mu1pt", "mumubb Muon1 PT (GeV)", 25,0,125);
      h_mumubb_mu1pt->Sumw2();
      h_mumubb_mu2pt=fs->make<TH1F>("h_mumubb_mu2pt", "emubb Muon2 PT (GeV)", 25,0,125);
      h_mumubb_mu2pt->Sumw2();

      h_muhadbb_DeltaRMuHad=fs->make<TH1F>("h_muhadbb_DeltaRMuHad", "muhadbb DeltaR(Mu,Tau)", 60,0,6);
      h_muhadbb_DeltaRMuHad->Sumw2();
      h_muhadbb_DeltaRMuj1=fs->make<TH1F>("h_muhadbb_DeltaRMuj1", "muhadbb DeltaR(Mu,J1)", 60,0,6);
      h_muhadbb_DeltaRMuj1->Sumw2();

      h_mumubb_DeltaRMuMu=fs->make<TH1F>("h_mumubb_DeltaRMuMu", "mumubb DeltaR(Mu,Mu)", 60,0,6);
      h_mumubb_DeltaRMuMu->Sumw2();
      h_mumubb_DeltaRMuj1=fs->make<TH1F>("h_mumubb_DeltaRMuj1", "mumubb DeltaR(Mu,J1)", 60,0,6);
      h_mumubb_DeltaRMuj1->Sumw2();

      h_emubb_DeltaRMuE=fs->make<TH1F>("h_emubb_DeltaRMuE", "emubb DeltaR(e,Mu)", 60,0,6);
      h_emubb_DeltaRMuE->Sumw2();
      h_emubb_DeltaRMuj1=fs->make<TH1F>("h_emubb_DeltaRMuj1", "emubb DeltaR(Mu,J1)", 60,0,6);
      h_emubb_DeltaRMuj1->Sumw2();

      h_eebb_DeltaREE=fs->make<TH1F>("h_eebb_DeltaREE", "eebb DeltaR(e,e)", 60,0,6);
      h_eebb_DeltaREE->Sumw2();
      h_eebb_DeltaREj1=fs->make<TH1F>("h_eebb_DeltaREj1", "eebb DeltaR(e,J1)", 60,0,6);
      h_eebb_DeltaREj1->Sumw2();   

      h_ehadbb_DeltaREHad=fs->make<TH1F>("h_ehadbb_DeltaREHad", "ehadbb DeltaR(e,tau)", 60,0,6);
      h_ehadbb_DeltaREHad->Sumw2();
      h_ehadbb_DeltaREJ1=fs->make<TH1F>("h_ehadbb_DeltaREJ1", "ehadbb DeltaR(e,J1)", 60,0,6);
      h_ehadbb_DeltaREJ1->Sumw2();  

      h_DeltaR_GenTauGenTau=fs->make<TH1F>("h_DeltaR_GenTauGenTau", "DeltaR(GenTau,GenTau)", 60,0,6);
      h_DeltaR_GenTauGenTau->Sumw2();  

      h_DeltaR_GenBGenE_SS=fs->make<TH1F>("h_DeltaR_GenBGenE_SS", "DeltaR(GenE,GenB) Same Sign", 60,0,6);
      h_DeltaR_GenBGenE_SS->Sumw2();  

      h_DeltaR_GenBGenE_OS=fs->make<TH1F>("h_DeltaR_GenBGenE_OS", "DeltaR(GenE,GenB) Opposite Sign", 60,0,6);
      h_DeltaR_GenBGenE_OS->Sumw2();  

      h_DeltaR_GenTauGenB=fs->make<TH1F>("h_DeltaR_GenTauGenB", "DeltaR(GenB,GenTau)", 60,0,6);
      h_DeltaR_GenTauGenB->Sumw2();  

      h_DeltaR_GenBGenB=fs->make<TH1F>("h_DeltaR_GenBGenB", "DeltaR(GenB,GenB)", 60,0,6);
      h_DeltaR_GenBGenB->Sumw2();  

      h_DeltaR_GenBGenMu_SS=fs->make<TH1F>("h_DeltaR_GenBGenMu_SS", "DeltaR(GenMu,GenB) Opposite Sign", 60,0,6);
      h_DeltaR_GenBGenMu_SS->Sumw2();  

      h_DeltaR_GenBGenMu_OS=fs->make<TH1F>("h_DeltaR_GenBGenMu_OS", "DeltaR(GenMu,GenB) Opposite Sign", 60,0,6);
      h_DeltaR_GenBGenMu_OS->Sumw2(); 

      h_nEStatus1_Gen=fs->make<TH1F>("h_nEStatus1_Gen", "N Electrons Status 1  (Gen)", 8,0,8);
      h_nEStatus1_Gen->Sumw2(); 

      h_nMuStatus1_Gen=fs->make<TH1F>("h_nMuStatus1_Gen", "N Muons Status 1 (Gen)", 8,0,8);
      h_nMuStatus1_Gen->Sumw2(); 

      h_nEStatus1_NotfromTau_Gen=fs->make<TH1F>("h_nEStatus1_NotfromTau_Gen", "N Electrons Status 1 (!fromtau) (Gen)", 8,0,8);
      h_nEStatus1_NotfromTau_Gen->Sumw2(); 

      h_nMuStatus1_NotfromTau_Gen=fs->make<TH1F>("h_nMuStatus1_NotfromTau_Gen", "N Muons Status (!fromtau) 1 (Gen)", 8,0,8);
      h_nMuStatus1_NotfromTau_Gen->Sumw2(); 

      h_nEHPFinal_Gen=fs->make<TH1F>("h_nEHPFinal_Gen", "N Electrons HP FinalState  (Gen)", 8,0,8);
      h_nEHPFinal_Gen->Sumw2(); 

      h_nMuHPFinal_Gen=fs->make<TH1F>("h_nMuHPFinal_Gen", "N Muons HP Finalstate (Gen)", 8,0,8);
      h_nMuHPFinal_Gen->Sumw2(); 

      h_nEHPFinal_NotfromTau_Gen=fs->make<TH1F>("h_nEHPFinal_NotfromTau_Gen", "N Electrons HP Finalstate (!fromtau) (Gen)", 8,0,8);
      h_nEHPFinal_NotfromTau_Gen->Sumw2(); 

      h_nMuHPFinal_NotfromTau_Gen=fs->make<TH1F>("h_nMuHPFinal_NotfromTau_Gen", "N Muons HP Finalstate (!fromtau) 1 (Gen)", 8,0,8);
      h_nMuHPFinal_NotfromTau_Gen->Sumw2(); 

      
      h_BDaughter_GenParticleID=fs->make<TH1F>("h_BDaughter_GenParticleID", "Gen BdaughterPDG ID", 37,0,37);
      h_BDaughter_GenParticleID->Sumw2();

      h_BDaughter_Status1_GenParticleID=fs->make<TH1F>("h_BDaughter_Status1_GenParticleID", "Gen Bdaughter(satus1) PDG ID", 37,0,37);
      h_BDaughter_Status1_GenParticleID->Sumw2();

      h_nBDaughters=fs->make<TH1F>("h_nBDaughters", "N Daughters", 8,0,8);
      h_nBDaughters->Sumw2(); 

      h_ENullMother=fs->make<TH1F>("h_ENullMother", "N E with Null Mother", 8,0,8);
      h_ENullMother->Sumw2(); 

      h_EnStepsToFinal=fs->make<TH1F>("h_EnStepsToFinal", "N E steps to B", 16,-8,8);
      h_EnStepsToFinal->Sumw2(); 

      h_EfirstMother=fs->make<TH1F>("h_EfirstMother", "E First Mother PDG ID", 37,0,37);
      h_EfirstMother->Sumw2();

      h_ELastMother=fs->make<TH1F>("h_ELastMother", "E Last Mother PDG ID", 37,0,37);
      h_ELastMother->Sumw2();

      h_EMotherB=fs->make<TH1F>("h_EMotherB", "N E with B Mother", 8,0,8);
      h_EMotherB->Sumw2(); 

      h_TauDaughtersIDs=fs->make<TH1F>("h_TauDaughtersIDs", "Tau Daughters", 38,-1,37);
      h_TauDaughtersIDs->Sumw2(); 

      h_nTauDaughters=fs->make<TH1F>("h_nTauDaughters", "Tau Daughters ", 5,0,5);
      h_nTauDaughters->Sumw2(); 

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

   edm::TriggerResults triggerResults = *triggerBits.product();
   cout<<to_string(triggerResults.size());
   const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerBits);  

   //for (unsigned int i = 0; i<trigNames.size();i++){cout<< to_string(i)<<" "<< trigNames.triggerName(i)<<" "<< triggerResults.accept(i)<<"\n";}
   double rho = h2AA::getObject<double>(eventrho_,iEvent);
   edm::Handle<GenEventInfoProduct>  genInfoHandle;
   try{iEvent.getByToken(genInfoProduct_, genInfoHandle);}
   catch(...){;}
   double genWeight=genInfoHandle->weight();

   h_MET->Fill(mets[0].pt(),genWeight);
   //////Gen Particle Sorting
   std::vector<reco::GenParticle*> genTauElectrons;
   std::vector<reco::GenParticle*> genTauMuons;
   std::vector<reco::GenParticle*> genTaus;
   std::vector<reco::GenParticle*> genBs;
   std::vector<reco::GenParticle*> genS1Electrons;
   std::vector<reco::GenParticle*> genHPFinalStateElectrons;
   std::vector<reco::GenParticle*> genHPFinalStateElectrons_NonTau;
   std::vector<reco::GenParticle*> genS1Muons;
   std::vector<reco::GenParticle*> genS1Electrons_NonTau;
   std::vector<reco::GenParticle*> genHPFinalStateMuons;
   std::vector<reco::GenParticle*> genHPFinalStateMuons_NonTau;
   std::vector<reco::GenParticle*> genS1Muons_NonTau;
   reco::GenParticle *genNuTau=NULL; //////These are just the tau neutrinos We know there is a maximum of 2 per event either one anti, one normal, both or none.
   reco::GenParticle *genNuTauBar=NULL; 
   for(auto& gen : genParticles){
      if (gen.status()==1){
         h_GenParticleID->Fill(abs(gen.pdgId())+.5,genWeight);
      }
      if(abs(gen.pdgId())==11 &&
         gen.isDirectHardProcessTauDecayProductFinalState()) {
         h_GenTauElectronPt->Fill(gen.pt(),genWeight);
         genTauElectrons.push_back(&gen);

      }
      if(abs(gen.pdgId())==11&&gen.status()==1) {
         genHPFinalStateElectrons.push_back(&gen);
         if (!gen.isDirectHardProcessTauDecayProductFinalState()){genHPFinalStateElectrons_NonTau.push_back(&gen);}
      }
      if(abs(gen.pdgId())==13&&gen.status()==1) {
         genHPFinalStateMuons.push_back(&gen);
         if (!gen.isDirectHardProcessTauDecayProductFinalState()){genHPFinalStateMuons_NonTau.push_back(&gen);}
      }
      if(abs(gen.pdgId())==11&&gen.status()==1) {
         genS1Electrons.push_back(&gen);
         if (!gen.isDirectHardProcessTauDecayProductFinalState()){genS1Electrons_NonTau.push_back(&gen);}
      }
      if(abs(gen.pdgId())==13 &&gen.status()==1) {
         genS1Muons.push_back(&gen);
         if (!gen.isDirectHardProcessTauDecayProductFinalState()){genS1Muons_NonTau.push_back(&gen);}
      }
      if (abs(gen.pdgId())==13 && gen.isDirectHardProcessTauDecayProductFinalState()){
         h_GenTauMuonPt->Fill(gen.pt(),genWeight);
         genTauMuons.push_back(&gen);
      }
      if (abs(gen.pdgId())==15 && gen.isHardProcess()){
         h_GenTauPt->Fill(gen.pt(),genWeight);
         genTaus.push_back(&gen);
         vector<int>::size_type n = gen.numberOfDaughters();
         h_nTauDaughters->Fill(n,genWeight);
         for(vector<int>::size_type j = 0; j < n; ++ j) {
            auto *d = gen.daughter( j );//Try type auto or some other

            if (d!=NULL){

               int dauId = abs(d->pdgId());
               h_TauDaughtersIDs->Fill(dauId,genWeight);
            }
            else{h_TauDaughtersIDs->Fill(-1,genWeight);}
         }
      }

      if (gen.pdgId()==16 && gen.isDirectHardProcessTauDecayProductFinalState()){
         genNuTau = &gen;

      }
      if (gen.pdgId()==-16 && gen.isDirectHardProcessTauDecayProductFinalState()){
         genNuTauBar = &gen;
      } 

      if (abs(gen.pdgId())==5&&gen.isHardProcess()){
         h_GenBJetPt->Fill(gen.pt(),genWeight);
         genBs.push_back(&gen);
         vector<int>::size_type n = gen.numberOfDaughters();
         h_nBDaughters->Fill(n,genWeight);
         for(vector<int>::size_type j = 0; j < n; ++ j) {
            auto *d = gen.daughter( j );//Try type auto or some other

            if (d!=NULL){

               int dauId = abs(d->pdgId());
               h_BDaughter_GenParticleID->Fill(dauId,genWeight);
               int particle_status=d->status();
               if (particle_status==1){
                  h_BDaughter_Status1_GenParticleID->Fill(dauId,genWeight);
               }
            }
         }
      
      }

   }

   int ENullMother = 0;
   int EMotherB=0;
   for(auto& ele : genHPFinalStateElectrons_NonTau){
      int n = 0;
      int particleID = ele->pdgId();
      auto *mommy = ele->mother();
      bool endinB=false;
      
      while(mommy!=NULL){
         particleID = mommy->pdgId();//get the particle id of the nth mother
         if (n==0){h_EfirstMother->Fill(abs(particleID),genWeight);}//fills if we have a nonnull first mother
         mommy=mommy->mother();
         n++;
         if (abs(particleID)==5){
            endinB = true;
            break;
         }
      }
      if (n==0){ENullMother++;}
      if (n==1&&endinB){EMotherB++;}
      if (!endinB){n=n*-1;}
      h_EnStepsToFinal->Fill(n,genWeight);
      h_ELastMother->Fill(particleID,genWeight);
   }
   h_ENullMother->Fill(ENullMother,genWeight);
   h_EMotherB->Fill(EMotherB,genWeight);
   
   std::cout<<"GenParticle Sorting\n";
   //Determin which is the neutrino associated with 
   reco::GenParticle *genTauHad=NULL;
   reco::GenParticle *genNuTauHad=NULL;
   bool DiTauHad=false;
   h_nEStatus1_Gen->Fill(genS1Electrons.size(),genWeight);
   h_nMuStatus1_Gen->Fill(genS1Muons.size(),genWeight);
   h_nEStatus1_NotfromTau_Gen->Fill(genS1Electrons_NonTau.size(),genWeight);
   h_nMuStatus1_NotfromTau_Gen->Fill(genS1Muons_NonTau.size(),genWeight);

   h_nEHPFinal_Gen->Fill(genHPFinalStateElectrons.size(),genWeight);
   h_nMuHPFinal_Gen->Fill(genHPFinalStateMuons.size(),genWeight);
   h_nEHPFinal_NotfromTau_Gen->Fill(genHPFinalStateElectrons_NonTau.size(),genWeight);
   h_nMuHPFinal_NotfromTau_Gen->Fill(genHPFinalStateMuons_NonTau.size(),genWeight);


   if(genBs.size()>1){
      TLorentzVector B1;
      TLorentzVector B2;
      B1.SetPtEtaPhiM(genBs[0]->pt(), genBs[0]->eta(), genBs[0]->phi(), genBs[0]->mass());
      B2.SetPtEtaPhiM(genBs[1]->pt(), genBs[1]->eta(), genBs[1]->phi(), genBs[1]->mass());
      h_DeltaR_GenBGenB->Fill(B1.DeltaR(B2),genWeight);
   }

   if(genTaus.size()>1){
      TLorentzVector T1;
      TLorentzVector T2;
      T1.SetPtEtaPhiM(genTaus[0]->pt(), genTaus[0]->eta(), genTaus[0]->phi(), genTaus[0]->mass());
      T1.SetPtEtaPhiM(genTaus[1]->pt(), genTaus[1]->eta(), genTaus[1]->phi(), genTaus[1]->mass());
      h_DeltaR_GenTauGenTau->Fill(T1.DeltaR(T2),genWeight);
   } 
      if(genTaus.size()>0&&genBs.size()>0){
      TLorentzVector T;
      TLorentzVector B;
      T.SetPtEtaPhiM(genTaus[0]->pt(), genTaus[0]->eta(), genTaus[0]->phi(), genTaus[0]->mass());
      B.SetPtEtaPhiM(genBs[0]->pt(), genBs[0]->eta(), genBs[0]->phi(), genBs[0]->mass());
      h_DeltaR_GenTauGenB->Fill(T.DeltaR(B),genWeight);
   }  

   std::cout<<"Tau DR plotting\n";

   for (auto& gen : genTaus){
      if ((genTauElectrons.size()==1&& genTauMuons.size()==0)){
         if ((*genTauElectrons[0]).charge()<0&&(*gen).charge()>0){
               //std::cout<<"e-\n";

         genTauHad = gen;
         genNuTauHad = genNuTauBar;
         }
         if ((*genTauElectrons[0]).charge()>0&&(*gen).charge()<0){
                     //std::cout<<"e+\n";
         genTauHad = gen;
         genNuTauHad = genNuTau;
         }
      }
      else if((genTauElectrons.size()==0&& genTauMuons.size()==1)){
         if ((*genTauMuons[0]).charge()<0&&(*gen).charge()>0){
         //std::cout<<"m-\n";
         genTauHad = gen;
         genNuTauHad = genNuTauBar;
         }
         if ((*genTauMuons[0]).charge()>0&&(*gen).charge()<0){
         //std::cout<<"m+\n";
         genTauHad= gen;
         genNuTauHad = genNuTau;
         }

      }
      else{
         DiTauHad=true;
      }
   }
   
   TLorentzVector GenTau1;
   TLorentzVector GenTau2;
   GenTau1.SetPtEtaPhiM(genTaus[0]->pt(), genTaus[0]->eta(), genTaus[0]->phi(), genTaus[0]->mass());
   GenTau2.SetPtEtaPhiM(genTaus[1]->pt(), genTaus[1]->eta(), genTaus[1]->phi(), genTaus[1]->mass());
   h_DeltaR_GenTauGenTau->Fill(GenTau1.DeltaR(GenTau2),genWeight);
   cout<<"GenTaus";
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
            for(unsigned int iGenE=0; iGenE<genTauElectrons.size(); iGenE++){
               reco::GenParticle gen = *genTauElectrons[iGenE];
               TLorentzVector e;
               TLorentzVector genE;
               e.SetPtEtaPhiM(ele.pt(), ele.eta(), ele.phi(), ele.mass());
               genE.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), gen.mass());
               float dr = genE.DeltaR(e);
               if(dr <.1){
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
   cout<<"Electron selection";
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
      for (unsigned int iGen=0; iGen<genTauMuons.size();iGen++){
         reco::GenParticle gen = *genTauMuons[iGen];
         TLorentzVector mu;
         TLorentzVector genMu;
         mu.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), muon.mass());
         genMu.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), gen.mass());
         float dr = genMu.DeltaR(mu);
         cout<<"DR= "<<to_string(dr)<<"\n";
         
         if(dr <.1){
            matched = true;
         }
      }
      if (matched==false){continue;}
      h_nMuons->Fill(3.5,genWeight);

      h_MuonPt->Fill(muon.pt(),genWeight);
      selected_muons.push_back(muon);
      h_MuonIsolation->Fill(h2AA::muonIsolation(muon),genWeight);
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
   h_nemubb->Fill(0.5);
   if (selected_jets.size()>1){
      h_nemubb->Fill(1.5);
      if(selected_muons.size()>0){
         h_nemubb->Fill(2.5);
         if(selected_electrons.size()>0){
            TLorentzVector ele;
            TLorentzVector mu;
            TLorentzVector jet;
            ele.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass());
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass());
            jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass());
            h_emubb_DeltaRMuE->Fill(mu.DeltaR(ele),genWeight);
            h_emubb_DeltaRMuj1->Fill(mu.DeltaR(jet),genWeight);
            h_nemubb->Fill(3.5);
            h_emubb_ept->Fill(selected_electrons[0].pt(),genWeight);
            h_emubb_mupt->Fill(selected_muons[0].pt(),genWeight);
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
               if (passTrig){h_eMu_singleMu->Fill(Trigger.c_str(), genWeight);}
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
   h_nGenMu->Fill(genTauMuons.size(), genWeight);
   cout<<"emubb";
   h_neebb->Fill(0.5, genWeight);
   if (selected_jets.size()>1){//If We have at least 2 jets
      h_neebb->Fill(1.5, genWeight);
      if(selected_electrons.size()>1){//if we have at least 2 electrons
         h_neebb->Fill(2.5, genWeight);
         h_eebb_e1pt->Fill(selected_electrons[0].pt(),genWeight);
         h_eebb_e2pt->Fill(selected_electrons[1].pt(),genWeight);
         TLorentzVector ele1;
         TLorentzVector ele2;
         TLorentzVector jet;
         ele1.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass());
         ele2.SetPtEtaPhiM(selected_electrons[1].pt(), selected_electrons[1].eta(), selected_electrons[1].phi(), selected_electrons[1].mass());
         jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass());
         h_eebb_DeltaREE->Fill(ele1.DeltaR(ele2),genWeight);
         h_eebb_DeltaREj1->Fill(ele1.DeltaR(jet),genWeight);
         for (unsigned int iTrigger = 0; iTrigger<DoubleElectron.size();iTrigger++){//Loop over all possible double electron triggers
            std::string Trigger = DoubleElectron[iTrigger];
            cout<<Trigger<<"\n";
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));//check to see if a given trigger passes slection
            if (passTrig){h_ee_DiEle->Fill(Trigger.c_str(), genWeight);}//if it does, fill our hist with the trigger name
         }
         for (unsigned int iTrigger = 0; iTrigger<singleEle.size();iTrigger++){//do the same again for single E trigger
            std::string Trigger = singleEle[iTrigger];
            bool passTrig=triggerBits->accept(trigNames.triggerIndex(Trigger));
            if (passTrig){h_ee_singleE->Fill(Trigger.c_str(), genWeight);}
         } 
      }
   }
   cout<<"eebb";
   h_nmumubb->Fill(0.5, genWeight);
   if (selected_jets.size()>1){
      h_nmumubb->Fill(1.5, genWeight);
      if(selected_muons.size()>1){
         h_nmumubb->Fill(2.5, genWeight);
         h_mumubb_mu1pt->Fill(selected_muons[0].pt(),genWeight);
         h_mumubb_mu2pt->Fill(selected_muons[1].pt(),genWeight);
         TLorentzVector mu1;
         TLorentzVector mu2;
         TLorentzVector jet;
         mu1.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass());
         mu2.SetPtEtaPhiM(selected_muons[1].pt(), selected_muons[1].eta(), selected_muons[1].phi(), selected_muons[1].mass());
         jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass());
         h_mumubb_DeltaRMuMu->Fill(mu1.DeltaR(mu2),genWeight);
         h_mumubb_DeltaRMuj1->Fill(mu1.DeltaR(jet),genWeight);
         for (unsigned int iTrigger = 0; iTrigger<DoubleMuon.size();iTrigger++){
            std::string Trigger = DoubleMuon[iTrigger];
            cout<<Trigger<<"\n";
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
   cout<<"mumubb";
   h_nehadbb->Fill(0.5, genWeight);
   if (selected_jets.size()>1){
      h_nehadbb->Fill(1.5, genWeight);
      cout<<"|Event Test1|";
      if(selected_taus.size()>0){
         h_nehadbb->Fill(2.5, genWeight);
         cout<<"|Event Test2|";
         if(selected_electrons.size()>0){
            h_nehadbb->Fill(3.5, genWeight);
            h_ehadbb_hadpt->Fill(selected_taus[0].pt(),genWeight);
            h_ehadbb_ept->Fill(selected_electrons[0].pt(),genWeight);
            cout<<"|Event Test3|";
            TLorentzVector ele;
            TLorentzVector tau;
            TLorentzVector jet;
            ele.SetPtEtaPhiM(selected_electrons[0].pt(), selected_electrons[0].eta(), selected_electrons[0].phi(), selected_electrons[0].mass());
            tau.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass());
            jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass());
            h_ehadbb_DeltaREHad->Fill(ele.DeltaR(tau),genWeight);
            h_ehadbb_DeltaREJ1->Fill(ele.DeltaR(jet),genWeight);
            for (unsigned int iTrigger = 0; iTrigger<singleEle.size();iTrigger++){
               cout<<"|ETrigger|";
               std::string Trigger = singleEle[iTrigger];
               cout<<Trigger<<"\n";
               int Index = trigNames.triggerIndex(Trigger.c_str());
               cout<<"INDEX = "<< std::to_string(Index); 
               //bool passTrig=triggerBits->accept(Index)
               //cout<<trigNames.triggerName(Index)<< "AND "<<trigNames[Index-1];
               bool passTrig = triggerResults.accept(Index);
               cout << "PASTRIGGERTEST";
               if (passTrig){h_ehad_singleE->Fill(Trigger.c_str(), genWeight);}
            }
            for (unsigned int iTrigger = 0; iTrigger<SingleTau.size();iTrigger++){
               std::string Trigger = SingleTau[iTrigger];
               cout<<Trigger<<"\n";

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
   cout<<"ehadbb";
   h_nmuhadbb->Fill(0.5, genWeight);
   if (selected_jets.size()>1){
      h_nmuhadbb->Fill(1.5, genWeight);
      if(selected_taus.size()>0){
         h_nmuhadbb->Fill(2.5, genWeight);
         if(selected_muons.size()>0){
            h_nmuhadbb->Fill(3.5, genWeight);
            h_muhadbb_mupt->Fill(selected_muons[0].pt(),genWeight);
            h_muhadbb_hadpt->Fill(selected_taus[0].pt(),genWeight);
            TLorentzVector mu;
            TLorentzVector tau;
            TLorentzVector jet;
            mu.SetPtEtaPhiM(selected_muons[0].pt(), selected_muons[0].eta(), selected_muons[0].phi(), selected_muons[0].mass());
            tau.SetPtEtaPhiM(selected_taus[0].pt(), selected_taus[0].eta(), selected_taus[0].phi(), selected_taus[0].mass());
            jet.SetPtEtaPhiM(selected_jets[0].pt(), selected_jets[0].eta(), selected_jets[0].phi(), selected_jets[0].mass());
            h_muhadbb_DeltaRMuHad->Fill(mu.DeltaR(tau),genWeight);
            h_muhadbb_DeltaRMuj1->Fill(mu.DeltaR(jet),genWeight);
            for (unsigned int iTrigger = 0; iTrigger<singleMu.size();iTrigger++){
               std::string Trigger = singleMu[iTrigger];
               cout<<Trigger<<"\n";
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
