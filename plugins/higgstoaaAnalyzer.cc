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
#include <utility>
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
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

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
      const edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;





      
      TH1F *h_std_electron_selection;
      TH1F *h_bJet_selection;
      TH1F *h_muon_selection;
      TH1F *h_tau_selection;



      TH1F *h_nbJet;
      TH1F *h_nbJet_m;
      TH1F *h_nbJet_f;

      TH1F *h_nTau;
      TH1F *h_nTau_m;
      TH1F *h_nTau_f;

      TH1F *h_nStdEle;
      TH1F *h_nStdEle_taum;
      TH1F *h_nStdEle_bm;
      TH1F *h_StdEle_taum_pt;
      TH1F *h_StdEle_bm_pt;
      TH1F *h_nMu;
      TH1F *h_nMu_taum;
      TH1F *h_nMu_bm;
      TH1F *h_Mu_taum_pt;
      TH1F *h_Mu_bm_pt;
      TH1F *h_Triggers;
      TH1F *h_BPHTriggers;

//B Daughter studies
      TH1F *h_deepFlavor_probb;
      TH1F *h_deepFlavor_probbb;
      TH1F *h_deepFlavor_problepb;
      TH1F *h_deepFlavor_sum;
      TH2F *h_deepFlavor_probbvsprobbb;

      
      TH1F *h_deepFlavor_probb_m;
      TH1F *h_deepFlavor_probbb_m;
      TH1F *h_deepFlavor_problepb_m;
      TH1F *h_deepFlavor_sum_m;
      TH2F *h_deepFlavor_probbvsprobbb_m;



      TH1F *h_deepFlavor_probb_f;
      TH1F *h_deepFlavor_probbb_f;
      TH1F *h_deepFlavor_problepb_f;
      TH1F *h_deepFlavor_sum_f;
      TH2F *h_deepFlavor_probbvsprobbb_f;



      TH1F *h_deepFlavor_subLeading_probb;
      TH1F *h_deepFlavor_subLeading_probbb;
      TH1F *h_deepFlavor_subLeading_problepb;
      TH1F *h_deepFlavor_subLeading_sum;
      TH2F *h_deepFlavor_subLeading_probbvsprobbb;
      
      TH1F *h_deepFlavor_subLeading_probb_m;
      TH1F *h_deepFlavor_subLeading_probbb_m;
      TH1F *h_deepFlavor_subLeading_problepb_m;
      TH1F *h_deepFlavor_subLeading_sum_m;
      TH2F *h_deepFlavor_subLeading_probbvsprobbb_m;

      TH1F *h_deepFlavor_subLeading_probb_f;
      TH1F *h_deepFlavor_subLeading_probbb_f;
      TH1F *h_deepFlavor_subLeading_problepb_f;
      TH1F *h_deepFlavor_subLeading_sum_f;
      TH2F *h_deepFlavor_subLeading_probbvsprobbb_f;


      TH1F *h_deepFlavor_leading_probb;
      TH1F *h_deepFlavor_leading_probbb;
      TH1F *h_deepFlavor_leading_problepb;
      TH1F *h_deepFlavor_leading_sum;
      TH2F *h_deepFlavor_leading_probbvsprobbb;
      
      TH1F *h_deepFlavor_leading_probb_m;
      TH1F *h_deepFlavor_leading_probbb_m;
      TH1F *h_deepFlavor_leading_problepb_m;
      TH1F *h_deepFlavor_leading_sum_m;
      TH2F *h_deepFlavor_leading_probbvsprobbb_m;

      TH1F *h_deepFlavor_leading_probb_f;
      TH1F *h_deepFlavor_leading_probbb_f;
      TH1F *h_deepFlavor_leading_problepb_f;
      TH1F *h_deepFlavor_leading_sum_f;
      TH2F *h_deepFlavor_leading_probbvsprobbb_f;



      TH1F *h_deepFlavor_Merged_probb;
      TH1F *h_deepFlavor_Merged_probbb;
      TH1F *h_deepFlavor_Merged_problepb;
      TH1F *h_deepFlavor_Merged_sum;
      TH2F *h_deepFlavor_Merged_probbvsprobbb;
      
      TH1F *h_deepFlavor_Merged_probb_m;
      TH1F *h_deepFlavor_Merged_probbb_m;
      TH1F *h_deepFlavor_Merged_problepb_m;
      TH1F *h_deepFlavor_Merged_sum_m;
      TH2F *h_deepFlavor_Merged_probbvsprobbb_m;



      TH1F *h_deepFlavor_Leptonic_probb;
      TH1F *h_deepFlavor_Leptonic_probbb;
      TH1F *h_deepFlavor_Leptonic_problepb;
      TH1F *h_deepFlavor_Leptonic_sum;
      TH2F *h_deepFlavor_Leptonic_probbvsprobbb;

      TH1F *h_deepFlavor_probb_m_Trigger;
      TH1F *h_deepFlavor_probbb_m_Trigger;
      TH1F *h_deepFlavor_problepb_m_Trigger;
      TH1F *h_deepFlavor_sum_m_Trigger;
      TH2F *h_deepFlavor_probbvsprobbb_m_Trigger;
////Event Selection
      TH1F *h_nEvent_tauE_tauE;
      TH1F *h_nEvent_tauMu_tauMu;
      TH1F *h_nEvent_tauE_tauMu;
      TH1F *h_nEvent_tauHad_tauE;
      TH1F *h_nEvent_tauHad_tauMu;
      TH1F *h_nEvent_gen_tauHad_tauE;
      TH1F *h_nEvent_gen_tauHad_tauMu;




////Trigger plots
      TH1F *h_tauE_tauE_Mvis_OS;
      TH1F *h_tauE_tauE_bMass_OS;
      TH1F *h_tauE_tauE_bDRl_OS;
      TH1F *h_tauE_tauE_DiTauPt_OS;
      TH1F *h_tauE_tauE_lepDR_OS;
      TH1F *h_tauE_tauE_MET_OS;
      TH2F *h_tauE_tauE_nMatchedEle_OS;
      TH1F *h_tauE_tauE_Trigger_OS;
      TH1F *h_tauE_tauE_BPH_Trigger_OS;
      TH1F *h_tauE_tauE_bPt_OS;
      TH1F *h_tauE_tauE_nE_OS;
      TH1F *h_tauE_tauE_nMu_OS;
      TH1F *h_tauE_tauE_nTau_OS;
      TH1F *h_tauE_tauE_nB_OS;


      TH1F *h_tauE_tauE_probb_OS;
      TH1F *h_tauE_tauE_probbb_OS;
      TH1F *h_tauE_tauE_problepb_OS;
      TH1F *h_tauE_tauE_sum_OS;

////Trigger plots
      TH1F *h_tauE_tauE_Mvis_SS;
      TH1F *h_tauE_tauE_bMass_SS;
      TH1F *h_tauE_tauE_bDRl_SS;
      TH1F *h_tauE_tauE_DiTauPt_SS;
      TH1F *h_tauE_tauE_lepDR_SS;
      TH1F *h_tauE_tauE_MET_SS;
      TH2F *h_tauE_tauE_nMatchedEle_SS;
      TH1F *h_tauE_tauE_Trigger_SS;
      TH1F *h_tauE_tauE_BPH_Trigger_SS;
      TH1F *h_tauE_tauE_bPt_SS;
      TH1F *h_tauE_tauE_nE_SS;
      TH1F *h_tauE_tauE_nMu_SS;
      TH1F *h_tauE_tauE_nTau_SS;
      TH1F *h_tauE_tauE_nB_SS;


      TH1F *h_tauE_tauE_probb_SS;
      TH1F *h_tauE_tauE_probbb_SS;
      TH1F *h_tauE_tauE_problepb_SS;
      TH1F *h_tauE_tauE_sum_SS;

//Tau E Tau Mu
      TH1F *h_tauMu_tauE_Mvis_OS;
      TH1F *h_tauMu_tauE_bMass_OS;
      TH1F *h_tauMu_tauE_bDRl_OS;
      TH1F *h_tauMu_tauE_DiTauPt_OS;
      TH1F *h_tauMu_tauE_lepDR_OS;
      TH1F *h_tauMu_tauE_muIso_OS;
      TH1F *h_tauMu_tauE_MET_OS;
      TH2F *h_tauMu_tauE_nMatchedEle_OS;
      TH2F *h_tauMu_tauE_nMatchedMu_OS;
      TH1F *h_tauMu_tauE_Trigger_OS;
      TH1F *h_tauMu_tauE_BPH_Trigger_OS;
      TH1F *h_tauMu_tauE_bPt_OS;
      TH1F *h_tauMu_tauE_nE_OS;
      TH1F *h_tauMu_tauE_nMu_OS;
      TH1F *h_tauMu_tauE_nTau_OS;
      TH1F *h_tauMu_tauE_nB_OS;

      TH1F *h_tauMu_tauE_probb_OS;
      TH1F *h_tauMu_tauE_probbb_OS;
      TH1F *h_tauMu_tauE_problepb_OS;
      TH1F *h_tauMu_tauE_sum_OS;

//Tau E Tau Mu
      TH1F *h_tauMu_tauE_Mvis_SS;
      TH1F *h_tauMu_tauE_bMass_SS;
      TH1F *h_tauMu_tauE_bDRl_SS;
      TH1F *h_tauMu_tauE_DiTauPt_SS;
      TH1F *h_tauMu_tauE_lepDR_SS;
      TH1F *h_tauMu_tauE_muIso_SS;
      TH1F *h_tauMu_tauE_MET_SS;
      TH2F *h_tauMu_tauE_nMatchedEle_SS;
      TH2F *h_tauMu_tauE_nMatchedMu_SS;
      TH1F *h_tauMu_tauE_Trigger_SS;
      TH1F *h_tauMu_tauE_BPH_Trigger_SS;
      TH1F *h_tauMu_tauE_bPt_SS;
      TH1F *h_tauMu_tauE_nE_SS;
      TH1F *h_tauMu_tauE_nMu_SS;
      TH1F *h_tauMu_tauE_nTau_SS;
      TH1F *h_tauMu_tauE_nB_SS;

      TH1F *h_tauMu_tauE_probb_SS;
      TH1F *h_tauMu_tauE_probbb_SS;
      TH1F *h_tauMu_tauE_problepb_SS;
      TH1F *h_tauMu_tauE_sum_SS;


      ///tauMu Tau Mu

      TH1F *h_tauMu_tauMu_Mvis_OS;
      TH1F *h_tauMu_tauMu_bMass_OS;
      TH1F *h_tauMu_tauMu_bDRl_OS;
      TH1F *h_tauMu_tauMu_DiTauPt_OS;
      TH1F *h_tauMu_tauMu_lepDR_OS;
      TH1F *h_tauMu_tauMu_muIso_OS;
      TH1F *h_tauMu_tauMu_MET_OS;      
      TH2F *h_tauMu_tauMu_nMatchedMu_OS;
      TH1F *h_tauMu_tauMu_Trigger_OS;
      TH1F *h_tauMu_tauMu_BPH_Trigger_OS;
      TH1F *h_tauMu_tauMu_bPt_OS;
      TH1F *h_tauMu_tauMu_nE_OS;
      TH1F *h_tauMu_tauMu_nMu_OS;
      TH1F *h_tauMu_tauMu_nTau_OS;
      TH1F *h_tauMu_tauMu_nB_OS;

   
      TH1F *h_tauMu_tauMu_probb_OS;
      TH1F *h_tauMu_tauMu_probbb_OS;
      TH1F *h_tauMu_tauMu_problepb_OS;
      TH1F *h_tauMu_tauMu_sum_OS;
      ///tauMu Tau Mu

      TH1F *h_tauMu_tauMu_Mvis_SS;
      TH1F *h_tauMu_tauMu_bMass_SS;
      TH1F *h_tauMu_tauMu_bDRl_SS;
      TH1F *h_tauMu_tauMu_DiTauPt_SS;
      TH1F *h_tauMu_tauMu_lepDR_SS;
      TH1F *h_tauMu_tauMu_muIso_SS;
      TH1F *h_tauMu_tauMu_MET_SS;      
      TH2F *h_tauMu_tauMu_nMatchedMu_SS;
      TH1F *h_tauMu_tauMu_Trigger_SS;
      TH1F *h_tauMu_tauMu_BPH_Trigger_SS;
      TH1F *h_tauMu_tauMu_bPt_SS;
      TH1F *h_tauMu_tauMu_nE_SS;
      TH1F *h_tauMu_tauMu_nMu_SS;
      TH1F *h_tauMu_tauMu_nTau_SS;
      TH1F *h_tauMu_tauMu_nB_SS;
   
      TH1F *h_tauMu_tauMu_probb_SS;
      TH1F *h_tauMu_tauMu_probbb_SS;
      TH1F *h_tauMu_tauMu_problepb_SS;
      TH1F *h_tauMu_tauMu_sum_SS;
//Tau E Tau Mu
      TH1F *h_tauHad_tauE_Mvis_OS;
      TH1F *h_tauHad_tauE_bMass_OS;
      TH1F *h_tauHad_tauE_bDRl_OS;
      TH1F *h_tauHad_tauE_DiTauPt_OS;
      TH1F *h_tauHad_tauE_lepDR_OS;
      TH1F *h_tauHad_tauE_muIso_OS;
      TH1F *h_tauHad_tauE_MET_OS;
      TH2F *h_tauHad_tauE_nMatchedTau_OS;
      TH2F *h_tauHad_tauE_nMatchedEle_OS;
      TH1F *h_tauHad_tauE_Trigger_OS;
      TH1F *h_tauHad_tauE_BPH_Trigger_OS;
      TH1F *h_tauHad_tauE_bPt_OS;
      TH1F *h_tauHad_tauE_nE_OS;
      TH1F *h_tauHad_tauE_nMu_OS;
      TH1F *h_tauHad_tauE_nTau_OS;
      TH1F *h_tauHad_tauE_nB_OS;

      TH1F *h_tauHad_tauE_probb_OS;
      TH1F *h_tauHad_tauE_probbb_OS;
      TH1F *h_tauHad_tauE_problepb_OS;
      TH1F *h_tauHad_tauE_sum_OS;

      TH1F *h_tauHad_tauE_Mvis_SS;
      TH1F *h_tauHad_tauE_bMass_SS;
      TH1F *h_tauHad_tauE_bDRl_SS;
      TH1F *h_tauHad_tauE_DiTauPt_SS;
      TH1F *h_tauHad_tauE_lepDR_SS;
      TH1F *h_tauHad_tauE_muIso_SS;
      TH1F *h_tauHad_tauE_MET_SS;
      TH2F *h_tauHad_tauE_nMatchedTau_SS;
      TH2F *h_tauHad_tauE_nMatchedEle_SS;
      TH1F *h_tauHad_tauE_Trigger_SS;
      TH1F *h_tauHad_tauE_BPH_Trigger_SS;
      TH1F *h_tauHad_tauE_bPt_SS;
      TH1F *h_tauHad_tauE_nE_SS;
      TH1F *h_tauHad_tauE_nMu_SS;
      TH1F *h_tauHad_tauE_nTau_SS;
      TH1F *h_tauHad_tauE_nB_SS;

      TH1F *h_tauHad_tauE_probb_SS;
      TH1F *h_tauHad_tauE_probbb_SS;
      TH1F *h_tauHad_tauE_problepb_SS;
      TH1F *h_tauHad_tauE_sum_SS;





//Tau E Tau Mu
      TH1F *h_tauHad_tauMu_Mvis_OS;
      TH1F *h_tauHad_tauMu_bMass_OS;
      TH1F *h_tauHad_tauMu_bDRl_OS;
      TH1F *h_tauHad_tauMu_DiTauPt_OS;
      TH1F *h_tauHad_tauMu_lepDR_OS;
      TH1F *h_tauHad_tauMu_muIso_OS;
      TH1F *h_tauHad_tauMu_MET_OS;
      TH1F *h_tauHad_tauMu_Trigger_OS;
      TH1F *h_tauHad_tauMu_BPH_Trigger_OS;
      TH1F *h_tauHad_tauMu_bPt_OS;
      TH1F *h_tauHad_tauMu_nE_OS;
      TH1F *h_tauHad_tauMu_nMu_OS;
      TH1F *h_tauHad_tauMu_nTau_OS;
      TH1F *h_tauHad_tauMu_nB_OS;

      TH2F *h_tauHad_tauMu_nMatchedMu_OS;
      TH1F *h_tauHad_tauMu_probb_OS;
      TH1F *h_tauHad_tauMu_probbb_OS;
      TH1F *h_tauHad_tauMu_problepb_OS;
      TH1F *h_tauHad_tauMu_sum_OS;
      
      TH1F *h_tauHad_tauMu_Mvis_SS;
      TH1F *h_tauHad_tauMu_bMass_SS;
      TH1F *h_tauHad_tauMu_bDRl_SS;
      TH1F *h_tauHad_tauMu_DiTauPt_SS;
      TH1F *h_tauHad_tauMu_lepDR_SS;
      TH1F *h_tauHad_tauMu_muIso_SS;
      TH1F *h_tauHad_tauMu_MET_SS;
      TH1F *h_tauHad_tauMu_Trigger_SS;
      TH1F *h_tauHad_tauMu_BPH_Trigger_SS;
      TH1F *h_tauHad_tauMu_bPt_SS;
      TH1F *h_tauHad_tauMu_nE_SS;
      TH1F *h_tauHad_tauMu_nMu_SS;
      TH1F *h_tauHad_tauMu_nTau_SS;
      TH1F *h_tauHad_tauMu_nB_SS;

      TH2F *h_tauHad_tauMu_nMatchedMu_SS;
      TH1F *h_tauHad_tauMu_probb_SS;
      TH1F *h_tauHad_tauMu_probbb_SS;
      TH1F *h_tauHad_tauMu_problepb_SS;
      TH1F *h_tauHad_tauMu_sum_SS;


//Tau E Tau Mu
      TH1F *h_gen_tauHad_tauE_Mvis;
      TH1F *h_gen_tauHad_tauE_bMass;
      TH1F *h_gen_tauHad_tauE_bDRl;
      TH1F *h_gen_tauHad_tauE_DiTauPt;
      TH1F *h_gen_tauHad_tauE_lepDR;
      TH1F *h_gen_tauHad_tauE_muIso;
      TH1F *h_gen_tauHad_tauE_MET;

      TH1F *h_gen_tauHad_tauE_probb;
      TH1F *h_gen_tauHad_tauE_probbb;
      TH1F *h_gen_tauHad_tauE_problepb;
      TH1F *h_gen_tauHad_tauE_sum;



//Tau E Tau Mu
      TH1F *h_gen_tauHad_tauMu_Mvis;
      TH1F *h_gen_tauHad_tauMu_bMass;
      TH1F *h_gen_tauHad_tauMu_bDRl;
      TH1F *h_gen_tauHad_tauMu_DiTauPt;
      TH1F *h_gen_tauHad_tauMu_lepDR;
      TH1F *h_gen_tauHad_tauMu_muIso;
      TH1F *h_gen_tauHad_tauMu_MET;

      TH1F *h_gen_tauHad_tauMu_probb;
      TH1F *h_gen_tauHad_tauMu_probbb;
      TH1F *h_gen_tauHad_tauMu_problepb;
      TH1F *h_gen_tauHad_tauMu_sum;

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
   //taus_{consumes<std::vector<pat::Tau>> (edm::InputTag(std::string("slimmedTausMuonCleaned")))},
   taus_{consumes<std::vector<pat::Tau>> (edm::InputTag(std::string("slimmedTaus")))},
   std_electrons_{consumes< std::vector<pat::Electron> >(edm::InputTag(std::string("slimmedElectrons")))},
   lowPt_electrons_{consumes< std::vector<pat::Electron> >(edm::InputTag(std::string("slimmedLowPtElectrons")))},
   triggerBits_{consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT"))},
   jetFlavorMatching_{consumes<reco::JetFlavourInfoMatchingCollection>(edm::InputTag(std::string("slimmedGenJetsFlavourInfos")))},
   triggerPrescales_{consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger","l1max","PAT"))}
{
   //now do what ever initialization is needed
      edm::Service<TFileService> fs;

      h_std_electron_selection=fs->make<TH1F>("h_std_electron_selection", "StdElectron Cuts", 3,0,3);
      h_std_electron_selection->SetCanExtend(TH1::kAllAxes);
      h_std_electron_selection->Sumw2();
      h_muon_selection=fs->make<TH1F>("h_muon_selection", "Muon Cuts", 3,0,3);
      h_muon_selection->SetCanExtend(TH1::kAllAxes);
      h_muon_selection->Sumw2();
      h_tau_selection=fs->make<TH1F>("h_tau_selection", "Tau Cuts", 3,0,3);
      h_tau_selection->SetCanExtend(TH1::kAllAxes);
      h_tau_selection->Sumw2();
      h_bJet_selection=fs->make<TH1F>("h_bJet_selection", "b Jet Cuts", 3,0,3);
      h_bJet_selection->SetCanExtend(TH1::kAllAxes);
      h_bJet_selection->Sumw2();
     

      h_nStdEle=fs->make<TH1F>("h_nStdEle", "nTau Ele per Event (all)", 4,0,4);
      h_nStdEle->Sumw2();
      h_nStdEle_taum=fs->make<TH1F>("h_nStdEle_taum", "nTau Ele per Event (Tau Matched)", 4,0,4);
      h_nStdEle_taum->Sumw2(); 
      h_nStdEle_bm=fs->make<TH1F>("h_nStdEle_bm", "nTau Ele per Event (b Matched)", 4,0,4);
      h_nStdEle_bm->Sumw2(); 
      h_StdEle_taum_pt=fs->make<TH1F>("h_StdEle_taum_pt", "Tau matched Electron Pt", 40,0,200);
      h_StdEle_taum_pt->Sumw2();
      h_StdEle_bm_pt=fs->make<TH1F>("h_StdEle_bm_pt", "b matched Electron Pt",40,0,200);
      h_StdEle_bm_pt->Sumw2();   

      h_nMu=fs->make<TH1F>("h_nMu", "nTau Mu per Event (all)", 4,0,4);
      h_nMu->Sumw2();
      h_nMu_taum=fs->make<TH1F>("h_nMu_taum", "nTau Mu per Event (Tau Matched)", 4,0,4);
      h_nMu_taum->Sumw2();   
      h_nMu_bm=fs->make<TH1F>("h_nMu_bm", "nTau Mu per Event (b Matched)", 4,0,4);
      h_nMu_bm->Sumw2();  
      h_Mu_taum_pt=fs->make<TH1F>("h_Mu_taum_pt", "Tau matched Muon Pt", 40,0,200);
      h_Mu_taum_pt->Sumw2();
      h_Mu_bm_pt=fs->make<TH1F>("h_Mu_bm_pt", "b matched Muon Pt",40,0,200);
      h_Mu_bm_pt->Sumw2();  

      h_nbJet=fs->make<TH1F>("h_nbJet", "bJet per Event (all)", 4,0,4);
      h_nbJet->Sumw2();
      h_nbJet_m=fs->make<TH1F>("h_nbJet_m", "bJet per Event (Matched)", 4,0,4);
      h_nbJet_m->Sumw2();  
      h_nbJet_f=fs->make<TH1F>("h_nbJet_f", "bJet per Event (!Matched)", 4,0,4);
      h_nbJet_f->Sumw2();  

      
      h_Triggers=fs->make<TH1F>("h_Triggers", "Triggers (no cuts applied)", 3,0,3);
      h_Triggers->SetCanExtend(TH1::kAllAxes);
      h_Triggers->Sumw2();

// b Discriminator
      h_deepFlavor_probb=fs->make<TH1F>("h_deepFlavor_probb", "pfDeepFlavourJetTags:probb All", 20,0,1);
      h_deepFlavor_probb->Sumw2();
      h_deepFlavor_probbb=fs->make<TH1F>("h_deepFlavor_probbb", "pfDeepFlavourJetTags:probbb All", 20,0,1);
      h_deepFlavor_probbb->Sumw2();
      h_deepFlavor_problepb=fs->make<TH1F>("h_deepFlavor_problepb", "pfDeepFlavourJetTags:problepb All", 20,0,1);
      h_deepFlavor_problepb->Sumw2();
      h_deepFlavor_sum=fs->make<TH1F>("h_deepFlavor_sum", "Total Btag All", 20,0,1);
      h_deepFlavor_sum->Sumw2();
      h_deepFlavor_probbvsprobbb=fs->make<TH2F>("h_deepFlavor_probbvsprobbb", "prob b vs prob bb All", 20,0,1,20,0,1);
      h_deepFlavor_probbvsprobbb->Sumw2();

      h_deepFlavor_probb_m=fs->make<TH1F>("h_deepFlavor_probb_m", "pfDeepFlavourJetTags:probb All (matched)", 20,0,1);
      h_deepFlavor_probb_m->Sumw2();
      h_deepFlavor_probbb_m=fs->make<TH1F>("h_deepFlavor_probbb_m", "pfDeepFlavourJetTags:probbb All (matched)", 20,0,1);
      h_deepFlavor_probbb_m->Sumw2();
      h_deepFlavor_problepb_m=fs->make<TH1F>("h_deepFlavor_probbb_m", "pfDeepFlavourJetTags:problepb All (matched)", 20,0,1);
      h_deepFlavor_problepb_m->Sumw2();
      h_deepFlavor_sum_m=fs->make<TH1F>("h_deepFlavor_sum_m", "Total Btag All (matched)", 20,0,1);
      h_deepFlavor_sum_m->Sumw2();
      h_deepFlavor_probbvsprobbb_m=fs->make<TH2F>("h_deepFlavor_probbvsprobbb_m", "prob b vs prob bb All (matched)", 20,0,1,20,0,1);
      h_deepFlavor_probbvsprobbb_m->Sumw2();


      h_deepFlavor_probb_f=fs->make<TH1F>("h_deepFlavor_probb_f", "pfDeepFlavourJetTags:probb All (fake)", 20,0,1);
      h_deepFlavor_probb_f->Sumw2();
      h_deepFlavor_probbb_f=fs->make<TH1F>("h_deepFlavor_probbb_f", "pfDeepFlavourJetTags:probbb All (fake)", 20,0,1);
      h_deepFlavor_probbb_f->Sumw2();
      h_deepFlavor_problepb_f=fs->make<TH1F>("h_deepFlavor_problepb_f", "pfDeepFlavourJetTags:problepb All (fake)", 20,0,1);
      h_deepFlavor_problepb_f->Sumw2();
      h_deepFlavor_sum_f=fs->make<TH1F>("h_deepFlavor_sum_f", "Total Btag All (fake)", 20,0,1);
      h_deepFlavor_sum_f->Sumw2();
      h_deepFlavor_probbvsprobbb_f=fs->make<TH2F>("h_deepFlavor_probbvsprobbb_f", "prob b vs prob bb All(fake)", 20,0,1,20,0,1);
      h_deepFlavor_probbvsprobbb_f->Sumw2();

      h_deepFlavor_leading_probb=fs->make<TH1F>("h_deepFlavor_leading_probb", "pfDeepFlavourJetTags:probb Leading", 20,0,1);
      h_deepFlavor_leading_probb->Sumw2();
      h_deepFlavor_leading_probbb=fs->make<TH1F>("h_deepFlavor_leading_probbb", "pfDeepFlavourJetTags:probbb Leading", 20,0,1);
      h_deepFlavor_leading_probbb->Sumw2();
      h_deepFlavor_leading_problepb=fs->make<TH1F>("h_deepFlavor_leading_problepb", "pfDeepFlavourJetTags:problepb Leading", 20,0,1);
      h_deepFlavor_leading_problepb->Sumw2();
      h_deepFlavor_leading_sum=fs->make<TH1F>("h_deepFlavor_leading_sum", "Total Btag Leading", 20,0,1);
      h_deepFlavor_leading_sum->Sumw2();
      h_deepFlavor_leading_probbvsprobbb=fs->make<TH2F>("h_deepFlavor_leading_probbvsprobbb", "prob b vs prob bb Leading", 20,0,1,20,0,1);
      h_deepFlavor_leading_probbvsprobbb->Sumw2();

      h_deepFlavor_leading_probb_m=fs->make<TH1F>("h_deepFlavor_leading_probb_m", "pfDeepFlavourJetTags:probb Leading (matched)", 20,0,1);
      h_deepFlavor_leading_probb_m->Sumw2();
      h_deepFlavor_leading_probbb_m=fs->make<TH1F>("h_deepFlavor_leading_probbb_m", "pfDeepFlavourJetTags:probbb Leading (matched)", 20,0,1);
      h_deepFlavor_leading_probbb_m->Sumw2();
      h_deepFlavor_leading_problepb_m=fs->make<TH1F>("h_deepFlavor_leading_problepb_m", "pfDeepFlavourJetTags:problepb Leading (matched)", 20,0,1);
      h_deepFlavor_leading_problepb_m->Sumw2();
      h_deepFlavor_leading_sum_m=fs->make<TH1F>("h_deepFlavor_leading_sum_m", "Total Btag Leading (matched)", 20,0,1);
      h_deepFlavor_leading_sum_m->Sumw2();
      h_deepFlavor_leading_probbvsprobbb_m=fs->make<TH2F>("h_deepFlavor_leading_probbvsprobbb_m", "prob b vs prob bb Leading (matched)", 20,0,1,20,0,1);
      h_deepFlavor_leading_probbvsprobbb_m->Sumw2();


      h_deepFlavor_leading_probb_f=fs->make<TH1F>("h_deepFlavor_leading_probb_f", "pfDeepFlavourJetTags:probb Leading (fake)", 20,0,1);
      h_deepFlavor_leading_probb_f->Sumw2();
      h_deepFlavor_leading_probbb_f=fs->make<TH1F>("h_deepFlavor_leading_probbb_f", "pfDeepFlavourJetTags:probbb Leading (fake)", 20,0,1);
      h_deepFlavor_leading_probbb_f->Sumw2();
      h_deepFlavor_leading_problepb_f=fs->make<TH1F>("h_deepFlavor_leading_problepb_f", "pfDeepFlavourJetTags:problepb Leading (fake)", 20,0,1);
      h_deepFlavor_leading_problepb_f->Sumw2();
      h_deepFlavor_leading_sum_f=fs->make<TH1F>("h_deepFlavor_leading_sum_f", "Total Btag Leading (fake)", 20,0,1);
      h_deepFlavor_leading_sum_f->Sumw2();
      h_deepFlavor_leading_probbvsprobbb_f=fs->make<TH2F>("h_deepFlavor_leading_probbvsprobbb_f", "prob b vs prob bb Leading(fake) ", 20,0,1,20,0,1);
      h_deepFlavor_leading_probbvsprobbb_f->Sumw2();

      h_deepFlavor_subLeading_probb=fs->make<TH1F>("h_deepFlavor_subLeading_probb", "pfDeepFlavourJetTags:probb SubLeading", 20,0,1);
      h_deepFlavor_subLeading_probb->Sumw2();
      h_deepFlavor_subLeading_probbb=fs->make<TH1F>("h_deepFlavor_subLeading_probbb", "pfDeepFlavourJetTags:probbb SubLeading", 20,0,1);
      h_deepFlavor_subLeading_probbb->Sumw2();
      h_deepFlavor_subLeading_problepb=fs->make<TH1F>("h_deepFlavor_subLeading_problepb", "pfDeepFlavourJetTags:problepb SubLeading", 20,0,1);
      h_deepFlavor_subLeading_problepb->Sumw2();
      h_deepFlavor_subLeading_sum=fs->make<TH1F>("h_deepFlavor_subLeading_sum", "Total Btag SubLeading", 20,0,1);
      h_deepFlavor_subLeading_sum->Sumw2();
      h_deepFlavor_subLeading_probbvsprobbb=fs->make<TH2F>("h_deepFlavor_subLeading_probbvsprobbb", "prob b vs prob bb Subleading", 20,0,1,20,0,1);
      h_deepFlavor_subLeading_probbvsprobbb->Sumw2();

      h_deepFlavor_subLeading_probb_m=fs->make<TH1F>("h_deepFlavor_subLeading_probb_m", "pfDeepFlavourJetTags:probb SubLeading (matched)", 20,0,1);
      h_deepFlavor_subLeading_probb_m->Sumw2();
      h_deepFlavor_subLeading_probbb_m=fs->make<TH1F>("h_deepFlavor_subLeading_probbb_m", "pfDeepFlavourJetTags:probbb SubLeading (matched)", 20,0,1);
      h_deepFlavor_subLeading_probbb_m->Sumw2();
      h_deepFlavor_subLeading_problepb_m=fs->make<TH1F>("h_deepFlavor_subLeading_problepb_m", "pfDeepFlavourJetTags:problepb SubLeading (matched)", 20,0,1);
      h_deepFlavor_subLeading_problepb_m->Sumw2();
      h_deepFlavor_subLeading_sum_m=fs->make<TH1F>("h_deepFlavor_subLeading_sum_m", "Total Btag SubLeading (matched)", 20,0,1);
      h_deepFlavor_subLeading_sum_m->Sumw2();
      h_deepFlavor_subLeading_probbvsprobbb_m=fs->make<TH2F>("h_deepFlavor_subLeading_probbvsprobbb_m", "prob b vs prob bb Subleading (matched)", 20,0,1,20,0,1);
      h_deepFlavor_subLeading_probbvsprobbb_m->Sumw2();

      h_deepFlavor_subLeading_probb_f=fs->make<TH1F>("h_deepFlavor_subLeading_probb_f", "pfDeepFlavourJetTags:probb SubLeading (fake)", 20,0,1);
      h_deepFlavor_subLeading_probb_f->Sumw2();
      h_deepFlavor_subLeading_probbb_f=fs->make<TH1F>("h_deepFlavor_subLeading_probbb_f", "pfDeepFlavourJetTags:probbb SubLeading (fake)", 20,0,1);
      h_deepFlavor_subLeading_probbb_f->Sumw2();
      h_deepFlavor_subLeading_problepb_f=fs->make<TH1F>("h_deepFlavor_subLeading_problepb_f", "pfDeepFlavourJetTags:problepb SubLeading (fake)", 20,0,1);
      h_deepFlavor_subLeading_problepb_f->Sumw2();
      h_deepFlavor_subLeading_sum_f=fs->make<TH1F>("h_deepFlavor_subLeading_sum_f", "Total Btag SubLeading (fake)", 20,0,1);
      h_deepFlavor_subLeading_sum_f->Sumw2();
      h_deepFlavor_subLeading_probbvsprobbb_f=fs->make<TH2F>("h_deepFlavor_subLeading_probbvsprobbb_f", "prob b vs prob bb Subleading (fake)", 20,0,1,20,0,1);
      h_deepFlavor_subLeading_probbvsprobbb_f->Sumw2();

      h_deepFlavor_Merged_probb=fs->make<TH1F>("h_deepFlavor_Merged_probb", "pfDeepFlavourJetTags:probb nJets = 1", 20,0,1);
      h_deepFlavor_Merged_probb->Sumw2();
      h_deepFlavor_Merged_probbb=fs->make<TH1F>("h_deepFlavor_Merged_probbb", "pfDeepFlavourJetTags:probbb nJets = 1", 20,0,1);
      h_deepFlavor_Merged_probbb->Sumw2();
      h_deepFlavor_Merged_problepb=fs->make<TH1F>("h_deepFlavor_Merged_problepb", "pfDeepFlavourJetTags:problepb nJets = 1", 20,0,1);
      h_deepFlavor_Merged_problepb->Sumw2();
      h_deepFlavor_Merged_sum=fs->make<TH1F>("h_deepFlavor_Merged_sum", "Total Btag nJets = 1", 20,0,1);
      h_deepFlavor_Merged_sum->Sumw2();
      h_deepFlavor_Merged_probbvsprobbb=fs->make<TH2F>("h_deepFlavor_Merged_probbvsprobbb", "prob b vs prob bb nJets = 1", 20,0,1,20,0,1);
      h_deepFlavor_Merged_probbvsprobbb->Sumw2();

      h_deepFlavor_Merged_probb_m=fs->make<TH1F>("h_deepFlavor_Merged_probb_m", "pfDeepFlavourJetTags:probb nJets = 1 (matched)", 20,0,1);
      h_deepFlavor_Merged_probb_m->Sumw2();
      h_deepFlavor_Merged_probbb_m=fs->make<TH1F>("h_deepFlavor_Merged_probbb_m", "pfDeepFlavourJetTags:probbb nJets = 1 (matched)", 20,0,1);
      h_deepFlavor_Merged_probbb_m->Sumw2();
      h_deepFlavor_Merged_problepb_m=fs->make<TH1F>("h_deepFlavor_Merged_probbb_m", "pfDeepFlavourJetTags:problepb nJets = 1 (matched)", 20,0,1);
      h_deepFlavor_Merged_problepb_m->Sumw2();
      h_deepFlavor_Merged_sum_m=fs->make<TH1F>("h_deepFlavor_Merged_sum_m", "Total Btag nJets = 1 (matched)", 20,0,1);
      h_deepFlavor_Merged_sum_m->Sumw2();
      h_deepFlavor_Merged_probbvsprobbb_m=fs->make<TH2F>("h_deepFlavor_Merged_probbvsprobbb_m", "prob b vs prob bb nJets = 1 (matched)", 20,0,1,20,0,1);
      h_deepFlavor_Merged_probbvsprobbb_m->Sumw2();


      h_deepFlavor_Leptonic_probb=fs->make<TH1F>("h_deepFlavor_Leptonic_probb", "pfDeepFlavourJetTags:probb Leptonic B", 20,0,1);
      h_deepFlavor_Leptonic_probb->Sumw2();
      h_deepFlavor_Leptonic_probbb=fs->make<TH1F>("h_deepFlavor_Leptonic_probbb", "pfDeepFlavourJetTags:probbb Leptonic B", 20,0,1);
      h_deepFlavor_Leptonic_probbb->Sumw2();
      h_deepFlavor_Leptonic_problepb=fs->make<TH1F>("h_deepFlavor_Leptonic_problepb", "pfDeepFlavourJetTags:problepb Leptonic B", 20,0,1);
      h_deepFlavor_Leptonic_problepb->Sumw2();
      h_deepFlavor_Leptonic_sum=fs->make<TH1F>("h_deepFlavor_Leptonic_sum", "Total Btag Leptonic B", 20,0,1);
      h_deepFlavor_Leptonic_sum->Sumw2();
      h_deepFlavor_Leptonic_probbvsprobbb=fs->make<TH2F>("h_deepFlavor_Leptonic_probbvsprobbb", "prob b vs prob bb Leptonic B", 20,0,1,20,0,1);
      h_deepFlavor_Leptonic_probbvsprobbb->Sumw2();

      h_deepFlavor_probb_m_Trigger=fs->make<TH1F>("h_deepFlavor_probb_m_Trigger", "pfDeepFlavourJetTags:probb (BTagMu_AK4DiJet20_Mu5) (matched)", 20,0,1);
      h_deepFlavor_probb_m_Trigger->Sumw2();
      h_deepFlavor_probbb_m_Trigger=fs->make<TH1F>("h_deepFlavor_probbb_m_Trigger", "pfDeepFlavourJetTags:probbb All (BTagMu_AK4DiJet20_Mu5)(matched)", 20,0,1);
      h_deepFlavor_probbb_m_Trigger->Sumw2();
      h_deepFlavor_problepb_m_Trigger=fs->make<TH1F>("h_deepFlavor_probbb_m_Trigger", "pfDeepFlavourJetTags:problepb All (BTagMu_AK4DiJet20_Mu5)(matched)", 20,0,1);
      h_deepFlavor_problepb_m_Trigger->Sumw2();
      h_deepFlavor_sum_m_Trigger=fs->make<TH1F>("h_deepFlavor_sum_m_Trigger", "Total Btag All (BTagMu_AK4DiJet20_Mu5)(matched)", 20,0,1);
      h_deepFlavor_sum_m_Trigger->Sumw2();
      h_deepFlavor_probbvsprobbb_m_Trigger=fs->make<TH2F>("h_deepFlavor_probbvsprobbb_m_Trigger", "prob b vs prob bb All (BTagMu_AK4DiJet20_Mu5)(matched)", 20,0,1,20,0,1);
      h_deepFlavor_probbvsprobbb_m_Trigger->Sumw2();
///Event Selection Plots
      h_nEvent_tauE_tauE=fs->make<TH1F>("h_nEvent_tauE_tauE", "NEvent tauE tauE",3,0,3);
      h_nEvent_tauE_tauE->SetCanExtend(TH1::kAllAxes);
      h_nEvent_tauE_tauE->Sumw2();

      h_nEvent_tauE_tauMu=fs->make<TH1F>("h_nEvent_tauE_tauMu", "NEvent tauE tauMu",4,0,4);
      h_nEvent_tauE_tauMu->SetCanExtend(TH1::kAllAxes);
      h_nEvent_tauE_tauMu->Sumw2();

      h_nEvent_tauMu_tauMu=fs->make<TH1F>("h_nEvent_tauMu_tauMu", "NEvent tauMu tauMu",3,0,3);
      h_nEvent_tauMu_tauMu->SetCanExtend(TH1::kAllAxes);
      h_nEvent_tauMu_tauMu->Sumw2();

      h_nEvent_tauHad_tauE=fs->make<TH1F>("h_nEvent_tauHad_tauE", "NEvent TauHad tauE",3,0,3);
      h_nEvent_tauHad_tauE->SetCanExtend(TH1::kAllAxes);
      h_nEvent_tauHad_tauE->Sumw2();

      h_nEvent_tauHad_tauMu=fs->make<TH1F>("h_nEvent_tauHad_tauMu", "NEvent TauHad tauMu",4,0,4);
      h_nEvent_tauHad_tauMu->SetCanExtend(TH1::kAllAxes);
      h_nEvent_tauHad_tauMu->Sumw2();


      h_nEvent_gen_tauHad_tauE=fs->make<TH1F>("h_nEvent_gen_tauHad_tauE", "NEvent Gen TauHad tauE",3,0,3);
      h_nEvent_gen_tauHad_tauE->SetCanExtend(TH1::kAllAxes);
      h_nEvent_gen_tauHad_tauE->Sumw2();

      h_nEvent_gen_tauHad_tauMu=fs->make<TH1F>("h_nEvent_gen_tauHad_tauMu", "NEvent Gen TauHad tauMu",4,0,4);
      h_nEvent_gen_tauHad_tauMu->SetCanExtend(TH1::kAllAxes);
      h_nEvent_gen_tauHad_tauMu->Sumw2();
///Event PT Plots







//bDiscriminator








//tauE tauE
      h_tauE_tauE_Mvis_OS=fs->make<TH1F>("h_tauE_tauE_Mvis_OS", "E||E Mvis (OS)",20,3,13);
      h_tauE_tauE_Mvis_OS->Sumw2();
      h_tauE_tauE_bMass_OS=fs->make<TH1F>("h_tauE_tauE_bMass_OS", "E||E bMass (OS)",40,0,20);
      h_tauE_tauE_bMass_OS->Sumw2();
      h_tauE_tauE_bDRl_OS=fs->make<TH1F>("h_tauE_tauE_bDRl_OS", "E||E bDRl (OS)",40,0,5);
      h_tauE_tauE_bDRl_OS->Sumw2();
      h_tauE_tauE_lepDR_OS=fs->make<TH1F>("h_tauE_tauE_lepDR_OS", "E||E lepDR (OS)",20,0,1.2);
      h_tauE_tauE_lepDR_OS->Sumw2();
      h_tauE_tauE_DiTauPt_OS=fs->make<TH1F>("h_tauE_tauE_DiTauPt_OS", "E||E DiTauPt(OS)",50,0,200);
      h_tauE_tauE_DiTauPt_OS->Sumw2();
      h_tauE_tauE_MET_OS=fs->make<TH1F>("h_tauE_tauE_MET_OS", "E||E MET (OS)",40,0,250);
      h_tauE_tauE_MET_OS->Sumw2();
      h_tauE_tauE_Trigger_OS=fs->make<TH1F>("h_tauE_tauE_Trigger_OS", "E||E Triggers (OS)",18,0,18);
      h_tauE_tauE_Trigger_OS->SetCanExtend(TH1::kAllAxes);
      h_tauE_tauE_Trigger_OS->Sumw2();
      h_tauE_tauE_BPH_Trigger_OS=fs->make<TH1F>("h_tauE_tauE_BPH_Trigger_OS", "E||E BPH Triggers (OS)",25,0,25);
      h_tauE_tauE_BPH_Trigger_OS->SetCanExtend(TH1::kAllAxes);
      h_tauE_tauE_BPH_Trigger_OS->Sumw2();
      h_tauE_tauE_bPt_OS=fs->make<TH1F>("h_tauE_tauE_bPt_OS", "E||E bPt (OS)",20,0,120);
      h_tauE_tauE_bPt_OS->Sumw2();
      h_tauE_tauE_nE_OS=fs->make<TH1F>("h_tauE_tauE_nE_OS", "E||E nE (OS)",5,0,5);
      h_tauE_tauE_nE_OS->Sumw2();
      h_tauE_tauE_nTau_OS=fs->make<TH1F>("h_tauE_tauE_nTau_OS", "E||E nTau (OS)",5,0,5);
      h_tauE_tauE_nTau_OS->Sumw2();
      h_tauE_tauE_nB_OS=fs->make<TH1F>("h_tauE_tauE_nB_OS", "E||E nB (OS)",5,0,5);
      h_tauE_tauE_nB_OS->Sumw2();
      h_tauE_tauE_nMu_OS=fs->make<TH1F>("h_tauE_tauE_nMu_OS", "E||E nMu (OS)",5,0,5);
      h_tauE_tauE_nMu_OS->Sumw2();

      h_tauE_tauE_nMatchedEle_OS=fs->make<TH2F>("h_tauE_tauE_nMatchedEle_OS", "E||E n Matched Electron(bx,ty) (OS)",4,0,4,4,0,4);
      h_tauE_tauE_nMatchedEle_OS->Sumw2(); 
      h_tauE_tauE_probb_OS=fs->make<TH1F>("h_tauE_tauE_probb_OS", "pfDeepFlavourJetTags:probb OS", 20,0,1);
      h_tauE_tauE_probb_OS->Sumw2();
      h_tauE_tauE_probbb_OS=fs->make<TH1F>("h_tauE_tauE_probbb_OS", "pfDeepFlavourJetTags:probbb OS", 20,0,1);
      h_tauE_tauE_probbb_OS->Sumw2();
      h_tauE_tauE_problepb_OS=fs->make<TH1F>("h_tauE_tauE_problepb_OS", "pfDeepFlavourJetTags:problepb OS", 20,0,1);
      h_tauE_tauE_problepb_OS->Sumw2();
      h_tauE_tauE_sum_OS=fs->make<TH1F>("h_tauE_tauE_sum_OS", "Total Btag OS", 20,0,1);
      h_tauE_tauE_sum_OS->Sumw2();

      h_tauE_tauE_Mvis_SS=fs->make<TH1F>("h_tauE_tauE_Mvis_SS", "E||E Mvis (SS)",20,3,13);
      h_tauE_tauE_Mvis_SS->Sumw2();
      h_tauE_tauE_bMass_SS=fs->make<TH1F>("h_tauE_tauE_bMass_SS", "E||E bMass (SS)",40,0,20);
      h_tauE_tauE_bMass_SS->Sumw2();
      h_tauE_tauE_bDRl_SS=fs->make<TH1F>("h_tauE_tauE_bDRl_SS", "E||E bDRl (SS)",40,0,5);
      h_tauE_tauE_bDRl_SS->Sumw2();
      h_tauE_tauE_lepDR_SS=fs->make<TH1F>("h_tauE_tauE_lepDR_SS", "E||E lepDR (SS)",20,0,1.2);
      h_tauE_tauE_lepDR_SS->Sumw2();
      h_tauE_tauE_DiTauPt_SS=fs->make<TH1F>("h_tauE_tauE_DiTauPt_SS", "E||E DiTauPt (SS)",50,0,200);
      h_tauE_tauE_DiTauPt_SS->Sumw2();
      h_tauE_tauE_MET_SS=fs->make<TH1F>("h_tauE_tauE_MET_SS", "E||E MET (SS)",40,0,250);
      h_tauE_tauE_MET_SS->Sumw2();
      h_tauE_tauE_Trigger_SS=fs->make<TH1F>("h_tauE_tauE_Trigger_SS", "E||E Triggers (SS)",18,0,18);
      h_tauE_tauE_Trigger_SS->SetCanExtend(TH1::kAllAxes);
      h_tauE_tauE_Trigger_SS->Sumw2();
      h_tauE_tauE_BPH_Trigger_SS=fs->make<TH1F>("h_tauE_tauE_BPH_Trigger_SS", "E||E BPH Triggers (SS)",25,0,25);
      h_tauE_tauE_BPH_Trigger_SS->SetCanExtend(TH1::kAllAxes);
      h_tauE_tauE_BPH_Trigger_SS->Sumw2();
      h_tauE_tauE_bPt_SS=fs->make<TH1F>("h_tauE_tauE_bPt_SS", "E||E bPt (SS)",20,0,120);
      h_tauE_tauE_bPt_SS->Sumw2();
      h_tauE_tauE_nE_SS=fs->make<TH1F>("h_tauE_tauE_nE_SS", "E||E nE (SS)",5,0,5);
      h_tauE_tauE_nE_SS->Sumw2();
      h_tauE_tauE_nTau_SS=fs->make<TH1F>("h_tauE_tauE_nTau_SS", "E||E nTau (SS)",5,0,5);
      h_tauE_tauE_nTau_SS->Sumw2();
      h_tauE_tauE_nB_SS=fs->make<TH1F>("h_tauE_tauE_nB_SS", "E||E nB (SS)",5,0,5);
      h_tauE_tauE_nB_SS->Sumw2();
      h_tauE_tauE_nMu_SS=fs->make<TH1F>("h_tauE_tauE_nMu_SS", "E||E nMu (SS)",5,0,5);
      h_tauE_tauE_nMu_SS->Sumw2();

      h_tauE_tauE_nMatchedEle_SS=fs->make<TH2F>("h_tauE_tauE_nMatchedEle_SS", "E||E n Matched Electron(bx,ty) (SS)",4,0,4,4,0,4);
      h_tauE_tauE_nMatchedEle_SS->Sumw2(); 
      h_tauE_tauE_probb_SS=fs->make<TH1F>("h_tauE_tauE_probb_SS", "E||E pfDeepFlavourJetTags:probb SS", 20,0,1);
      h_tauE_tauE_probb_SS->Sumw2();
      h_tauE_tauE_probbb_SS=fs->make<TH1F>("h_tauE_tauE_probbb_SS", "E||E pfDeepFlavourJetTags:probbb SS", 20,0,1);
      h_tauE_tauE_probbb_SS->Sumw2();
      h_tauE_tauE_problepb_SS=fs->make<TH1F>("h_tauE_tauE_problepb_SS", "E||E pfDeepFlavourJetTags:problepb SS", 20,0,1);
      h_tauE_tauE_problepb_SS->Sumw2();
      h_tauE_tauE_sum_SS=fs->make<TH1F>("h_tauE_tauE_sum_SS", "E||E Total Btag SS", 20,0,1);
      h_tauE_tauE_sum_SS->Sumw2();


//Mvis
      h_tauMu_tauE_Mvis_OS=fs->make<TH1F>("h_tauMu_tauE_Mvis_OS", "Mu||E Mvis OS",20,3,13);
      h_tauMu_tauE_Mvis_OS->Sumw2();
      h_tauMu_tauE_bMass_OS=fs->make<TH1F>("h_tauMu_tauE_bMass_OS", "Mu||E bMass OS",45,0,20);
      h_tauMu_tauE_bMass_OS->Sumw2();
      h_tauMu_tauE_bDRl_OS=fs->make<TH1F>("h_tauMu_tauE_bDRl_OS", "Mu||E bDRl OS",45,0,5);
      h_tauMu_tauE_bDRl_OS->Sumw2();
      h_tauMu_tauE_lepDR_OS=fs->make<TH1F>("h_tauMu_tauE_lepDR_OS", "Mu||E lepDR OS",20,0,1.2);
      h_tauMu_tauE_lepDR_OS->Sumw2();
      h_tauMu_tauE_muIso_OS=fs->make<TH1F>("h_tauMu_tauE_muIso_OS", "Mu||E muIso OS",20,0,1);
      h_tauMu_tauE_muIso_OS->Sumw2();
      h_tauMu_tauE_MET_OS=fs->make<TH1F>("h_tauMu_tauE_MET_OS", "Mu||E MET OS",40,0,250);
      h_tauMu_tauE_MET_OS->Sumw2();
      h_tauMu_tauE_DiTauPt_OS=fs->make<TH1F>("h_tauMu_tauE_DiTauPt_OS", "Mu||E DiTauPt OS",50,0,200);
      h_tauMu_tauE_DiTauPt_OS->Sumw2();
      h_tauMu_tauE_nMatchedEle_OS=fs->make<TH2F>("h_tauMu_tauE_nMatchedEle_OS", "Mu||E n Matched Electron(bx,ty) OS",4,0,4,4,0,4);
      h_tauMu_tauE_nMatchedEle_OS->Sumw2(); 
      h_tauMu_tauE_nMatchedMu_OS=fs->make<TH2F>("h_tauMu_tauE_nMatchedMu_OS", "Mu||E n Matched Muons(bx,ty) OS",4,0,4,4,0,4);
      h_tauMu_tauE_nMatchedMu_OS->Sumw2(); 
      h_tauMu_tauE_Trigger_OS=fs->make<TH1F>("h_tauMu_tauE_Trigger_OS", "Mu||E Triggers (OS)",18,0,18);
      h_tauMu_tauE_Trigger_OS->SetCanExtend(TH1::kAllAxes);
      h_tauMu_tauE_Trigger_OS->Sumw2(); 
      h_tauMu_tauE_BPH_Trigger_OS=fs->make<TH1F>("h_tauMu_tauE_BPH_Trigger_OS", "Mu||E BPH Triggers (OS)",25,0,25);
      h_tauMu_tauE_BPH_Trigger_OS->SetCanExtend(TH1::kAllAxes);
      h_tauMu_tauE_BPH_Trigger_OS->Sumw2();
      h_tauMu_tauE_bPt_OS=fs->make<TH1F>("h_tauMu_tauE_bPt_OS", "Mu||E bPt (OS)",20,0,120);
      h_tauMu_tauE_bPt_OS->Sumw2();
      h_tauMu_tauE_nE_OS=fs->make<TH1F>("h_tauMu_tauE_nE_OS", "Mu||E nE (OS)",5,0,5);
      h_tauMu_tauE_nE_OS->Sumw2();
      h_tauMu_tauE_nTau_OS=fs->make<TH1F>("h_tauMu_tauE_nTau_OS", "Mu||E nTau (OS)",5,0,5);
      h_tauMu_tauE_nTau_OS->Sumw2();
      h_tauMu_tauE_nB_OS=fs->make<TH1F>("h_tauMu_tauE_nB_OS", "Mu||E nB (OS)",5,0,5);
      h_tauMu_tauE_nB_OS->Sumw2();
      h_tauMu_tauE_nMu_OS=fs->make<TH1F>("h_tauMu_tauE_nMu_OS", "Mu||E nMu (OS)",5,0,5);
      h_tauMu_tauE_nMu_OS->Sumw2();
 
      h_tauMu_tauE_probb_OS=fs->make<TH1F>("h_tauMu_tauE_probb_OS", "Mu||E pfDeepFlavourJetTags:probb OS", 20,0,1);
      h_tauMu_tauE_probb_OS->Sumw2();
      h_tauMu_tauE_probbb_OS=fs->make<TH1F>("h_tauMu_tauE_probbb_OS", "Mu||E pfDeepFlavourJetTags:probbb OS", 20,0,1);
      h_tauMu_tauE_probbb_OS->Sumw2();
      h_tauMu_tauE_problepb_OS=fs->make<TH1F>("h_tauMu_tauE_problepb_OS", " Mu||E pfDeepFlavourJetTags:problepb OS", 20,0,1);
      h_tauMu_tauE_problepb_OS->Sumw2();
      h_tauMu_tauE_sum_OS=fs->make<TH1F>("h_tauMu_tauE_sum_OS", "Mu||E Total Btag OS", 20,0,1);
      h_tauMu_tauE_sum_OS->Sumw2();

      h_tauMu_tauE_Mvis_SS=fs->make<TH1F>("h_tauMu_tauE_Mvis_SS", "Mu||E Mvis SS",20,3,13);
      h_tauMu_tauE_Mvis_SS->Sumw2();
      h_tauMu_tauE_bMass_SS=fs->make<TH1F>("h_tauMu_tauE_bMass_SS", "Mu||E bMass SS",45,0,20);
      h_tauMu_tauE_bMass_SS->Sumw2();
      h_tauMu_tauE_bDRl_SS=fs->make<TH1F>("h_tauMu_tauE_bDRl_SS", "Mu||E bDRl SS",45,0,5);
      h_tauMu_tauE_bDRl_SS->Sumw2();
      h_tauMu_tauE_lepDR_SS=fs->make<TH1F>("h_tauMu_tauE_lepDR_SS", "Mu||E lepDR SS",20,0,1.2);
      h_tauMu_tauE_lepDR_SS->Sumw2();
      h_tauMu_tauE_muIso_SS=fs->make<TH1F>("h_tauMu_tauE_muIso_SS", "Mu||E muIso SS",20,0,1);
      h_tauMu_tauE_muIso_SS->Sumw2();
      h_tauMu_tauE_MET_SS=fs->make<TH1F>("h_tauMu_tauE_MET_SS", "Mu||E MET SS",40,0,250);
      h_tauMu_tauE_MET_SS->Sumw2();
      h_tauMu_tauE_DiTauPt_SS=fs->make<TH1F>("h_tauMu_tauE_DiTauPt_SS", "Mu||E DiTauPt SS",50,0,200);
      h_tauMu_tauE_DiTauPt_SS->Sumw2();
      h_tauMu_tauE_nMatchedEle_SS=fs->make<TH2F>("h_tauMu_tauE_nMatchedEle_SS", "Mu||E n Matched Electron(bx,ty) SS",4,0,4,4,0,4);
      h_tauMu_tauE_nMatchedEle_SS->Sumw2(); 
      h_tauMu_tauE_nMatchedMu_SS=fs->make<TH2F>("h_tauMu_tauE_nMatchedMu_SS", "Mu||E n Matched Muons(bx,ty) SS",4,0,4,4,0,4);
      h_tauMu_tauE_nMatchedMu_SS->Sumw2(); 
      h_tauMu_tauE_Trigger_SS=fs->make<TH1F>("h_tauMu_tauE_Trigger_SS", "Mu||E Triggers (SS)",18,0,18);
      h_tauMu_tauE_Trigger_SS->SetCanExtend(TH1::kAllAxes);
      h_tauMu_tauE_Trigger_SS->Sumw2();
      h_tauMu_tauE_BPH_Trigger_SS=fs->make<TH1F>("h_tauMu_tauE_BPH_Trigger_SS", "Mu||E BPH Triggers (SS)",25,0,25);
      h_tauMu_tauE_BPH_Trigger_SS->SetCanExtend(TH1::kAllAxes);
      h_tauMu_tauE_BPH_Trigger_SS->Sumw2();
      h_tauMu_tauE_bPt_SS=fs->make<TH1F>("h_tauMu_tauE_bPt_SS", "Mu||E bPt (SS)",20,0,120);
      h_tauMu_tauE_bPt_SS->Sumw2();
      h_tauMu_tauE_nE_SS=fs->make<TH1F>("h_tauMu_tauE_nE_SS", "Mu||E nE (SS)",5,0,5);
      h_tauMu_tauE_nE_SS->Sumw2();
      h_tauMu_tauE_nTau_SS=fs->make<TH1F>("h_tauMu_tauE_nTau_SS", "Mu||E nTau (SS)",5,0,5);
      h_tauMu_tauE_nTau_SS->Sumw2();
      h_tauMu_tauE_nB_SS=fs->make<TH1F>("h_tauMu_tauE_nB_SS", "Mu||E nB (SS)",5,0,5);
      h_tauMu_tauE_nB_SS->Sumw2();
      h_tauMu_tauE_nMu_SS=fs->make<TH1F>("h_tauMu_tauE_nMu_SS", "Mu||E nMu (SS)",5,0,5);
      h_tauMu_tauE_nMu_SS->Sumw2();
 
      h_tauMu_tauE_probb_SS=fs->make<TH1F>("h_tauMu_tauE_probb_SS", "Mu||E pfDeepFlavourJetTags:probb SS", 20,0,1);
      h_tauMu_tauE_probb_SS->Sumw2();
      h_tauMu_tauE_probbb_SS=fs->make<TH1F>("h_tauMu_tauE_probbb_SS", "Mu||E pfDeepFlavourJetTags:probbb SS", 20,0,1);
      h_tauMu_tauE_probbb_SS->Sumw2();
      h_tauMu_tauE_problepb_SS=fs->make<TH1F>("h_tauMu_tauE_problepb_SS", " Mu||E pfDeepFlavourJetTags:problepb SS", 20,0,1);
      h_tauMu_tauE_problepb_SS->Sumw2();
      h_tauMu_tauE_sum_SS=fs->make<TH1F>("h_tauMu_tauE_sum_SS", "Mu||E Total Btag SS", 20,0,1);
      h_tauMu_tauE_sum_SS->Sumw2();
//Mvis
      h_tauMu_tauMu_Mvis_OS=fs->make<TH1F>("h_tauMu_tauMu_Mvis_OS", "Mu||Mu Mvis OS",20,3,13);
      h_tauMu_tauMu_Mvis_OS->Sumw2();
      h_tauMu_tauMu_bMass_OS=fs->make<TH1F>("h_tauMu_tauMu_bMass_OS", "Mu||Mu bMass OS",45,0,20);
      h_tauMu_tauMu_bMass_OS->Sumw2();
      h_tauMu_tauMu_bDRl_OS=fs->make<TH1F>("h_tauMu_tauMu_bDRl_OS", "Mu||Mu bDRmu OS",45,0,5);
      h_tauMu_tauMu_bDRl_OS->Sumw2();
      h_tauMu_tauMu_DiTauPt_OS=fs->make<TH1F>("h_tauMu_tauMu_DiTauPt_OS", "Mu||Mu DiTauPt OS",50,0,200);
      h_tauMu_tauMu_DiTauPt_OS->Sumw2();
      h_tauMu_tauMu_lepDR_OS=fs->make<TH1F>("h_tauMu_tauMu_lepDR_OS", "Mu||Mu lepDR OS",20,0,1.2);
      h_tauMu_tauMu_lepDR_OS->Sumw2();
      h_tauMu_tauMu_muIso_OS=fs->make<TH1F>("h_tauMu_tauMu_muIso_OS", "Mu||Mu muIso OS",20,0,1);
      h_tauMu_tauMu_muIso_OS->Sumw2();
      h_tauMu_tauMu_MET_OS=fs->make<TH1F>("h_tauMu_tauMu_MET_OS", "Mu||Mu MET OS",40,0,250);
      h_tauMu_tauMu_MET_OS->Sumw2();
      h_tauMu_tauMu_nMatchedMu_OS=fs->make<TH2F>("h_tauMu_tauMu_nMatchedMu_OS", "Mu||Mu n Matched Muons(bx,ty) OS",4,0,4,4,0,4);
      h_tauMu_tauMu_nMatchedMu_OS->Sumw2();
      h_tauMu_tauMu_Trigger_OS=fs->make<TH1F>("h_tauMu_tauMu_Trigger_OS", "Mu||Mu Triggers (OS)",18,0,18);
      h_tauMu_tauMu_Trigger_OS->SetCanExtend(TH1::kAllAxes);
      h_tauMu_tauMu_Trigger_OS->Sumw2();
      h_tauMu_tauMu_BPH_Trigger_OS=fs->make<TH1F>("h_tauMu_tauMu_BPH_Trigger_OS", "Mu||Mu BPH Triggers (OS)",25,0,25);
      h_tauMu_tauMu_BPH_Trigger_OS->SetCanExtend(TH1::kAllAxes);
      h_tauMu_tauMu_BPH_Trigger_OS->Sumw2();
      h_tauMu_tauMu_bPt_OS=fs->make<TH1F>("h_tauMu_tauMu_bPt_OS", "Mu||Mu bPt (SS)",20,0,120);
      h_tauMu_tauMu_bPt_OS->Sumw2();
      h_tauMu_tauMu_nE_OS=fs->make<TH1F>("h_tauMu_tauMu_nE_OS", "Mu||Mu nE (OS)",5,0,5);
      h_tauMu_tauMu_nE_OS->Sumw2();
      h_tauMu_tauMu_nTau_OS=fs->make<TH1F>("h_tauMu_tauMu_nTau_OS", "Mu||Mu nTau (OS)",5,0,5);
      h_tauMu_tauMu_nTau_OS->Sumw2();
      h_tauMu_tauMu_nB_OS=fs->make<TH1F>("h_tauMu_tauMu_nB_OS", "Mu||Mu nB (OS)",5,0,5);
      h_tauMu_tauMu_nB_OS->Sumw2();
      h_tauMu_tauMu_nMu_OS=fs->make<TH1F>("h_tauMu_tauMu_nMu_OS", "Mu||Mu nMu (OS)",5,0,5);
      h_tauMu_tauMu_nMu_OS->Sumw2();

      h_tauMu_tauMu_probb_OS=fs->make<TH1F>("h_tauMu_tauMu_probb_OS", "Mu||Mu pfDeepFlavourJetTags:probb OS", 20,0,1);
      h_tauMu_tauMu_probb_OS->Sumw2();
      h_tauMu_tauMu_probbb_OS=fs->make<TH1F>("h_tauMu_tauMu_probbb_OS", "Mu||Mu pfDeepFlavourJetTags:probbb OS", 20,0,1);
      h_tauMu_tauMu_probbb_OS->Sumw2();
      h_tauMu_tauMu_problepb_OS=fs->make<TH1F>("h_tauMu_tauMu_problepb_OS", " Mu||Mu pfDeepFlavourJetTags:problepb OS", 20,0,1);
      h_tauMu_tauMu_problepb_OS->Sumw2();
      h_tauMu_tauMu_sum_OS=fs->make<TH1F>("h_tauMu_tauMu_sum_OS", "Mu||Mu Total Btag OS", 20,0,1);
      h_tauMu_tauMu_sum_OS->Sumw2(); 

      h_tauMu_tauMu_Mvis_SS=fs->make<TH1F>("h_tauMu_tauMu_Mvis_SS", "Mu||Mu Mvis SS",20,3,13);
      h_tauMu_tauMu_Mvis_SS->Sumw2();
      h_tauMu_tauMu_bMass_SS=fs->make<TH1F>("h_tauMu_tauMu_bMass_SS", "Mu||Mu bMass SS",45,0,20);
      h_tauMu_tauMu_bMass_SS->Sumw2();
      h_tauMu_tauMu_bDRl_SS=fs->make<TH1F>("h_tauMu_tauMu_bDRl_SS", "Mu||Mu bDRl SS",45,0,5);
      h_tauMu_tauMu_bDRl_SS->Sumw2();
      h_tauMu_tauMu_DiTauPt_SS=fs->make<TH1F>("h_tauMu_tauMu_DiTauPt_SS", "Mu||Mu DiTauPt SS",50,0,200);
      h_tauMu_tauMu_DiTauPt_SS->Sumw2();
      h_tauMu_tauMu_lepDR_SS=fs->make<TH1F>("h_tauMu_tauMu_lepDR_SS", "Mu||Mu lepDR SS",20,0,1.2);
      h_tauMu_tauMu_lepDR_SS->Sumw2();
      h_tauMu_tauMu_muIso_SS=fs->make<TH1F>("h_tauMu_tauMu_muIso_SS", "Mu||Mu muIso SS",20,0,1);
      h_tauMu_tauMu_muIso_SS->Sumw2();
      h_tauMu_tauMu_MET_SS=fs->make<TH1F>("h_tauMu_tauMu_MET_SS", "Mu||Mu MET SS",40,0,250);
      h_tauMu_tauMu_MET_SS->Sumw2();
      h_tauMu_tauMu_nMatchedMu_SS=fs->make<TH2F>("h_tauMu_tauMu_nMatchedMu_SS", "Mu||Mu n Matched Muons(bx,ty) SS",4,0,4,4,0,4);
      h_tauMu_tauMu_nMatchedMu_SS->Sumw2();
      h_tauMu_tauMu_Trigger_SS=fs->make<TH1F>("h_tauMu_tauMu_Trigger_SS", "Mu||Mu Triggers (SS)",18,0,18);
      h_tauMu_tauMu_Trigger_SS->SetCanExtend(TH1::kAllAxes);
      h_tauMu_tauMu_Trigger_SS->Sumw2();
      h_tauMu_tauMu_BPH_Trigger_SS=fs->make<TH1F>("h_tauMu_tauMu_BPH_Trigger_SS", "Mu||Mu BPH Triggers (SS)",25,0,25);
      h_tauMu_tauMu_BPH_Trigger_SS->SetCanExtend(TH1::kAllAxes);
      h_tauMu_tauMu_BPH_Trigger_SS->Sumw2();
      h_tauMu_tauMu_bPt_SS=fs->make<TH1F>("h_tauMu_tauMu_bPt_SS", "Mu||Mu bPt (SS)",20,0,120);
      h_tauMu_tauMu_bPt_SS->Sumw2();
      h_tauMu_tauMu_nE_SS=fs->make<TH1F>("h_tauMu_tauMu_nE_SS", "Mu||Mu nE (SS)",5,0,5);
      h_tauMu_tauMu_nE_SS->Sumw2();
      h_tauMu_tauMu_nTau_SS=fs->make<TH1F>("h_tauMu_tauMu_nTau_SS", "Mu||Mu nTau (SS)",5,0,5);
      h_tauMu_tauMu_nTau_SS->Sumw2();
      h_tauMu_tauMu_nB_SS=fs->make<TH1F>("h_tauMu_tauMu_nB_SS", "Mu||Mu nB (SS)",5,0,5);
      h_tauMu_tauMu_nB_SS->Sumw2();
      h_tauMu_tauMu_nMu_SS=fs->make<TH1F>("h_tauMu_tauMu_nMu_SS", "Mu||Mu nMu (SS)",5,0,5);
      h_tauMu_tauMu_nMu_SS->Sumw2();

      h_tauMu_tauMu_probb_SS=fs->make<TH1F>("h_tauMu_tauMu_probb_SS", "Mu||Mu pfDeepFlavourJetTags:probb SS", 20,0,1);
      h_tauMu_tauMu_probb_SS->Sumw2();
      h_tauMu_tauMu_probbb_SS=fs->make<TH1F>("h_tauMu_tauMu_probbb_SS", "Mu||Mu pfDeepFlavourJetTags:probbb SS", 20,0,1);
      h_tauMu_tauMu_probbb_SS->Sumw2();
      h_tauMu_tauMu_problepb_SS=fs->make<TH1F>("h_tauMu_tauMu_problepb_SS", " Mu||Mu pfDeepFlavourJetTags:problepb SS", 20,0,1);
      h_tauMu_tauMu_problepb_SS->Sumw2();
      h_tauMu_tauMu_sum_SS=fs->make<TH1F>("h_tauMu_tauMu_sum_SS", "Mu||Mu Total Btag SS", 20,0,1);
      h_tauMu_tauMu_sum_SS->Sumw2(); 
//tau Mu Tau Had
//tauHad tauE
      h_tauHad_tauE_Mvis_OS=fs->make<TH1F>("h_tauHad_tauE_Mvis_OS", "tauHad||E Mvis (OS)",20,3,13);
      h_tauHad_tauE_Mvis_OS->Sumw2();
      h_tauHad_tauE_bMass_OS=fs->make<TH1F>("h_tauHad_tauE_bMass_OS", "tauHad||E bMass (OS)",40,0,20);
      h_tauHad_tauE_bMass_OS->Sumw2();
      h_tauHad_tauE_bDRl_OS=fs->make<TH1F>("h_tauHad_tauE_bDRl_OS", "tauHad||E bDRl (OS)",40,0,5);
      h_tauHad_tauE_bDRl_OS->Sumw2();
      h_tauHad_tauE_lepDR_OS=fs->make<TH1F>("h_tauHad_tauE_lepDR_OS", "tauHad||E lepDR (OS)",20,0,1.2);
      h_tauHad_tauE_lepDR_OS->Sumw2();
      h_tauHad_tauE_DiTauPt_OS=fs->make<TH1F>("h_tauHad_tauE_DiTauPt_OS", "tauHad||E DiTauPt (OS)",50,0,200);
      h_tauHad_tauE_DiTauPt_OS->Sumw2();
      h_tauHad_tauE_MET_OS=fs->make<TH1F>("h_tauHad_tauE_MET_OS", "tauHad||E MET (OS)",40,0,250);
      h_tauHad_tauE_MET_OS->Sumw2();
      h_tauHad_tauE_Trigger_OS=fs->make<TH1F>("h_tauHad_tauE_Trigger_OS", "Had||E Triggers (OS)",18,0,18);
      h_tauHad_tauE_Trigger_OS->SetCanExtend(TH1::kAllAxes);
      h_tauHad_tauE_Trigger_OS->Sumw2();
      h_tauHad_tauE_BPH_Trigger_OS=fs->make<TH1F>("h_tauHad_tauE_BPH_Trigger_OS", "Had||E BPH Triggers (OS)",25,0,25);
      h_tauHad_tauE_BPH_Trigger_OS->SetCanExtend(TH1::kAllAxes);
      h_tauHad_tauE_BPH_Trigger_OS->Sumw2();
      h_tauHad_tauE_bPt_OS=fs->make<TH1F>("h_tauHad_tauE_bPt_OS", "Had||E bPt (OS)",20,0,120);
      h_tauHad_tauE_bPt_OS->Sumw2();
      h_tauHad_tauE_nE_OS=fs->make<TH1F>("h_tauHad_tauE_nE_OS", "Had||E nE (OS)",5,0,5);
      h_tauHad_tauE_nE_OS->Sumw2();
      h_tauHad_tauE_nTau_OS=fs->make<TH1F>("h_tauHad_tauE_nTau_OS", "Had||E nTau (OS)",5,0,5);
      h_tauHad_tauE_nTau_OS->Sumw2();
      h_tauHad_tauE_nB_OS=fs->make<TH1F>("h_tauHad_tauE_nB_OS", "Had||E nB (OS)",5,0,5);
      h_tauHad_tauE_nB_OS->Sumw2();
      h_tauHad_tauE_nMu_OS=fs->make<TH1F>("h_tauHad_tauE_nMu_OS", "Had||E nMu (OS)",5,0,5);
      h_tauHad_tauE_nMu_OS->Sumw2();

      h_tauHad_tauE_nMatchedEle_OS=fs->make<TH2F>("h_tauHad_tauE_nMatchedEle_OS", "tauHad||E n Matched Electron(bx,ty) (OS)",4,0,4,4,0,4);
      h_tauHad_tauE_nMatchedEle_OS->Sumw2(); 


      h_tauHad_tauE_probb_OS=fs->make<TH1F>("h_tauHad_tauE_probb_OS", "tauHad||E pfDeepFlavourJetTags:probb (OS)", 20,0,1);
      h_tauHad_tauE_probb_OS->Sumw2();
      h_tauHad_tauE_probbb_OS=fs->make<TH1F>("h_tauHad_tauE_probbb_OS", "tauHad||E pfDeepFlavourJetTags:probbb (OS)", 20,0,1);
      h_tauHad_tauE_probbb_OS->Sumw2();
      h_tauHad_tauE_problepb_OS=fs->make<TH1F>("h_tauHad_tauE_problepb_OS", " tauHad||E pfDeepFlavourJetTags:problepb (OS)", 20,0,1);
      h_tauHad_tauE_problepb_OS->Sumw2();
      h_tauHad_tauE_sum_OS=fs->make<TH1F>("h_tauHad_tauE_sum_OS", "tauHad||E Total Btag (OS)", 20,0,1);
      h_tauHad_tauE_sum_OS->Sumw2();       
      
      
      h_tauHad_tauE_Mvis_SS=fs->make<TH1F>("h_tauHad_tauE_Mvis_SS", "tauHad||E Mvis (SS)",45,0,20);
      h_tauHad_tauE_Mvis_SS->Sumw2();
      h_tauHad_tauE_bMass_SS=fs->make<TH1F>("h_tauHad_tauE_bMass_SS", "tauHad||E bMass (SS)",40,0,20);
      h_tauHad_tauE_bMass_SS->Sumw2();
      h_tauHad_tauE_bDRl_SS=fs->make<TH1F>("h_tauHad_tauE_bDRl_SS", "tauHad||E bDRl (SS)",40,0,5);
      h_tauHad_tauE_bDRl_SS->Sumw2();
      h_tauHad_tauE_lepDR_SS=fs->make<TH1F>("h_tauHad_tauE_lepDR_SS", "tauHad||E lepDR (SS)",20,0,1.2);
      h_tauHad_tauE_lepDR_SS->Sumw2();
      h_tauHad_tauE_DiTauPt_SS=fs->make<TH1F>("h_tauHad_tauE_DiTauPt_SS", "tauHad||E DiTauPt (SS)",50,0,200);
      h_tauHad_tauE_DiTauPt_SS->Sumw2();
      h_tauHad_tauE_MET_SS=fs->make<TH1F>("h_tauHad_tauE_MET_SS", "tauHad||E MET (SS)",40,0,250);
      h_tauHad_tauE_MET_SS->Sumw2();
      h_tauHad_tauE_Trigger_SS=fs->make<TH1F>("h_tauHad_tauE_Trigger_SS", "Had||E Triggers (SS)",18,0,18);
      h_tauHad_tauE_Trigger_SS->SetCanExtend(TH1::kAllAxes);
      h_tauHad_tauE_Trigger_SS->Sumw2();
      h_tauHad_tauE_BPH_Trigger_SS=fs->make<TH1F>("h_tauHad_tauE_BPH_Trigger_SS", "Had||E BPH Triggers (SS)",25,0,25);
      h_tauHad_tauE_BPH_Trigger_SS->SetCanExtend(TH1::kAllAxes);
      h_tauHad_tauE_BPH_Trigger_SS->Sumw2();
      h_tauHad_tauE_bPt_SS=fs->make<TH1F>("h_tauHad_tauE_bPt_SS", "Had||E bPt (SS)",20,0,120);
      h_tauHad_tauE_bPt_SS->Sumw2();
      
      h_tauHad_tauE_nE_SS=fs->make<TH1F>("h_tauHad_tauE_nE_SS", "Had||E nE (SS)",5,0,5);
      h_tauHad_tauE_nE_SS->Sumw2();
      h_tauHad_tauE_nTau_SS=fs->make<TH1F>("h_tauHad_tauE_nTau_SS", "Had||E nTau (SS)",5,0,5);
      h_tauHad_tauE_nTau_SS->Sumw2();
      h_tauHad_tauE_nB_SS=fs->make<TH1F>("h_tauHad_tauE_nB_SS", "Had||E nB (SS)",5,0,5);
      h_tauHad_tauE_nB_SS->Sumw2();
      h_tauHad_tauE_nMu_SS=fs->make<TH1F>("h_tauHad_tauE_nMu_SS", "Had||E nMu (SS)",5,0,5);
      h_tauHad_tauE_nMu_SS->Sumw2();

      h_tauHad_tauE_nMatchedEle_SS=fs->make<TH2F>("h_tauHad_tauE_nMatchedEle_SS", "tauHad||E n Matched Electron(bx,ty) (SS)",4,0,4,4,0,4);
      h_tauHad_tauE_nMatchedEle_SS->Sumw2(); 


      h_tauHad_tauE_probb_SS=fs->make<TH1F>("h_tauHad_tauE_probb_SS", "tauHad||E pfDeepFlavourJetTags:probb (SS)", 20,0,1);
      h_tauHad_tauE_probb_SS->Sumw2();
      h_tauHad_tauE_probbb_SS=fs->make<TH1F>("h_tauHad_tauE_probbb_SS", "tauHad||E pfDeepFlavourJetTags:probbb (SS)", 20,0,1);
      h_tauHad_tauE_probbb_SS->Sumw2();
      h_tauHad_tauE_problepb_SS=fs->make<TH1F>("h_tauHad_tauE_problepb_SS", " tauHad||E pfDeepFlavourJetTags:problepb (SS)", 20,0,1);
      h_tauHad_tauE_problepb_SS->Sumw2();
      h_tauHad_tauE_sum_SS=fs->make<TH1F>("h_tauHad_tauE_sum_SS", "tauHad||E Total Btag (SS)", 20,0,1);
      h_tauHad_tauE_sum_SS->Sumw2();   

//tauHad tauMu
      h_tauHad_tauMu_Mvis_OS=fs->make<TH1F>("h_tauHad_tauMu_Mvis_OS", "tauHad||Mu Mvis (OS)",20,3,13);
      h_tauHad_tauMu_Mvis_OS->Sumw2();
      h_tauHad_tauMu_bMass_OS=fs->make<TH1F>("h_tauHad_tauMu_bMass_OS", "tauHad||Mu bMass (OS)",40,0,20);
      h_tauHad_tauMu_bMass_OS->Sumw2();
      h_tauHad_tauMu_bDRl_OS=fs->make<TH1F>("h_tauHad_tauMu_bDRl_OS", "tauHad||Mu bDRl (OS)",40,0,5);
      h_tauHad_tauMu_bDRl_OS->Sumw2();
      h_tauHad_tauMu_lepDR_OS=fs->make<TH1F>("h_tauHad_tauMu_lepDR_OS", "tauHad||Mu lepDR (OS)",20,0,1.2);
      h_tauHad_tauMu_lepDR_OS->Sumw2();
      h_tauHad_tauMu_DiTauPt_OS=fs->make<TH1F>("h_tauHad_tauMu_DiTauPt_OS", "tauHad||Mu DiTauPt (OS)",50,0,200);
      h_tauHad_tauMu_DiTauPt_OS->Sumw2();
      h_tauHad_tauMu_MET_OS=fs->make<TH1F>("h_tauHad_tauMu_MET_OS", "tauHad||Mu MET (OS)",40,0,250);
      h_tauHad_tauMu_MET_OS->Sumw2();
      h_tauHad_tauMu_Trigger_OS=fs->make<TH1F>("h_tauHad_tauMu_Trigger_OS", "Had||Mu Triggers (OS)",18,0,18);
      h_tauHad_tauMu_Trigger_OS->SetCanExtend(TH1::kAllAxes);
      h_tauHad_tauMu_Trigger_OS->Sumw2();
      h_tauHad_tauMu_BPH_Trigger_OS=fs->make<TH1F>("h_tauHad_tauMu_BPH_Trigger_OS", "Had||Mu BPH Triggers (OS)",25,0,25);
      h_tauHad_tauMu_BPH_Trigger_OS->SetCanExtend(TH1::kAllAxes);
      h_tauHad_tauMu_BPH_Trigger_OS->Sumw2();
      h_tauHad_tauMu_bPt_OS=fs->make<TH1F>("h_tauHad_tauMu_bPt_OS", "Had||Mu bPt (OS)",20,0,120);
      h_tauHad_tauMu_bPt_OS->Sumw2();
      h_tauHad_tauMu_nE_OS=fs->make<TH1F>("h_tauHad_tauMu_nE_OS", "Had||Mu nE (OS)",5,0,5);
      h_tauHad_tauMu_nE_OS->Sumw2();
      h_tauHad_tauMu_nTau_OS=fs->make<TH1F>("h_tauHad_tauMu_nTau_OS", "Had||Mu nTau (OS)",5,0,5);
      h_tauHad_tauMu_nTau_OS->Sumw2();
      h_tauHad_tauMu_nB_OS=fs->make<TH1F>("h_tauHad_tauMu_nB_OS", "Had||Mu nB (OS)",5,0,5);
      h_tauHad_tauMu_nB_OS->Sumw2();
      h_tauHad_tauMu_nMu_OS=fs->make<TH1F>("h_tauHad_tauMu_nMu_OS", "Had||Mu nMu (OS)",5,0,5);
      h_tauHad_tauMu_nMu_OS->Sumw2();



      h_tauHad_tauMu_nMatchedMu_OS=fs->make<TH2F>("h_tauHad_tauMu_nMatchedMu_OS", "tauHad||Mu n Matched Electron(bx,ty) (OS)",4,0,4,4,0,4);
      h_tauHad_tauMu_nMatchedMu_OS->Sumw2(); 
      h_tauHad_tauMu_probb_OS=fs->make<TH1F>("h_tauHad_tauMu_probb_OS", "tauHad||Mu pfDeepFlavourJetTags:probb (OS)", 20,0,1);
      h_tauHad_tauMu_probb_OS->Sumw2();
      h_tauHad_tauMu_probbb_OS=fs->make<TH1F>("h_tauHad_tauMu_probbb_OS", "tauHad||Mu pfDeepFlavourJetTags:probbb (OS)", 20,0,1);
      h_tauHad_tauMu_probbb_OS->Sumw2();
      h_tauHad_tauMu_problepb_OS=fs->make<TH1F>("h_tauHad_tauMu_problepb_OS", " tauHad||Mu pfDeepFlavourJetTags:problepb (OS)", 20,0,1);
      h_tauHad_tauMu_problepb_OS->Sumw2();
      h_tauHad_tauMu_sum_OS=fs->make<TH1F>("h_tauHad_tauMu_sum_OS", "tauHad||Mu Total Btag (OS)", 20,0,1);
      h_tauHad_tauMu_sum_OS->Sumw2(); 


      h_tauHad_tauMu_Mvis_SS=fs->make<TH1F>("h_tauHad_tauMu_Mvis_SS", "tauHad||Mu Mvis (SS)",20,3,13);
      h_tauHad_tauMu_Mvis_SS->Sumw2();
      h_tauHad_tauMu_bMass_SS=fs->make<TH1F>("h_tauHad_tauMu_bMass_SS", "tauHad||Mu bMass (SS)",40,0,20);
      h_tauHad_tauMu_bMass_SS->Sumw2();
      h_tauHad_tauMu_bDRl_SS=fs->make<TH1F>("h_tauHad_tauMu_bDRl_SS", "tauHad||Mu bDRl (SS)",40,0,5);
      h_tauHad_tauMu_bDRl_SS->Sumw2();
      h_tauHad_tauMu_lepDR_SS=fs->make<TH1F>("h_tauHad_tauMu_lepDR_SS", "tauHad||Mu lepDR (SS)",20,0,1.2);
      h_tauHad_tauMu_lepDR_SS->Sumw2();
      h_tauHad_tauMu_DiTauPt_SS=fs->make<TH1F>("h_tauHad_tauMu_DiTauPt_SS", "tauHad||Mu DiTauPt (SS)",50,0,200);
      h_tauHad_tauMu_DiTauPt_SS->Sumw2();
      h_tauHad_tauMu_MET_SS=fs->make<TH1F>("h_tauHad_tauMu_MET_SS", "tauHad||Mu MET (SS)",40,0,250);
      h_tauHad_tauMu_MET_SS->Sumw2();
      h_tauHad_tauMu_Trigger_SS=fs->make<TH1F>("h_tauHad_tauMu_Trigger_SS", "Had||Mu Triggers (SS)",18,0,18);
      h_tauHad_tauMu_Trigger_SS->SetCanExtend(TH1::kAllAxes);
      h_tauHad_tauMu_Trigger_SS->Sumw2();
      h_tauHad_tauMu_BPH_Trigger_SS=fs->make<TH1F>("h_tauHad_tauMu_BPH_Trigger_SS", "Had||Mu BPH Triggers (SS)",25,0,25);
      h_tauHad_tauMu_BPH_Trigger_SS->SetCanExtend(TH1::kAllAxes);
      h_tauHad_tauMu_BPH_Trigger_SS->Sumw2();
      h_tauHad_tauMu_bPt_SS=fs->make<TH1F>("h_tauHad_tauMu_bPt_SS", "Had||Mu bPt (SS)",20,0,120);
      h_tauHad_tauMu_bPt_SS->Sumw2();
      h_tauHad_tauMu_nE_SS=fs->make<TH1F>("h_tauHad_tauMu_nE_SS", "Had||Mu nE (SS)",5,0,5);
      h_tauHad_tauMu_nE_SS->Sumw2();
      h_tauHad_tauMu_nTau_SS=fs->make<TH1F>("h_tauHad_tauMu_nTau_SS", "Had||Mu nTau (SS)",5,0,5);
      h_tauHad_tauMu_nTau_SS->Sumw2();
      h_tauHad_tauMu_nB_SS=fs->make<TH1F>("h_tauHad_tauMu_nB_SS", "Had||Mu nB (SS)",5,0,5);
      h_tauHad_tauMu_nB_SS->Sumw2();
      h_tauHad_tauMu_nMu_SS=fs->make<TH1F>("h_tauHad_tauMu_nMu_SS", "Had||Mu nMu (SS)",5,0,5);
      h_tauHad_tauMu_nMu_SS->Sumw2();



      h_tauHad_tauMu_nMatchedMu_SS=fs->make<TH2F>("h_tauHad_tauMu_nMatchedMu_SS", "tauHad||Mu n Matched Electron(bx,ty) (SS)",4,0,4,4,0,4);
      h_tauHad_tauMu_nMatchedMu_SS->Sumw2(); 
      h_tauHad_tauMu_probb_SS=fs->make<TH1F>("h_tauHad_tauMu_probb_SS", "tauHad||Mu pfDeepFlavourJetTags:probb (SS)", 20,0,1);
      h_tauHad_tauMu_probb_SS->Sumw2();
      h_tauHad_tauMu_probbb_SS=fs->make<TH1F>("h_tauHad_tauMu_probbb_SS", "tauHad||Mu pfDeepFlavourJetTags:probbb (SS)", 20,0,1);
      h_tauHad_tauMu_probbb_SS->Sumw2();
      h_tauHad_tauMu_problepb_SS=fs->make<TH1F>("h_tauHad_tauMu_problepb_SS", " tauHad||Mu pfDeepFlavourJetTags:problepb (SS)", 20,0,1);
      h_tauHad_tauMu_problepb_SS->Sumw2();
      h_tauHad_tauMu_sum_SS=fs->make<TH1F>("h_tauHad_tauMu_sum_SS", "tauHad||Mu Total Btag (SS)", 20,0,1);
      h_tauHad_tauMu_sum_SS->Sumw2(); 
///GEN Tau Had
      h_gen_tauHad_tauMu_Mvis=fs->make<TH1F>("h_gen_tauHad_tauMu_Mvis", "Gen tauHad||Mu Mvis",20,3,13);
      h_gen_tauHad_tauMu_Mvis->Sumw2();
      h_gen_tauHad_tauMu_bMass=fs->make<TH1F>("h_gen_tauHad_tauMu_bMass", "Gen tauHad||Mu bMass",40,0,20);
      h_gen_tauHad_tauMu_bMass->Sumw2();
      h_gen_tauHad_tauMu_bDRl=fs->make<TH1F>("h_gen_tauHad_tauMu_bDRl", "Gen tauHad||Mu bDRl",40,0,5);
      h_gen_tauHad_tauMu_bDRl->Sumw2();
      h_gen_tauHad_tauMu_lepDR=fs->make<TH1F>("h_gen_tauHad_tauMu_lepDR", "Gen tauHad||Mu lepDR",20,0,1.2);
      h_gen_tauHad_tauMu_lepDR->Sumw2();
      h_gen_tauHad_tauMu_DiTauPt=fs->make<TH1F>("h_gen_tauHad_tauMu_DiTauPt", "Gen tauHad||Mu DiTauPt",50,0,200);
      h_gen_tauHad_tauMu_DiTauPt->Sumw2();
      h_gen_tauHad_tauMu_MET=fs->make<TH1F>("h_gen_tauHad_tauMu_MET", "Gen tauHad||Mu MET",40,0,250);
      h_gen_tauHad_tauMu_MET->Sumw2();

 

      h_gen_tauHad_tauMu_probb=fs->make<TH1F>("h_gen_tauHad_tauMu_probb", "Gen tauHad||Mu pfDeepFlavourJetTags:probb", 20,0,1);
      h_gen_tauHad_tauMu_probb->Sumw2();
      h_gen_tauHad_tauMu_probbb=fs->make<TH1F>("h_gen_tauHad_tauMu_probbb", "Gen tauHad||Mu pfDeepFlavourJetTags:probbb", 20,0,1);
      h_gen_tauHad_tauMu_probbb->Sumw2();
      h_gen_tauHad_tauMu_problepb=fs->make<TH1F>("h_gen_tauHad_tauMu_problepb", " Gen tauHad||Mu pfDeepFlavourJetTags:problepb", 20,0,1);
      h_gen_tauHad_tauMu_problepb->Sumw2();
      h_gen_tauHad_tauMu_sum=fs->make<TH1F>("h_gen_tauHad_tauMu_sum", "Gen tauHad||Mu Total Btag", 20,0,1);
      h_gen_tauHad_tauMu_sum->Sumw2(); 

      h_gen_tauHad_tauE_Mvis=fs->make<TH1F>("h_gen_tauHad_tauE_Mvis", "Gen tauHad||E Mvis",20,3,13);
      h_gen_tauHad_tauE_Mvis->Sumw2();
      h_gen_tauHad_tauE_bMass=fs->make<TH1F>("h_gen_tauHad_tauE_bMass", "Gen tauHad||E bMass",40,0,20);
      h_gen_tauHad_tauE_bMass->Sumw2();
      h_gen_tauHad_tauE_bDRl=fs->make<TH1F>("h_gen_tauHad_tauE_bDRl", "Gen tauHad||E bDRl",40,0,5);
      h_gen_tauHad_tauE_bDRl->Sumw2();
      h_gen_tauHad_tauE_lepDR=fs->make<TH1F>("h_gen_tauHad_tauE_lepDR", "Gen tauHad||E lepDR",20,0,1.2);
      h_gen_tauHad_tauE_lepDR->Sumw2();
      h_gen_tauHad_tauE_DiTauPt=fs->make<TH1F>("h_gen_tauHad_tauE_DiTauPt", "Gen tauHad||E DiTauPt",50,0,200);
      h_gen_tauHad_tauE_DiTauPt->Sumw2();
      h_gen_tauHad_tauE_MET=fs->make<TH1F>("h_gen_tauHad_tauE_MET", "Gen tauHad||E MET",40,0,250);
      h_gen_tauHad_tauE_MET->Sumw2();



      h_gen_tauHad_tauE_probb=fs->make<TH1F>("h_gen_tauHad_tauE_probb", "Gen tauHad||E pfDeepFlavourJetTags:probb", 20,0,1);
      h_gen_tauHad_tauE_probb->Sumw2();
      h_gen_tauHad_tauE_probbb=fs->make<TH1F>("h_gen_tauHad_tauE_probbb", "Gen tauHad||E pfDeepFlavourJetTags:probbb", 20,0,1);
      h_gen_tauHad_tauE_probbb->Sumw2();
      h_gen_tauHad_tauE_problepb=fs->make<TH1F>("h_gen_tauHad_tauE_problepb", " Gen tauHad||E pfDeepFlavourJetTags:problepb", 20,0,1);
      h_gen_tauHad_tauE_problepb->Sumw2();
      h_gen_tauHad_tauE_sum=fs->make<TH1F>("h_gen_tauHad_tauE_sum", "Gen tauHad||E Total Btag", 20,0,1);
      h_gen_tauHad_tauE_sum->Sumw2(); 

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
   //std::cout<<to_string(triggerResults.size());
   edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
   iEvent.getByToken(triggerPrescales_, triggerPrescales);

   const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerBits);  
   //for (unsigned int i = 0, n = trigNames.size(); i < n; ++i) {
      //std::cout << "Trigger " << trigNames.triggerName(i) <<
      //   ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
      //   ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
      //   << std::endl;
   //}
   edm::Handle<GenEventInfoProduct>  genInfoHandle;
   try{iEvent.getByToken(genInfoProduct_, genInfoHandle);}
   catch(...){;}
   double genWeight=genInfoHandle->weight();
   pat::MET met=(*mets.product())[0];
   float MET = met.pt();


   std::vector<std::pair<std::string, float>> NewTriggers;
   //NewTriggers.push_back(make_pair("HLT_Mu24_v3",59.96e3));
   NewTriggers.push_back(make_pair("HLT_IsoMu24_v13",59.96e3));
   NewTriggers.push_back(make_pair("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v5",59.96e3));

   NewTriggers.push_back(make_pair("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v15",59.96e3));
   NewTriggers.push_back(make_pair("HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_v3",59.96e3));
   NewTriggers.push_back(make_pair("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v14",2.98e3));
   NewTriggers.push_back(make_pair("HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_v2",0.21e3));
   NewTriggers.push_back(make_pair("HLT_DoublePFJets40_CaloBTagDeepCSV_p71_v2",0.0065e3));
   NewTriggers.push_back(make_pair("HLT_DoubleMu4_3_Bs_v14",59.96e3));
   NewTriggers.push_back(make_pair("HLT_Mu3_PFJet40_v16",0.0046e3));
   NewTriggers.push_back(make_pair("HLT_IsoMu20_v15",.27e3));
   NewTriggers.push_back(make_pair("HLT_Mu8_TrkIsoVVL_v12",.0086e3));
   NewTriggers.push_back(make_pair("HLT_Ele15_WPLoose_Gsf_v3",1));
   NewTriggers.push_back(make_pair("HLT_BTagMu_AK4DiJet40_Mu5_v13",0.13e3));
   NewTriggers.push_back(make_pair("HLT_BTagMu_AK8DiJet170_Mu5_noalgo_v9",16.64e3));
   NewTriggers.push_back(make_pair("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v8",0.0011e3));
   NewTriggers.push_back(make_pair("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v16",.0064e3));
   NewTriggers.push_back(make_pair("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v4",59.96e3));
   NewTriggers.push_back(make_pair("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v1",59.96e3));
   
   std::vector<std::pair<std::string, float>> BParkingTriggers;
   BParkingTriggers.push_back(make_pair("HLT_Mu7_IP4_part0_v2",1.40e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu7_IP4_part1_v2",1.40e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu7_IP4_part2_v2",1.40e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu7_IP4_part3_v2",1.40e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu7_IP4_part4_v2",1.40e3));


   BParkingTriggers.push_back(make_pair("HLT_Mu8_IP3_part0_v3",.29e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu8_IP3_part1_v3",.29e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu8_IP3_part2_v3",.29e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu8_IP3_part3_v3",.29e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu8_IP3_part4_v3",.29e3));

   BParkingTriggers.push_back(make_pair("HLT_Mu8_IP6_part0_v2",1.66e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu8_IP6_part1_v2",1.66e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu8_IP6_part2_v2",1.66e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu8_IP6_part3_v2",1.66e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu8_IP6_part4_v2",1.66e3));

   BParkingTriggers.push_back(make_pair("HLT_Mu9_IP6_part0_v3",6.51e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu9_IP6_part1_v3",6.51e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu9_IP6_part2_v3",6.51e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu9_IP6_part3_v3",6.51e3));
   BParkingTriggers.push_back(make_pair("HLT_Mu9_IP6_part4_v3",6.51e3));

   //////Gen Particle Sorting
   std::vector<const reco::GenParticle*> gentauElectrons;
   std::vector<const reco::GenParticle*> gentauMuons;
   std::vector<const reco::GenParticle*> genTaus;
   std::vector<const reco::GenParticle*> genBs;
   std::vector<const reco::GenParticle*> genS1Electrons_NonTau;
   std::vector<const reco::GenParticle*> genS1Muons_NonTau;
   std::vector<const reco::GenParticle*> genNuTau; //////These are just the tau neutrinos We know there is a maximum of 2 per event either one anti, one normal, both or none.
   std::vector<const reco::GenParticle*> genNuTauBar; 
   for(const auto& gen : *genParticles){

      if(abs(gen.pdgId())==11 &&
         gen.isDirectHardProcessTauDecayProductFinalState()) {
         gentauElectrons.push_back(&gen);

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
         gentauMuons.push_back(&gen);
      }
      if (abs(gen.pdgId())==15 && gen.isHardProcess()){
         genTaus.push_back(&gen);
      }

      if (gen.pdgId()==16 && gen.isDirectHardProcessTauDecayProductFinalState()){
         genNuTau.push_back(&gen);

      }
      if (gen.pdgId()==-16 && gen.isDirectHardProcessTauDecayProductFinalState()){
         genNuTauBar.push_back(&gen);
      } 
      if (abs(gen.pdgId())==5&&gen.isHardProcess()){
         genBs.push_back(&gen);

      
      }

   }

   ////////////////////BJet ELECTRONS///////////////////////////////////////////
   std::vector<const reco::GenParticle*> genEfromB_S1;

   for(auto& ele : genS1Electrons_NonTau){
      int particleID = ele->pdgId();
      auto *mom = ele->mother();
      bool isfromb=false;
      bool isfroma=false;
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
      }
      if (isfromb&&isfroma){
         genEfromB_S1.push_back(ele);
      }

   }
   /////////////////////////////////B Muon Sorting/////////////////////////

   std::vector<const reco::GenParticle*> genMufromB_S1;

   for(auto& muon : genS1Muons_NonTau){
      int particleID = muon->pdgId();
      auto *mom = muon->mother();
      bool isfromb=false;
      bool isfroma=false;
      TLorentzVector mu;
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
         }
         if (abs(particleID)==36){isfroma=true;}
      }
      if (isfromb&&isfroma){
         genMufromB_S1.push_back(muon);
      }

   }
   //Determin which is the neutrino associated with 


   bool DiTauHad=false;
   std::vector<const reco::GenParticle*> genTauHad;
   std::vector<const reco::GenParticle*> genNuTauHad;
   for (auto& gen : genTaus){
      if ((gentauElectrons.size()==1&& gentauMuons.size()==0)){
         if ((*gentauElectrons[0]).charge()<0&&(*gen).charge()>0){

            genTauHad.push_back(gen);
            genNuTauHad.push_back(genNuTauBar[0]);
         }
         if ((*gentauElectrons[0]).charge()>0&&(*gen).charge()<0){
            genTauHad.push_back(gen);
            genNuTauHad.push_back(genNuTau[0]);
         }
      }
      else if((gentauElectrons.size()==0&& gentauMuons.size()==1)){
            if ((*gentauMuons[0]).charge()<0&&(*gen).charge()>0){
            genTauHad.push_back(gen);
            genNuTauHad.push_back(genNuTauBar[0]);
         }
         if ((*gentauMuons[0]).charge()>0&&(*gen).charge()<0){
            genTauHad.push_back(gen);
            genNuTauHad.push_back(genNuTauBar[0]);
         }

      }
      
      else{
         DiTauHad=true;
      }
   }
   

   std::vector<pat::Electron> selected_std_tau_electrons;
   std::vector<pat::Electron> selected_std_b_electrons;
   std::vector<pat::Electron> selected_std_electrons;
   std::vector<int> matchedtauEleGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level
   std::vector<int> matchedbEleGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level

   for(auto& ele : *std_electrons){
      h_std_electron_selection->Fill("all Ele",genWeight);
      if (ele.pt()>7 && abs(ele.eta())<2.5){
         h_std_electron_selection->Fill("Pt>7 ete<2.5",genWeight);
         if(ele.electronID("cutBasedElectronID-Fall17-94X-V2-loose")==1){
            h_std_electron_selection->Fill("Loose ID",genWeight);
            selected_std_electrons.push_back(ele);
            bool tauEmatched= false;
            float tauDr_min=.1;
            float tauMatchedGen_index=-1;
            for(unsigned int iGenE=0; iGenE<gentauElectrons.size(); iGenE++){
               if (std::count(matchedtauEleGenIndexs.begin(), matchedtauEleGenIndexs.end(), iGenE)==0){
                  reco::GenParticle gen = *gentauElectrons[iGenE];
                  TLorentzVector e;
                  TLorentzVector genE;
                  e.SetPtEtaPhiM(ele.pt(), ele.eta(), ele.phi(), ele.mass());
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
               h_std_electron_selection->Fill("tau Matched*",genWeight);
               selected_std_tau_electrons.push_back(ele);
               h_StdEle_taum_pt->Fill(ele.pt(),genWeight);
               matchedtauEleGenIndexs.push_back(tauMatchedGen_index);
            }

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
               h_std_electron_selection->Fill("b Matched**",genWeight);
               h_StdEle_bm_pt->Fill(ele.pt(),genWeight);
               selected_std_b_electrons.push_back(ele);
               matchedbEleGenIndexs.push_back(bMatchedGen_index);
           
            }
         }
      }
   }
   //std::cout<<"LowPt Electron selection\n";
   std::vector<pat::Jet> selected_bjets;
   std::vector<pat::Jet> selected_bjets_m;
   std::vector<pat::Jet> selected_bjets_f;

   std::vector<int> matchedBIDs;
   std::vector<bool> jetLocks(jetFlavorMatching->size(),false);
   for (auto& jet:*jets){
      h_bJet_selection->Fill("All Jets",genWeight);
      if (jet.pt()>10 && abs(jet.eta())<2.5 ){
         h_bJet_selection->Fill("pt>10, |eta|<2.5",genWeight);

         float NHF  = jet.neutralHadronEnergyFraction();
         float NEMF = jet.neutralEmEnergyFraction();
         float CHF  = jet.chargedHadronEnergyFraction();
         float MUF  = jet.muonEnergyFraction();
         float CEMF = jet.chargedEmEnergyFraction();
         float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
         float CHM = jet.chargedMultiplicity();
         if (CEMF>0.8){continue;}
         h_bJet_selection->Fill("CEMF<.8",genWeight);
         if (CHM<0){continue;}
         h_bJet_selection->Fill("CHM>0",genWeight);
         if (CHF <0){continue;}
         h_bJet_selection->Fill("CHF>0",genWeight);

         if (NumConst<1){continue;}
         h_bJet_selection->Fill("charge +neutral multi>1",genWeight);

         if (NEMF>0.9){continue;}
         h_bJet_selection->Fill("NEMF<.9",genWeight);
         if (MUF>0.8){continue;}
         h_bJet_selection->Fill("MUF<.8",genWeight);
         if (NHF>.9){continue;}

         h_bJet_selection->Fill("NHF<.9",genWeight);
         bool match=false;
         float dRMin=9999.;
         int matchedIdx=-1;

         h_deepFlavor_sum->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb")  + jet.bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
         h_deepFlavor_probb->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
         h_deepFlavor_probbvsprobbb->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probb"),jet.bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);

         h_deepFlavor_probbb->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
         h_deepFlavor_problepb->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
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
            jetLocks.at(matchedIdx) = true;
            h_deepFlavor_sum_m->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb")  + jet.bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
            h_deepFlavor_probb_m->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
            h_deepFlavor_probbb_m->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
            h_deepFlavor_problepb_m->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
            h_deepFlavor_probbvsprobbb_m->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probb"),jet.bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);

         }
         else{
            h_deepFlavor_probbvsprobbb_f->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probb"),jet.bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);

            h_deepFlavor_sum_f->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb")  + jet.bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
            h_deepFlavor_probb_f->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
            h_deepFlavor_probbb_f->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
            h_deepFlavor_problepb_f->Fill(jet.bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
         }
         if (jet.bDiscriminator("pfDeepFlavourJetTags:probb") + jet.bDiscriminator("pfDeepFlavourJetTags:probbb")  + jet.bDiscriminator("pfDeepFlavourJetTags:problepb")<.7476){continue;}

         h_bJet_selection->Fill("DeepJet Sum>.7476",genWeight);
         if (match==true){
            selected_bjets_m.push_back(jet);
         }
         else{
            selected_bjets_f.push_back(jet);
         }
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
         //h_bJet_selection->Fill(10.5,genWeight);
         selected_bjets.push_back(jet);

      }
   }
   std::vector<pat::Muon> selected_tau_gen_matched_Muons;
   std::vector<pat::Muon> selected_b_gen_matched_Muons;
   std::vector<pat::Muon> selected_Muons;
   std::vector<pat::Muon> selected_b_Muons;   
   std::vector<int> matchedtauMuGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level
   std::vector<int> matchedbMuGenIndexs;////We keep track of all matched Gen Electrons to make sure we dont match multiple candidates to the same gen level

   for (auto& muon: *muons){
      h_muon_selection->Fill("All Muons",genWeight);
      if (muon::isLooseMuon(muon)==false){continue;}////////////Loose Muon ID
      h_muon_selection->Fill("IsLooseMuon",genWeight);
      if((muon.pt())<3|| abs(muon.eta())>2.4){continue;}
      h_muon_selection->Fill("Pt>3,|eta|<2.4",genWeight);
 
      if(h2AA::muonIsolation(muon)>.25){continue;}
      h_muon_selection->Fill("Isolation<.25",genWeight);
      bool tauMumatched= false;
      float tauDr_min=.1;
      float tauMatchedGen_index=-1;
      bool isFromB=false;
      if (selected_bjets.size()>0){
         TLorentzVector mu;
         TLorentzVector b;
         mu.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), muon.mass());
         b.SetPtEtaPhiM(selected_bjets[0].pt(), selected_bjets[0].eta(), selected_bjets[0].phi(), selected_bjets[0].mass());
         float dr = b.DeltaR(mu);
         if (dr<.4){
            isFromB=true;
         }
      }
      if (isFromB){
         selected_b_Muons.push_back(muon);
         h_muon_selection->Fill("bDRMu<.4",genWeight);
         h_deepFlavor_Leptonic_sum->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
         h_deepFlavor_Leptonic_probb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
         h_deepFlavor_Leptonic_probbb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
         h_deepFlavor_Leptonic_problepb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
         h_deepFlavor_Leptonic_probbvsprobbb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb"),selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
      }
      else {
         selected_Muons.push_back(muon);
         h_muon_selection->Fill("bDRMu>.4",genWeight);
      }


      for (unsigned int iGen=0; iGen<gentauMuons.size();iGen++){
         if (std::count(matchedtauMuGenIndexs.begin(), matchedtauMuGenIndexs.end(), iGen)==0){

            reco::GenParticle gen = *gentauMuons[iGen];
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
      if (tauMumatched==true){
         h_muon_selection->Fill("tauMatched *",genWeight);
         matchedtauMuGenIndexs.push_back(tauMatchedGen_index);
         h_Mu_taum_pt->Fill(muon.pt(),genWeight);

         selected_tau_gen_matched_Muons.push_back(muon);
      }
      bool bMumatched= false;
      float bDr_min=.1;
      float bMatchedGen_index=-1;
      for (unsigned int iGen=0; iGen<genMufromB_S1.size();iGen++){
         if (std::count(matchedbMuGenIndexs.begin(), matchedbMuGenIndexs.end(), iGen)==0){
            reco::GenParticle gen = *genMufromB_S1[iGen];
            TLorentzVector mu;
            TLorentzVector genMu;
            mu.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), muon.mass());
            genMu.SetPtEtaPhiM(gen.pt(), gen.eta(), gen.phi(), gen.mass());
            float dr = genMu.DeltaR(mu);
            if(dr <bDr_min){
               bMumatched = true;
               bDr_min=dr;
               bMatchedGen_index=iGen;///Keep track of the gen particle that is closest to the reco candidate
            }
         }
      }
      
      if (bMumatched==true){
         h_muon_selection->Fill("b matched *",genWeight);
         matchedbMuGenIndexs.push_back(bMatchedGen_index);
         h_Mu_bm_pt->Fill(muon.pt(),genWeight);

         selected_b_gen_matched_Muons.push_back(muon);
      }

   }
   std::sort(selected_std_electrons.begin(), selected_std_electrons.end(), h2AA::sortByPt<pat::Electron>);
   std::sort(selected_std_tau_electrons.begin(), selected_std_tau_electrons.end(), h2AA::sortByPt<pat::Electron>);
   std::sort(selected_std_b_electrons.begin(), selected_std_b_electrons.end(), h2AA::sortByPt<pat::Electron>);
   std::sort(selected_Muons.begin(), selected_Muons.end(), h2AA::sortByPt<pat::Muon>);
   std::sort(selected_b_gen_matched_Muons.begin(), selected_b_gen_matched_Muons.end(), h2AA::sortByPt<pat::Muon>);
   std::sort(selected_tau_gen_matched_Muons.begin(), selected_tau_gen_matched_Muons.end(), h2AA::sortByPt<pat::Muon>);
   std::sort(genEfromB_S1.begin(), genEfromB_S1.end(), h2AA::sortGenByPt);
   std::sort(gentauElectrons.begin(), gentauElectrons.end(), h2AA::sortGenByPt);
   std::sort(gentauMuons.begin(), gentauMuons.end(), h2AA::sortGenByPt);
   //std::cout<<"Mu selection\n";
   //std::cout<<to_string(NewTriggers.size())<<"\n";
   for (unsigned int index = 0; index < NewTriggers.size(); ++index){
      //std::cout<<to_string(index)<<"\n";

      if (triggerBits->accept(trigNames.triggerIndex(NewTriggers[index].first))){
         ////std::cout<<NewTriggers[index]<<"\n";
         h_Triggers->Fill(index,NewTriggers[index].second*genWeight);
      }
      
   }
   h_nbJet->Fill(selected_bjets.size(), genWeight);
   h_nbJet_m->Fill(selected_bjets_m.size(), genWeight);
   h_nbJet_f->Fill(selected_bjets_f.size(), genWeight);
   h_nStdEle->Fill(selected_std_electrons.size(), genWeight);
   h_nStdEle_bm->Fill(selected_std_b_electrons.size(), genWeight);
   h_nStdEle_taum->Fill(selected_std_tau_electrons.size(), genWeight);
   h_nMu->Fill(selected_Muons.size(), genWeight);
   h_nMu_bm->Fill(selected_b_gen_matched_Muons.size(), genWeight);
   h_nMu_taum->Fill(selected_tau_gen_matched_Muons.size(), genWeight);

   //std::cout<<"jet selection\n";
   std::sort(selected_bjets.begin(), selected_bjets.end(), h2AA::sortByPt<pat::Jet>);
   std::sort(selected_bjets_m.begin(), selected_bjets_m.end(), h2AA::sortByPt<pat::Jet>);
   std::sort(selected_bjets_f.begin(), selected_bjets_f.end(), h2AA::sortByPt<pat::Jet>);
   //std::cout<<"bsorting\n";
   std::vector<pat::Tau> selected_taus;
   for (auto& tau:*taus){

      h_tau_selection->Fill("Total Candidates", genWeight);
      if (tau.pt()<20 || abs(tau.eta())>2.3){continue;}
      h_tau_selection->Fill("pt>10,eta<2.3", genWeight);

      if (tau.tauID("decayModeFindingNewDMs")==false){continue;}
      h_tau_selection->Fill("DM Finding", genWeight);

      if (tau.tauID("againstElectronLooseMVA6") ==false&& selected_std_electrons.size()==0){continue;}
      h_tau_selection->Fill("againstElectronLooseMVA6*", genWeight);
      if (tau.tauID("byLooseIsolationMVArun2v1PWoldDMwLT")==false){continue;}
      h_tau_selection->Fill("byLooseIsolationMVArun2v1PWoldDMwLT", genWeight);
      if (tau.tauID("againstMuonLoose3")==false&&selected_Muons.size()==0){continue;}
      h_tau_selection->Fill("againstMuonLoose3**", genWeight);

      TLorentzVector t;
      t.SetPtEtaPhiM(tau.pt(), tau.eta(), tau.phi(), tau.mass());
      bool muFake= false;
      for (auto& mu: selected_Muons){
         TLorentzVector m;
         m.SetPtEtaPhiM(mu.pt(), mu.eta(), mu.phi(), mu.mass());
         if (m.DeltaR(t)<.1){
            muFake=true;
            break;
         }
      }
      if (muFake){continue;}
      h_tau_selection->Fill("t.DR(all mu)>.1 ", genWeight);
      bool eleFake= false;
      for (auto& ele: selected_std_electrons){
         TLorentzVector e;
         e.SetPtEtaPhiM(ele.pt(), ele.eta(), ele.phi(), ele.mass());
         if (e.DeltaR(t)<.1){
            eleFake=true;
            break;
         }
      }
      if (eleFake){continue;}
      h_tau_selection->Fill("t.DR(all e)>.1", genWeight);
      if (selected_bjets.size()>0){
         TLorentzVector b;
         b.SetPtEtaPhiM(selected_bjets[0].pt(),selected_bjets[0].eta(),selected_bjets[0].phi(),selected_bjets[0].mass());
         if(b.DeltaR(t)<.6){continue;}
      }
      h_tau_selection->Fill("t.DR(b[0])>.6", genWeight);
      selected_taus.push_back(tau);

      if(genTauHad.size()==0||genNuTauHad.size()==0){continue;}
      if(DiTauHad==true){continue;}
      reco::GenParticle thad =*genTauHad[0];
      reco::GenParticle nt = *genNuTauHad[0];
      TLorentzVector genT;
      TLorentzVector genNu;
      genT.SetPtEtaPhiM(thad.pt(), thad.eta(), thad.phi(), thad.mass());
      genNu.SetPtEtaPhiM(nt.pt(), nt.eta(), nt.phi(), nt.mass());
      if(t.DeltaR(genT-genNu)>.4){continue;}
      h_tau_selection->Fill("Gen Matching>.4***", genWeight);

   }
   std::sort(selected_taus.begin(), selected_taus.end(), h2AA::sortByPt<pat::Tau>);
   if (selected_bjets.size()>0){
      
      h_deepFlavor_leading_sum->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
      h_deepFlavor_leading_probb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
      h_deepFlavor_leading_probbb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
      h_deepFlavor_leading_problepb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
      h_deepFlavor_leading_probbvsprobbb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb"),selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);

      if (selected_bjets.size()>1){
         h_deepFlavor_subLeading_sum->Fill(selected_bjets[1].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[1].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[1].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
         h_deepFlavor_subLeading_probb->Fill(selected_bjets[1].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
         h_deepFlavor_subLeading_probbb->Fill(selected_bjets[1].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
         h_deepFlavor_subLeading_problepb->Fill(selected_bjets[1].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
         h_deepFlavor_subLeading_probbvsprobbb->Fill(selected_bjets[1].bDiscriminator("pfDeepFlavourJetTags:probb"),selected_bjets[1].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
      }
   }
   //std::cout<<"jet selection\n";
   if (selected_bjets.size()==1){
      
      h_deepFlavor_Merged_sum->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
      h_deepFlavor_Merged_probb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
      h_deepFlavor_Merged_probbb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
      h_deepFlavor_Merged_problepb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
      h_deepFlavor_Merged_probbvsprobbb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb"),selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
   }
   if (selected_bjets_m.size()==1){
      
      h_deepFlavor_Merged_sum_m->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
      h_deepFlavor_Merged_probb_m->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
      h_deepFlavor_Merged_probbb_m->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
      h_deepFlavor_Merged_problepb_m->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
      h_deepFlavor_Merged_probbvsprobbb_m->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb"),selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
   }
   if (selected_bjets_m.size()>0){
      
      h_deepFlavor_leading_sum_m->Fill(selected_bjets_m[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets_m[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets_m[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
      h_deepFlavor_leading_probb_m->Fill(selected_bjets_m[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
      h_deepFlavor_leading_probbb_m->Fill(selected_bjets_m[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
      h_deepFlavor_leading_problepb_m->Fill(selected_bjets_m[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
      h_deepFlavor_leading_probbvsprobbb_m->Fill(selected_bjets_m[0].bDiscriminator("pfDeepFlavourJetTags:probb"),selected_bjets_m[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);

      if (selected_bjets_m.size()>1){
         h_deepFlavor_subLeading_sum_m->Fill(selected_bjets_m[1].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets_m[1].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets_m[1].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
         h_deepFlavor_subLeading_probb_m->Fill(selected_bjets_m[1].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
         h_deepFlavor_subLeading_probbb_m->Fill(selected_bjets_m[1].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
         h_deepFlavor_subLeading_problepb_m->Fill(selected_bjets_m[1].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
         h_deepFlavor_subLeading_probbvsprobbb_m->Fill(selected_bjets_m[1].bDiscriminator("pfDeepFlavourJetTags:probb"),selected_bjets_m[1].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);

      }
   }
   if (selected_bjets_f.size()>0){
      
      h_deepFlavor_leading_sum_f->Fill(selected_bjets_f[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets_f[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets_f[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
      h_deepFlavor_leading_probb_f->Fill(selected_bjets_f[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
      h_deepFlavor_leading_probbb_f->Fill(selected_bjets_f[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
      h_deepFlavor_leading_problepb_f->Fill(selected_bjets_f[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
      h_deepFlavor_leading_probbvsprobbb_f->Fill(selected_bjets_f[0].bDiscriminator("pfDeepFlavourJetTags:probb"),selected_bjets_f[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);

      if (selected_bjets_f.size()>1){
         h_deepFlavor_subLeading_sum_f->Fill(selected_bjets_f[1].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets_f[1].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets_f[1].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
         h_deepFlavor_subLeading_probb_f->Fill(selected_bjets_f[1].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
         h_deepFlavor_subLeading_probbb_f->Fill(selected_bjets_f[1].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
         h_deepFlavor_subLeading_problepb_f->Fill(selected_bjets_f[1].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
         h_deepFlavor_subLeading_probbvsprobbb_f->Fill(selected_bjets_f[1].bDiscriminator("pfDeepFlavourJetTags:probb"),selected_bjets_f[1].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);

      }
   }


   for (auto &muon:selected_Muons){
      TLorentzVector mu;
      
      mu.SetPtEtaPhiM(muon.pt(),muon.eta(),muon.phi(),muon.mass());
      for(auto &jet:selected_bjets){
         TLorentzVector b;
         b.SetPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),jet.mass());
         if (b.DeltaR(mu)<.4){

         }
      }
   }

   //std::cout<<"flavor plotting\n";

//////////////////////Event Selection////////////////////////
///////////////////
   std::sort(genMufromB_S1.begin(), genMufromB_S1.end(), h2AA::sortGenByPt);
   std::sort(genEfromB_S1.begin(), genEfromB_S1.end(), h2AA::sortGenByPt);
   std::sort(gentauElectrons.begin(), gentauElectrons.end(), h2AA::sortGenByPt);
   std::sort(gentauMuons.begin(), gentauMuons.end(), h2AA::sortGenByPt);
   //std::string EGLoose = "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v7";
   //std::string singleMuLoose="HLT_IsoMu27_v13";
   //std::string doubleMuonLoose="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v3";
   //std::string MuEGLoose="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v12";
   //std::string DoubleElectronLoose="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v17";
   std::string EGLoose = "HLT_Ele15_WPLoose_Gsf_v3";
   std::string singleMuLoose="HLT_IsoMu20_v15";
   std::string doubleMuonLoose="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v14";
   std::string MuEGLoose="HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v15";
   std::string DoubleElectronLoose="HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_v3";
   bool passTrig_EGLoose=triggerBits->accept(trigNames.triggerIndex(EGLoose));
   bool passTrig_singleMuLoose=triggerBits->accept(trigNames.triggerIndex(singleMuLoose));
   bool passTrig_doubleMuonLoose=triggerBits->accept(trigNames.triggerIndex(doubleMuonLoose));
   bool passTrig_DoubleElectronLoose=triggerBits->accept(trigNames.triggerIndex(DoubleElectronLoose));
   bool passTrig_MuEGLoose=triggerBits->accept(trigNames.triggerIndex(MuEGLoose));
   bool passTrig_Mu7IP4= (
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu7_IP4_part0_v2"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu7_IP4_part1_v2"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu7_IP4_part2_v2"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu7_IP4_part3_v2"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu7_IP4_part4_v2"))
   );
   bool passTrig_Mu8IP3= (
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu8_IP3_part0_v3"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu8_IP3_part1_v3"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu8_IP3_part2_v3"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu8_IP3_part3_v3"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu8_IP3_part4_v3"))
   );
   bool passTrig_Mu8IP6= (
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu8_IP6_part0_v2"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu8_IP6_part1_v2"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu8_IP6_part2_v2"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu8_IP6_part3_v2"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu8_IP6_part4_v2"))
   );
   bool passTrig_Mu9IP6= (
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu9_IP6_part0_v3"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu9_IP6_part1_v3"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu9_IP6_part2_v3"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu9_IP6_part3_v3"))||
      triggerBits->accept(trigNames.triggerIndex("HLT_Mu9_IP6_part4_v3"))
   );
   if (passTrig_doubleMuonLoose==false){
      //std::cout<<"double mu false\n";
   }
   else{
      //std::cout<<"double mu true\n";

   }
//tauE tauE
   h_nEvent_tauE_tauE->Fill("Total Events",genWeight);
   if (selected_std_electrons.size()>1){
      h_nEvent_tauE_tauE->Fill("nStdtauE>1",genWeight);
      if (selected_bjets.size()>0){
         h_nEvent_tauE_tauE->Fill("nbJets>0",genWeight);
         TLorentzVector e1;
         e1.SetPtEtaPhiM(selected_std_electrons[0].pt(),selected_std_electrons[0].eta(),selected_std_electrons[0].phi(),selected_std_electrons[0].mass());
         TLorentzVector e2;
         e2.SetPtEtaPhiM(selected_std_electrons[1].pt(),selected_std_electrons[1].eta(),selected_std_electrons[1].phi(),selected_std_electrons[1].mass());
         TLorentzVector b;
         if (selected_bjets.size()>1){
            TLorentzVector b1;
            b1.SetPtEtaPhiM(selected_bjets[0].pt(),selected_bjets[0].eta(),selected_bjets[0].phi(),selected_bjets[0].mass());
            TLorentzVector b2;
            b2.SetPtEtaPhiM(selected_bjets[1].pt(),selected_bjets[1].eta(),selected_bjets[1].phi(),selected_bjets[1].mass());
            b = b1+b2;
         }
         else{
            b.SetPtEtaPhiM(selected_bjets[0].pt(),selected_bjets[0].eta(),selected_bjets[0].phi(),selected_bjets[0].mass());
         }  
         float mvis = (e1+e2).M();
         float DiTauPt=(e2+e1).Pt();
         float bMass=b.M();
         float bDRl=b.DeltaR(e1);
         float lepDR=e1.DeltaR(e2);
         float bPt=b.Pt();
         if (bMass<18){
            h_nEvent_tauE_tauE->Fill("bMass<18",genWeight);
            if (mvis<12&&mvis>3){
               h_nEvent_tauE_tauE->Fill("3<mvis<12",genWeight);
               //if (passTrig_DoubleElectronLoose){
                  //h_nEvent_tauE_tauE->Fill("Double Electron", genWeight);
               if (lepDR<1){
                  h_nEvent_tauE_tauE->Fill("(lepDR<1)",genWeight);
                  if(selected_std_electrons[0].charge()!=selected_std_electrons[1].charge()){
                     
                     for (unsigned int index = 0; index < NewTriggers.size(); ++index){
                        if (triggerBits->accept(trigNames.triggerIndex(NewTriggers[index].first))){
                           h_tauE_tauE_Trigger_OS->Fill(index,NewTriggers[index].second*genWeight);
                        }
                        
                     }
                     for (unsigned int index = 0; index < BParkingTriggers.size(); ++index){
                        if (triggerBits->accept(trigNames.triggerIndex(BParkingTriggers[index].first))){
                           h_tauE_tauE_BPH_Trigger_OS->Fill(index,BParkingTriggers[index].second*genWeight);
                        }
                        
                     }
                     
                     if (passTrig_Mu7IP4){
                        h_tauE_tauE_BPH_Trigger_OS->Fill(21.,7e3*genWeight);
                     }
                     if (passTrig_Mu8IP3){
                        h_tauE_tauE_BPH_Trigger_OS->Fill(22.,1.45e3*genWeight);
                     }
                     if (passTrig_Mu8IP6){
                        h_tauE_tauE_BPH_Trigger_OS->Fill(23.,8.3e3*genWeight);
                     }
                     if (passTrig_Mu9IP6){
                        h_tauE_tauE_BPH_Trigger_OS->Fill(24.,33.7e3*genWeight);
                     }
                     h_nEvent_tauE_tauE->Fill("OS",genWeight);

                     h_tauE_tauE_bPt_OS->Fill(bPt, genWeight);

                     h_tauE_tauE_Mvis_OS->Fill(mvis, genWeight);
                     h_tauE_tauE_bMass_OS->Fill(bMass, genWeight);
                     h_tauE_tauE_DiTauPt_OS->Fill(DiTauPt, genWeight);
                     h_tauE_tauE_bDRl_OS->Fill(bDRl, genWeight);
                     h_tauE_tauE_lepDR_OS->Fill(lepDR, genWeight);
                     h_tauE_tauE_MET_OS->Fill(MET, genWeight);
                     h_tauE_tauE_nE_OS->Fill(selected_std_electrons.size(), genWeight);
                     h_tauE_tauE_nTau_OS->Fill(selected_taus.size(), genWeight);
                     h_tauE_tauE_nMu_OS->Fill(selected_Muons.size(), genWeight);
                     h_tauE_tauE_nB_OS->Fill(selected_bjets.size(), genWeight);
                     h_tauE_tauE_nMatchedEle_OS->Fill(selected_std_b_electrons.size(),selected_std_tau_electrons.size(),genWeight);
                     h_tauE_tauE_sum_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
                     h_tauE_tauE_probb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
                     h_tauE_tauE_probbb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
                     h_tauE_tauE_problepb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
                  }
                  else{
                     h_nEvent_tauE_tauE->Fill("SS",genWeight);
                     h_tauE_tauE_Mvis_SS->Fill(mvis, genWeight);
                     h_tauE_tauE_bMass_SS->Fill(bMass, genWeight);
                     h_tauE_tauE_DiTauPt_SS->Fill(DiTauPt, genWeight);
                     h_tauE_tauE_bDRl_SS->Fill(bDRl, genWeight);
                     h_tauE_tauE_lepDR_SS->Fill(lepDR, genWeight);
                     h_tauE_tauE_MET_SS->Fill(MET, genWeight);
                     h_tauE_tauE_bPt_SS->Fill(bPt, genWeight);
                     h_tauE_tauE_nE_SS->Fill(selected_std_electrons.size(), genWeight);
                     h_tauE_tauE_nTau_SS->Fill(selected_taus.size(), genWeight);
                     h_tauE_tauE_nMu_SS->Fill(selected_Muons.size(), genWeight);
                     h_tauE_tauE_nB_SS->Fill(selected_bjets.size(), genWeight);
                     for (unsigned int index = 0; index < NewTriggers.size(); ++index){
                        if (triggerBits->accept(trigNames.triggerIndex(NewTriggers[index].first))){
                           h_tauE_tauE_Trigger_SS->Fill(index,NewTriggers[index].second*genWeight);
                        }
                        
                     }
                     for (unsigned int index = 0; index < BParkingTriggers.size(); ++index){
                        if (triggerBits->accept(trigNames.triggerIndex(BParkingTriggers[index].first))){
                           h_tauE_tauE_BPH_Trigger_SS->Fill(index,BParkingTriggers[index].second*genWeight);
                        }
                        
                     }
                     if (passTrig_Mu7IP4){
                        h_tauE_tauE_BPH_Trigger_SS->Fill(21.,7.e3*genWeight);
                     }
                     if (passTrig_Mu8IP3){
                        h_tauE_tauE_BPH_Trigger_SS->Fill(22.,1.45e3*genWeight);
                     }
                     if (passTrig_Mu8IP6){
                        h_tauE_tauE_BPH_Trigger_SS->Fill(23.,8.3e3*genWeight);
                     }
                     if (passTrig_Mu9IP6){
                        h_tauE_tauE_BPH_Trigger_SS->Fill(24.,33.7e3*genWeight);
                     }
                     h_tauE_tauE_nMatchedEle_SS->Fill(selected_std_b_electrons.size(),selected_std_tau_electrons.size(),genWeight);
                     h_tauE_tauE_sum_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
                     h_tauE_tauE_probb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
                     h_tauE_tauE_probbb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
                     h_tauE_tauE_problepb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
                  }
               }
            }
         }
      }
   }
   //std::cout<<"Tau E tauE\n";


   h_nEvent_tauE_tauMu->Fill("Total Events",genWeight);
   if (selected_std_electrons.size()>0){
      h_nEvent_tauE_tauMu->Fill("Std Tau Ele>0",genWeight);
      if (selected_Muons.size()>0){
         h_nEvent_tauE_tauMu->Fill("Tau Muons>0",genWeight);
         if (selected_bjets.size()>0){
            h_nEvent_tauE_tauMu->Fill("bjets>0",genWeight);
            TLorentzVector b;
            if (selected_bjets.size()>1){
               TLorentzVector b1;
               b1.SetPtEtaPhiM(selected_bjets[0].pt(),selected_bjets[0].eta(),selected_bjets[0].phi(),selected_bjets[0].mass());
               TLorentzVector b2;
               b2.SetPtEtaPhiM(selected_bjets[1].pt(),selected_bjets[1].eta(),selected_bjets[1].phi(),selected_bjets[1].mass());
               b = b1+b2;
            }
            else{
               b.SetPtEtaPhiM(selected_bjets[0].pt(),selected_bjets[0].eta(),selected_bjets[0].phi(),selected_bjets[0].mass());
            } 
            TLorentzVector e;
            e.SetPtEtaPhiM(selected_std_electrons[0].pt(),selected_std_electrons[0].eta(),selected_std_electrons[0].phi(),selected_std_electrons[0].mass());
            TLorentzVector mu;
            mu.SetPtEtaPhiM(selected_Muons[0].pt(),selected_Muons[0].eta(),selected_Muons[0].phi(),selected_Muons[0].mass());
            float mvis = (e+mu).M();
            float DiTauPt=(e+mu).Pt();
            float bMass=b.M();
            float lepDR = e.DeltaR(mu);
            float bDRl=b.DeltaR(mu);
            float bPt=b.Pt();
            if (bMass<18){
               h_nEvent_tauE_tauMu->Fill("bMass<18",genWeight);
               if (mvis<12&&mvis>3){
                  h_nEvent_tauE_tauMu->Fill("3<mvis<12",genWeight);
                  //if (passTrig_DoubleElectronLoose){
                     //h_nEvent_tauE_tauE->Fill("Double Electron", genWeight);
                  if (lepDR<1){
                     h_nEvent_tauE_tauMu->Fill("(lepDR<1)",genWeight);
                  //if (passTrig_MuEGLoose){
                     //h_nEvent_tauE_tauMu->Fill("MuEG Trigger",genWeight);
                     if(selected_std_electrons[0].charge()!=selected_Muons[0].charge()){
                        h_nEvent_tauE_tauMu->Fill("OS",genWeight);
                        h_tauMu_tauE_Mvis_OS->Fill(mvis, genWeight);
                        h_tauMu_tauE_bMass_OS->Fill(bMass, genWeight);
                        h_tauMu_tauE_DiTauPt_OS->Fill(DiTauPt, genWeight);
                        h_tauMu_tauE_bDRl_OS->Fill(bDRl, genWeight);
                        h_tauMu_tauE_lepDR_OS->Fill(lepDR, genWeight);
                        h_tauMu_tauE_MET_OS->Fill(MET, genWeight);
                        h_tauMu_tauE_bPt_OS->Fill(bPt, genWeight);
                        h_tauMu_tauE_nE_OS->Fill(selected_std_electrons.size(), genWeight);
                        h_tauMu_tauE_nTau_OS->Fill(selected_taus.size(), genWeight);
                        h_tauMu_tauE_nMu_OS->Fill(selected_Muons.size(), genWeight);
                        h_tauMu_tauE_nB_OS->Fill(selected_bjets.size(), genWeight);

                        for (unsigned int index = 0; index < NewTriggers.size(); ++index){
                           if (triggerBits->accept(trigNames.triggerIndex(NewTriggers[index].first))){
                              h_tauMu_tauE_Trigger_OS->Fill(index,NewTriggers[index].second*genWeight);
                           }
                           
                        }                        
                        for (unsigned int index = 0; index < BParkingTriggers.size(); ++index){
                           if (triggerBits->accept(trigNames.triggerIndex(BParkingTriggers[index].first))){
                              h_tauMu_tauE_BPH_Trigger_OS->Fill(index,BParkingTriggers[index].second*genWeight);
                           }
                           
                        }
                        if (passTrig_Mu7IP4){
                           h_tauMu_tauE_BPH_Trigger_OS->Fill(21.,7e3*genWeight);
                        }
                        if (passTrig_Mu8IP3){
                           h_tauMu_tauE_BPH_Trigger_OS->Fill(22.,1.45e3*genWeight);
                        }
                        if (passTrig_Mu8IP6){
                           h_tauMu_tauE_BPH_Trigger_OS->Fill(23.,8.3e3*genWeight);
                        }
                        if (passTrig_Mu9IP6){
                           h_tauMu_tauE_BPH_Trigger_OS->Fill(24.,33.7e3*genWeight);
                        }
                        h_tauMu_tauE_nMatchedEle_OS->Fill(selected_std_b_electrons.size(),selected_std_tau_electrons.size(),genWeight);
                        h_tauMu_tauE_nMatchedMu_OS->Fill(selected_b_gen_matched_Muons.size(),selected_tau_gen_matched_Muons.size(),genWeight);
                        h_tauMu_tauE_sum_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
                        h_tauMu_tauE_probb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
                        h_tauMu_tauE_probbb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
                        h_tauMu_tauE_problepb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
                     }
                     else{
                        h_nEvent_tauE_tauMu->Fill("SS",genWeight);
                        h_tauMu_tauE_Mvis_SS->Fill(mvis, genWeight);
                        h_tauMu_tauE_bMass_SS->Fill(bMass, genWeight);
                        h_tauMu_tauE_DiTauPt_SS->Fill(DiTauPt, genWeight);
                        h_tauMu_tauE_bDRl_SS->Fill(bDRl, genWeight);
                        h_tauMu_tauE_lepDR_SS->Fill(lepDR, genWeight);
                        h_tauMu_tauE_MET_SS->Fill(MET, genWeight);
                        h_tauMu_tauE_bPt_SS->Fill(bPt, genWeight);
                        h_tauMu_tauE_nE_SS->Fill(selected_std_electrons.size(), genWeight);
                        h_tauMu_tauE_nTau_SS->Fill(selected_taus.size(), genWeight);
                        h_tauMu_tauE_nMu_SS->Fill(selected_Muons.size(), genWeight);
                        h_tauMu_tauE_nB_SS->Fill(selected_bjets.size(), genWeight);
                        h_tauMu_tauE_nMatchedEle_SS->Fill(selected_std_b_electrons.size(), selected_std_tau_electrons.size(),genWeight);
                        h_tauMu_tauE_nMatchedMu_SS->Fill(selected_b_gen_matched_Muons.size(),selected_tau_gen_matched_Muons.size(),genWeight);
                        for (unsigned int index = 0; index < NewTriggers.size(); ++index){
                           if (triggerBits->accept(trigNames.triggerIndex(NewTriggers[index].first))){
                              h_tauMu_tauE_Trigger_SS->Fill(index,NewTriggers[index].second*genWeight);
                           }
                           
                        }                        
                        for (unsigned int index = 0; index < BParkingTriggers.size(); ++index){
                           if (triggerBits->accept(trigNames.triggerIndex(BParkingTriggers[index].first))){
                              h_tauMu_tauE_BPH_Trigger_SS->Fill(index,BParkingTriggers[index].second*genWeight);
                           }
                           
                        }
                        if (passTrig_Mu7IP4){
                           h_tauMu_tauE_BPH_Trigger_SS->Fill(21.,7.e3*genWeight);
                        }
                        if (passTrig_Mu8IP3){
                           h_tauMu_tauE_BPH_Trigger_SS->Fill(22.,1.45e3*genWeight);
                        }
                        if (passTrig_Mu8IP6){
                           h_tauMu_tauE_BPH_Trigger_SS->Fill(23.,8.3e3*genWeight);
                        }
                        if (passTrig_Mu9IP6){
                           h_tauMu_tauE_BPH_Trigger_SS->Fill(24.,33.7e3*genWeight);
                        }
                        h_tauMu_tauE_sum_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
                        h_tauMu_tauE_probb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
                        h_tauMu_tauE_probbb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
                        h_tauMu_tauE_problepb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
                     }
                  //}
                  //if (passTrig_doubleMuonLoose){
                  //   h_nEvent_tauE_tauE->Fill("doubleMu Trig*", genWeight);
                  }
               }
            }
         }
      }
   }

   //std::cout<<"Tau E tauMu\n";


   h_nEvent_tauMu_tauMu->Fill("Total Events",genWeight);
   if (selected_Muons.size()>1){
      //std::cout<<"size>1\n";

      h_nEvent_tauMu_tauMu->Fill("Std Tau Mu>1",genWeight);
      if (selected_bjets.size()>0){
         TLorentzVector b;
         if (selected_bjets.size()>1){
            TLorentzVector b1;
            b1.SetPtEtaPhiM(selected_bjets[0].pt(),selected_bjets[0].eta(),selected_bjets[0].phi(),selected_bjets[0].mass());
            TLorentzVector b2;
            b2.SetPtEtaPhiM(selected_bjets[1].pt(),selected_bjets[1].eta(),selected_bjets[1].phi(),selected_bjets[1].mass());
            b = b1+b2;
         }
         else{
            b.SetPtEtaPhiM(selected_bjets[0].pt(),selected_bjets[0].eta(),selected_bjets[0].phi(),selected_bjets[0].mass());
         } 
         h_nEvent_tauMu_tauMu->Fill("Njets>0",genWeight);
         TLorentzVector mu1;
         mu1.SetPtEtaPhiM(selected_Muons[0].pt(),selected_Muons[0].eta(),selected_Muons[0].phi(),selected_Muons[0].mass());
         TLorentzVector mu2;
         mu2.SetPtEtaPhiM(selected_Muons[1].pt(),selected_Muons[1].eta(),selected_Muons[1].phi(),selected_Muons[1].mass());
         float mvis = (mu1+mu2).M();
         float DiTauPt=(mu2+mu1).Pt();
         float bMass=b.M();
         float bDRl=b.DeltaR(mu1);
         float lepDR=mu1.DeltaR(mu2);
         float bPt=b.Pt();
         if (bMass<18){
            h_nEvent_tauMu_tauMu->Fill("bMass<18",genWeight);
            if (mvis<12&&mvis>3){
               h_nEvent_tauMu_tauMu->Fill("3<mvis<12",genWeight);
               //if (passTrig_DoubleElectronLoose){
                  //h_nEvent_tauE_tauE->Fill("Double Electron", genWeight);
               if (lepDR<1){
                  h_nEvent_tauMu_tauMu->Fill("(lepDR<1)",genWeight);
                  if(selected_Muons[1].charge()!=selected_Muons[0].charge()){
                     //std::cout<<"Charge1\n";
                     for (unsigned int index = 0; index < NewTriggers.size(); ++index){
                        if (triggerBits->accept(trigNames.triggerIndex(NewTriggers[index].first))){
                           h_tauMu_tauMu_Trigger_OS->Fill(index,NewTriggers[index].second*genWeight);
                        }
                        
                     }
                     for (unsigned int index = 0; index < BParkingTriggers.size(); ++index){
                        if (triggerBits->accept(trigNames.triggerIndex(BParkingTriggers[index].first))){
                           h_tauMu_tauMu_BPH_Trigger_OS->Fill(index,BParkingTriggers[index].second*genWeight);
                        }
                        
                     }
                     if (passTrig_Mu7IP4){
                        h_tauMu_tauMu_BPH_Trigger_OS->Fill(21.,7e3*genWeight);
                     }
                     if (passTrig_Mu8IP3){
                        h_tauMu_tauMu_BPH_Trigger_OS->Fill(22.,1.45e3*genWeight);
                     }
                     if (passTrig_Mu8IP6){
                        h_tauMu_tauMu_BPH_Trigger_OS->Fill(23.,8.3e3*genWeight);
                     }
                     if (passTrig_Mu9IP6){
                        h_tauMu_tauMu_BPH_Trigger_OS->Fill(24.,33.7e3*genWeight);
                     }
                     h_nEvent_tauMu_tauMu->Fill("OS",genWeight);
                     h_tauMu_tauMu_Mvis_OS->Fill(mvis, genWeight);
                     h_tauMu_tauMu_bMass_OS->Fill(bMass, genWeight);
                     h_tauMu_tauMu_DiTauPt_OS->Fill(DiTauPt, genWeight);
                     h_tauMu_tauMu_bDRl_OS->Fill(bDRl, genWeight);
                     h_tauMu_tauMu_lepDR_OS->Fill(lepDR, genWeight);
                     h_tauMu_tauMu_MET_OS->Fill(MET, genWeight);
                     h_tauMu_tauMu_bPt_OS->Fill(bPt, genWeight);
                     
                     h_tauMu_tauMu_nE_OS->Fill(selected_std_electrons.size(), genWeight);
                     h_tauMu_tauMu_nTau_OS->Fill(selected_taus.size(), genWeight);
                     h_tauMu_tauMu_nMu_OS->Fill(selected_Muons.size(), genWeight);
                     h_tauMu_tauMu_nB_OS->Fill(selected_bjets.size(), genWeight);
                     h_tauMu_tauMu_nMatchedMu_OS->Fill(selected_b_gen_matched_Muons.size(),selected_tau_gen_matched_Muons.size(),genWeight);
                     h_tauMu_tauMu_sum_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
                     h_tauMu_tauMu_probb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
                     h_tauMu_tauMu_probbb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
                     h_tauMu_tauMu_problepb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
                  }
                  else{
                     //std::cout<<"Charge2\n";
                     h_nEvent_tauMu_tauMu->Fill("SS",genWeight);
                     h_tauMu_tauMu_Mvis_SS->Fill(mvis, genWeight);
                     h_tauMu_tauMu_bMass_SS->Fill(bMass, genWeight);
                     h_tauMu_tauMu_DiTauPt_SS->Fill(DiTauPt, genWeight);
                     h_tauMu_tauMu_bDRl_SS->Fill(bDRl, genWeight);
                     h_tauMu_tauMu_lepDR_SS->Fill(lepDR, genWeight);
                     h_tauMu_tauMu_MET_SS->Fill(MET, genWeight);
                     h_tauMu_tauMu_bPt_SS->Fill(bPt, genWeight);
                     h_tauMu_tauMu_nE_SS->Fill(selected_std_electrons.size(), genWeight);
                     h_tauMu_tauMu_nTau_SS->Fill(selected_taus.size(), genWeight);
                     h_tauMu_tauMu_nMu_SS->Fill(selected_Muons.size(), genWeight);
                     h_tauMu_tauMu_nB_SS->Fill(selected_bjets.size(), genWeight);
                     //std::cout<<"Charge1\n";
                     for (unsigned int index = 0; index < NewTriggers.size(); ++index){
                        if (triggerBits->accept(trigNames.triggerIndex(NewTriggers[index].first))){
                           h_tauMu_tauMu_Trigger_SS->Fill(index,NewTriggers[index].second*genWeight);
                        }
                        
                     }
                     for (unsigned int index = 0; index < BParkingTriggers.size(); ++index){
                        if (triggerBits->accept(trigNames.triggerIndex(BParkingTriggers[index].first))){
                           h_tauMu_tauMu_BPH_Trigger_SS->Fill(index,BParkingTriggers[index].second*genWeight);
                        }
                        
                     }
                     if (passTrig_Mu7IP4){
                        h_tauMu_tauMu_BPH_Trigger_SS->Fill(21.,7.e3*genWeight);
                     }
                     if (passTrig_Mu8IP3){
                        h_tauMu_tauMu_BPH_Trigger_SS->Fill(22.,1.45e3*genWeight);
                     }
                     if (passTrig_Mu8IP6){
                        h_tauMu_tauMu_BPH_Trigger_SS->Fill(23.,8.3e3*genWeight);
                     }
                     if (passTrig_Mu9IP6){
                        h_tauMu_tauMu_BPH_Trigger_SS->Fill(24.,33.7e3*genWeight);
                     }
                     h_tauMu_tauMu_nMatchedMu_SS->Fill(selected_b_gen_matched_Muons.size(),selected_tau_gen_matched_Muons.size(),genWeight);
                     h_tauMu_tauMu_sum_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
                     h_tauMu_tauMu_probb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
                     h_tauMu_tauMu_probbb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
                     h_tauMu_tauMu_problepb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
                  }
               //}
               }
               
               //std::cout<<"Posttrigger1\n";
            }
         }
         //std::cout<<"postrigger2\n";
      }
     //std::cout<<"postrigger3\n"; 
   }
   

   //std::cout<<"Tau Mu tauMu\n";
   
   h_nEvent_tauHad_tauE->Fill("Total Events",genWeight);
   if (selected_taus.size()>0){
      h_nEvent_tauHad_tauE->Fill("Tau Had>0",genWeight);
      if (selected_std_electrons.size()>0){
         h_nEvent_tauHad_tauE->Fill("Tau Electrons>0",genWeight);
         if (selected_bjets.size()>0){
            TLorentzVector b;
            if (selected_bjets.size()>1){
               TLorentzVector b1;
               b1.SetPtEtaPhiM(selected_bjets[0].pt(),selected_bjets[0].eta(),selected_bjets[0].phi(),selected_bjets[0].mass());
               TLorentzVector b2;
               b2.SetPtEtaPhiM(selected_bjets[1].pt(),selected_bjets[1].eta(),selected_bjets[1].phi(),selected_bjets[1].mass());
               b = b1+b2;
            }
            else{
               b.SetPtEtaPhiM(selected_bjets[0].pt(),selected_bjets[0].eta(),selected_bjets[0].phi(),selected_bjets[0].mass());
            } 
            h_nEvent_tauHad_tauE->Fill("selected bjets>0",genWeight);

            TLorentzVector t;
            t.SetPtEtaPhiM(selected_taus[0].pt(),selected_taus[0].eta(),selected_taus[0].phi(),selected_taus[0].mass());
            TLorentzVector e;
            e.SetPtEtaPhiM(selected_std_electrons[0].pt(),selected_std_electrons[0].eta(),selected_std_electrons[0].phi(),selected_std_electrons[0].mass());
            float mvis = (t+e).M();
            float DiTauPt=(t+e).Pt();
            float bMass=b.M();
            float bDRl=b.DeltaR(e);
            float lepDR=e.DeltaR(t);
            float bPt=b.Pt();
            if (bMass<18){
               h_nEvent_tauHad_tauE->Fill("bMass<18",genWeight);
               if (mvis<12&&mvis>3){
                  h_nEvent_tauHad_tauE->Fill("3<mvis<12",genWeight);
                  //if (passTrig_DoubleElectronLoose){
                     //h_nEvent_tauE_tauE->Fill("Double Electron", genWeight);
                  if (lepDR<1){
                     h_nEvent_tauHad_tauE->Fill("(lepDR<1)",genWeight);
                     if(selected_taus[0].charge()!=selected_std_electrons[0].charge()){
                        h_nEvent_tauHad_tauE->Fill("OS",genWeight);
                        h_tauHad_tauE_Mvis_OS->Fill(mvis,genWeight);
                        h_tauHad_tauE_bDRl_OS->Fill(bDRl,genWeight);
                        h_tauHad_tauE_bMass_OS->Fill(bMass,genWeight);
                        h_tauHad_tauE_DiTauPt_OS->Fill(DiTauPt,genWeight);
                        h_tauHad_tauE_lepDR_OS->Fill(lepDR,genWeight);
                        h_tauHad_tauE_bPt_OS->Fill(bPt, genWeight);
                        h_tauHad_tauE_MET_OS->Fill(MET,genWeight);
                        h_tauHad_tauE_nE_OS->Fill(selected_std_electrons.size(), genWeight);
                        h_tauHad_tauE_nTau_OS->Fill(selected_taus.size(), genWeight);
                        h_tauHad_tauE_nMu_OS->Fill(selected_Muons.size(), genWeight);
                        h_tauHad_tauE_nB_OS->Fill(selected_bjets.size(), genWeight);
                        h_tauHad_tauE_nMatchedEle_OS->Fill(selected_std_b_electrons.size(),selected_std_tau_electrons.size(),genWeight);
                        h_tauHad_tauE_sum_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
                        h_tauHad_tauE_probb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
                        h_tauHad_tauE_probbb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
                        h_tauHad_tauE_problepb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
                        for (unsigned int index = 0; index < NewTriggers.size(); ++index){
                           if (triggerBits->accept(trigNames.triggerIndex(NewTriggers[index].first))){
                              h_tauHad_tauE_Trigger_OS->Fill(index,NewTriggers[index].second*genWeight);
                           }
                           
                        }
                        for (unsigned int index = 0; index < BParkingTriggers.size(); ++index){
                           if (triggerBits->accept(trigNames.triggerIndex(BParkingTriggers[index].first))){
                              h_tauHad_tauE_BPH_Trigger_OS->Fill(index,BParkingTriggers[index].second*genWeight);
                           }
                           
                        }
                        if (passTrig_Mu7IP4){
                           h_tauHad_tauE_BPH_Trigger_OS->Fill(21.,7e3*genWeight);
                        }
                        if (passTrig_Mu8IP3){
                           h_tauHad_tauE_BPH_Trigger_OS->Fill(22.,1.45e3*genWeight);
                        }
                        if (passTrig_Mu8IP6){
                           h_tauHad_tauE_BPH_Trigger_OS->Fill(23.,8.3e3*genWeight);
                        }
                        if (passTrig_Mu9IP6){
                           h_tauHad_tauE_BPH_Trigger_OS->Fill(24.,33.7e3*genWeight);
                        }
                     }
                     else{
                        h_nEvent_tauHad_tauE->Fill("SS",genWeight);
                        h_tauHad_tauE_Mvis_SS->Fill(mvis,genWeight);
                        h_tauHad_tauE_bDRl_SS->Fill(bDRl,genWeight);
                        h_tauHad_tauE_bMass_SS->Fill(bMass,genWeight);
                        h_tauHad_tauE_DiTauPt_SS->Fill(DiTauPt,genWeight);
                        h_tauHad_tauE_lepDR_SS->Fill(lepDR,genWeight);
                        h_tauHad_tauE_MET_SS->Fill(MET,genWeight);
                        h_tauHad_tauE_bPt_SS->Fill(bPt, genWeight);
                        h_tauHad_tauE_nE_SS->Fill(selected_std_electrons.size(), genWeight);
                        h_tauHad_tauE_nTau_SS->Fill(selected_taus.size(), genWeight);
                        h_tauHad_tauE_nMu_SS->Fill(selected_Muons.size(), genWeight);
                        h_tauHad_tauE_nB_SS->Fill(selected_bjets.size(), genWeight);
                        h_tauHad_tauE_nMatchedEle_SS->Fill(selected_std_b_electrons.size(),selected_std_tau_electrons.size(),genWeight);
                        h_tauHad_tauE_sum_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
                        h_tauHad_tauE_probb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
                        h_tauHad_tauE_probbb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
                        h_tauHad_tauE_problepb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
                        for (unsigned int index = 0; index < NewTriggers.size(); ++index){
                           if (triggerBits->accept(trigNames.triggerIndex(NewTriggers[index].first))){
                              h_tauHad_tauE_Trigger_SS->Fill(index,NewTriggers[index].second*genWeight);
                           }
                           
                        }
                        for (unsigned int index = 0; index < BParkingTriggers.size(); ++index){
                           if (triggerBits->accept(trigNames.triggerIndex(BParkingTriggers[index].first))){
                              h_tauHad_tauE_BPH_Trigger_SS->Fill(index,BParkingTriggers[index].second*genWeight);
                           }
                           
                        }
                        if (passTrig_Mu7IP4){
                           h_tauHad_tauE_BPH_Trigger_SS->Fill(21.,7.e3*genWeight);
                        }
                        if (passTrig_Mu8IP3){
                           h_tauHad_tauE_BPH_Trigger_SS->Fill(22.,1.45e3*genWeight);
                        }
                        if (passTrig_Mu8IP6){
                           h_tauHad_tauE_BPH_Trigger_SS->Fill(23.,8.3e3*genWeight);
                        }
                        if (passTrig_Mu9IP6){
                           h_tauHad_tauE_BPH_Trigger_SS->Fill(24.,33.7e3*genWeight);
                        }
                     }
                  }
               }
            }   
         }
      }
   }

   
   //std::cout<<"Tau E tauHad\n";

   h_nEvent_tauHad_tauMu->Fill("Total Events",genWeight);
   if (selected_taus.size()>0){
      //std::cout<<"Muhad Size\n";
      h_nEvent_tauHad_tauMu->Fill("Selected Tau Had>0",genWeight);
      if (selected_Muons.size()>0){
         h_nEvent_tauHad_tauMu->Fill("Tau Muons>0",genWeight);
         if (selected_bjets.size()>0){
            TLorentzVector b;
            if (selected_bjets.size()>1){
               TLorentzVector b1;
               b1.SetPtEtaPhiM(selected_bjets[0].pt(),selected_bjets[0].eta(),selected_bjets[0].phi(),selected_bjets[0].mass());
               TLorentzVector b2;
               b2.SetPtEtaPhiM(selected_bjets[1].pt(),selected_bjets[1].eta(),selected_bjets[1].phi(),selected_bjets[1].mass());
               b = b1+b2;
            }
            else{
               b.SetPtEtaPhiM(selected_bjets[0].pt(),selected_bjets[0].eta(),selected_bjets[0].phi(),selected_bjets[0].mass());
            } 
            h_nEvent_tauHad_tauMu->Fill("selected bjets>0",genWeight);

            TLorentzVector t;
            t.SetPtEtaPhiM(selected_taus[0].pt(),selected_taus[0].eta(),selected_taus[0].phi(),selected_taus[0].mass());
            TLorentzVector mu;
            mu.SetPtEtaPhiM(selected_Muons[0].pt(),selected_Muons[0].eta(),selected_Muons[0].phi(),selected_Muons[0].mass());
            float mvis = (t+mu).M();
            float DiTauPt=(t+mu).Pt();
            float bMass=b.M();
            float bDRl=b.DeltaR(mu);
            float lepDR=mu.DeltaR(t);
            float bPt=b.Pt();
            if (bMass<18){
               h_nEvent_tauHad_tauMu->Fill("bMass<18",genWeight);
               if (mvis<12&&mvis>3){
                  h_nEvent_tauHad_tauMu->Fill("3<mvis<12",genWeight);
                  //if (passTrig_DoubleElectronLoose){
                     //h_nEvent_tauE_tauE->Fill("Double Electron", genWeight);
                  if (lepDR<1){
                     h_nEvent_tauHad_tauMu->Fill("(lepDR<1)",genWeight);
                     if(selected_taus[0].charge()!=selected_Muons[0].charge()){
                        h_nEvent_tauHad_tauMu->Fill("OS",genWeight);
                        h_tauHad_tauMu_Mvis_OS->Fill(mvis,genWeight);
                        h_tauHad_tauMu_bDRl_OS->Fill(bDRl,genWeight);
                        h_tauHad_tauMu_bMass_OS->Fill(bMass,genWeight);
                        h_tauHad_tauMu_DiTauPt_OS->Fill(DiTauPt,genWeight);
                        h_tauHad_tauMu_lepDR_OS->Fill(lepDR,genWeight);
                        h_tauHad_tauMu_MET_OS->Fill(MET,genWeight);
                        h_tauHad_tauMu_bPt_OS->Fill(bPt, genWeight);
                        h_tauHad_tauMu_nMatchedMu_OS->Fill(selected_b_gen_matched_Muons.size(),selected_tau_gen_matched_Muons.size(),genWeight);
                        h_tauHad_tauMu_sum_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
                        h_tauHad_tauMu_probb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
                        h_tauHad_tauMu_probbb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
                        h_tauHad_tauMu_problepb_OS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
                        h_tauHad_tauMu_nE_OS->Fill(selected_std_electrons.size(), genWeight);
                        h_tauHad_tauMu_nTau_OS->Fill(selected_taus.size(), genWeight);
                        h_tauHad_tauMu_nMu_OS->Fill(selected_Muons.size(), genWeight);
                        h_tauHad_tauMu_nB_OS->Fill(selected_bjets.size(), genWeight);
                        for (unsigned int index = 0; index < NewTriggers.size(); ++index){
                           //std::cout<<NewTriggers[index].first<<"\n";
                           if (triggerBits->accept(trigNames.triggerIndex(NewTriggers[index].first))){
                              //std::cout<<"passcheck\n";
                              h_tauHad_tauMu_Trigger_OS->Fill(index,NewTriggers[index].second*genWeight);
                           }
                           
                        }
                        //std::cout<<"MuHad Normal Triggers OS\n";
                        for (unsigned int index = 0; index < BParkingTriggers.size(); ++index){
                           if (triggerBits->accept(trigNames.triggerIndex(BParkingTriggers[index].first))){
                              h_tauHad_tauMu_BPH_Trigger_OS->Fill(index,BParkingTriggers[index].second*genWeight);
                           }
                           
                        }
                        
                        if (passTrig_Mu7IP4){
                           h_tauHad_tauMu_BPH_Trigger_OS->Fill(21.,7e3*genWeight);
                        }
                        if (passTrig_Mu8IP3){
                           h_tauHad_tauMu_BPH_Trigger_OS->Fill(22.,1.45e3*genWeight);
                        }
                        if (passTrig_Mu8IP6){
                           h_tauHad_tauMu_BPH_Trigger_OS->Fill(23.,8.3e3*genWeight);
                        }
                        if (passTrig_Mu9IP6){
                           h_tauHad_tauMu_BPH_Trigger_OS->Fill(24.,33.7e3*genWeight);
                        }
                        //std::cout<<"MuHad BPH Triggers OS\n";
                     }
                     else{
                        h_nEvent_tauHad_tauMu->Fill("SS",genWeight);
                        h_tauHad_tauMu_Mvis_SS->Fill(mvis,genWeight);
                        h_tauHad_tauMu_bDRl_SS->Fill(bDRl,genWeight);
                        h_tauHad_tauMu_bMass_SS->Fill(bMass,genWeight);
                        h_tauHad_tauMu_DiTauPt_SS->Fill(DiTauPt,genWeight);
                        h_tauHad_tauMu_lepDR_SS->Fill(lepDR,genWeight);
                        h_tauHad_tauMu_MET_SS->Fill(MET,genWeight);
                        h_tauHad_tauMu_bPt_SS->Fill(bPt, genWeight);
                        h_tauHad_tauMu_nE_SS->Fill(selected_std_electrons.size(), genWeight);
                        h_tauHad_tauMu_nTau_SS->Fill(selected_taus.size(), genWeight);
                        h_tauHad_tauMu_nMu_SS->Fill(selected_Muons.size(), genWeight);
                        h_tauHad_tauMu_nB_SS->Fill(selected_bjets.size(), genWeight);
                        h_tauHad_tauMu_nMatchedMu_SS->Fill(selected_b_gen_matched_Muons.size(),selected_tau_gen_matched_Muons.size(),genWeight);
                        h_tauHad_tauMu_sum_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
                        h_tauHad_tauMu_probb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
                        h_tauHad_tauMu_probbb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
                        h_tauHad_tauMu_problepb_SS->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
                        for (unsigned int index = 0; index < NewTriggers.size(); ++index){
                           //std::cout<<NewTriggers[index].first<<"\n";
                           if (triggerBits->accept(trigNames.triggerIndex(NewTriggers[index].first))){
                              //std::cout<<"passcheck\n";
                              h_tauHad_tauMu_Trigger_SS->Fill(index,NewTriggers[index].second*genWeight);
                           }
                           
                        }   
                        //std::cout<<"MuHad Norm Triggers SS\n";               
                        for (unsigned int index = 0; index < BParkingTriggers.size(); ++index){
                           if (triggerBits->accept(trigNames.triggerIndex(BParkingTriggers[index].first))){
                              h_tauHad_tauMu_BPH_Trigger_SS->Fill(index,BParkingTriggers[index].second*genWeight);
                           }
                           
                        }
                        if (passTrig_Mu7IP4){
                           h_tauHad_tauMu_BPH_Trigger_SS->Fill(21.,7.e3*genWeight);
                        }
                        if (passTrig_Mu8IP3){
                           h_tauHad_tauMu_BPH_Trigger_SS->Fill(22.,1.45e3*genWeight);
                        }
                        if (passTrig_Mu8IP6){
                           h_tauHad_tauMu_BPH_Trigger_SS->Fill(23.,8.3e3*genWeight);
                        }
                        if (passTrig_Mu9IP6){
                           h_tauHad_tauMu_BPH_Trigger_SS->Fill(24.,33.7e3*genWeight);
                        }
                        //std::cout<<"MuHad BPH Triggers SS\n"; 
                     }
                  }
               }
            }
      
         }
      }
   }

   //std::cout<<"Tau Mu tauHad\n";

   h_nEvent_gen_tauHad_tauMu->Fill("Total Events",genWeight);
   if (genTauHad.size()>0&&genNuTauHad.size()>0&&DiTauHad==false){
      h_nEvent_gen_tauHad_tauMu->Fill("Gen Tau Had==1",genWeight);
      //std::cout<<"Ehad Size\n";
      if (gentauMuons.size()>0){
         h_nEvent_gen_tauHad_tauMu->Fill("Gen Tau Muons>0",genWeight);
         //std::cout<<"GenTauMu\n";
         if (selected_bjets.size()>0){
            //std::cout<<"Tau bsize\n";
            TLorentzVector b;
            b.SetPtEtaPhiM(selected_bjets[0].pt(),selected_bjets[0].eta(),selected_bjets[0].phi(),selected_bjets[0].mass());
            h_nEvent_gen_tauHad_tauMu->Fill("selected bjets>0",genWeight);
            TLorentzVector nu;
            nu.SetPtEtaPhiM(genNuTauHad[0]->pt(),genNuTauHad[0]->eta(),genNuTauHad[0]->phi(),genNuTauHad[0]->mass());
            TLorentzVector t;
            t.SetPtEtaPhiM(genTauHad[0]->pt(),genTauHad[0]->eta(),genTauHad[0]->phi(),genTauHad[0]->mass());
            TLorentzVector mu;
            mu.SetPtEtaPhiM(gentauMuons[0]->pt(),gentauMuons[0]->eta(),gentauMuons[0]->phi(),gentauMuons[0]->mass());
            float mvis = (t+mu-nu).M();
            float DiTauPt=(t+mu-nu).Pt();
            float bMass=b.M();
            float bDRl=b.DeltaR(mu);
            float lepDR=mu.DeltaR(t-nu);

            if (passTrig_singleMuLoose){
               //std::cout<<"Trigger\n";
               h_nEvent_gen_tauHad_tauMu->Fill("Single Mu Trigger",genWeight);
               h_gen_tauHad_tauMu_Mvis->Fill(mvis,genWeight);
               h_gen_tauHad_tauMu_bDRl->Fill(bDRl,genWeight);
               h_gen_tauHad_tauMu_bMass->Fill(bMass,genWeight);
               h_gen_tauHad_tauMu_DiTauPt->Fill(DiTauPt,genWeight);
               h_gen_tauHad_tauMu_lepDR->Fill(lepDR,genWeight);
               h_gen_tauHad_tauMu_MET->Fill(MET,genWeight);
               h_gen_tauHad_tauMu_sum->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
               h_gen_tauHad_tauMu_probb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
               h_gen_tauHad_tauMu_probbb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
               h_gen_tauHad_tauMu_problepb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
            
            }
                  
         }
      }
   }
   //std::cout<<"Tau Mu gen_tauHad\n";
   
   h_nEvent_gen_tauHad_tauE->Fill("Total Events",genWeight);
   if (genTauHad.size()>0&&genNuTauHad.size()>0&&DiTauHad==false){
      h_nEvent_gen_tauHad_tauE->Fill("Gen Tau Had==1",genWeight);
      if (gentauElectrons.size()>0){
         h_nEvent_gen_tauHad_tauE->Fill("Gen Tau Electrons>0",genWeight);
         if (selected_bjets.size()>0){
            TLorentzVector b;
            b.SetPtEtaPhiM(selected_bjets[0].pt(),selected_bjets[0].eta(),selected_bjets[0].phi(),selected_bjets[0].mass());
            h_nEvent_gen_tauHad_tauE->Fill("selected bjets>0",genWeight);
            TLorentzVector nu;
            nu.SetPtEtaPhiM(genNuTauHad[0]->pt(),genNuTauHad[0]->eta(),genNuTauHad[0]->phi(),genNuTauHad[0]->mass());

            TLorentzVector t;
            t.SetPtEtaPhiM(genTauHad[0]->pt(),genTauHad[0]->eta(),genTauHad[0]->phi(),genTauHad[0]->mass());
            TLorentzVector e;
            e.SetPtEtaPhiM(gentauElectrons[0]->pt(),gentauElectrons[0]->eta(),gentauElectrons[0]->phi(),gentauElectrons[0]->mass());
            float mvis = (t+e-nu).M();
            float DiTauPt=(t+e-nu).Pt();
            float bMass=b.M();
            float bDRl=b.DeltaR(e);
            float lepDR=e.DeltaR(t-nu);

            if (passTrig_EGLoose){
               h_nEvent_gen_tauHad_tauE->Fill("Single Ele Trigger",genWeight);
               h_gen_tauHad_tauE_Mvis->Fill(mvis,genWeight);
               h_gen_tauHad_tauE_bDRl->Fill(bDRl,genWeight);
               h_gen_tauHad_tauE_bMass->Fill(bMass,genWeight);
               h_gen_tauHad_tauE_DiTauPt->Fill(DiTauPt,genWeight);
               h_gen_tauHad_tauE_lepDR->Fill(lepDR,genWeight);
               h_gen_tauHad_tauE_MET->Fill(MET,genWeight);
               h_gen_tauHad_tauE_sum->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb")  + selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb"),genWeight);
               h_gen_tauHad_tauE_probb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probb") ,genWeight);
               h_gen_tauHad_tauE_probbb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:probbb") ,genWeight);
               h_gen_tauHad_tauE_problepb->Fill(selected_bjets[0].bDiscriminator("pfDeepFlavourJetTags:problepb") ,genWeight);
            
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
