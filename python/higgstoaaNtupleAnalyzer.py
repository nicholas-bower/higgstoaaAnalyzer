import ROOT,sys,os
import numpy as np
ROOT.gInterpreter.Declare('#include "../interface/JetInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../interface/MuonInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../interface/ElectronInfoDS.h"')
ROOT.gInterpreter.Declare('#include "../interface/TauInfoDS.h"')
location = os.getcwd() 
query = "ls *.root"
infiles=os.popen(query).read().split()
inputFileList=[]
for line in infiles:
    inputFileList.append("file:"+location+"/"+line)
    print ("file:"+location+"/"+line)


#out=ROOT.TFile.Open("testEMu.root",'recreate')
histos={}
out=ROOT.TFile.Open("Outfile.root",'recreate')
histos["h_nEvents"]=ROOT.TH1F("h_nEvents", "nEvents", 2,0,2)
histos["h_std_electron_selection"]=ROOT.TH1F("h_std_electron_selection", "StdElectron Cuts", 3,0,3)
histos["h_std_electron_selection"].SetCanExtend( ROOT.TH1.kAllAxes)
histos['h_muon_selection']=ROOT.TH1F("h_muon_selection", "Muon Cuts", 3,0,3)
histos['h_muon_selection'].SetCanExtend(ROOT.TH1.kAllAxes)
histos['h_eCtau_selection']=ROOT.TH1F("h_eCtau_selection", "e Clean Tau Cuts", 3,0,3)
histos['h_eCtau_selection'].SetCanExtend(ROOT.TH1.kAllAxes)
histos['h_muCtau_selection']=ROOT.TH1F("h_muCtau_selection", "Mu Clean Tau Cuts", 3,0,3)
histos['h_muCtau_selection'].SetCanExtend(ROOT.TH1.kAllAxes)
histos['h_bJet_selection']=ROOT.TH1F("h_bJet_selection", "b Jet Cuts", 3,0,3)
histos['h_bJet_selection'].SetCanExtend(ROOT.TH1.kAllAxes)

histos['h_nStdEle']=ROOT.TH1F("h_nStdEle", "nTau Ele per Event (all)", 4,0,4)
histos['h_nStdEle_taum']=ROOT.TH1F("h_nStdEle_taum", "nTau Ele per Event (Tau Matched)", 4,0,4)
histos['h_nStdEle_bm']=ROOT.TH1F("h_nStdEle_bm", "nTau Ele per Event (b Matched)", 4,0,4)
histos['h_StdEle_taum_pt']=ROOT.TH1F("h_StdEle_taum_pt", "Tau matched Electron Pt", 40,0,200)
histos['h_StdEle_bm_pt']=ROOT.TH1F("h_StdEle_bm_pt", "b matched Electron Pt",40,0,200)

histos['h_nMu']=ROOT.TH1F("h_nMu", "nTau Mu per Event (all)", 4,0,4)
histos['h_nMu_taum']=ROOT.TH1F("h_nMu_taum", "nTau Mu per Event (Tau Matched)", 4,0,4)
histos['h_nMu_bm']=ROOT.TH1F("h_nMu_bm", "nTau Mu per Event (b Matched)", 4,0,4)
histos['h_Mu_taum_pt']=ROOT.TH1F("h_Mu_taum_pt", "Tau matched Muon Pt", 40,0,200)
histos['h_Mu_bm_pt']=ROOT.TH1F("h_Mu_bm_pt", "b matched Muon Pt",40,0,200)

histos['h_nbJet']=ROOT.TH1F("h_nbJet", "bJet per Event (all)", 4,0,4)
histos['h_nbJet_m']=ROOT.TH1F("h_nbJet_m", "bJet per Event (Matched)", 4,0,4)
histos['h_nbJet_f']=ROOT.TH1F("h_nbJet_f", "bJet per Event (!Matched)", 4,0,4)

histos['h_deepFlavor_Merged_probb']=ROOT.TH1F("h_deepFlavor_Merged_probb", "pfDeepFlavourJetTags:probb nJets = 1", 20,0,1)
histos['h_deepFlavor_Merged_probbb']=ROOT.TH1F("h_deepFlavor_Merged_probbb", "pfDeepFlavourJetTags:probbb nJets = 1", 20,0,1)
histos['h_deepFlavor_Merged_problepb']=ROOT.TH1F("h_deepFlavor_Merged_problepb", "pfDeepFlavourJetTags:problepb nJets = 1", 20,0,1)
histos['h_deepFlavor_Merged_sum']=ROOT.TH1F("h_deepFlavor_Merged_sum", "Total DeepFlavor nJets = 1", 20,0,1)
histos['h_deepFlavor_Merged_probbvsprobbb']=ROOT.TH2F("h_deepFlavor_Merged_probbvsprobbb", "prob b vs prob bb nJets = 1", 20,0,1,20,0,1)

histos['h_deepFlavor_Split_Lead_probb']=ROOT.TH1F("h_deepFlavor_Split_Lead_probb", "pfDeepFlavourJetTags:probb nJets>1 Leading", 20,0,1)
histos['h_deepFlavor_Split_Lead_probbb']=ROOT.TH1F("h_deepFlavor_Split_Lead_probbb", "pfDeepFlavourJetTags:probbb nJets>1 Leading", 20,0,1)
histos['h_deepFlavor_Split_Lead_problepb']=ROOT.TH1F("h_deepFlavor_Split_Lead_problepb", "pfDeepFlavourJetTags:problepb nJets>1 Leading", 20,0,1)
histos['h_deepFlavor_Split_Lead_sum']=ROOT.TH1F("h_deepFlavor_Split_Lead_sum", "Total DeepFlavor nJets>1 Leading", 20,0,1)
histos['h_deepFlavor_Split_Lead_probbvsprobbb']=ROOT.TH2F("h_deepFlavor_Split_Lead_probbvsprobbb", "prob b vs prob bb nJets>1 Leading", 20,0,1,20,0,1)

histos['h_deepFlavor_Split_Sublead_probb']=ROOT.TH1F("h_deepFlavor_Split_Sublead_probb", "pfDeepFlavourJetTags:probb nJets>1 Subleading", 20,0,1)
histos['h_deepFlavor_Split_Sublead_probbb']=ROOT.TH1F("h_deepFlavor_Split_Sublead_probbb", "pfDeepFlavourJetTags:probbb nJets>1 Subleading", 20,0,1)
histos['h_deepFlavor_Split_Sublead_problepb']=ROOT.TH1F("h_deepFlavor_Split_Sublead_problepb", "pfDeepFlavourJetTags:problepb nJets>1 Subleading", 20,0,1)
histos['h_deepFlavor_Split_Sublead_sum']=ROOT.TH1F("h_deepFlavor_Split_Sublead_sum", "Total DeepFlavor nJets>1 Subleading", 20,0,1)
histos['h_deepFlavor_Split_Sublead_probbvsprobbb']=ROOT.TH2F("h_deepFlavor_Split_Sublead_probbvsprobbb", "prob b vs prob bb nJets>1 Subleading", 20,0,1,20,0,1)



histos['h_CSV_Merged_probb']=ROOT.TH1F("h_CSV_Merged_probb", "pfCSVJetTags:probb nJets = 1", 20,0,1)
histos['h_CSV_Merged_probbb']=ROOT.TH1F("h_CSV_Merged_probbb", "pfCSVJetTags:probbb nJets = 1", 20,0,1)
histos['h_CSV_Merged_sum']=ROOT.TH1F("h_CSV_Merged_sum", "Total CSV nJets = 1", 20,0,1)
histos['h_CSV_Merged_probbvsprobbb']=ROOT.TH2F("h_CSV_Merged_probbvsprobbb", "CSV prob b vs prob bb nJets = 1", 20,0,1,20,0,1)

histos['h_CSV_Split_Lead_probb']=ROOT.TH1F("h_CSV_Split_Lead_probb", "pfCSVJetTags:probb nJets>1 Leading", 20,0,1)
histos['h_CSV_Split_Lead_probbb']=ROOT.TH1F("h_CSV_Split_Lead_probbb", "pfCSVJetTags:probbb nJets>1 Leading", 20,0,1)
histos['h_CSV_Split_Lead_sum']=ROOT.TH1F("h_CSV_Split_Lead_sum", "Total CSV nJets>1 Leading", 20,0,1)
histos['h_CSV_Split_Lead_probbvsprobbb']=ROOT.TH2F("h_CSV_Split_Lead_probbvsprobbb", "CSV prob b vs prob bb nJets>1 Leading", 20,0,1,20,0,1)

histos['h_CSV_Split_Sublead_probb']=ROOT.TH1F("h_CSV_Split_Sublead_probb", "pfCSVJetTags:probb nJets>1 Subleading", 20,0,1)
histos['h_CSV_Split_Sublead_probbb']=ROOT.TH1F("h_CSV_Split_Sublead_probbb", "pfCSVJetTags:probbb nJets>1 Subleading", 20,0,1)
histos['h_CSV_Split_Sublead_sum']=ROOT.TH1F("h_CSV_Split_Sublead_sum", "Total CSV nJets>1 Subleading", 20,0,1)
histos['h_CSV_Split_Sublead_probbvsprobbb']=ROOT.TH2F("h_CSV_Split_Sublead_probbvsprobbb", "CSV prob b vs prob bb nJets>1 Subleading", 20,0,1,20,0,1)



histos['h_nEvent_tauHad_tauE']=ROOT.TH1F("h_nEvent_tauHad_tauE", "NEvent TauHad tauE",3,0,3)
histos['h_nEvent_tauHad_tauE'].SetCanExtend(ROOT.TH1.kAllAxes)

histos['h_tauHad_tauE_Mvis']=ROOT.TH1F("h_tauHad_tauE_Mvis", "tauHad||E Mvis",45,0,20)
histos['h_tauHad_tauE_bMass']=ROOT.TH1F("h_tauHad_tauE_bMass", "tauHad||E bMass",40,0,20)
histos['h_tauHad_tauE_bDRl']=ROOT.TH1F("h_tauHad_tauE_bDRl", "tauHad||E bDRtau",40,0,5)
histos['h_tauHad_tauE_lepDR']=ROOT.TH1F("h_tauHad_tauE_lepDR", "tauHad||E lepDR",40,0,4)
histos['h_tauHad_tauE_DiTauPt']=ROOT.TH1F("h_tauHad_tauE_DiTauPt", "tauHad||E DiTauPt",50,0,200)
histos['h_tauHad_tauE_MET']=ROOT.TH1F("h_tauHad_tauE_MET", "tauHad||E MET",40,0,250)


histos['h_nEvent_tauHad_tauMu']=ROOT.TH1F("h_nEvent_tauHad_tauMu", "NEvent TauHad tauMu",4,0,4)
histos['h_nEvent_tauHad_tauMu'].SetCanExtend(ROOT.TH1.kAllAxes)

histos['h_tauHad_tauMu_Mvis']=ROOT.TH1F("h_tauHad_tauMu_Mvis", "tauHad||Mu Mvis",45,0,20)
histos['h_tauHad_tauMu_bMass']=ROOT.TH1F("h_tauHad_tauMu_bMass", "tauHad||Mu bMass",40,0,20)
histos['h_tauHad_tauMu_bDRl']=ROOT.TH1F("h_tauHad_tauMu_bDRl", "tauHad||Mu bDRtau",40,0,5)
histos['h_tauHad_tauMu_lepDR']=ROOT.TH1F("h_tauHad_tauMu_lepDR", "tauHad||Mu lepDR",40,0,4)
histos['h_tauHad_tauMu_DiTauPt']=ROOT.TH1F("h_tauHad_tauMu_DiTauPt", "tauHad||Mu DiTauPt",50,0,200)
histos['h_tauHad_tauMu_MET']=ROOT.TH1F("h_tauHad_tauMu_MET", "tauHad||Mu MET",40,0,250)

fchain = ROOT.TChain('tcpNtuples/analysisTree')
for h in histos.keys():
    histos[h].Sumw2()

#gchain = ROOT.TChain('tcpGenNtuples/genTree')

pi = np.pi

h = {}





#prefix = "root://cmseos.fnal.gov//store/user/mwulansa/TCPNtuple/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/Ntuple_DYJetsToLL_M-50_Summer20UL17_v5-2/211117_035810/0000/"

#prefix = "root://cmseos.fnal.gov//store/user/mwulansa/TCPNtuple/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/Ntuple_DYJetsToLL_M-50_Summer20UL17_v5-2/211117_035810/0000/"


for inputFileName in inputFileList:
    inputFileName=inputFileName.replace("\n","")
    print(inputFileName.replace("\n",""))

    fchain.Add(inputFileName)
    fchain.AddFriend('tcpTrigNtuples/triggerTree', inputFileName)
    fchain.AddFriend('lumiSummary/lumiTree', inputFileName)

#    print(fchain.GetEntries())
#    gchain.Add(inputFileName)

jets = ROOT.JetInfoDS()
muons = ROOT.MuonInfoDS()
electrons = ROOT.ElectronInfoDS()
muCtaus = ROOT.TauInfoDS()
eCtaus = ROOT.TauInfoDS()


fchain.SetBranchAddress("Jets", ROOT.AddressOf(jets))
fchain.SetBranchAddress("TausMCleaned", ROOT.AddressOf(muCtaus))
fchain.SetBranchAddress("TausECleaned", ROOT.AddressOf(eCtaus))
fchain.SetBranchAddress("Muons", ROOT.AddressOf(muons))
fchain.SetBranchAddress("Electrons", ROOT.AddressOf(electrons))

for iev in range(fchain.GetEntries()): # Be careful!!!                                                                                                   
    fchain.GetEntry(iev)

    mets = fchain.GetBranch("Mets")
    met_pt = mets.GetLeaf('pt').GetValue()
    met_phi = mets.GetLeaf('phi').GetValue()

    weight = fchain.GetBranch("lumiInfo")
    genWeight = weight.GetLeaf('weight').GetValue()

    histos['h_nEvents'].Fill(0.5, 1)
    histos['h_nEvents'].Fill(1.5, genWeight)


    isMu = fchain.GetLeaf('isMu').GetValue()
    isIsoMu = fchain.GetLeaf('isIsoMu').GetValue()
    isIsoMuTau = fchain.GetLeaf('isIsoMuTau').GetValue()
    isIsoEle = fchain.GetLeaf('isIsoEle').GetValue()
    isEleTau = fchain.GetLeaf('isEleTau').GetValue()
    isMuonEG = fchain.GetLeaf('isMuonEG').GetValue()
    isDoubleMu = fchain.GetLeaf('isDoubleMu').GetValue()
    isDoubleEG = fchain.GetLeaf('isDoubleEG').GetValue()

    

    selected_muclean_taus = []
    selected_eclean_taus = []
    selected_electrons = []
    selected_bs = []
    selected_mus = []
    selected_b_mus = []

    if jets.size()>0:
        for i in range(jets.size()):
            jet = jets.at(i)
            if jet.pt>10 and abs(jet.eta)<2.5:
                histos['h_bJet_selection'].Fill("pt>10, |eta|<2.5",genWeight)
                if jet.id >= 1:
                    histos['h_bJet_selection'].Fill("Loose Jet ID",genWeight)
                    if (jet.flavorprobb+jet.flavorprobbb+jet.flavorproblepb) > 0.7476:
                        histos['h_bJet_selection'].Fill("total Deep Flavor >.7476",genWeight)
                        selected_bs+=[jet]

    selected_bs.sort(key=lambda x: x.pt, reverse=True)
    if len(selected_bs)==1:
        histos["h_deepFlavor_Merged_probb"].Fill(selected_bs[0].flavorprobb,genWeight)
        histos["h_deepFlavor_Merged_probbb"].Fill(selected_bs[0].flavorprobbb,genWeight)
        histos["h_deepFlavor_Merged_problepb"].Fill(selected_bs[0].flavorproblepb,genWeight)
        histos["h_deepFlavor_Merged_sum"].Fill(selected_bs[0].flavorprobb+selected_bs[0].flavorprobbb+selected_bs[0].flavorproblepb,genWeight)
        histos["h_deepFlavor_Merged_probbvsprobbb"].Fill(selected_bs[0].flavorprobb, selected_bs[0].flavorprobbb, genWeight)
    if len(selected_bs)>1:
        histos["h_deepFlavor_Split_Lead_probb"].Fill(selected_bs[0].flavorprobb,genWeight)
        histos["h_deepFlavor_Split_Lead_probbb"].Fill(selected_bs[0].flavorprobbb,genWeight)
        histos["h_deepFlavor_Split_Lead_problepb"].Fill(selected_bs[0].flavorproblepb,genWeight)
        histos["h_deepFlavor_Split_Lead_sum"].Fill(selected_bs[0].flavorprobb+selected_bs[0].flavorproblepb+selected_bs[0].flavorprobbb,genWeight)
        histos["h_deepFlavor_Split_Lead_probbvsprobbb"].Fill(selected_bs[0].flavorprobb, selected_bs[0].flavorprobbb, genWeight)
        histos["h_deepFlavor_Split_Sublead_probb"].Fill(selected_bs[1].flavorprobb,genWeight)
        histos["h_deepFlavor_Split_Sublead_probbb"].Fill(selected_bs[1].flavorprobbb,genWeight)
        histos["h_deepFlavor_Split_Sublead_problepb"].Fill(selected_bs[1].flavorproblepb,genWeight)
        histos["h_deepFlavor_Split_Sublead_sum"].Fill(selected_bs[1].flavorprobb+selected_bs[1].flavorproblepb+selected_bs[1].flavorprobbb,genWeight)
        histos["h_deepFlavor_Split_Sublead_probbvsprobbb"].Fill(selected_bs[1].flavorprobb, selected_bs[1].flavorprobbb, genWeight)

    if len(selected_bs)==1:
        histos["h_CSV_Merged_probb"].Fill(selected_bs[0].csvprobb,genWeight)
        histos["h_CSV_Merged_probbb"].Fill(selected_bs[0].csvprobbb,genWeight)
        histos["h_CSV_Merged_sum"].Fill(selected_bs[0].csvprobb+selected_bs[0].csvprobbb,genWeight)
        histos["h_CSV_Merged_probbvsprobbb"].Fill(selected_bs[0].csvprobb, selected_bs[0].csvprobbb, genWeight)
    if len(selected_bs)>1:
        histos["h_CSV_Split_Lead_probb"].Fill(selected_bs[0].csvprobb,genWeight)
        histos["h_CSV_Split_Lead_probbb"].Fill(selected_bs[0].csvprobbb,genWeight)
        histos["h_CSV_Split_Lead_sum"].Fill(selected_bs[0].csvprobb+selected_bs[0].csvprobbb,genWeight)
        histos["h_CSV_Split_Lead_probbvsprobbb"].Fill(selected_bs[0].csvprobb, selected_bs[0].csvprobbb, genWeight)
        histos["h_CSV_Split_Sublead_probb"].Fill(selected_bs[1].csvprobb,genWeight)
        histos["h_CSV_Split_Sublead_probbb"].Fill(selected_bs[1].csvprobbb,genWeight)
        histos["h_CSV_Split_Sublead_sum"].Fill(selected_bs[1].csvprobb+selected_bs[1].csvprobbb,genWeight)
        histos["h_CSV_Split_Sublead_probbvsprobbb"].Fill(selected_bs[1].csvprobb, selected_bs[1].csvprobbb, genWeight)
    histos['h_nbJet'].Fill(len(selected_bs),genWeight)
    if muons.size()>0:
        for i in range(muons.size()):
            muon = muons.at(i)
            histos['h_muon_selection'].Fill("Total Muons",genWeight)
            histos['h_muon_selection'].Fill("Total Muons",genWeight)
            if muon.pt>3 and abs(muon.eta)<2.4:
                histos['h_muon_selection'].Fill("pt>10, |eta|<2.5",genWeight)

                if muon.id >= 1:
                    histos['h_muon_selection'].Fill("Loose ID",genWeight)
                    if muon.iso<.25:
                        histos['h_muon_selection'].Fill("Iso<.25",genWeight)
                        selected_mus+=[muon]
                        isFromB=False
                        if (len(selected_bs)>0):
                            mu = ROOT.TLorentzVector()
                            b=ROOT.TLorentzVector()
                            mu.SetPtEtaPhiM(muon.pt, muon.eta, muon.phi, muon.mass)
                            b.SetPtEtaPhiM(selected_bs[0].pt, selected_bs[0].eta, selected_bs[0].phi, selected_bs[0].mass)
                            dr = b.DeltaR(mu)
                            if (dr<.4):
                                isFromB=True
                        if isFromB:
                            histos['h_muon_selection'].Fill("DR(Mu,b)<.4",genWeight)
                            selected_b_mus+=[muon]
                        else:
                            histos['h_muon_selection'].Fill("DR(Mu,b)>.4",genWeight)
                            selected_bs+=[muon]
                                
                            

    if electrons.size()>0:
        for i in range(electrons.size()):
            electron = electrons.at(i)
            histos['h_std_electron_selection'].Fill("Total Electrons",genWeight)
            if electron.pt>7 and abs(electron.eta)<2.4:
                histos['h_std_electron_selection'].Fill("pt>7, eta<2.4",genWeight)
                if electron.id >= 1 :
                    histos['h_std_electron_selection'].Fill("LooseID",genWeight)
                    selected_electrons+=[electron]
    if eCtaus.size()>0:
        for i in range(eCtaus.size()):
            tau=eCtaus.at(i)
            histos['h_eCtau_selection'].Fill("Total Taus",genWeight)
            if (tau.pt>20 and abs(tau.eta)<2.4):
                histos['h_eCtau_selection'].Fill("Pt>20, eta<2.4",genWeight)
                if (tau.mvaid>=3):
                    histos['h_eCtau_selection'].Fill("Loose MVA ID",genWeight)
                    JetSeperation=False
                    isLepFake=False

                    t=ROOT.TLorentzVector()
                    t.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass)
                    if len(selected_bs)>0:
                        b=ROOT.TLorentzVector()
                        b.SetPtEtaPhiM(selected_bs[0].pt,selected_bs[0].eta,selected_bs[0].phi,selected_bs[0].mass)
                        bdr=b.DeltaR(t)
                        if bdr>.5: JetSeperation=True
                    for electron in selected_electrons:
                        ele = ROOT.TLorentzVector()
                        ele.SetPtEtaPhiM(electron.pt, electron.eta, electron.phi, electron.mass)
                        edr = t.DeltaR(ele)
                        if(edr<.1):
                            isLepFake=True
                    for muon in selected_mus:
                        mu = ROOT.TLorentzVector()
                        mu.SetPtEtaPhiM(muon.pt, muon.eta, muon.phi, muon.mass)
                        edr = t.DeltaR(mu)
                        if(dr<.1):
                            eisLepFake=True
                    if JetSeperation:
                        histos['h_eCtau_selection'].Fill("DR(tau,Bjets)>.5",genWeight)
                        if not isLepFake:
                            histos['h_eCtau_selection'].Fill("DR(tau,Lep)>.1",genWeight)
                            selected_eclean_taus+=[tau]
    if muCtaus.size()>0:
        for i in range(muCtaus.size()):
            tau=muCtaus.at(i)
            histos['h_muCtau_selection'].Fill("Total Taus",genWeight)
            if (tau.pt>20 and abs(tau.eta)<2.4):
                histos['h_muCtau_selection'].Fill("Pt>20, eta<2.4",genWeight)
                if (tau.mvaid>=3):
                    histos['h_muCtau_selection'].Fill("Loose MVA ID",genWeight)
                    JetSeperation=False
                    isLepFake=False
                    t=ROOT.TLorentzVector()
                    t.SetPtEtaPhiM(tau.pt, tau.eta, tau.phi, tau.mass)
                    if len(selected_bs)>0:
                        b=ROOT.TLorentzVector()
                        b.SetPtEtaPhiM(selected_bs[0].pt,selected_bs[0].eta,selected_bs[0].phi,selected_bs[0].mass)

                        bdr=b.DeltaR(t)
                        if bdr>.5: JetSeperation=True
                    for electron in selected_electrons:
                        ele = ROOT.TLorentzVector()
                        ele.SetPtEtaPhiM(electron.pt, electron.eta, electron.phi, electron.mass)
                        edr = t.DeltaR(ele)

                        if(edr<.1):
                            isLepFake=True
                    for muon in selected_mus:
                        mu = ROOT.TLorentzVector()
                        mu.SetPtEtaPhiM(muon.pt, muon.eta, muon.phi, muon.mass)
                        mudr = t.DeltaR(mu)
                        if(mudr<.1):
                            isLepFake=True
                    if JetSeperation:
                        histos['h_muCtau_selection'].Fill("DR(tau,All Bjets)>.5",genWeight)
                        if not isLepFake:
                            histos['h_muCtau_selection'].Fill("DR(tau,Lep)>.1",genWeight)
                            selected_muclean_taus+=[tau]

    selected_bs.sort(key=lambda x: x.pt, reverse=True)
    selected_electrons.sort(key=lambda x: x.pt, reverse=True)
    selected_mus.sort(key=lambda x: x.pt, reverse=True)
    selected_muclean_taus.sort(key=lambda x: x.pt, reverse=True)
    selected_eclean_taus.sort(key=lambda x: x.pt, reverse=True)

#       if isHT == 1 or isSingleJet == 1 or isMu == 1 or \\
#          isIsoMu == 1 or isIsoMuTau == 1 or \\
#          isIsoEle == 1 or isEleTau == 1 or \\
#          isMuonEG == 1: 
    histos['h_nEvent_tauHad_tauE'].Fill("TotalEvents",genWeight)
    if len(selected_electrons)>0:
        histos['h_nEvent_tauHad_tauE'].Fill("NEle>0",genWeight)
        if len(selected_eclean_taus)>0:
            histos['h_nEvent_tauHad_tauE'].Fill("NEleCleanTau>0",genWeight)
            if len(selected_bs)>0:
                histos['h_nEvent_tauHad_tauE'].Fill("Nb>0",genWeight)
                b=ROOT.TLorentzVector()
                #if (len(selected_bs)>1):

                #    b1=ROOT.TLorentzVector()
                #    b1.SetPtEtaPhiM(selected_bs[0].pt,selected_bs[0].eta,selected_bs[0].phi,selected_bs[0].mass)
                #    b2=ROOT.TLorentzVector()
                #    b2.SetPtEtaPhiM(selected_bs[1].pt,selected_bs[1].eta,selected_bs[1].phi,selected_bs[1].mass)
                #    b = b1+b2
                #else:
                #    b.SetPtEtaPhiM(selected_bs[0].pt,selected_bs[0].eta,selected_bs[0].phi,selected_bs[0].mass)
                b.SetPtEtaPhiM(selected_bs[0].pt,selected_bs[0].eta,selected_bs[0].phi,selected_bs[0].mass)
                e=ROOT.TLorentzVector()
                e.SetPtEtaPhiM(selected_electrons[0].pt,selected_electrons[0].eta,selected_electrons[0].phi,selected_electrons[0].mass)
                t=ROOT.TLorentzVector()
                t.SetPtEtaPhiM(selected_eclean_taus[0].pt,selected_eclean_taus[0].eta,selected_eclean_taus[0].phi,selected_eclean_taus[0].mass)
                mvis = (t+e).M()
                DiTauPt=(t+e).Pt()
                bMass=b.M()
                bDRl=b.DeltaR(t)
                lepDR=e.DeltaR(t)
                if isEleTau:
                    histos['h_nEvent_tauHad_tauE'].Fill("ETau Trigger*",genWeight)
                if isMuonEG:
                    histos['h_nEvent_tauHad_tauE'].Fill("MuEG Trigger*",genWeight)
                if isIsoMu:
                    histos['h_nEvent_tauHad_tauE'].Fill("IsoMuTrigger*",genWeight)
                if isIsoEle:
                    histos['h_nEvent_tauHad_tauE'].Fill("IsoEle Trigger*",genWeight)
                    histos['h_tauHad_tauE_Mvis'].Fill(mvis,genWeight)
                    histos['h_tauHad_tauE_bMass'].Fill(bMass,genWeight)
                    histos['h_tauHad_tauE_bDRl'].Fill(bDRl,genWeight)
                    histos['h_tauHad_tauE_lepDR'].Fill(lepDR,genWeight)
                    histos['h_tauHad_tauE_DiTauPt'].Fill(DiTauPt,genWeight)
                    histos['h_tauHad_tauE_MET'].Fill(met_pt,genWeight)
    histos['h_nEvent_tauHad_tauMu'].Fill("TotalEvents",genWeight)
    if len(selected_mus)>0:
        histos['h_nEvent_tauHad_tauMu'].Fill("nMu>0",genWeight)
        if len(selected_muclean_taus)>0:
            histos['h_nEvent_tauHad_tauMu'].Fill("NMuCleanTau>0",genWeight)
            if len(selected_bs)>0:
                histos['h_nEvent_tauHad_tauMu'].Fill("N b>0",genWeight)
                b=ROOT.TLorentzVector()
                #if (len(selected_bs)>1):

                #    b1=ROOT.TLorentzVector()
                #    b1.SetPtEtaPhiM(selected_bs[0].pt,selected_bs[0].eta,selected_bs[0].phi,selected_bs[0].mass)
                #    b2=ROOT.TLorentzVector()
                #    b2.SetPtEtaPhiM(selected_bs[1].pt,selected_bs[1].eta,selected_bs[1].phi,selected_bs[1].mass)
                #    b = b1+b2
                #else:
                #    b.SetPtEtaPhiM(selected_bs[0].pt,selected_bs[0].eta,selected_bs[0].phi,selected_bs[0].mass)
                b.SetPtEtaPhiM(selected_bs[0].pt,selected_bs[0].eta,selected_bs[0].phi,selected_bs[0].mass)
                mu=ROOT.TLorentzVector()
                mu.SetPtEtaPhiM(selected_mus[0].pt,selected_mus[0].eta,selected_mus[0].phi,selected_mus[0].mass)
                t=ROOT.TLorentzVector()
                t.SetPtEtaPhiM(selected_muclean_taus[0].pt,selected_muclean_taus[0].eta,selected_muclean_taus[0].phi,selected_muclean_taus[0].mass)
                mvis = (t+mu).M()
                DiTauPt=(t+mu).Pt()
                bMass=b.M()
                bDRl=b.DeltaR(t)
                lepDR=mu.DeltaR(t)
                if isIsoMuTau:
                    histos['h_nEvent_tauHad_tauMu'].Fill("MuTau TRigger*",genWeight)
                if isIsoMu:
                    histos['h_nEvent_tauHad_tauMu'].Fill("IsoMu Trigger**",genWeight)
                    histos['h_tauHad_tauMu_Mvis'].Fill(mvis,genWeight)
                    histos['h_tauHad_tauMu_bMass'].Fill(bMass,genWeight)
                    histos['h_tauHad_tauMu_bDRl'].Fill(bDRl,genWeight)
                    histos['h_tauHad_tauMu_lepDR'].Fill(lepDR,genWeight)
                    histos['h_tauHad_tauMu_DiTauPt'].Fill(DiTauPt,genWeight)
                    histos['h_tauHad_tauMu_MET'].Fill(met_pt,genWeight)                


out.cd()

for key in histos.keys():
    histos[key].Write()

out.Close()
