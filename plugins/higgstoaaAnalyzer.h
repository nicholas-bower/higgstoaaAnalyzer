using namespace std;
namespace h2AA{
    template<class T>
    T getObject(const edm::EDGetTokenT<T> token_, const edm::Event& event){
        edm::Handle<T> handle;
        event.getByToken(token_, handle);
        T physicsObject = *handle.product();
        return physicsObject;
    }
    
    bool checkID(reco::GsfElectron e, int iD, double rho, float Esc) {
        float hoverECut=0.0;
        float SigmaIeIeCut=0;
        float EinvPinvCut = 0;
        unsigned int missingHitsCut=0;
        float dEtaInSeedCut=0;
        float dPhiCut=0;
        float ePt = e.pt();
        bool result=false;
        float dEtaInSeed = e.deltaEtaSuperClusterTrackAtVtx()-e.superCluster()->eta()+e.superCluster()->seed()->eta();
        float GsfEleEInverseMinusPInverse = abs((1.0-e.eSuperClusterOverP())/e.ecalEnergy());
        constexpr reco::HitPattern::HitCategory missingHitType = reco::HitPattern::MISSING_INNER_HITS;
        const unsigned int mhits = e.gsfTrack()->hitPattern().numberOfAllHits(missingHitType);
        if (e.isEB()){
            if (iD==0){
            hoverECut = .05+1.16/Esc+.0324*rho/Esc;
            SigmaIeIeCut =.0126;
            dEtaInSeedCut = .00463;
            dPhiCut = .148;
            EinvPinvCut = .209;
            missingHitsCut = 2;
            }
            else if (iD==1){
            hoverECut = .05+1.16/Esc+.0324*rho/Esc;
            SigmaIeIeCut =.0112;
            dEtaInSeedCut = .00377;
            dPhiCut = .0884;
            EinvPinvCut = .193;
            missingHitsCut = 1;
            }
            else if (iD==3){
            hoverECut = .05+1.16/Esc+.0324*rho/Esc;
            SigmaIeIeCut =.0104;
            dEtaInSeedCut = .00255;
            dPhiCut = .022;
            EinvPinvCut = .159;
            missingHitsCut = 1;
            }
            else if (iD==2){
            hoverECut = .05+1.16/Esc+.0324*rho/Esc;
            SigmaIeIeCut =.0106;
            dEtaInSeedCut = .0032;
            dPhiCut = .0547;
            EinvPinvCut = .184;
            missingHitsCut = 1;
            } 
        }
        if (e.isEE()){
            if (iD==0){
            hoverECut = .05+2.54/Esc+.183*rho/Esc;
            SigmaIeIeCut =.0457;
            dEtaInSeedCut = .00814;
            dPhiCut = .19;
            EinvPinvCut = .132;
            missingHitsCut = 3;
            }
            else if (iD==1){
            hoverECut = .05+2.54/Esc+.183*rho/Esc;
            SigmaIeIeCut =.0425;
            dEtaInSeedCut = .00674;
            dPhiCut = .169;
            EinvPinvCut = .111;
            missingHitsCut = 1;
            }
            else if (iD==2){
            hoverECut = .05+2.54/Esc+.183*rho/Esc;
            SigmaIeIeCut =.0387;
            dEtaInSeedCut = .00632;
            dPhiCut = .0394;
            EinvPinvCut = .0721;
            missingHitsCut = 1;
            }
            else if (iD==3){
            hoverECut = .05+2.54/Esc+.183*rho/Esc;
            SigmaIeIeCut =.0353;
            dEtaInSeedCut = .00501;
            dPhiCut = .0236;
            EinvPinvCut = .0197;
            missingHitsCut = 1;
            } 
        }
        if(e.full5x5_sigmaIetaIeta()<SigmaIeIeCut &&
            e.hadronicOverEm()<hoverECut &&
            dEtaInSeed<dEtaInSeedCut &&
            e.deltaPhiSuperClusterTrackAtVtx()<dPhiCut &&
            GsfEleEInverseMinusPInverse<EinvPinvCut&&
            mhits <= missingHitsCut&& ePt>7
            ){
            result =  true;
        }
        else{result =  false;}
        return result;
    }

    float muonIsolation(reco::Muon m){
        float  iso = (m.pfIsolationR04().sumPhotonEt+m.pfIsolationR04().sumNeutralHadronEt-0.5*m.pfIsolationR04().sumPUPt)/m.pt();
        if  (iso<0){iso=0;}
        iso = iso+ m.pfIsolationR04().sumChargedHadronPt/m.pt();
        return iso;
    }
    template<class T>
    bool sortByPt(T i, T j){return i.pt()> j.pt();}

}