#ifndef elecfake_h
#define elecfake_h

#include "AnalyzerCore.h"

//==============================================================
// Electron fake rate 측정 (Run 2 2016, prescaled single electron trigger)
//
// 방법: QCD-enriched 영역에서 loose electron 중 tight 를 통과하는 비율을 잰다
//   FR(ptcorr, |scEta|) = N(tight) / N(loose)
//
// Userflags:
//   MeasFakeEle8  - Ele8 트리거 경로만 측정
//   MeasFakeEle12 - Ele12 트리거 경로만 측정
//   MeasFakeEle23 - Ele23 트리거 경로만 측정
//   (플래그 없음)  - 세 경로 모두 측정
//   away jet pT 변화 (30 / 60 GeV) 시스테마틱은 항상 채운다
//==============================================================
class elecfake : public AnalyzerCore {
public:
    elecfake();
    ~elecfake();

    void initializeAnalyzer();
    void executeEvent();

    // rootcling 이 vector<elecfake::TriggerPath> dictionary 를 만들기 때문에 public 필요
    struct TriggerPath {
        TString name;      // 히스토그램 prefix: Ele8, Ele12, Ele23
        TString trigger;   // HLT 경로 이름
        float ptCut;       // offline electron pT 컷
        float ptCorrCut;   // cone-corrected pT 영역 컷
    };

    // Userflags
    bool MeasFakeEle8, MeasFakeEle12, MeasFakeEle23;

    // 트리거 경로 / 시스테마틱 / binning
    vector<TriggerPath> paths;
    RVec<TString> systs;                 // Central + AwayJetPt30 + AwayJetPt60
    RVec<float> ptCorrBins, absEtaBins;  // FR 2D 히스토그램 binning

    // loose WP 의 raw MVANoIso 컷 (etaRegion IB / OB / EC, Run 별로 다름)
    float looseMvaIB, looseMvaOB, looseMvaEC;

    // Analysis cuts (electron ID 표 + 측정영역 정의)
    struct AnalysisCuts {
        // electron 공통 (baseline): trigger emulation + 아래 컷
        float electron_pt_min       = 15.0;   // 측정 대상 electron pT 컷
        float electron_veto_pt_min  = 10.0;   // 공식 측정의 veto electron pT 컷
        float electron_eta_max      =  2.5;   // object 선택은 |eta| (공식 SelectElectrons 와 동일)
        float electron_dz_max       =  0.1;   // cm
        int   electron_losthits_max =  1;     // expected missing inner hits
        // tight WP (MVANoIso POG wp90 는 PassElectronWP 에서)
        float tight_sip3d_max   = 4.0;
        float tight_miniiso_max = 0.1;
        // loose WP
        float loose_sip3d_max   = 8.0;
        float loose_miniiso_max = 0.4;
        // muon veto (fake.cc 의 loose muon WP 와 동일)
        float muon_pt_min      = 10.0;
        float muon_eta_max     =  2.4;
        float muon_dz_max      =  0.1;
        float muon_tkiso_max   =  0.4;
        float muon_sip3d_max   =  5.0;
        float muon_miniiso_max =  0.6;
        // jet
        float jet_pt_min   = 25.0;
        float jet_eta_max  =  2.4;   // 2016: 2.4, 그 외: 2.5 (initializeAnalyzer 에서 설정)
        float jet_lep_dr   =  0.4;   // lepton 과 겹치는 jet 제거
        float awayjet_pt   = 40.0;   // Central away jet pT 컷
        float awayjet_dr   =  0.7;   // dR(el, away jet)
        // 측정 영역 (W/Z prompt 오염 억제)
        float met_max = 25.0;
        float mt_max  = 25.0;
        // Z-enriched 영역 (MC normalization 용)
        float z_mass    = 91.2;
        float z_window  = 15.0;
        float zjet_pt   = 40.0;
    } cuts;

    // 이벤트마다 executeEvent 에서 새로 채우는 physics objects
    RVec<Electron> vetoElectrons;                    // pT > 10 (event veto + jet cleaning)
    RVec<Electron> looseElectrons, tightElectrons;   // 측정 대상 pT > 15, tight ⊂ loose ⊂ veto
    RVec<Muon> vetoMuons;
    RVec<Jet> rawJets;                   // lepton cleaning 이전 (flavor 매칭용)
    RVec<Jet> jets;
    Particle METv;
    RVec<Gen> gens;                      // prompt / fake 구분용 (MC only)
    RVec<GenJet> genJets;                // pileup jet ID SF 용 (MC only)

    // Helper functions
    bool PassElectronWP(const Electron &el, const TString &wp) const;
    bool PassMuonVeto(const Muon &mu) const;
    bool PassVetoMapJet(const Jet &jet, const RVec<Muon> &muons) const;
    void measureFakeRate(const TriggerPath &path, float weight);
    void fillZEnriched(const TriggerPath &path, const Event &ev, float weight);
    void fillElectronKinematics(const TString &prefix, const Electron &el,
                                float ptCorr, float MT, int nJets, float weight);
    TString LeptonTypeToString(int leptonType) const;

    // fake electron 의 source jet parton flavor 분류 (MC only, flavor-weighted FR 용)
    //   parton/hadron flavour -> {b,c,s,d,u,g,unmatched} 문자열 (fake.cc 와 동일)
    TString FlavorTag(int partonFlavour, int hadronFlavour) const;
    //   electron 을 가장 가까운 (raw reco / gen) jet 에 dR<0.4 로 매칭해 flavor 반환
    //   주의: 측정 jet 컬렉션(jets)은 lepton cleaning 으로 source jet 이 제거되므로
    //         반드시 rawJets / genJets (cleaning 이전) 에 매칭한다
    TString RecoJetFlavor(const Electron &el) const;
    TString GenJetFlavor(const Electron &el) const;
};

#endif
