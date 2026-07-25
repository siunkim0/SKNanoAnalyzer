#ifndef chig_h
#define chig_h

#include "AnalyzerCore.h"

//==============================================================
// charged Higgs -> tb search (heavy H+, pp > H+ > tb s-channel)
//
// signal: pp -> H+ -> t b,  t -> W+ b  (+ charge conjugate)
//   최종 상태 = b-jet 2개 + W,  m(H+) = 200 GeV ~ 3 TeV
//
// 두 채널:
//   Lep : W -> l nu (e/mu). 1 tight lepton + MET + >=2 jets
//         nu pz 는 mW constraint 로 복원 (판별식 < 0 이면 실수부)
//   Had : W -> qq'. lepton veto + >=4 jets,
//         W = b-candidate 를 제외한 jet pair 중 m(jj) 가 80.4 에
//         가장 가까운 것 (mass window)
//
// background modeling (arXiv:2605.02848, HIG-20-012 §6 방식):
//   b-candidate = DeepJet score 상위 2개 jet. 이걸로 dataset 분류:
//     2b   : 둘 다 Medium 통과            -> signal dataset
//     1b1L : 하나 Medium, 하나 Loose-only  -> background model dataset
//   (score 정렬이라 분류가 모호하지 않다: s1 >= s2)
//   analysis region 은 |m_top - 172.5| 로: SR < 25 < VR < 45 < CR < 90
//   -> CR 에서 GBReweighter 로 1b1L -> 2b 학습, SR(1b1L) reweight 가
//      SR(2b) 의 background m(tb) template
//
// 공통 재구성:
//   top = W + b (2개 b-candidate 중 m(Wb) 이 172.5 에 가까운 쪽)
//   H+  = top + (남은 b-candidate)  ->  m(tb) 가 signal peak
//
// 출력: 히스토그램 (<ch>/<cat>/<region>/...) + flat TTree "tree"
//       (BDT reweighting 학습용 per-event 변수)
//==============================================================
class chig : public AnalyzerCore {
public:
    chig();
    ~chig();

    void initializeAnalyzer();
    void executeEvent();

    // Trigger (era 별 unprescaled)
    RVec<TString> muTriggers;
    RVec<TString> elTriggers;
    RVec<TString> hadTriggers;

    // Analysis cuts
    struct AnalysisCuts {
        // tight lepton (trigger-safe pT 는 era 별로 initializeAnalyzer 에서 설정)
        float mu_pt_min      = 26.0;
        float mu_eta_max     =  2.4;
        float el_pt_min      = 35.0;
        float el_eta_max     =  2.5;
        // loose lepton (veto 용)
        float veto_lep_pt_min = 10.0;
        // jet
        float jet_pt_min  = 30.0;
        float jet_eta_max =  2.4;   // b-tagging acceptance
        float jet_lep_dr  =  0.4;
        // event selection
        float met_min     = 30.0;   // Lep channel
        int   njet_lep_min = 2;
        int   njet_had_min = 4;
        // Had: PFHT1050 plateau 위에서만 (trigger eff ~1, 그리고 Delphes 신호의
        // 가짜 HLT bit 문제도 이 offline cut 이 해결한다)
        float had_ht_min  = 1100.0;
        // mass reconstruction
        float w_mass      = 80.4;
        float top_mass    = 172.5;
        float w_window    = 20.0;   // Had: |m(jj) - 80.4| < 20
        // analysis regions: |m_top - 172.5| 경계
        float sr_dm       = 25.0;
        float vr_dm       = 45.0;
        float cr_dm       = 90.0;
    } cuts;

    // Delphes-mimic 신호 샘플 (private MadGraph 생산) 여부:
    // Pileup_nTrueInt=0, L1PreFiringWeight/Jet_puId 브랜치 없음 →
    // PU/prefire weight 와 jet PU ID 컷을 건너뛴다
    bool isDelphes = false;

    // ID / tagger 설정
    TString muonTightID, muonVetoID;
    TString electronTightID, electronVetoID;
    JetTagging::JetFlavTagger bTagger = JetTagging::JetFlavTagger::DeepJet;
    float bTagWPMedium = -1.;
    float bTagWPLoose = -1.;

    // 이벤트마다 executeEvent 에서 새로 채우는 physics objects
    RVec<Muon> tightMuons, vetoMuons;
    RVec<Electron> tightElectrons, vetoElectrons;
    RVec<Jet> jets;      // 전체 selected jets (pT 순)
    RVec<Jet> bCands;    // DeepJet score 상위 2개 (b-candidates, score 순)
    RVec<Jet> otherJets; // 나머지 (Had 채널 W 후보용, pT 순)
    Particle METv;

    // Helper functions
    // b-candidate 분류: "2b" / "1b1L" / "" (reject). bCands/otherJets 를 채운다
    TString ClassifyBCands();
    // |m_top - 172.5| 로 region 결정: "SR" / "VR" / "CR" / "" (reject)
    TString GetRegion(float mtop) const;
    // mW constraint 로 nu pz 후보 반환 (실수해 2개 or 실수부 1개)
    RVec<float> SolveNeutrinoPz(const Particle &lep, const Particle &met) const;
    // 히스토그램 + BDT 학습용 tree 를 함께 채운다
    //   channel: 0=mu, 1=el, 2=had / cat: 2="2b", 1="1b1L" / region: 0=SR,1=VR,2=CR
    void fillOutputs(const TString &prefix, int channel, int cat, int region,
                     const Particle &W, const Particle &top, const Particle &Hc,
                     const Jet &bTop, const Jet &bHc,
                     float lepPt, float mT, float weight);
};

#endif
