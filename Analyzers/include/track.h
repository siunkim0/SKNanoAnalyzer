#ifndef track_h
#define track_h

#include "AnalyzerCore.h"

//==============================================================
// Jet 안의 track/constituent 개수 & quark vs gluon 스터디 (v6)
//
//   - Tight jet ID, pT > 30 GeV, |eta| < 2.4 인 jet 선택
//   - HLT_IsoMu24 + MET>20 (W->munu enriched, muon 채널)
//   - jet 마다: nConstituents(total) + chMultiplicity/neMultiplicity
//     (charged/neutral, Run3 v12+ 전용 → 2022+2022EE 에서 사용)
//   - MC는 partonFlavour truth로 u/d/s/c/b/gluon 6분할 (v6, 이전 v4 = uds 통합)
//   - NoSel/: offline muon impact parameter (dXY/dZ/IP3D/SIP3D) — cut 전혀 없음,
//     QCD MuEnriched 스터디용
//   - pT축 TeV까지 확장 → ⟨N⟩ vs pT log-law 외삽 + flavour 비교
//==============================================================
class track : public AnalyzerCore {
public:
    track();
    ~track();

    void initializeAnalyzer();
    void executeEvent();

    // Analysis cuts
    struct AnalysisCuts {
        float jet_pt_min  = 30.0;
        float jet_eta_max =  2.4;
        float met_pt_min  = 20.0;   // W->munu enriched
    } cuts;

    // era 별 single muon trigger
    TString IsoMuTriggerName;
    float   TriggerSafePtCut;
};

#endif
