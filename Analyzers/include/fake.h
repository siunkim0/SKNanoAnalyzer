#ifndef fake_h
#define fake_h

#include "AnalyzerCore.h"

//==============================================================
// Muon fake rate 측정 (Run 2 2016, prescaled single muon trigger)
//
// 방법: QCD-enriched 영역에서 loose muon 중 tight 를 통과하는 비율을 잰다
//   FR(ptcorr, |eta|) = N(tight) / N(loose)
//
// Userflags:
//   MeasFakeMu8  - Mu8 트리거 경로만 측정
//   MeasFakeMu17 - Mu17 트리거 경로만 측정
//   (플래그 없음) - 두 경로 모두 측정
//   away jet pT 변화 (30 / 60 GeV) 시스테마틱은 항상 채운다
//   MakeTree     - flavor BDT 학습용 flat ntuple ("mu" 트리) 를 추가로 쓴다
//==============================================================
class fake : public AnalyzerCore {
public:
    fake();
    ~fake();

    void initializeAnalyzer();
    void executeEvent();

    // rootcling 이 vector<fake::TriggerPath> dictionary 를 만들기 때문에 public 필요
    struct TriggerPath {
        TString name;      // 히스토그램 prefix: Mu8, Mu17
        TString trigger;   // HLT 경로 이름
        float ptCut;       // offline muon pT 컷
        float ptCorrCut;   // cone-corrected pT 영역 컷
    };

    // Userflags
    bool MeasFakeMu8, MeasFakeMu17;
    bool MakeTree;   // flavor BDT 학습용 flat ntuple 출력

    // 트리거 경로 / 시스테마틱 / binning
    vector<TriggerPath> paths;
    //==== electron veto 의 raw MVANoIso 컷: etaRegion (IB, OB, EC) 별.
    //==== Run 2 / Run 3 에서 MVA 가 retraining 되어 값이 다르다
    //==== (elecfake.h / closnff.h 와 같은 이름, 같은 값을 쓴다)
    float looseMvaIB, looseMvaOB, looseMvaEC;

    RVec<TString> systs;                 // Central + AwayJetPt30 + AwayJetPt60
    RVec<float> ptCorrBins, absEtaBins;  // FR 2D 히스토그램 binning

    // Analysis cuts (analysis note 의 WP 표 + 측정영역 정의)
    struct AnalysisCuts {
        // muon 공통: POG medium ID + 아래 컷
        float muon_pt_min    = 10.0;
        float muon_eta_max   =  2.4;
        float muon_dz_max    =  0.1;   // cm
        float muon_tkiso_max =  0.4;   // rel. tracker iso R03 (trigger emulation)
        // tight WP
        float tight_sip3d_max   = 3.0;
        float tight_miniiso_max = 0.1;
        // loose WP (Run 2)
        // loose WP: **Run 별로 다르다** (공식 MeasFakeRateV4 는
        //   Run2 -> Muon::Pass_HcToWALooseRun2 (SIP3D<5, miniiso<0.6)
        //   Run3 -> Muon::Pass_HcToWALooseRun3 (SIP3D<8, miniiso<0.4)
        // 를 쓴다). 아래는 Run2 기본값이고 Run3 값은 initializeAnalyzer 에서
        // 덮어쓴다. tight WP 은 Run 무관하게 같다.
        float loose_sip3d_max   = 5.0;
        float loose_miniiso_max = 0.6;
        // electron veto (공식 측정의 veto electron 과 동일: pT > 10)
        // loose_el_sip3d_max 는 Run 별로 다르다 (Run2 8, Run3 6) — initializeAnalyzer
        float loose_el_sip3d_max   = 8.0;
        float loose_el_miniiso_max = 0.4;   // Run 무관
        float electron_pt_min  = 10.0;
        float electron_eta_max =  2.5;
        // jet
        float jet_pt_min   = 25.0;
        float jet_eta_max  =  2.4;   // 2016: 2.4, 그 외: 2.5 (initializeAnalyzer 에서 설정)
        float jet_lep_dr   =  0.4;   // lepton 과 겹치는 jet 제거
        float awayjet_pt   = 40.0;   // Central away jet pT 컷
        float awayjet_dr   =  0.7;   // dR(mu, away jet)
        // 측정 영역 (W/Z prompt 오염 억제)
        float met_max = 25.0;
        float mt_max  = 25.0;
        // Z-enriched 영역 (MC normalization 용)
        float z_mass    = 91.2;
        float z_window  = 15.0;
        float zjet_pt   = 40.0;
    } cuts;

    // 이벤트마다 executeEvent 에서 새로 채우는 physics objects
    RVec<Muon> looseMuons, tightMuons;   // tight ⊂ loose
    RVec<Electron> looseElectrons;
    RVec<Jet> rawJets;                   // lepton cleaning 이전 (flavor 매칭용)
    RVec<Jet> jets;
    Particle METv;
    RVec<Gen> gens;                      // prompt / fake 구분용 (MC only)
    RVec<GenJet> genJets;                // pileup jet ID SF 용 (MC only)

    // Helper functions
    bool PassMuonWP(const Muon &mu, const TString &wp) const;
    bool PassElectronVeto(const Electron &el) const;
    bool PassVetoMapJet(const Jet &jet, const RVec<Muon> &muons) const;
    void measureFakeRate(const TriggerPath &path, float weight);
    void fillZEnriched(const TriggerPath &path, const Event &ev, float weight);
    void fillMuonKinematics(const TString &prefix, const Muon &mu,
                            float ptCorr, float MT, int nJets, float weight);
    TString LeptonTypeToString(int leptonType) const;

    // fake muon 의 source jet parton flavor 분류 (MC only, flavor-weighted FR 용)
    //   parton/hadron flavour -> {b,c,s,d,u,g,pileup,unmatched} 문자열
    //   isPileup: 대응 gen jet 이 없는 reco jet (genJetIdx<0) -> flavour 보다 우선
    TString FlavorTag(int partonFlavour, int hadronFlavour, bool isPileup) const;
    //   reco: muon 의 NanoAOD jetIdx 가 가리키는 jet, gen: 최근접 gen jet (dR<0.4)
    //   주의: 측정 jet 컬렉션(jets)은 lepton cleaning 으로 source jet 이 제거되므로
    //         반드시 rawJets / genJets (cleaning 이전) 에서 찾는다
    TString RecoJetFlavor(const Muon &mu) const;
    TString GenJetFlavor(const Muon &mu) const;

    //==== flavor BDT 학습용 flat ntuple (userflag MakeTree)
    //   FlavorTag 문자열 -> 정수 코드 (SetBranch 는 문자열을 지원하지 않는다)
    int FlavorCode(const TString &flav) const;
    //   loose muon 1개짜리 측정 영역 이벤트마다 정확히 한 줄.
    //   trigger path loop 밖에서 부르므로 Mu8/Mu17 중복 기록이 없다
    //   (경로 통과 여부는 pass_mu8 / pass_mu17 컬럼으로 남긴다).
    //   ev 는 GetTriggerLumi 가 non-const 라 const ref 로 받을 수 없다
    void fillTreeRow(Event &ev, float evtSF);
};

#endif
