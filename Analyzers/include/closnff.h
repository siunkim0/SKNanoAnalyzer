#ifndef closnff_h
#define closnff_h

#include "AnalyzerCore.h"

//==============================================================
// Fake rate closure test, AN-25-154 Fig. 50 의 3Mu + 1E2Mu 채널 선택
//
// aa/SKNanoAnalyzer 의 ClosFakeRate.cc (choij1589, B2G-25-013 branch) 를
// 그대로 옮기되, fake rate 를 analyzer 안에서 correctionlib 으로 계산하지 않고
// **fail-side muon 의 feature 를 flat tree 로 뱉는다**. neural fake factor 는
// 2D 히스토그램으로 표현할 수 없어 per-muon feature vector 가 필요하기 때문이다.
// binned map 과 neural FR 을 같은 이벤트에 offline 에서 나란히 적용한다.
//
// 선택 (ClosFakeRate::selectEvent 와 동일). 두 채널을 한 번에 돌린다:
// lepton multiplicity 가 서로 배타적이라 한 이벤트가 양쪽에 들어갈 수 없으므로
// job 을 두 번 돌릴 이유가 없다. 어느 쪽인지는 channel 브랜치로 남긴다.
//
//   channel = 0 (3Mu):
//     loose muon 3개 == veto muon 3개, veto electron 0개
//     DoubleMuon trigger, pT > 20 / 10 / 10, |eta| < 2.4
//     |q1+q2+q3| == 1, OS pair 두 개 모두 M > 12 GeV
//
//   channel = 1 (1E2Mu):
//     loose muon 2개 == veto muon 2개, loose electron 1개 == veto electron 1개
//     EMu trigger, (mu1 > 25 && e > 15) || (mu1 > 10 && e > 25)  <- raw pT 로 자른다
//     q(mu1) + q(mu2) == 0, M(mu mu) > 12 GeV
//     electron 은 pT > 15 (veto electron 만 pT > 10)
//
//   공통: jet >= 2, b-jet (DeepJet Medium) >= 1
//         MET 컷 없음, MT 컷 없음, Z veto 없음
//
//   SR (region=0) : loose lepton 이 모두 tight       -> Observed
//   SB (region=1) : 하나 이상이 loose-not-tight      -> Expected 의 재료
//   fail 은 muon 이든 electron 이든 똑같이 세고 부호도 같다
//   (TriLeptonBase::GetFakeWeight 는 두 flavour 를 한 곱셈에 넣는다).
//
// Observed / Expected 는 offline 에서 만든다:
//   Observed(X) = sum_{SR}  w                        1[X in bin]
//   Expected(X) = sum_{SB}  w (-1)^(N_fail+1) prod_fail f/(1-f)
// 두 번째 줄의 부호 규칙은 TriLeptonBase::GetFakeWeight 와 같다.
//
// **WP 은 FR 을 측정한 analyzer 와 반드시 같아야 한다.**
//   muon     : fake.cc::PassMuonWP
//   electron : elecfake.cc::PassElectronWP
// 둘 다 자체 구현이다 — 우리 fork 의 HcToWA ID 에는 muon TkRelIso 컷이 no-op
// 이 되는 버그와 electron raw-MVA 컷의 부호가 뒤집힌 버그가 있다.
// 측정과 적용의 WP 이 어긋나면 closure 자체가 의미를 잃는다.
//
// (이전 버전의 PassElectronVeto 는 Electron::Pass_HcToWABaseline + SIP3D/miniiso
//  + raw-MVA 로 짜여 있었는데, 그 baseline 이 elecfake 의 공통 컷과 글자 그대로
//  같아서 loose WP 과 동치였다. 그래서 3Mu 의 electron veto 는 이 교체로
//  달라지지 않는다 — 3Mu 결과가 움직이면 그건 버그다.)
//
// Userflags: 없음 (tree 는 항상 쓴다)
// 출력 tree "clos": 선택된 이벤트당 정확히 한 줄
//==============================================================
class closnff : public AnalyzerCore {
public:
    closnff();
    ~closnff();

    void initializeAnalyzer();
    void executeEvent();

    //==== era 별 HLT 경로 (TriLeptonBase.cc 와 동일)
    RVec<TString> DblMuTriggers;   // 3Mu 채널
    RVec<TString> EMuTriggers;     // 1E2Mu 채널

    //==== loose electron WP 의 raw MVANoIso 컷 (elecfake.cc 와 같은 값)
    float looseMvaIB, looseMvaOB, looseMvaEC;

    //==== Analysis cuts
    //   muon WP 값은 fake.h::AnalysisCuts 와 같은 숫자여야 한다 (위 주석 참고).
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
        // electron (elecfake.cc::AnalysisCuts 와 같은 숫자)
        float electron_veto_pt_min = 10.0;   // veto / jet cleaning
        float electron_pt_min      = 15.0;   // 측정 대상 (1E2Mu 의 fake 후보)
        float electron_eta_max     =  2.5;
        float electron_dz_max      =  0.1;
        int   electron_losthits_max = 1;
        float tight_el_sip3d_max   = 4.0;
        float tight_el_miniiso_max = 0.1;
        float loose_el_sip3d_max   = 8.0;
        float loose_el_miniiso_max = 0.4;
        // trilepton 선택 (3Mu)
        float lead_mu_pt   = 20.0;
        float sub_mu_pt    = 10.0;
        // 1E2Mu 의 비대칭 leg (둘 중 하나만 만족하면 된다)
        float emu_lead_mu_pt = 25.0;
        float emu_sub_el_pt  = 15.0;
        float emu_lead_el_pt = 25.0;
        float emu_sub_mu_pt  = 10.0;
        float pair_mass_min = 12.0;
        int   n_jet_min     = 2;
        int   n_bjet_min    = 1;
        // jet
        float jet_pt_min   = 20.0;     // ClosFakeRate 는 20 GeV (fake.cc 는 25)
        float jet_eta_max  =  2.4;     // 2016: 2.4, 그 외: 2.5
        float jet_lep_dr   =  0.4;
        float z_mass       = 91.2;
    } cuts;

    //==== 이벤트마다 새로 채우는 physics objects
    RVec<Muon> looseMuons, tightMuons;
    //   veto 는 pT > 10, loose/tight 는 pT > 15 (ClosFakeRate 와 같은 두 단계).
    //   jet cleaning 과 채널 판정은 veto 를, fake 후보는 loose 를 쓴다.
    RVec<Electron> vetoElectrons, looseElectrons, tightElectrons;
    RVec<Jet> rawJets, jets, bjets;
    Particle METv;
    RVec<Gen> gens;
    RVec<GenJet> genJets;

    //==== 측정 analyzer 와 글자 그대로 같은 helper (WP 이 어긋나면 안 된다)
    bool PassMuonWP(const Muon &mu, const TString &wp) const;        // fake.cc
    bool PassElectronWP(const Electron &el, const TString &wp) const; // elecfake.cc
    bool PassVetoMapJet(const Jet &jet, const RVec<Muon> &muons) const;

    //==== ClosFakeRate::configureChargeOf 와 동일: (ss1, ss2, os) 로 정렬
    std::tuple<Muon, Muon, Muon> configureChargeOf(const RVec<Muon> &muons) const;

    //==== lepton 한 개의 FR 입력 feature 를 mu{slot}_* / el{slot}_* 브랜치로 쓴다.
    //   source-jet / SV 블록은 각각 fake.cc / elecfake.cc 의 fillTreeRow 에서
    //   그대로 옮겼다 — 모델이 학습한 입력과 정의가 한 글자도 달라선 안 된다.
    //
    //   **포인터를 받는다.** SetBranch 는 deque 에 push 하고 FillTrees 가 비우는
    //   구조라 어떤 row 에서든 브랜치를 하나라도 건너뛰면 다음 row 의 주소가
    //   어긋난다. 채널마다 안 쓰는 slot 이 생기므로(3Mu 는 el1, 1E2Mu 는 mu3)
    //   nullptr 을 넘기면 같은 브랜치 목록을 sentinel 로 채우게 해서 그 불변식을
    //   호출자의 규율이 아니라 구조로 보장한다.
    void fillMuonSlot(int slot, const Muon *mu);
    void fillElectronSlot(int slot, const Electron *el);
};

#endif
