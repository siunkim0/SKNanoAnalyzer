#ifndef muid_h
#define muid_h

#include "AnalyzerCore.h"

//==============================================================
// Muon loose WP 의 "fake muon 에 대한 효율"을 컷 하나씩 분해한다.
//
// 왜 필요한가: trilepton 영역의 `== 3 loose mu` 단계가 전체 이벤트의
// 0.00107% 만 남긴다 (DY Run2, 3.385e9 -> 36,285). 이 숫자가
//   (a) DY 이벤트에 애초에 세 번째 muon **후보**가 거의 없어서인지,
//   (b) loose ID 가 fake muon 에게 그만큼 가혹해서인지
// 를 구분할 수 없었다. 디스크의 어떤 트리로도 답할 수 없다 —
// fake/MakeTree 도 closnff 도 **이미 loose 를 통과한** muon 만 저장하기
// 때문이다. 그래서 GetAllMuons() 에서 시작하는 새 pass 가 필요하다.
//
// 측정: 컷을 전혀 걸지 않은 모든 muon 을 GetLeptonType(mu, gens) <= 0
// 으로 fake 라 부르고 (fake.cc:473-478 과 같은 관례), P(fake 가 loose 통과)
// 와 **다섯 컷 중 무엇이 떨어뜨렸는지**를 센다.
//
// loose WP (fake::PassMuonWP, fake.cc:134-149 / closnff.cc:147-160 과 동일):
//   1 isPOGMediumId()
//   2 |dZ()|          < 0.1
//   3 TkRelIso()      < 0.4
//   4 |SIP3D()|       < 5   (Run3: 8)
//   5 MiniPFRelIso()  < 0.6 (Run3: 0.4)
// dXY 컷은 **없다**.
//
// **acceptance (pT > 10, |eta| < 2.4) 는 맨 뒤에 붙인다.** closnff.cc:236 은
// ID 앞에 걸지만, 여기서 보고 싶은 것은 "GetAllMuons() 바로 다음에 loose ID 를
// 걸면 fake 가 얼마나 남는가" 다. pT/eta 는 ID 가 아니라 운동학이므로 앞에
// 두면 ID 효율이 acceptance 효율에 섞여버린다. 순차 단계는 그래서
//   0 컷 없음 -> 1..5 loose ID -> 6 acceptance
// 이고, N-1 / marginal 의 분모도 **컷을 전혀 걸지 않은 모든 muon** 이다.
// nLooseMu (closnff 대조용) 만은 closnff 정의 그대로 acceptance + ID 다.
//
// **WP 숫자 중복에 대하여.** 이 값은 이제 다섯 군데에 있다 (fake.h,
// closnff.h, 두 initializeAnalyzer, 그리고 여기). fake.cc/closnff.cc 는
// 서로 어긋나면 closure 자체가 무의미해지므로 그 둘은 반드시 같아야 한다.
// muid 는 **측정만** 하고 closure 에 아무것도 먹이지 않으므로 여기서
// 어긋나도 closure 는 깨지지 않는다 — 대신 이 스터디의 결론이 틀려진다.
// 기준은 언제나 fake.cc:134 다.
//
// 출력은 히스토그램뿐이다 (tree 없음). 디렉토리 두 개:
//   All/  : 모든 muon
//   Fake/ : GetLeptonType <= 0 인 muon 만
// prompt 는 All - Fake 로 얻으므로 따로 채우지 않는다.
//
// 이벤트 preamble 은 closnff::executeEvent 와 **일부러 글자 그대로 같다**
// (weight / noise filter / jet veto map). 그래야 여기서 잰 per-muon 효율을
// 기존 CutFlow_3Mu 에 그대로 꽂을 수 있다.
//==============================================================
class muid : public AnalyzerCore {
public:
    muid();
    ~muid();

    void initializeAnalyzer();
    void executeEvent();

    //==== 컷 값. muon 부분은 closnff.h::AnalysisCuts 와 같은 숫자여야 한다.
    struct AnalysisCuts {
        float muon_pt_min    = 10.0;
        float muon_eta_max   =  2.4;
        float muon_dz_max    =  0.1;   // cm
        float muon_tkiso_max =  0.4;   // rel. tracker iso R03
        // loose WP (Run 2 기본값, Run 3 는 initializeAnalyzer 에서 덮어쓴다)
        float loose_sip3d_max   = 5.0;
        float loose_miniiso_max = 0.6;
    } cuts;

    //==== source jet flavour -> fake.cc / closnff.cc 와 같은 5-class 정수 코드.
    //     0=b 1=c 2=uds 3=g 4=pileup, -1=jet 없음/미확인.
    //     공유하지 않고 복사한 이유는 PassMuonWP 과 같다.
    int FlavorCode(const Jet *jet) const;

    //==== muon 이 clustering 된 PF jet 을 NanoAOD jetIdx 로 찾는다
    //     (dR 최근접 매칭은 muon 이 실제로 안 들어있는 이웃 jet 을 집는다).
    const Jet *MatchedJet(const Muon &mu) const;

    //==== 이벤트마다 새로 채운다
    RVec<Jet> rawJets;
    RVec<Gen> gens;
};

#endif
