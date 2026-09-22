#include "muid.h"

//==== Constructor and Destructor
muid::muid() {}
muid::~muid() {}

//==============================================================
// Initialize
//==============================================================
void muid::initializeAnalyzer() {

    //==== loose WP 은 Run 별로 다르다 (closnff.cc:64 / fake.cc:57 과 같은 값).
    //     여기가 어긋나면 이 스터디의 결론이 측정 영역과 다른 WP 을 재게 된다.
    if (Run == 3) { cuts.loose_sip3d_max = 8.0; cuts.loose_miniiso_max = 0.4; }

    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);
}

//==============================================================
// source jet (fake.cc:504-515 / closnff.cc:615-623 과 같은 방식)
//==============================================================
const Jet *muid::MatchedJet(const Muon &mu) const {
    const short jetIdx = mu.JetIdx();
    if (jetIdx < 0) return nullptr;            // 연결된 jet 자체가 없는 muon
    for (const auto &jet : rawJets) {
        if (jet.OriginalIndex() != jetIdx) continue;
        return &jet;
    }
    return nullptr;                            // jetIdx 는 있으나 컬렉션에 없음 (<15 GeV)
}

//==============================================================
// jet flavour code (closnff.cc:591-602 과 글자 그대로 같다)
//==============================================================
int muid::FlavorCode(const Jet *jet) const {
    if (IsDATA) return -2;
    if (!jet)   return -1;
    if (jet->genJetIdx() < 0) return 4;        // pileup, flavour 보다 먼저
    const int hf = jet->hadronFlavour();
    if (hf == 5) return 0;
    if (hf == 4) return 1;
    const int ap = abs(jet->partonFlavour());
    if (ap == 21) return 3;
    if (ap == 1 || ap == 2 || ap == 3) return 2;
    return -1;
}

//==============================================================
// Event loop
//==============================================================
void muid::executeEvent() {

    //==== truth 스터디다. data 에는 gens 가 없으므로 아예 돌지 않는다.
    if (IsDATA) return;

    Event ev = GetEvent();
    rawJets = GetAllJets();

    //==== weight 은 closnff.cc:215 의 cfW 와 **같은 정의**여야 한다.
    //     그래야 여기 CutFlow bin 0-2 가 CutFlow_3Mu bin 0-2 와 자릿수까지 맞고,
    //     per-muon 효율을 기존 trilepton cutflow 에 그대로 꽂을 수 있다.
    const float w = MCweight() * ev.GetTriggerLumi("Full") * GetL1PrefireWeight()
                    * myCorr->GetPUWeight(ev.nTrueInt());

    FillHist("CutFlow", 0.5, w, 3, 0., 3.);
    if (!PassNoiseFilter(rawJets, ev)) return;
    FillHist("CutFlow", 1.5, w, 3, 0., 3.);

    RVec<Muon> rawMuons = GetAllMuons();

    if (!PassVetoMap(rawJets, rawMuons, "jetvetomap")) return;
    FillHist("CutFlow", 2.5, w, 3, 0., 3.);

    gens = GetAllGens();

    //==== 축 정의 ------------------------------------------------------
    //   step 축은 7 bin: 0 모든 muon (컷 없음) / 1 +Medium / 2 +dZ /
    //                    3 +TkRelIso / 4 +SIP3D / 5 +miniIso (= loose ID) /
    //                    6 +acceptance
    //   **acceptance 를 맨 뒤에 둔다.** 알고 싶은 것은 "GetAllMuons() 바로
    //   다음에 loose ID 를 걸면 fake 가 얼마나 남는가" 이고, pT/eta 는 ID 가
    //   아니라 운동학이다. 앞에 두면 ID 효율이 acceptance 효율에 섞인다.
    static const RVec<float> STEPBINS = {0., 1., 2., 3., 4., 5., 6., 7.};

    //   pT 축은 fake.cc:45 의 ptCorr binning 에 10 GeV 아래를 붙였다.
    //   step 0 (컷 없음) 에는 pT < 10 muon 이 들어있으므로 그것들이 조용히
    //   underflow 로 사라지지 않게 하려는 것이다.
    static const RVec<float> PTBINS =
        {0., 5., 10., 12., 14., 17., 20., 30., 50., 100., 200.};
    //   |eta| 도 같은 이유로 2.4 위에 한 칸을 더 둔다 (fake.cc:46 + overflow bin).
    static const RVec<float> ETABINS = {0., 0.9, 1.6, 2.4, 5.0};

    //   isolation 은 균일 binning 에서 b-jet muon 의 절반이 overflow 로 빠진다
    //   (track.cc:128-137). 0 을 첫 bin 으로 두고 1e-3 ~ 1e2 log 등간격.
    static const RVec<float> ISOBINS = [] {
        RVec<float> b{0.f};
        for (int k = 0; k <= 100; k++) b.push_back(std::pow(10.f, -3.f + 0.05f * k));
        return b;
    }();

    int nMu = 0, nFake = 0, nMuInAcc = 0, nFakeInAcc = 0, nLoose = 0, nLooseFake = 0;

    for (const auto &mu : rawMuons) {

        //==== fake 판정. GetLeptonType > 0 이 prompt 라는 관례는 fake.cc:473-478.
        //     ID 를 하나도 요구하지 않은 muon 에도 그대로 쓸 수 있다
        //     (track.cc:104 / TestGetLeptonType.cc:47 이 이미 그렇게 쓴다).
        const int  lt     = GetLeptonType(mu, gens);
        const bool isFake = (lt <= 0);

        //==== acceptance. closnff.cc:236 의 조건을 **부정 형태 그대로** 옮겼다
        //     (>= / <= 경계가 한 칸 어긋나면 nLooseMu 대조가 깨진다).
        //     순차 단계에서는 loose ID **뒤에** 적용한다 (위 축 정의 주석).
        const bool inAcc =
            !(mu.Pt() < cuts.muon_pt_min || fabs(mu.Eta()) > cuts.muon_eta_max);

        //==== loose WP 다섯 컷을 따로따로 (fake::PassMuonWP 를 펼친 것)
        const bool c[5] = { mu.isPOGMediumId(),
                            fabs(mu.dZ())     < cuts.muon_dz_max,
                            mu.TkRelIso()     < cuts.muon_tkiso_max,
                            fabs(mu.SIP3D())  < cuts.loose_sip3d_max,
                            mu.MiniPFRelIso() < cuts.loose_miniiso_max };
        const int npass = (int)c[0] + c[1] + c[2] + c[3] + c[4];
        const bool passID = (npass == 5);              // loose ID 만 (운동학 무관)
        const bool isLoose = inAcc && passID;          // closnff 의 looseMuons 정의

        nMu++;                        if (isFake) nFake++;
        if (inAcc)   { nMuInAcc++;    if (isFake) nFakeInAcc++; }
        if (isLoose) { nLoose++;      if (isFake) nLooseFake++; }

        //==== 순차 단계 (0..6). 처음 실패한 컷에서 멈춘다.
        //     0 은 컷이 전혀 없는 모든 muon, 1-5 가 loose ID 다섯 컷,
        //     6 이 그 뒤에 붙는 acceptance.
        int nseq = 0;
        for (int i = 0; i < 5; i++) { if (!c[i]) break; nseq = 1 + i; }
        if (passID && inAcc) nseq = 6;

        //==== source jet flavour 축: 0 nojet / 1 b / 2 c / 3 uds / 4 g /
        //     5 pileup / 6 undef. FlavorCode 는 "jet 없음" 과 "매칭됐지만
        //     flavour 미확인" 을 둘 다 -1 로 주므로 여기서만 갈라 둔다.
        const Jet *mj   = MatchedJet(mu);
        const int  fc   = FlavorCode(mj);
        const int  fbin = !mj ? 0 : (fc >= 0 ? fc + 1 : 6);

        const float absEta = fabs(mu.Eta());

        //==== All/ 은 모든 muon, Fake/ 는 fake 만. prompt 는 All - Fake 로 얻는다.
        for (int d = 0; d < 2; d++) {
            if (d == 1 && !isFake) continue;
            const TString p = (d == 0) ? "All/" : "Fake/";

            for (int b = 0; b <= nseq; b++) {
                FillHist(p + "Seq",           b + 0.5,           w, 7, 0., 7.);
                FillHist(p + "Seq_vs_Flav",   b + 0.5, fbin,     w, 7, 0., 7., 7, -0.5, 6.5);
                FillHist(p + "Seq_vs_LType",  b + 0.5, lt,       w, 7, 0., 7., 13, -6.5, 6.5);
                FillHist(p + "Seq_vs_Pt",     b + 0.5, mu.Pt(),  w, STEPBINS, PTBINS);
                FillHist(p + "Seq_vs_AbsEta", b + 0.5, absEta,   w, STEPBINS, ETABINS);
            }

            //==== N-1 과 marginal. bin 0 은 둘의 공통 분모 = **컷을 전혀 걸지
            //     않은 모든 muon** (acceptance 도 요구하지 않는다).
            //     Seq 와 달리 순서에 의존하지 않는다 — Seq 만 보면 손실이
            //     "먼저 적용된 컷" 쪽으로 몰려 보이기 때문에 둘 다 필요하다.
            //     bin 6 = 다섯 컷 전부 = loose ID (acceptance 없음).
            FillHist(p + "Nm1",  0.5, w, 7, 0., 7.);
            FillHist(p + "Marg", 0.5, w, 7, 0., 7.);
            for (int i = 0; i < 5; i++) {
                if (npass - (int)c[i] == 4) FillHist(p + "Nm1",  i + 1.5, w, 7, 0., 7.);
                if (c[i])                   FillHist(p + "Marg", i + 1.5, w, 7, 0., 7.);
            }
            if (npass == 5) {
                FillHist(p + "Nm1",  6.5, w, 7, 0., 7.);
                FillHist(p + "Marg", 6.5, w, 7, 0., 7.);
            }

            //==== N-1 분포: 나머지 네 컷을 통과한 muon 의 그 변수 분포.
            //     여기도 acceptance 를 요구하지 않는다 (위와 같은 이유).
            //     범위를 컷보다 훨씬 넓게 잡은 것은 의도적이다 — 떨어진 fake 가
            //     컷에서 얼마나 멀리 있는지가 이 플롯의 요점이다.
            if (npass - (int)c[0] == 4)
                FillHist(p + "Nm1_MediumId", (int)c[0],        w, 2, -0.5, 1.5);
            if (npass - (int)c[1] == 4)
                FillHist(p + "Nm1_dZ",       mu.dZ(),          w, 200, -1., 1.);
            if (npass - (int)c[2] == 4)
                FillHist(p + "Nm1_TkRelIso", mu.TkRelIso(),    w, ISOBINS);
            if (npass - (int)c[3] == 4)
                FillHist(p + "Nm1_SIP3D",    fabs(mu.SIP3D()), w, 200, 0., 200.);
            if (npass - (int)c[4] == 4)
                FillHist(p + "Nm1_MiniIso",  mu.MiniPFRelIso(), w, ISOBINS);

            //==== 같은 변수의 inclusive 분포 (다른 컷 요구 없음). N-1 과 나란히
            //     두면 "이 컷이 거르는 것" 과 "다른 컷이 이미 거른 것" 이 갈린다.
            FillHist(p + "Inc_MediumId", (int)c[0],         w, 2, -0.5, 1.5);
            FillHist(p + "Inc_dZ",       mu.dZ(),           w, 200, -1., 1.);
            FillHist(p + "Inc_TkRelIso", mu.TkRelIso(),     w, ISOBINS);
            FillHist(p + "Inc_SIP3D",    fabs(mu.SIP3D()),  w, 200, 0., 200.);
            FillHist(p + "Inc_MiniIso",  mu.MiniPFRelIso(), w, ISOBINS);
        }
    }

    //==== 이벤트당 한 번. 0.00107% 를 "후보가 있느냐" 와 "ID 를 통과하느냐" 로
    //     가르는 것이 이 네 히스토그램이다.
    //     nLooseMu 는 closnff.cc:290 의 CutFlow_nLooseMu 와 정의가 같아야 한다
    //     (같은 6 bin, 같은 min(n,5) 처리) — 대조용 게이트다.
    FillHist("nRawMu",       min(nMu,        5) + 0.5, w, 6, 0., 6.);
    FillHist("nFakeMu",      min(nFake,      5) + 0.5, w, 6, 0., 6.);
    FillHist("nRawMuInAcc",  min(nMuInAcc,   5) + 0.5, w, 6, 0., 6.);
    FillHist("nFakeMuInAcc", min(nFakeInAcc, 5) + 0.5, w, 6, 0., 6.);
    FillHist("nLooseMu",     min(nLoose,     5) + 0.5, w, 6, 0., 6.);
    FillHist("nLooseFakeMu", min(nLooseFake, 5) + 0.5, w, 6, 0., 6.);
}
