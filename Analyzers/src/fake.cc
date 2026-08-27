#include "fake.h"

//==== Constructor and Destructor
fake::fake() : MeasFakeMu8(false), MeasFakeMu17(false), MakeTree(false) {}
fake::~fake() {}

//==== Initialize variables
void fake::initializeAnalyzer() {

    //==== Userflags
    MeasFakeMu8  = HasFlag("MeasFakeMu8");   // Mu8 트리거 경로만 측정
    MeasFakeMu17 = HasFlag("MeasFakeMu17");  // Mu17 트리거 경로만 측정
    MakeTree     = HasFlag("MakeTree");      // flavor BDT 학습용 flat ntuple

    //==== 플래그가 없으면 두 경로 모두 측정
    if (!MeasFakeMu8 && !MeasFakeMu17) {
        MeasFakeMu8 = true;
        MeasFakeMu17 = true;
    }

    //==== prescaled single muon trigger 경로 정의 (analysis note Table 38)
    //==== {히스토그램 prefix, HLT 이름, muon pT 컷, cone-corrected pT 컷}
    if (MeasFakeMu8)  paths.push_back({"Mu8",  "HLT_Mu8_TrkIsoVVL",  10., 10.});
    if (MeasFakeMu17) paths.push_back({"Mu17", "HLT_Mu17_TrkIsoVVL", 20., 30.});

    //==== systematic 목록: Central = away jet pT > 40 GeV
    //==== away jet pT variation 은 2D 히스토그램 몇 개만 더 채우는 것이라
    //==== 부담이 없으므로 항상 켠다
    systs = {"Central", "AwayJetPt30", "AwayJetPt60"};

    //==== fake rate 2D 히스토그램 binning: cone-corrected pT x |eta|
    //==== 공식 측정처럼 [100,200] 을 따로 둔다 (FR 파일을 만들 때는 100 까지만 사용)
    ptCorrBins = {10., 12., 14., 17., 20., 30., 50., 100., 200.};
    absEtaBins = {0., 0.9, 1.6, 2.4};

    //==== veto electron 의 loose MVANoIso 컷 (Run 2 / Run 3 에서 값이 다르다)
    //==== elecfake.cc / closnff.cc 와 완전히 같은 값 — 어긋나면 측정과 적용의
    //==== electron veto 가 달라진다
    if (Run == 2) { looseMvaIB = 0.985; looseMvaOB = 0.96; looseMvaEC = 0.85; }
    else          { looseMvaIB = 0.8;   looseMvaOB = 0.5;  looseMvaEC = -0.8; }

    //==== muon loose WP 은 Run 별로 다르다 (공식 MeasFakeRateV4 의
    //==== HcToWALooseRun2 / HcToWALooseRun3). tight 은 Run 무관.
    //==== 이 값이 어긋나면 FR 측정 자체가 공식 map 과 달라진다 —
    //==== Run3 에서 우리 FR 이 공식 대비 ptcorr 에 따라 1.30 ~ 0.81 로
    //==== 어긋났던 원인이 정확히 이것이었다.
    if (Run == 3) { cuts.loose_sip3d_max = 8.0; cuts.loose_miniiso_max = 0.4; }

    //==== electron loose WP 도 Run 별로 다르다 (공식
    //==== Electron::Pass_HcToWALooseRun2 SIP3D<8 vs ...Run3 SIP3D<6;
    //==== miniiso 는 둘 다 0.4, MVA 컷은 위에서 이미 Run 별로 잡았다).
    //==== tight (HcToWATight) 는 Run 무관하다.
    if (Run == 3) cuts.loose_el_sip3d_max = 6.0;

    //==== era 별 jet |eta| 컷
    if (DataEra == "2016preVFP" || DataEra == "2016postVFP") {
        cuts.jet_eta_max = 2.4;
    } else {
        cuts.jet_eta_max = 2.5;
    }

    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);

    //==== flavor BDT 학습용 flat ntuple (loose muon 한 개당 한 줄)
    if (MakeTree) NewTree("mu");
}

//==============================================================
// Muon WP (analysis note 의 WP 표)
// 공통       : POG medium ID, |dz| < 0.1 cm, rel. tracker iso (R03) < 0.4
// tight      : |SIP3D| < 3, rel. mini-iso < 0.1
// loose(Run2): |SIP3D| < 5, rel. mini-iso < 0.6
//==============================================================
bool fake::PassMuonWP(const Muon &mu, const TString &wp) const {
    //==== 공통 컷
    if (!mu.isPOGMediumId()) return false;
    if (!(fabs(mu.dZ()) < cuts.muon_dz_max)) return false;
    if (!(mu.TkRelIso() < cuts.muon_tkiso_max)) return false;   // trigger emulation

    //==== tight / loose 별 컷
    if (wp == "tight") {
        if (!(fabs(mu.SIP3D()) < cuts.tight_sip3d_max)) return false;
        if (!(mu.MiniPFRelIso() < cuts.tight_miniiso_max)) return false;
    } else {
        if (!(fabs(mu.SIP3D()) < cuts.loose_sip3d_max)) return false;
        if (!(mu.MiniPFRelIso() < cuts.loose_miniiso_max)) return false;
    }
    return true;
}

//==============================================================
// Electron veto WP: 공식 HcToWALooseRun2 (Run 2 전용)
// framework Electron::Pass_HcToWALooseRun2() 는 raw MVA 컷 부호가
// 반전된 버그가 있어 직접 구현한다 (elecfake.cc 의 loose WP 와 동일)
//==============================================================
bool fake::PassElectronVeto(const Electron &el) const {
    if (!el.Pass_HcToWABaseline()) return false;
    if (!(el.SIP3D() < cuts.loose_el_sip3d_max)) return false;
    if (!(el.MiniPFRelIso() < cuts.loose_el_miniiso_max)) return false;
    //==== raw MVANoIso 컷 (IB, OB, EC); wp90 통과도 인정 (tight ⊂ loose)
    //==== 값은 initializeAnalyzer 에서 Run 2 / Run 3 별로 정한다
    float mvaCut = looseMvaEC;
    if (el.etaRegion() == Electron::ETAREGION::IB)      mvaCut = looseMvaIB;
    else if (el.etaRegion() == Electron::ETAREGION::OB) mvaCut = looseMvaOB;
    if (!(el.isMVANoIsoWP90() || el.MvaNoIso() > mvaCut)) return false;
    return true;
}

//==============================================================
// jet 단위 veto map: 공식 framework 의 로직을 직접 구현
// (우리 fork 의 AnalyzerCore::PassVetoMap 은 muon dR < 0.2 이거나
//  loose PU ID 를 통과 못 한 jet 자체를 veto 해 버려 동작이 다르다)
// EM fraction > 0.9 또는 muon dR < 0.2 인 jet 은 검사 대상에서 제외
//==============================================================
bool fake::PassVetoMapJet(const Jet &jet, const RVec<Muon> &muons) const {
    if (jet.chEmEF() + jet.neEmEF() > 0.9) return true;
    for (const auto &mu : muons)
        if (jet.DeltaR(mu) < 0.2) return true;
    return !myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), "jetvetomap");
}

void fake::executeEvent() {

    //==== 이벤트와 raw physics object 읽기
    Event ev = GetEvent();
    rawJets = GetAllJets();   // flavor 매칭용으로 lepton cleaning 이전 컬렉션을 보관

    //==== noise (MET) filter 통과 요구
    if (!PassNoiseFilter(rawJets, ev)) return;

    RVec<Muon> rawMuons = GetAllMuons();

    //==== jet veto map (공식 측정과 동일하게 event 단위로 적용)
    if (!PassVetoMap(rawJets, rawMuons, "jetvetomap")) return;

    RVec<Electron> rawElectrons = GetAllElectrons();
    sort(rawMuons.begin(), rawMuons.end(), PtComparing);

    //==== Step 1: loose / tight muon 선택 (tight 는 loose 의 부분집합)
    looseMuons.clear();
    tightMuons.clear();
    for (const auto &mu : rawMuons) {
        if (mu.Pt() < cuts.muon_pt_min || fabs(mu.Eta()) > cuts.muon_eta_max) continue;
        if (!PassMuonWP(mu, "loose")) continue;
        looseMuons.push_back(mu);
        if (PassMuonWP(mu, "tight")) tightMuons.push_back(mu);
    }

    //==== Step 2: veto 용 loose electron 선택 (framework loose ID 의
    //====         raw MVA 컷 반전 버그 때문에 PassElectronVeto 로 직접 선택)
    looseElectrons.clear();
    for (const auto &el : rawElectrons) {
        if (!(el.Pt() > cuts.electron_pt_min)) continue;
        if (!(fabs(el.Eta()) < cuts.electron_eta_max)) continue;
        if (PassElectronVeto(el)) looseElectrons.push_back(el);
    }

    //==== Step 3: jet 선택 (공식 측정과 같은 순서)
    //====         tight ID → loose lepton 과 dR < 0.4 로 겹치는 jet 제거
    //====         → jet 단위 veto map → loose PU ID (Run2)
    jets = SelectJets(rawJets, "tight", cuts.jet_pt_min, cuts.jet_eta_max);
    jets = JetsVetoLeptonInside(jets, looseElectrons, looseMuons, cuts.jet_lep_dr);
    RVec<Jet> jetsNoPuId;   // PU ID SF 계산용 (PU ID 적용 전)
    if (Run == 2) {
        RVec<Jet> vetoMapped;
        for (const auto &jet : jets)
            if (PassVetoMapJet(jet, rawMuons)) vetoMapped.push_back(jet);
        jetsNoPuId = vetoMapped;
        jets = SelectJets(vetoMapped, "loosePuId", cuts.jet_pt_min, cuts.jet_eta_max);
    }
    sort(jets.begin(), jets.end(), PtComparing);

    //==== Step 4: PUPPI MET + Type-I correction
    METv = ApplyTypeICorrection(ev.GetMETVector(Event::MET_Type::PUPPI),
                                rawJets, rawElectrons, rawMuons);

    //==== Step 5: gen 정보 (MC 에서 prompt / fake 구분 + PU ID SF 용)
    gens = IsDATA ? RVec<Gen>() : GetAllGens();
    genJets = IsDATA ? RVec<GenJet>() : GetAllGenJets();

    //==== 트리거 경로 공통 event weight 보정 (MC only, 공식 측정과 동일한 세트)
    float evtSF = 1.;
    if (!IsDATA) {
        //==== top pT reweight (TT 샘플만 해당)
        if (MCSample.Contains("TTLL") || MCSample.Contains("TTLJ"))
            evtSF *= myCorr->GetTopPtReweight(gens);

        //==== pileup jet ID SF: PU ID 적용 전, away jet pT 컷(40)을 넘는 jet 대상
        if (Run == 2) {
            RVec<Jet> sfJets;
            for (const auto &jet : jetsNoPuId)
                if (jet.Pt() > cuts.awayjet_pt) sfJets.push_back(jet);
            unordered_map<int, int> matchedIdx =
                GenJetMatching(sfJets, genJets, fixedGridRhoFastjetAll, 0.4, 10.);
            evtSF *= myCorr->GetPileupJetIDSF(sfJets, matchedIdx, "loose");
        }
    }

    //==== Step 6: 트리거 경로(Mu8 / Mu17)마다 측정
    for (const auto &path : paths) {

        //==== 해당 트리거 통과 요구
        if (!ev.PassTrigger(path.trigger)) continue;

        //==== weight 계산 (MC only)
        float weight = 1.;
        if (!IsDATA) {
            //==== prescale 은 trigger lumi 에 반영되어 있지 않으므로 여기의
            //==== normalization 은 임시값; 최종 MC normalization 은
            //==== ZEnriched 영역에서 offline 으로 뽑는다
            weight = MCweight() * ev.GetTriggerLumi(path.trigger);
            weight *= GetL1PrefireWeight();
            weight *= myCorr->GetPUWeight(ev.nTrueInt());
            weight *= evtSF;   // top pT reweight + pileup jet ID SF
        }

        measureFakeRate(path, weight);    // fake rate 측정 영역
        fillZEnriched(path, ev, weight);  // MC normalization 용 Z 영역
    }

    //==== Step 7: flavor BDT 학습용 ntuple (userflag MakeTree)
    //====         경로 loop 밖에서 이벤트당 한 번만 부른다 — loop 안에서 쓰면
    //====         두 트리거를 모두 통과한 muon 이 두 줄로 중복 기록된다
    fillTreeRow(ev, evtSF);
}

//==============================================================
// Fake rate 측정 영역
// selection: loose muon 딱 1개 + electron 0개
//            + away jet (pT > 40, dR > 0.7) + MET < 25 + MT < 25
// FR(ptcorr, |eta|) = tight 히스토그램 / loose 히스토그램 (offline 에서 계산)
//==============================================================
void fake::measureFakeRate(const TriggerPath &path, float weight) {

    //==== Step 1: single lepton 이벤트 요구 (loose muon 1개, electron 0개)
    if (looseMuons.size() != 1) return;
    if (looseElectrons.size() != 0) return;

    const Muon &mu = looseMuons.at(0);

    //==== muon RECO SF (공식 측정과 동일하게 측정 대상 muon 에 적용)
    if (!IsDATA) weight *= myCorr->GetMuonRECOSF(mu);

    //==== Step 2: 트리거별 muon pT 컷 (Mu8: 10, Mu17: 20)
    if (mu.Pt() < path.ptCut) return;

    //==== Step 3: cone-corrected pT = pT * (1 + max(0, miniIso - 0.1))
    //====         tight iso 기준(0.1)을 넘는 만큼 pT 를 보정해 준다
    const float ptCorr = mu.Pt() * (1. + max(0.f, mu.MiniPFRelIso() - cuts.tight_miniiso_max));
    if (ptCorr < path.ptCorrCut) return;   // Mu8: 10, Mu17: 30

    //==== tight WP 통과 여부, prompt/fake 구분 (MC only)
    const bool isTight = PassMuonWP(mu, "tight");
    const TString ltype = IsDATA ? "" : LeptonTypeToString(GetLeptonType(mu, gens));

    //==== transverse mass MT(mu, MET)
    const float MT = sqrt(2. * mu.Pt() * METv.Pt() * (1. - cos(mu.DeltaPhi(METv))));

    //==== loose 는 항상 채우고, tight 통과 시 tight 도 채운다 → FR = tight / loose
    RVec<TString> tags = {"loose"};
    if (isTight) tags.push_back("tight");

    //==== overflow 는 마지막 bin 에 넣는다
    const float xval = min(ptCorr, ptCorrBins.back() - 0.1f);
    const float yval = fabs(mu.Eta());

    //==== away jet pT 컷을 바꿔가며 반복 (Central = 40, syst = 30 / 60)
    for (const auto &syst : systs) {
        float awayJetPtCut = cuts.awayjet_pt;
        if (syst == "AwayJetPt30") awayJetPtCut = 30.;
        else if (syst == "AwayJetPt60") awayJetPtCut = 60.;

        //==== Step 4: away jet 요구 (pT > cut && dR(mu, jet) > 0.7)
        int nJets = 0;
        bool hasAwayJet = false;
        for (const auto &jet : jets) {
            if (jet.Pt() < awayJetPtCut) continue;
            nJets++;
            if (jet.DeltaR(mu) > cuts.awayjet_dr) hasAwayJet = true;
        }
        if (!hasAwayJet) continue;

        //==== Inclusive (validation) 영역: MET / MT 컷 없이 kinematics 채움
        //==== 공식(MeasFakeRateV4) QCD MC FR 은 이 영역에서 측정되므로
        //==== 비교용 2D FR 히스토그램도 같이 채운다
        if (syst == "Central") {
            for (const auto &tag : tags) {
                fillMuonKinematics(path.name + "/Inclusive/" + tag,
                                   mu, ptCorr, MT, nJets, weight);
                FillHist(path.name + "/Inclusive/" + tag + "/fake_ptcorr_abseta",
                         xval, yval, weight, ptCorrBins, absEtaBins);
                if (!IsDATA) {
                    fillMuonKinematics(path.name + "/Inclusive/" + tag + "/" + ltype,
                                       mu, ptCorr, MT, nJets, weight);
                    FillHist(path.name + "/Inclusive/" + tag + "/" + ltype + "/fake_ptcorr_abseta",
                             xval, yval, weight, ptCorrBins, absEtaBins);
                }
            }
        }

        //==== Step 5: 측정 영역 컷 (W/Z 의 prompt muon 오염 억제)
        if (!(METv.Pt() < cuts.met_max)) continue;
        if (!(MT < cuts.mt_max)) continue;

        //==== FR 계산용 2D 히스토그램 (cone-corrected pT x |eta|)
        const TString base = path.name + "/MeasReg/" + syst;
        for (const auto &tag : tags) {
            FillHist(base + "/" + tag + "/fake_ptcorr_abseta",
                     xval, yval, weight, ptCorrBins, absEtaBins);
            if (!IsDATA)
                FillHist(base + "/" + tag + "/" + ltype + "/fake_ptcorr_abseta",
                         xval, yval, weight, ptCorrBins, absEtaBins);
        }

        //==== per-parton-flavor 2D 히스토그램 (flavor-weighted FR 용, MC only)
        //==== source jet 을 gen / reco 두 방식으로 매칭해 b/c/s/d/u/g 로 분류
        //==== FR_flavor = tight_flavor / loose_flavor (offline)
        //==== fake muon 만 (GetLeptonType<=0): TTLJ 등 prompt 오염 표본에서 필수
        if (!IsDATA && syst == "Central" && GetLeptonType(mu, gens) <= 0) {
            const TString gf = GenJetFlavor(mu);
            const TString rf = RecoJetFlavor(mu);
            for (const auto &tag : tags) {
                FillHist(base + "/" + tag + "/genflav/" + gf + "/fake_ptcorr_abseta",
                         xval, yval, weight, ptCorrBins, absEtaBins);
                FillHist(base + "/" + tag + "/recoflav/" + rf + "/fake_ptcorr_abseta",
                         xval, yval, weight, ptCorrBins, absEtaBins);
            }
        }

        //==== 측정 영역 kinematics (Central only)
        if (syst == "Central") {
            for (const auto &tag : tags) {
                fillMuonKinematics(base + "/" + tag, mu, ptCorr, MT, nJets, weight);
                if (!IsDATA)
                    fillMuonKinematics(base + "/" + tag + "/" + ltype,
                                       mu, ptCorr, MT, nJets, weight);
            }
        }
    }
}

//==============================================================
// Z-enriched 영역: prescaled trigger 경로별 MC normalization 추출용
// selection: OS tight muon pair + |M(ll) - 91.2| < 15 + jet pT > 40 하나 이상
//==============================================================
void fake::fillZEnriched(const TriggerPath &path, const Event &ev, float weight) {

    //==== Step 1: loose 2개 == tight 2개 (딱 tight pair 만 있는 이벤트)
    if (!(looseMuons.size() == 2 && tightMuons.size() == 2)) return;
    if (!looseElectrons.empty()) return;

    const Muon &mu1 = tightMuons.at(0);
    const Muon &mu2 = tightMuons.at(1);

    //==== muon RECO SF (공식 측정과 동일하게 두 muon 모두에 적용)
    if (!IsDATA) weight *= myCorr->GetMuonRECOSF(tightMuons);

    //==== Step 2: opposite sign + pT 컷 (leading 은 트리거별 컷, subleading 은 10)
    if (mu1.Charge() + mu2.Charge() != 0) return;
    if (mu1.Pt() < path.ptCut || mu2.Pt() < cuts.muon_pt_min) return;

    //==== Step 3: jet pT > 40 하나 이상 (측정 영역과 같은 topology 맞추기)
    bool hasJet = false;
    for (const auto &jet : jets) {
        if (jet.Pt() > cuts.zjet_pt) { hasJet = true; break; }
    }
    if (!hasJet) return;

    //==== Step 4: Z mass window |M(ll) - 91.2| < 15
    const float mll = (mu1 + mu2).M();
    if (fabs(mll - cuts.z_mass) > cuts.z_window) return;

    //==== MC 는 두 muon 모두 prompt 인지에 따라 prompt / fake 로도 나눠 채운다
    RVec<TString> prefixes = {path.name + "/ZEnriched"};
    if (!IsDATA) {
        const bool isPrompt = (GetLeptonType(mu1, gens) > 0) && (GetLeptonType(mu2, gens) > 0);
        prefixes.push_back(prefixes[0] + (isPrompt ? "/prompt" : "/fake"));
    }

    for (const auto &prefix : prefixes) {
        FillHist(prefix + "/mll", mll, weight, 60, 76.2, 106.2);
        FillHist(prefix + "/mu1_pt", mu1.Pt(), weight, 200, 0., 200.);
        FillHist(prefix + "/mu2_pt", mu2.Pt(), weight, 200, 0., 200.);
        FillHist(prefix + "/mu1_eta", mu1.Eta(), weight, 48, -2.4, 2.4);
        FillHist(prefix + "/mu2_eta", mu2.Eta(), weight, 48, -2.4, 2.4);
        FillHist(prefix + "/nPV", ev.nPV(), weight, 70, 0., 70.);
    }
}

//==============================================================
// muon kinematics 히스토그램 묶음 (prefix 아래 7개)
//==============================================================
void fake::fillMuonKinematics(const TString &prefix, const Muon &mu,
                              float ptCorr, float MT, int nJets, float weight) {
    FillHist(prefix + "/pt", mu.Pt(), weight, 200, 0., 200.);
    FillHist(prefix + "/ptcorr", ptCorr, weight, 200, 0., 200.);
    FillHist(prefix + "/eta", mu.Eta(), weight, 48, -2.4, 2.4);
    FillHist(prefix + "/phi", mu.Phi(), weight, 64, -3.2, 3.2);
    FillHist(prefix + "/MET", METv.Pt(), weight, 100, 0., 100.);
    FillHist(prefix + "/MT", MT, weight, 100, 0., 100.);
    FillHist(prefix + "/nJets", nJets, weight, 10, 0., 10.);
}

//==============================================================
// GetLeptonType > 0 : prompt (EW/BSM prompt, tau daughter, internal conversion)
// GetLeptonType <= 0: hadron 기원 / external conversion / unmatched → fake
//==============================================================
TString fake::LeptonTypeToString(int leptonType) const {
    return (leptonType > 0) ? "prompt" : "fake";
}

//==============================================================
// Source jet parton flavor 분류
//   hadronFlavour (ghost B/C hadron matching) 로 b/c 를 먼저 잡고,
//   light (hadronFlavour==0) 는 |partonFlavour| 로 g/s/u/d 세분화한다.
//   partonFlavour: 21=gluon, 1=d, 2=u, 3=s (부호 = quark/antiquark)
//   매칭 실패(partonFlavour==0 등)는 unmatched
//   isPileup (대응하는 gen jet 이 없는 reco jet) 은 flavour 보다 먼저 본다:
//   공식 코드 MeasFakeRateV4::getMotherJetFlavour 가 genJetIdx<0 이면 hadronFlavour
//   를 보지 않고 pujet 을 돌려주므로 같은 우선순위를 쓴다.
//==============================================================
TString fake::FlavorTag(int partonFlavour, int hadronFlavour, bool isPileup) const {
    if (isPileup) return "pileup";
    if (hadronFlavour == 5) return "b";
    if (hadronFlavour == 4) return "c";
    const int ap = abs(partonFlavour);
    if (ap == 21) return "g";
    if (ap == 3)  return "s";
    if (ap == 2)  return "u";
    if (ap == 1)  return "d";
    return "unmatched";
}

//==== muon 이 clustering 된 PF jet 을 NanoAOD 의 jetIdx 포인터로 직접 찾는다.
//==== dR<0.4 최근접 매칭은 muon 이 실제로 들어있지 않은 이웃 jet 을 집을 수 있어서
//==== 공식 코드(MeasFakeRateV4::getMotherJetFlavour)와 같이 jetIdx 를 쓴다.
//==== rawJets 의 OriginalIndex 가 곧 NanoAOD Jet 인덱스라 그대로 대조하면 된다.
TString fake::RecoJetFlavor(const Muon &mu) const {
    const short jetIdx = mu.JetIdx();
    if (jetIdx < 0) return "unmatched";        // 연결된 jet 자체가 없는 muon
    for (const auto &jet : rawJets) {
        if (jet.OriginalIndex() != jetIdx) continue;
        return FlavorTag(jet.partonFlavour(), jet.hadronFlavour(),
                         jet.genJetIdx() < 0);
    }
    return "unmatched";                        // jetIdx 는 있으나 컬렉션에서 못 찾음
}

//==== muon 을 가장 가까운 gen jet 에 dR<0.4 매칭 (source jet 의 진짜 flavor)
TString fake::GenJetFlavor(const Muon &mu) const {
    float bestDR = 0.4f;
    int bestParton = 0, bestHadron = 0;
    bool matched = false;
    for (const auto &gjet : genJets) {
        const float dr = gjet.DeltaR(mu);
        if (dr < bestDR) {
            bestDR = dr;
            bestParton = gjet.partonFlavour();
            bestHadron = gjet.hadronFlavour();
            matched = true;
        }
    }
    return matched ? FlavorTag(bestParton, bestHadron, false) : "unmatched";
}

//==============================================================
// FlavorTag 문자열 -> 5-class 정수 코드 (SetBranch 는 문자열 branch 를 못 만든다)
//   0=b 1=c 2=uds(u+d+s) 3=g 4=pileup, -1=unmatched
//==============================================================
int fake::FlavorCode(const TString &flav) const {
    if (flav == "b") return 0;
    if (flav == "c") return 1;
    if (flav == "s" || flav == "u" || flav == "d") return 2;
    if (flav == "g") return 3;
    if (flav == "pileup") return 4;
    return -1;
}

//==============================================================
// flavor BDT 학습용 flat ntuple (userflag MakeTree, tree 이름 "mu")
//
// loose muon 이 정확히 1개인 측정 영역 이벤트마다 한 줄. FR 을 재는 모집단과
// 같은 population 이어야 하므로 selection 은 measureFakeRate 와 동일하되,
// away jet 컷만 가장 느슨한 30 GeV 로 열어 두고 hasAwayJet30/40/60 을
// 컬럼으로 남긴다 (offline 에서 Central=40 / syst=30,60 를 그대로 복원 가능).
//
// MET / MT 컷도 같은 이유로 **offline** 으로 옮겼다 (2026-08-11). 공식
// (MeasFakeRateV4) MC FR 은 Inclusive 영역 — away jet 만, MET/MT 컷 없음 — 에서
// 재므로, 컷을 트리에 박아 두면 그 맵을 재현할 수 없다. ev_met / ev_mt 컬럼이
// 그대로 있으므로 MeasReg (MET<25 && MT<25) 는 offline 에서 복원된다.
//
// 주의 1: 이 함수는 trigger path loop **밖**에서 불러야 한다. loop 안에서 쓰면
//         Mu8/Mu17 을 모두 통과한 muon 이 두 줄로 중복되어 train/test 폴드에
//         나뉘어 들어가고 검증 성능이 부풀려진다. 경로 정보는 pass_mu8 /
//         pass_mu17 / w_mu8 / w_mu17 컬럼으로 남긴다.
// 주의 2: AnalyzerCore::SetBranch 는 deque 에 값을 push 하고 그 원소의 주소를
//         TBranch 에 넘기는데 FillTrees 가 매번 clear()+shrink_to_fit() 한다.
//         따라서 **모든 branch 를 매 줄마다 빠짐없이 SetBranch 해야 한다.**
//         조건부로 건너뛰면 그 branch 는 해제된 메모리를 가리키게 된다.
//         Run2/Run3 전용 변수와 matched jet 이 없는 경우는 sentinel 로 채운다.
// 주의 3: DATA 에는 jet flavour 가 없어 Jet::partonFlavour/genJetIdx 가 ctor
//         기본값(-999)이고, FlavorTag 는 isPileup 을 먼저 보므로 모든 jet 이
//         "pileup" 으로 나온다. label 은 DATA 에서 -2 로 둔다.
//==============================================================
void fake::fillTreeRow(Event &ev, float evtSF) {

    if (!MakeTree) return;

    //==== 측정 영역과 동일: loose muon 1개 + electron 0개
    if (looseMuons.size() != 1) return;
    if (looseElectrons.size() != 0) return;

    const Muon &mu = looseMuons.at(0);

    //==== cone-corrected pT / MT (measureFakeRate 와 같은 정의)
    const float ptCorr = mu.Pt() * (1. + max(0.f, mu.MiniPFRelIso() - cuts.tight_miniiso_max));
    const float MT = sqrt(2. * mu.Pt() * METv.Pt() * (1. - cos(mu.DeltaPhi(METv))));

    //==== 트리거 경로별 통과 여부 (paths 는 userflag 로 잘릴 수 있어 직접 정의)
    const bool passMu8  = ev.PassTrigger("HLT_Mu8_TrkIsoVVL")
                          && mu.Pt() >= 10. && ptCorr >= 10.;
    const bool passMu17 = ev.PassTrigger("HLT_Mu17_TrkIsoVVL")
                          && mu.Pt() >= 20. && ptCorr >= 30.;
    if (!passMu8 && !passMu17) return;

    //==== away jet: 가장 느슨한 30 GeV 로 열어 두고 40/60 은 플래그로 남긴다
    int nJet25 = 0, nJet30 = 0, nJet40 = 0, nJet60 = 0;
    bool hasAway30 = false, hasAway40 = false, hasAway60 = false;
    float awayJetPt = -999., awayJetDR = -999., awayJetDPhi = -999., HT = 0.;
    for (const auto &jet : jets) {
        HT += jet.Pt();
        nJet25++;
        const bool away = (jet.DeltaR(mu) > cuts.awayjet_dr);
        if (jet.Pt() > 30.) { nJet30++; if (away) hasAway30 = true; }
        if (jet.Pt() > 40.) {
            nJet40++;
            if (away) {
                hasAway40 = true;
                if (awayJetPt < 0.) {   // jets 는 pT 내림차순 → 첫 번째가 leading
                    awayJetPt   = jet.Pt();
                    awayJetDR   = jet.DeltaR(mu);
                    awayJetDPhi = fabs(jet.DeltaPhi(mu));
                }
            }
        }
        if (jet.Pt() > 60.) { nJet60++; if (away) hasAway60 = true; }
    }
    if (!hasAway30) return;

    //==== 측정 영역의 MET / MT 컷 (W/Z prompt 오염 억제) 은 여기서 걸지 않는다.
    //==== ev_met / ev_mt 를 그대로 남겨 offline 에서 자르는 쪽이 넓다:
    //====   Inclusive (공식 MC FR)      = 컷 없음
    //====   MeasReg   (data FR / bkg)   = ev_met < 25 && ev_mt < 25

    //==== muon 의 source jet: lepton cleaning 이전 컬렉션에서 NanoAOD jetIdx 로 찾는다
    const short jetIdx = mu.JetIdx();
    const Jet *mj = nullptr;
    for (const auto &jet : rawJets) {
        if (jet.OriginalIndex() != jetIdx) continue;
        mj = &jet;
        break;
    }

    //==== matched jet 변수 (없으면 sentinel)
    const float FSENT = -999.;
    const int   ISENT = -99;
    float jetPt = FSENT, jetRawPt = FSENT, jetEta = FSENT, jetMass = FSENT;
    float jetArea = FSENT, jetDR = FSENT, jetPtRatio = FSENT, jetPtRel = FSENT;
    float jetChHEF = FSENT, jetNeHEF = FSENT, jetChEmEF = FSENT, jetNeEmEF = FSENT;
    float jetMuEF = FSENT, jetMuonSubtr = FSENT;
    float jetDeepB = FSENT, jetDeepCvB = FSENT, jetDeepCvL = FSENT, jetDeepQG = FSENT;
    float jetQGL = FSENT, jetPuIdDisc = FSENT;
    int   jetNConst = ISENT, jetNMuons = ISENT, jetNElectrons = ISENT;
    int   jetChMult = ISENT, jetNeMult = ISENT;   // Run3 전용 (아래 참조)
    bool  jetTightId = false, jetLoosePuId = false;
    if (mj) {
        jetPt        = mj->Pt();
        jetRawPt     = mj->GetRawPt();
        jetEta       = mj->Eta();
        jetMass      = mj->M();
        jetDR        = mj->DeltaR(mu);
        jetPtRatio   = (mj->Pt() > 0.) ? mu.Pt() / mj->Pt() : FSENT;
        //==== pTrel: (jet - muon) 축에 대한 muon 운동량의 수직 성분
        const TVector3 axis = mj->Vect() - mu.Vect();
        if (axis.Mag() > 1e-6) jetPtRel = mu.Vect().Perp(axis);
        jetChHEF     = mj->chHEF();
        jetNeHEF     = mj->neHEF();
        jetChEmEF    = mj->chEmEF();
        jetNeEmEF    = mj->neEmEF();
        jetMuEF      = mj->muEF();
        jetNConst    = mj->nConstituents();
        jetNMuons    = mj->nMuons();
        jetNElectrons= mj->nElectrons();
        jetDeepB     = mj->GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        const pair<float,float> cvals = mj->GetCTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        jetDeepCvB   = cvals.first;
        jetDeepCvL   = cvals.second;
        jetDeepQG    = mj->GetQvGTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        jetTightId   = mj->Pass_tightJetID();
        jetLoosePuId = mj->PassID("loosePuId");   // Run3 에서는 항상 참
        //==== charged / neutral multiplicity 는 Run3 (NanoAODv12+) 전용이다.
        //     Run2 v9 에는 Jet_chMultiplicity / Jet_neMultiplicity branch 자체가
        //     없고, AnalyzerCore 의 Run==3 가지에서만 SetHadronMultiplicities 가
        //     불린다. 게다가 Jet 의 생성자는 nConstituents 등과 달리 이 두 멤버를
        //     -999 로 초기화하지 않으므로 Run2 에서 읽으면 쓰레기 값이 나온다.
        //     따라서 반드시 Run==3 로 막고 나머지는 sentinel 로 남긴다
        //     (jet_qgl / jet_puIdDisc 의 Run2 전용 처리와 정확히 대칭).
        if (Run == 3) {
            jetChMult = int(mj->chMultiplicity());
            jetNeMult = int(mj->neMultiplicity());
        }
        if (jetIdx >= 0 && jetIdx < nJet) {
            jetArea      = Jet_area[jetIdx];
            jetMuonSubtr = Jet_muonSubtrFactor[jetIdx];
            if (Run == 2) {
                jetQGL      = mj->qgl();            // Run2 전용
                jetPuIdDisc = Jet_puIdDisc[jetIdx]; // Run2 전용
            }
        }
    }

    //==== source jet 안의 secondary vertex (IVF)
    //   b/c 분리의 고전적인 변수다: b hadron 은 M ~ 5 GeV, c hadron 은 ~ 2 GeV 이고
    //   b 쪽이 track 수도 많고 더 멀리 날아간다. DeepJet 이 내부적으로 SV 를 쓰지만
    //   출력은 3개 숫자로 압축돼 있어 raw SV 변수는 그 위에 정보를 더한다.
    //   jet 축 기준 dR < 0.4 로 매칭하고, 대표 SV 는 dlenSig 가 가장 큰 것
    //   (= 가장 유의미하게 변위된 vertex, HF 판별의 표준 선택) 을 쓴다.
    int   svN = 0, svNTracks = ISENT;
    float svMass = FSENT, svPt = FSENT, svPtRatio = FSENT, svDlenSig = FSENT;
    float svDxySig = FSENT, svPAngle = FSENT, svChi2 = FSENT;
    float svDRJet = FSENT, svDRMu = FSENT, svMassSum = FSENT;
    if (mj) {
        const float jetPhi = mj->Phi();
        float bestDlenSig = -1.;
        int   bestIdx = -1;
        float massSum = 0.;
        for (int i = 0; i < nSV; i++) {
            float dphi = SV_phi[i] - jetPhi;
            while (dphi >  M_PI) dphi -= 2. * M_PI;
            while (dphi < -M_PI) dphi += 2. * M_PI;
            const float deta = SV_eta[i] - jetEta;
            if (sqrt(deta * deta + dphi * dphi) > 0.4) continue;
            svN++;
            massSum += SV_mass[i];
            if (SV_dlenSig[i] > bestDlenSig) { bestDlenSig = SV_dlenSig[i]; bestIdx = i; }
        }
        svMassSum = massSum;
        if (bestIdx >= 0) {
            svMass     = SV_mass[bestIdx];
            svPt       = SV_pt[bestIdx];
            svPtRatio  = (jetPt > 0.) ? SV_pt[bestIdx] / jetPt : FSENT;
            svDlenSig  = SV_dlenSig[bestIdx];
            svDxySig   = SV_dxySig[bestIdx];
            svPAngle   = SV_pAngle[bestIdx];
            svChi2     = SV_chi2[bestIdx];
            svNTracks  = (int)SV_ntracks[bestIdx];
            float dphiJ = SV_phi[bestIdx] - jetPhi;
            while (dphiJ >  M_PI) dphiJ -= 2. * M_PI;
            while (dphiJ < -M_PI) dphiJ += 2. * M_PI;
            const float detaJ = SV_eta[bestIdx] - jetEta;
            svDRJet = sqrt(detaJ * detaJ + dphiJ * dphiJ);
            //==== muon 은 이 SV 에서 나온 것이므로 dR(SV, mu) 도 직접적인 handle
            float dphiM = SV_phi[bestIdx] - mu.Phi();
            while (dphiM >  M_PI) dphiM -= 2. * M_PI;
            while (dphiM < -M_PI) dphiM += 2. * M_PI;
            const float detaM = SV_eta[bestIdx] - mu.Eta();
            svDRMu = sqrt(detaM * detaM + dphiM * dphiM);
        }
    }

    //==== label (MC only). DATA 는 jet flavour 자체가 없으므로 -2
    int flavor = -2, flavorFine = -2, genFlavor = -2, leptonType = ISENT;
    if (!IsDATA) {
        const TString rf = RecoJetFlavor(mu);
        const TString gf = GenJetFlavor(mu);
        flavor    = FlavorCode(rf);
        genFlavor = FlavorCode(gf);
        //==== u/d/s 를 나중에 다시 쪼갤 수 있도록 세분 코드도 남긴다
        flavorFine = (rf == "b") ? 0 : (rf == "c") ? 1 : (rf == "s") ? 2
                   : (rf == "u") ? 3 : (rf == "d") ? 4 : (rf == "g") ? 5
                   : (rf == "pileup") ? 6 : -1;
        leptonType = GetLeptonType(mu, gens);
    }
    const bool isPrompt = (!IsDATA && leptonType > 0);

    //==== weight: 경로마다 trigger lumi 가 달라 하나로 합칠 수 없으므로 따로 남긴다
    float commonW = 1.;
    if (!IsDATA) {
        commonW = MCweight() * GetL1PrefireWeight()
                * myCorr->GetPUWeight(ev.nTrueInt()) * evtSF
                * myCorr->GetMuonRECOSF(mu);
    }
    const float wMu8  = passMu8  ? commonW * (IsDATA ? 1.f : ev.GetTriggerLumi("HLT_Mu8_TrkIsoVVL"))  : 0.f;
    const float wMu17 = passMu17 ? commonW * (IsDATA ? 1.f : ev.GetTriggerLumi("HLT_Mu17_TrkIsoVVL")) : 0.f;

    //==== 아래 SetBranch 블록은 반드시 무조건부 직선 코드여야 한다 (주의 2 참고)
    const TString t = "mu";

    //---- label / bookkeeping
    SetBranch(t, "flavor",      flavor);
    SetBranch(t, "flavor_fine", flavorFine);
    SetBranch(t, "genflavor",   genFlavor);
    SetBranch(t, "isTight",     PassMuonWP(mu, "tight"));
    SetBranch(t, "isPrompt",    isPrompt);
    SetBranch(t, "leptonType",  leptonType);
    SetBranch(t, "isData",      (bool)IsDATA);
    SetBranch(t, "pass_mu8",    passMu8);
    SetBranch(t, "pass_mu17",   passMu17);
    SetBranch(t, "w_mu8",       wMu8);
    SetBranch(t, "w_mu17",      wMu17);
    SetBranch(t, "mcweight",    IsDATA ? 1.f : (float)MCweight());
    SetBranch(t, "run",         ev.run());
    SetBranch(t, "lumi",        ev.lumi());
    SetBranch(t, "evt",         ev.event());

    //---- muon kinematics / IP / isolation (전부 data 에서도 얻을 수 있는 양)
    SetBranch(t, "mu_pt",          (float)mu.Pt());
    SetBranch(t, "mu_conept",      ptCorr);
    SetBranch(t, "mu_eta",         (float)mu.Eta());
    SetBranch(t, "mu_abseta",      (float)fabs(mu.Eta()));
    SetBranch(t, "mu_phi",         (float)mu.Phi());
    SetBranch(t, "mu_charge",      (int)mu.Charge());
    SetBranch(t, "mu_dxy",         mu.dXY());
    SetBranch(t, "mu_dxyerr",      mu.dXYerr());
    SetBranch(t, "mu_dxysig",      mu.dXY() / max(mu.dXYerr(), 1e-6f));
    SetBranch(t, "mu_dz",          mu.dZ());
    SetBranch(t, "mu_dzerr",       mu.dZerr());
    SetBranch(t, "mu_dzsig",       mu.dZ() / max(mu.dZerr(), 1e-6f));
    SetBranch(t, "mu_ip3d",        mu.IP3D());
    SetBranch(t, "mu_sip3d",       mu.SIP3D());
    SetBranch(t, "mu_miniiso",     mu.MiniPFRelIso());
    SetBranch(t, "mu_miniiso_chg", mu.MiniPFRelIsoChg());
    SetBranch(t, "mu_miniiso_neu", mu.MiniPFRelIso() - mu.MiniPFRelIsoChg());
    SetBranch(t, "mu_pfiso03",     mu.PfRelIso03());
    SetBranch(t, "mu_pfiso03_chg", mu.PfRelIso03Chg());
    SetBranch(t, "mu_pfiso04",     mu.PfRelIso04());
    SetBranch(t, "mu_tkreliso",    mu.TkRelIso());
    SetBranch(t, "mu_ptErr",       mu.PtErr());
    SetBranch(t, "mu_ptErrRel",    mu.PtErr() / max((float)mu.Pt(), 1e-6f));
    SetBranch(t, "mu_tunepRelPt",  mu.TunepRelPt());
    SetBranch(t, "mu_segmentComp", mu.SegmentComp());
    SetBranch(t, "mu_nStations",   mu.NStations());
    SetBranch(t, "mu_nTrkLayers",  mu.nTrackerLayers());
    SetBranch(t, "mu_mvaTTH",      mu.MvaTTH());
    SetBranch(t, "mu_mvaLowPt",    mu.MvaLowPt());
    SetBranch(t, "mu_softMva",     mu.SoftMva());
    SetBranch(t, "mu_tightCharge", mu.TightCharge());
    SetBranch(t, "mu_isGlobal",    mu.isGlobal());
    SetBranch(t, "mu_isTracker",   mu.isTracker());
    SetBranch(t, "mu_isStandalone",mu.isStandalone());
    SetBranch(t, "mu_isPFcand",    mu.isPFcand());
    SetBranch(t, "mu_highPurity",  mu.highPurity());
    SetBranch(t, "mu_pogMedium",   mu.isPOGMediumId());
    SetBranch(t, "mu_pogTight",    mu.isPOGTightId());

    //---- muon <-> jet 상관 변수 (jetIdx 가 없으면 NanoAOD 가 이미 fallback 값을 준다)
    SetBranch(t, "mu_jetPtRelv2",     mu.JetPtRelv2());
    SetBranch(t, "mu_jetRelIso",      mu.JetRelIso());
    SetBranch(t, "mu_jetNDauCharged", mu.JetNDauCharged());
    SetBranch(t, "mu_hasJet",         (bool)(mj != nullptr));

    //---- source jet 변수 (없으면 sentinel)
    SetBranch(t, "jet_pt",          jetPt);
    SetBranch(t, "jet_rawpt",       jetRawPt);
    SetBranch(t, "jet_eta",         jetEta);
    SetBranch(t, "jet_mass",        jetMass);
    SetBranch(t, "jet_area",        jetArea);
    SetBranch(t, "jet_dr_mu",       jetDR);
    SetBranch(t, "jet_ptratio",     jetPtRatio);
    SetBranch(t, "jet_ptrel",       jetPtRel);
    SetBranch(t, "jet_chHEF",       jetChHEF);
    SetBranch(t, "jet_neHEF",       jetNeHEF);
    SetBranch(t, "jet_chEmEF",      jetChEmEF);
    SetBranch(t, "jet_neEmEF",      jetNeEmEF);
    SetBranch(t, "jet_muEF",        jetMuEF);
    SetBranch(t, "jet_muonSubtr",   jetMuonSubtr);
    SetBranch(t, "jet_nconst",      jetNConst);
    SetBranch(t, "jet_chmult",      jetChMult);     // Run3 전용, Run2 는 sentinel
    SetBranch(t, "jet_nemult",      jetNeMult);     // Run3 전용, Run2 는 sentinel
    SetBranch(t, "jet_nmuons",      jetNMuons);
    SetBranch(t, "jet_nelectrons",  jetNElectrons);
    SetBranch(t, "jet_deepjet_b",   jetDeepB);
    SetBranch(t, "jet_deepjet_cvb", jetDeepCvB);
    SetBranch(t, "jet_deepjet_cvl", jetDeepCvL);
    SetBranch(t, "jet_deepjet_qg",  jetDeepQG);
    SetBranch(t, "jet_qgl",         jetQGL);        // Run2 전용, Run3 는 sentinel
    SetBranch(t, "jet_puIdDisc",    jetPuIdDisc);   // Run2 전용, Run3 는 sentinel
    SetBranch(t, "jet_tightId",     jetTightId);
    SetBranch(t, "jet_loosePuId",   jetLoosePuId);

    //---- source jet 안의 secondary vertex (매칭 SV 가 없으면 sentinel)
    SetBranch(t, "sv_n",         svN);
    SetBranch(t, "sv_mass",      svMass);
    SetBranch(t, "sv_masssum",   svMassSum);
    SetBranch(t, "sv_pt",        svPt);
    SetBranch(t, "sv_ptratio",   svPtRatio);
    SetBranch(t, "sv_dlensig",   svDlenSig);
    SetBranch(t, "sv_dxysig",    svDxySig);
    SetBranch(t, "sv_pangle",    svPAngle);
    SetBranch(t, "sv_chi2",      svChi2);
    SetBranch(t, "sv_ntracks",   svNTracks);
    SetBranch(t, "sv_dr_jet",    svDRJet);
    SetBranch(t, "sv_dr_mu",     svDRMu);

    //---- event / away jet context
    SetBranch(t, "ev_met",          (float)METv.Pt());
    SetBranch(t, "ev_mt",           MT);
    SetBranch(t, "ev_dphi_mu_met",  (float)fabs(mu.DeltaPhi(METv)));
    SetBranch(t, "ev_ht",           HT);
    SetBranch(t, "ev_njet25",       nJet25);
    SetBranch(t, "ev_njet30",       nJet30);
    SetBranch(t, "ev_njet40",       nJet40);
    SetBranch(t, "ev_njet60",       nJet60);
    SetBranch(t, "ev_hasAway30",    hasAway30);
    SetBranch(t, "ev_hasAway40",    hasAway40);
    SetBranch(t, "ev_hasAway60",    hasAway60);
    SetBranch(t, "ev_awayjet_pt",   awayJetPt);
    SetBranch(t, "ev_awayjet_dr",   awayJetDR);
    SetBranch(t, "ev_awayjet_dphi", awayJetDPhi);
    SetBranch(t, "ev_nPV",          ev.nPV());
    SetBranch(t, "ev_nPVsGood",     ev.nPVsGood());
    SetBranch(t, "ev_rho",          fixedGridRhoFastjetAll);
    SetBranch(t, "ev_nTrueInt",     IsDATA ? -1.f : ev.nTrueInt());

    FillTrees(t);
}
