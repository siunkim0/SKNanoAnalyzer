#include "elecfake.h"

//==== Constructor and Destructor
elecfake::elecfake() : MeasFakeEle8(false), MeasFakeEle12(false), MeasFakeEle23(false),
                       MakeTree(false) {}
elecfake::~elecfake() {}

//==== Initialize variables
void elecfake::initializeAnalyzer() {

    //==== Userflags
    MeasFakeEle8  = HasFlag("MeasFakeEle8");    // Ele8 트리거 경로만 측정
    MeasFakeEle12 = HasFlag("MeasFakeEle12");   // Ele12 트리거 경로만 측정
    MeasFakeEle23 = HasFlag("MeasFakeEle23");   // Ele23 트리거 경로만 측정
    MakeTree      = HasFlag("MakeTree");        // neural fake factor 학습용 flat ntuple

    //==== 플래그가 없으면 세 경로 모두 측정
    if (!MeasFakeEle8 && !MeasFakeEle12 && !MeasFakeEle23) {
        MeasFakeEle8 = true;
        MeasFakeEle12 = true;
        MeasFakeEle23 = true;
    }

    //==== prescaled single electron trigger 경로 정의 (analysis note Table 38)
    //==== {히스토그램 prefix, HLT 이름, electron pT 컷, cone-corrected pT 컷}
    if (MeasFakeEle8)  paths.push_back({"Ele8",  "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",  10., 10.});
    if (MeasFakeEle12) paths.push_back({"Ele12", "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", 15., 20.});
    if (MeasFakeEle23) paths.push_back({"Ele23", "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", 25., 35.});

    //==== systematic 목록: Central = away jet pT > 40 GeV
    //==== away jet pT variation 은 2D 히스토그램 몇 개만 더 채우는 것이라
    //==== 부담이 없으므로 항상 켠다
    systs = {"Central", "AwayJetPt30", "AwayJetPt60"};

    //==== fake rate 2D 히스토그램 binning: cone-corrected pT x |scEta|
    //==== pT 경계는 트리거 stitching 지점 (20, 35) 을 포함해야 한다
    //==== 공식 측정처럼 [100,200] 을 따로 둔다 (FR 파일을 만들 때는 100 까지만 사용)
    ptCorrBins = {15., 17., 20., 25., 35., 50., 100., 200.};
    absEtaBins = {0., 0.8, 1.479, 2.5};

    //==== loose WP 의 raw MVANoIso 컷: etaRegion (IB, OB, EC) 별,
    //==== Run 2 / Run 3 에서 MVA 가 retraining 되어 값이 다르다
    //==== 주의: Electron::Pass_HcToWALooseRun2() 는 이 컷의 부호가 반전된
    //====       버그가 있어 PassElectronWP 에서 직접 구현한다
    if (Run == 2) { looseMvaIB = 0.985; looseMvaOB = 0.96; looseMvaEC = 0.85; }
    else          { looseMvaIB = 0.8;   looseMvaOB = 0.5;  looseMvaEC = -0.8; }

    //==== electron loose WP 도 Run 별로 다르다 (공식
    //==== Electron::Pass_HcToWALooseRun2 SIP3D<8 vs ...Run3 SIP3D<6;
    //==== miniiso 는 둘 다 0.4, MVA 컷은 위에서 이미 Run 별로 잡았다).
    //==== tight (HcToWATight) 는 Run 무관하다.
    //==== veto muon 도 같은 이유로 Run 별 loose WP 을 쓴다
    //==== (Muon::Pass_HcToWALooseRun2 SIP3D<5, miniiso<0.6
    //====  -> ...Run3 SIP3D<8, miniiso<0.4)
    if (Run == 3) {
        cuts.loose_sip3d_max  = 6.0;
        cuts.muon_sip3d_max   = 8.0;
        cuts.muon_miniiso_max = 0.4;
    }

    //==== era 별 jet |eta| 컷
    if (DataEra == "2016preVFP" || DataEra == "2016postVFP") {
        cuts.jet_eta_max = 2.4;
    } else {
        cuts.jet_eta_max = 2.5;
    }

    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);

    if (MakeTree) NewTree("el");
}

//==============================================================
// Electron WP (electron ID 표)
// 공통 (baseline) : trigger emulation (CaloIdL_TrackIdL_IsoVL: sieie, dEtaIn,
//                   dPhiIn, H/E, ECAL/HCAL PF cluster iso (EA), tracker iso;
//                   GAP 영역 제외) + conversion veto + missing inner hits <= 1
//                   + |dz| < 0.1 cm
// tight           : MVANoIso POG wp90, |SIP3D| < 4, rel. mini-iso < 0.1
// loose           : raw MVANoIso > (IB, OB, EC) 컷 (또는 wp90 통과),
//                   |SIP3D| < 8, rel. mini-iso < 0.4
//==============================================================
bool elecfake::PassElectronWP(const Electron &el, const TString &wp) const {
    //==== 공통 컷
    if (!el.Pass_CaloIdL_TrackIdL_IsoVL()) return false;   // trigger emulation
    if (!el.ConvVeto()) return false;
    if (!(el.LostHits() <= cuts.electron_losthits_max)) return false;
    if (!(fabs(el.dZ()) < cuts.electron_dz_max)) return false;

    //==== tight / loose 별 컷
    if (wp == "tight") {
        if (!el.isMVANoIsoWP90()) return false;
        if (!(fabs(el.SIP3D()) < cuts.tight_sip3d_max)) return false;
        if (!(el.MiniPFRelIso() < cuts.tight_miniiso_max)) return false;
    } else {
        //==== tight ⊂ loose 를 보장하기 위해 wp90 통과도 loose MVA 컷으로 인정
        float mvaCut = looseMvaEC;
        if (el.etaRegion() == Electron::ETAREGION::IB)      mvaCut = looseMvaIB;
        else if (el.etaRegion() == Electron::ETAREGION::OB) mvaCut = looseMvaOB;
        if (!(el.isMVANoIsoWP90() || el.MvaNoIso() > mvaCut)) return false;
        if (!(fabs(el.SIP3D()) < cuts.loose_sip3d_max)) return false;
        if (!(el.MiniPFRelIso() < cuts.loose_miniiso_max)) return false;
    }
    return true;
}

//==============================================================
// Muon veto: fake.cc 의 loose muon WP 와 동일
// (POG medium ID, |dz| < 0.1, tkRelIso < 0.4, |SIP3D| < 5, mini-iso < 0.6)
//==============================================================
bool elecfake::PassMuonVeto(const Muon &mu) const {
    if (mu.Pt() < cuts.muon_pt_min || fabs(mu.Eta()) > cuts.muon_eta_max) return false;
    if (!mu.isPOGMediumId()) return false;
    if (!(fabs(mu.dZ()) < cuts.muon_dz_max)) return false;
    if (!(mu.TkRelIso() < cuts.muon_tkiso_max)) return false;
    if (!(fabs(mu.SIP3D()) < cuts.muon_sip3d_max)) return false;
    if (!(mu.MiniPFRelIso() < cuts.muon_miniiso_max)) return false;
    return true;
}

//==============================================================
// jet 단위 veto map: 공식 framework 의 로직을 직접 구현 (fake.cc 와 동일)
// EM fraction > 0.9 또는 muon dR < 0.2 인 jet 은 검사 대상에서 제외
//==============================================================
bool elecfake::PassVetoMapJet(const Jet &jet, const RVec<Muon> &muons) const {
    if (jet.chEmEF() + jet.neEmEF() > 0.9) return true;
    for (const auto &mu : muons)
        if (jet.DeltaR(mu) < 0.2) return true;
    return !myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), "jetvetomap");
}

void elecfake::executeEvent() {

    //==== 이벤트와 raw physics object 읽기
    Event ev = GetEvent();
    rawJets = GetAllJets();   // flavor 매칭용으로 lepton cleaning 이전 컬렉션을 보관

    //==== noise (MET) filter 통과 요구
    if (!PassNoiseFilter(rawJets, ev)) return;

    RVec<Muon> rawMuons = GetAllMuons();

    //==== jet veto map (공식 측정과 동일하게 event 단위로 적용; Run2 는 no-op)
    if (!PassVetoMap(rawJets, rawMuons, "jetvetomap")) return;

    RVec<Electron> rawElectrons = GetAllElectrons();
    sort(rawElectrons.begin(), rawElectrons.end(), PtComparing);

    //==== Step 1: veto / loose / tight electron 선택 (tight ⊂ loose ⊂ veto)
    //====         공식 측정처럼 veto 는 pT > 10, 측정(loose) 은 pT > 15;
    //====         object |eta| 컷은 공식(SelectElectrons)과 같이 일반 eta 를 쓰고
    //====         FR binning 만 supercluster eta 를 쓴다
    vetoElectrons.clear();
    looseElectrons.clear();
    tightElectrons.clear();
    for (const auto &el : rawElectrons) {
        if (fabs(el.Eta()) > cuts.electron_eta_max) continue;
        if (el.Pt() < cuts.electron_veto_pt_min) continue;
        if (!PassElectronWP(el, "loose")) continue;
        vetoElectrons.push_back(el);
        if (el.Pt() < cuts.electron_pt_min) continue;   // 측정 대상은 pT > 15
        looseElectrons.push_back(el);
        if (PassElectronWP(el, "tight")) tightElectrons.push_back(el);
    }

    //==== Step 2: muon veto 선택 (pT > 10, loose muon WP)
    vetoMuons.clear();
    for (const auto &mu : rawMuons) {
        if (PassMuonVeto(mu)) vetoMuons.push_back(mu);
    }

    //==== Step 3: jet 선택 (공식 측정과 같은 순서)
    //====         tight ID → veto lepton (pT>10 electron + muon) 과 dR<0.4 제거
    //====         → jet 단위 veto map → loose PU ID (Run2)
    jets = SelectJets(rawJets, "tight", cuts.jet_pt_min, cuts.jet_eta_max);
    jets = JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, cuts.jet_lep_dr);
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

    //==== Step 6: 트리거 경로(Ele8 / Ele12 / Ele23)마다 측정
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

    //==== Step 7: neural fake factor 학습용 ntuple (userflag MakeTree)
    //====         반드시 trigger path loop **밖**에서 부른다 (fillTreeRow 주의 1)
    fillTreeRow(ev, evtSF);
}

//==============================================================
// Fake rate 측정 영역
// selection: loose electron 딱 1개 + muon 0개
//            + away jet (pT > 40, dR > 0.7) + MET < 25 + MT < 25
// FR(ptcorr, |scEta|) = tight 히스토그램 / loose 히스토그램 (offline 에서 계산)
//==============================================================
void elecfake::measureFakeRate(const TriggerPath &path, float weight) {

    //==== Step 1: single lepton 이벤트 요구
    //====         측정 electron 1개 + soft(pT>10) 전자도 그 1개뿐 + muon 0개
    if (looseElectrons.size() != 1) return;
    if (vetoElectrons.size() != 1) return;
    if (vetoMuons.size() != 0) return;

    const Electron &el = looseElectrons.at(0);

    //==== electron RECO SF (공식 측정과 동일하게 측정 대상 전자에 적용)
    if (!IsDATA) weight *= myCorr->GetElectronRECOSF(looseElectrons);

    //==== Step 2: 트리거별 electron pT 컷 (Ele8: 10, Ele12: 15, Ele23: 25)
    if (el.Pt() < path.ptCut) return;

    //==== Step 3: cone-corrected pT = pT * (1 + max(0, miniIso - 0.1))
    //====         tight iso 기준(0.1)을 넘는 만큼 pT 를 보정해 준다
    const float ptCorr = el.Pt() * (1. + max(0.f, el.MiniPFRelIso() - cuts.tight_miniiso_max));
    if (ptCorr < path.ptCorrCut) return;   // Ele8: 10, Ele12: 20, Ele23: 35

    //==== tight WP 통과 여부, prompt/fake 구분 (MC only)
    const bool isTight = PassElectronWP(el, "tight");
    const TString ltype = IsDATA ? "" : LeptonTypeToString(GetLeptonType(el, gens));

    //==== transverse mass MT(el, MET)
    const float MT = sqrt(2. * el.Pt() * METv.Pt() * (1. - cos(el.DeltaPhi(METv))));

    //==== loose 는 항상 채우고, tight 통과 시 tight 도 채운다 → FR = tight / loose
    RVec<TString> tags = {"loose"};
    if (isTight) tags.push_back("tight");

    //==== overflow 는 마지막 bin 에 넣는다 (binning 은 |scEta|)
    const float xval = min(ptCorr, ptCorrBins.back() - 0.1f);
    const float yval = fabs(el.scEta());

    //==== away jet pT 컷을 바꿔가며 반복 (Central = 40, syst = 30 / 60)
    for (const auto &syst : systs) {
        float awayJetPtCut = cuts.awayjet_pt;
        if (syst == "AwayJetPt30") awayJetPtCut = 30.;
        else if (syst == "AwayJetPt60") awayJetPtCut = 60.;

        //==== Step 4: away jet 요구 (pT > cut && dR(el, jet) > 0.7)
        int nJets = 0;
        bool hasAwayJet = false;
        for (const auto &jet : jets) {
            if (jet.Pt() < awayJetPtCut) continue;
            nJets++;
            if (jet.DeltaR(el) > cuts.awayjet_dr) hasAwayJet = true;
        }
        if (!hasAwayJet) continue;

        //==== Inclusive (validation) 영역: MET / MT 컷 없이 kinematics + 2D FR
        //==== 공식(MeasFakeRateV4) QCD MC FR 은 이 영역에서 측정되므로
        //==== 비교용 2D FR 히스토그램도 같이 채운다
        if (syst == "Central") {
            for (const auto &tag : tags) {
                fillElectronKinematics(path.name + "/Inclusive/" + tag,
                                       el, ptCorr, MT, nJets, weight);
                FillHist(path.name + "/Inclusive/" + tag + "/fake_ptcorr_abseta",
                         xval, yval, weight, ptCorrBins, absEtaBins);
                if (!IsDATA) {
                    fillElectronKinematics(path.name + "/Inclusive/" + tag + "/" + ltype,
                                           el, ptCorr, MT, nJets, weight);
                    FillHist(path.name + "/Inclusive/" + tag + "/" + ltype + "/fake_ptcorr_abseta",
                             xval, yval, weight, ptCorrBins, absEtaBins);
                }
            }
        }

        //==== Step 5: 측정 영역 컷 (W/Z 의 prompt electron 오염 억제)
        if (!(METv.Pt() < cuts.met_max)) continue;
        if (!(MT < cuts.mt_max)) continue;

        //==== FR 계산용 2D 히스토그램 (cone-corrected pT x |scEta|)
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
        //==== fake electron 만 (GetLeptonType<=0): TTLJ 등 prompt 오염 표본에서 필수
        if (!IsDATA && syst == "Central" && GetLeptonType(el, gens) <= 0) {
            const TString gf = GenJetFlavor(el);
            const TString rf = RecoJetFlavor(el);
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
                fillElectronKinematics(base + "/" + tag, el, ptCorr, MT, nJets, weight);
                if (!IsDATA)
                    fillElectronKinematics(base + "/" + tag + "/" + ltype,
                                           el, ptCorr, MT, nJets, weight);
            }
        }
    }
}

//==============================================================
// Z-enriched 영역: prescaled trigger 경로별 MC normalization 추출용
// selection: OS tight electron pair + |M(ll) - 91.2| < 15 + jet pT > 40 하나 이상
//==============================================================
void elecfake::fillZEnriched(const TriggerPath &path, const Event &ev, float weight) {

    //==== Step 1: loose 2개 == tight 2개 (딱 tight pair 만 있는 이벤트)
    if (!(looseElectrons.size() == 2 && tightElectrons.size() == 2)) return;
    if (!vetoMuons.empty()) return;

    const Electron &el1 = tightElectrons.at(0);
    const Electron &el2 = tightElectrons.at(1);

    //==== electron RECO SF (공식 측정과 동일하게 두 전자 모두에 적용)
    if (!IsDATA) weight *= myCorr->GetElectronRECOSF(tightElectrons);

    //==== Step 2: opposite sign + pT 컷 (leading 은 트리거별 컷, subleading 은 10)
    if (el1.Charge() + el2.Charge() != 0) return;
    if (el1.Pt() < path.ptCut || el2.Pt() < cuts.electron_pt_min) return;

    //==== Step 3: jet pT > 40 하나 이상 (측정 영역과 같은 topology 맞추기)
    bool hasJet = false;
    for (const auto &jet : jets) {
        if (jet.Pt() > cuts.zjet_pt) { hasJet = true; break; }
    }
    if (!hasJet) return;

    //==== Step 4: Z mass window |M(ll) - 91.2| < 15
    const float mll = (el1 + el2).M();
    if (fabs(mll - cuts.z_mass) > cuts.z_window) return;

    //==== MC 는 두 electron 모두 prompt 인지에 따라 prompt / fake 로도 나눠 채운다
    RVec<TString> prefixes = {path.name + "/ZEnriched"};
    if (!IsDATA) {
        const bool isPrompt = (GetLeptonType(el1, gens) > 0) && (GetLeptonType(el2, gens) > 0);
        prefixes.push_back(prefixes[0] + (isPrompt ? "/prompt" : "/fake"));
    }

    for (const auto &prefix : prefixes) {
        FillHist(prefix + "/mll", mll, weight, 60, 76.2, 106.2);
        FillHist(prefix + "/el1_pt", el1.Pt(), weight, 200, 0., 200.);
        FillHist(prefix + "/el2_pt", el2.Pt(), weight, 200, 0., 200.);
        FillHist(prefix + "/el1_eta", el1.scEta(), weight, 50, -2.5, 2.5);
        FillHist(prefix + "/el2_eta", el2.scEta(), weight, 50, -2.5, 2.5);
        FillHist(prefix + "/nPV", ev.nPV(), weight, 70, 0., 70.);
    }
}

//==============================================================
// electron kinematics 히스토그램 묶음 (prefix 아래 7개)
//==============================================================
void elecfake::fillElectronKinematics(const TString &prefix, const Electron &el,
                                      float ptCorr, float MT, int nJets, float weight) {
    FillHist(prefix + "/pt", el.Pt(), weight, 200, 0., 200.);
    FillHist(prefix + "/ptcorr", ptCorr, weight, 200, 0., 200.);
    FillHist(prefix + "/eta", el.scEta(), weight, 50, -2.5, 2.5);
    FillHist(prefix + "/phi", el.Phi(), weight, 64, -3.2, 3.2);
    FillHist(prefix + "/MET", METv.Pt(), weight, 100, 0., 100.);
    FillHist(prefix + "/MT", MT, weight, 100, 0., 100.);
    FillHist(prefix + "/nJets", nJets, weight, 10, 0., 10.);
}

//==============================================================
// GetLeptonType > 0 : prompt (EW/BSM prompt, tau daughter, internal conversion)
// GetLeptonType <= 0: hadron 기원 / external conversion / unmatched → fake
//==============================================================
TString elecfake::LeptonTypeToString(int leptonType) const {
    return (leptonType > 0) ? "prompt" : "fake";
}

//==============================================================
// Source jet parton flavor 분류 (fake.cc 와 동일 로직)
//   hadronFlavour (ghost B/C hadron matching) 로 b/c 를 먼저 잡고,
//   light (hadronFlavour==0) 는 |partonFlavour| 로 g/s/u/d 세분화한다.
//   partonFlavour: 21=gluon, 1=d, 2=u, 3=s (부호 = quark/antiquark)
//   매칭 실패는 unmatched
//   isPileup (대응하는 gen jet 이 없는 reco jet) 은 flavour 보다 먼저 본다:
//   공식 코드 MeasFakeRateV4::getMotherJetFlavour 가 genJetIdx<0 이면 hadronFlavour
//   를 보지 않고 pujet 을 돌려주므로 같은 우선순위를 쓴다.
//==============================================================
TString elecfake::FlavorTag(int partonFlavour, int hadronFlavour, bool isPileup) const {
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

//==== electron 이 clustering 된 PF jet 을 NanoAOD 의 jetIdx 포인터로 직접 찾는다.
//==== dR<0.4 최근접 매칭은 electron 이 실제로 들어있지 않은 이웃 jet 을 집을 수 있어서
//==== 공식 코드(MeasFakeRateV4::getMotherJetFlavour)와 같이 jetIdx 를 쓴다.
//==== rawJets 의 OriginalIndex 가 곧 NanoAOD Jet 인덱스라 그대로 대조하면 된다.
TString elecfake::RecoJetFlavor(const Electron &el) const {
    const short jetIdx = el.JetIdx();
    if (jetIdx < 0) return "unmatched";        // 연결된 jet 자체가 없는 electron
    for (const auto &jet : rawJets) {
        if (jet.OriginalIndex() != jetIdx) continue;
        return FlavorTag(jet.partonFlavour(), jet.hadronFlavour(),
                         jet.genJetIdx() < 0);
    }
    return "unmatched";                        // jetIdx 는 있으나 컬렉션에서 못 찾음
}

//==== electron 을 가장 가까운 gen jet 에 dR<0.4 매칭 (source jet 의 진짜 flavor)
TString elecfake::GenJetFlavor(const Electron &el) const {
    float bestDR = 0.4f;
    int bestParton = 0, bestHadron = 0;
    bool matched = false;
    for (const auto &gjet : genJets) {
        const float dr = gjet.DeltaR(el);
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
//   0=b 1=c 2=uds(u+d+s) 3=g 4=pileup, -1=unmatched  (fake.cc 와 동일한 코드표)
//==============================================================
int elecfake::FlavorCode(const TString &flav) const {
    if (flav == "b") return 0;
    if (flav == "c") return 1;
    if (flav == "s" || flav == "u" || flav == "d") return 2;
    if (flav == "g") return 3;
    if (flav == "pileup") return 4;
    return -1;
}

//==============================================================
// neural fake factor 학습용 flat ntuple (userflag MakeTree, tree 이름 "el")
//
// fake.cc::fillTreeRow 의 electron 판. 규칙과 함정이 전부 같으므로 그쪽 주석을
// 먼저 읽을 것. electron 이라서 달라지는 점만 여기 적는다:
//
//   - 측정 영역은 measureFakeRate 와 같다: loose electron 1개 + soft(pT>10)
//     전자도 그 1개뿐 + muon 0개. veto 전자 조건을 빼먹으면 FR 분모의 모집단이
//     달라진다.
//   - 트리거가 셋 (Ele8/Ele12/Ele23) 이고 stitch 지점이 ptCorr 20 / 35 다.
//     경로별 pT · ptCorr 컷은 initializeAnalyzer 의 paths 와 같은 값을 직접 쓴다
//     (paths 는 userflag 로 잘릴 수 있다).
//   - binning 축은 |scEta| 다 (el_abssceta). object |eta| 컷만 일반 eta 를 쓴다.
//   - MET / MT 컷은 여기서 걸지 않는다. 공식 QCD MC FR 이 Inclusive (away jet
//     만) 에서 측정되므로 컷을 트리에 박으면 그 맵을 재현할 수 없다.
//     measureFakeRate 은 여전히 걸므로 히스토그램 쪽은 영향이 없다.
//   - Electron 은 MiniPFRelIsoChg / PtErr / SegmentComp 등의 accessor 가 없고,
//     Electron_jetPtRelv2 / jetRelIso / jetNDauCharged 는 loader 에는 있지만
//     Electron 객체로 옮겨지지 않는다. 학습에 쓰지 않는 양이라 DataFormats 를
//     건드리지 않고 그냥 뺐다.
//==============================================================
void elecfake::fillTreeRow(Event &ev, float evtSF) {

    if (!MakeTree) return;

    //==== 측정 영역과 동일: loose electron 1개 + veto electron 도 그 1개 + muon 0개
    if (looseElectrons.size() != 1) return;
    if (vetoElectrons.size() != 1) return;
    if (vetoMuons.size() != 0) return;

    const Electron &el = looseElectrons.at(0);

    //==== cone-corrected pT / MT (measureFakeRate 와 같은 정의)
    const float ptCorr = el.Pt() * (1. + max(0.f, el.MiniPFRelIso() - cuts.tight_miniiso_max));
    const float MT = sqrt(2. * el.Pt() * METv.Pt() * (1. - cos(el.DeltaPhi(METv))));

    //==== 트리거 경로별 통과 여부 (paths 는 userflag 로 잘릴 수 있어 직접 정의)
    const bool passEle8  = ev.PassTrigger("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30")
                           && el.Pt() >= 10. && ptCorr >= 10.;
    const bool passEle12 = ev.PassTrigger("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30")
                           && el.Pt() >= 15. && ptCorr >= 20.;
    const bool passEle23 = ev.PassTrigger("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30")
                           && el.Pt() >= 25. && ptCorr >= 35.;
    if (!passEle8 && !passEle12 && !passEle23) return;

    //==== away jet: 가장 느슨한 30 GeV 로 열어 두고 40/60 은 플래그로 남긴다
    int nJet25 = 0, nJet30 = 0, nJet40 = 0, nJet60 = 0;
    bool hasAway30 = false, hasAway40 = false, hasAway60 = false;
    float awayJetPt = -999., awayJetDR = -999., awayJetDPhi = -999., HT = 0.;
    for (const auto &jet : jets) {
        HT += jet.Pt();
        nJet25++;
        const bool away = (jet.DeltaR(el) > cuts.awayjet_dr);
        if (jet.Pt() > 30.) { nJet30++; if (away) hasAway30 = true; }
        if (jet.Pt() > 40.) {
            nJet40++;
            if (away) {
                hasAway40 = true;
                if (awayJetPt < 0.) {   // jets 는 pT 내림차순 → 첫 번째가 leading
                    awayJetPt   = jet.Pt();
                    awayJetDR   = jet.DeltaR(el);
                    awayJetDPhi = fabs(jet.DeltaPhi(el));
                }
            }
        }
        if (jet.Pt() > 60.) { nJet60++; if (away) hasAway60 = true; }
    }
    if (!hasAway30) return;

    //==== electron 의 source jet: lepton cleaning 이전 컬렉션에서 NanoAOD jetIdx 로 찾는다
    //====                        (RecoJetFlavor 와 같은 매칭, 여기서는 Jet 을 들고 있는다)
    const short jetIdx = el.JetIdx();
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
    bool  jetTightId = false, jetLoosePuId = false;
    if (mj) {
        jetPt        = mj->Pt();
        jetRawPt     = mj->GetRawPt();
        jetEta       = mj->Eta();
        jetMass      = mj->M();
        jetDR        = mj->DeltaR(el);
        jetPtRatio   = (mj->Pt() > 0.) ? el.Pt() / mj->Pt() : FSENT;
        //==== pTrel: (jet - electron) 축에 대한 electron 운동량의 수직 성분
        const TVector3 axis = mj->Vect() - el.Vect();
        if (axis.Mag() > 1e-6) jetPtRel = el.Vect().Perp(axis);
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
        if (jetIdx >= 0 && jetIdx < nJet) {
            jetArea      = Jet_area[jetIdx];
            jetMuonSubtr = Jet_muonSubtrFactor[jetIdx];
            if (Run == 2) {
                jetQGL      = mj->qgl();            // Run2 전용
                jetPuIdDisc = Jet_puIdDisc[jetIdx]; // Run2 전용
            }
        }
    }

    //==== source jet 안의 secondary vertex (IVF) — fake.cc 와 같은 jet 축 매칭
    //   b/c 분리의 고전적인 변수다: b hadron 은 M ~ 5 GeV, c hadron 은 ~ 2 GeV.
    //   jet 축 기준 dR < 0.4, 대표 SV 는 dlenSig 가 가장 큰 것.
    int   svN = 0, svNTracks = ISENT;
    float svMass = FSENT, svPt = FSENT, svPtRatio = FSENT, svDlenSig = FSENT;
    float svDxySig = FSENT, svPAngle = FSENT, svChi2 = FSENT;
    float svDRJet = FSENT, svDREl = FSENT, svMassSum = FSENT;
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
            //==== electron 은 이 SV 에서 나온 것이므로 dR(SV, el) 도 직접적인 handle
            float dphiE = SV_phi[bestIdx] - el.Phi();
            while (dphiE >  M_PI) dphiE -= 2. * M_PI;
            while (dphiE < -M_PI) dphiE += 2. * M_PI;
            const float detaE = SV_eta[bestIdx] - el.Eta();
            svDREl = sqrt(detaE * detaE + dphiE * dphiE);
        }
    }

    //==== label (MC only). DATA 는 jet flavour 자체가 없으므로 -2
    int flavor = -2, flavorFine = -2, genFlavor = -2, leptonType = ISENT;
    if (!IsDATA) {
        const TString rf = RecoJetFlavor(el);
        const TString gf = GenJetFlavor(el);
        flavor    = FlavorCode(rf);
        genFlavor = FlavorCode(gf);
        //==== u/d/s 를 나중에 다시 쪼갤 수 있도록 세분 코드도 남긴다
        flavorFine = (rf == "b") ? 0 : (rf == "c") ? 1 : (rf == "s") ? 2
                   : (rf == "u") ? 3 : (rf == "d") ? 4 : (rf == "g") ? 5
                   : (rf == "pileup") ? 6 : -1;
        leptonType = GetLeptonType(el, gens);
    }
    const bool isPrompt = (!IsDATA && leptonType > 0);

    //==== weight: 경로마다 trigger lumi 가 달라 하나로 합칠 수 없으므로 따로 남긴다
    //====         measureFakeRate 과 같은 세트 (+ electron RECO SF)
    float commonW = 1.;
    if (!IsDATA) {
        commonW = MCweight() * GetL1PrefireWeight()
                * myCorr->GetPUWeight(ev.nTrueInt()) * evtSF
                * myCorr->GetElectronRECOSF(looseElectrons);
    }
    const float wEle8  = passEle8  ? commonW * (IsDATA ? 1.f : ev.GetTriggerLumi("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30"))  : 0.f;
    const float wEle12 = passEle12 ? commonW * (IsDATA ? 1.f : ev.GetTriggerLumi("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30")) : 0.f;
    const float wEle23 = passEle23 ? commonW * (IsDATA ? 1.f : ev.GetTriggerLumi("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30")) : 0.f;

    //==== 아래 SetBranch 블록은 반드시 무조건부 직선 코드여야 한다
    //====   (SetBranch 는 deque 에 push 하고 그 원소의 주소를 TBranch 에 넘기는데
    //====    FillTrees 가 매번 clear() 하므로, 조건부로 건너뛴 branch 는 해제된
    //====    메모리를 가리키게 된다. 없는 값은 sentinel 로 채운다.)
    const TString t = "el";

    //---- label / bookkeeping
    SetBranch(t, "flavor",      flavor);
    SetBranch(t, "flavor_fine", flavorFine);
    SetBranch(t, "genflavor",   genFlavor);
    SetBranch(t, "isTight",     PassElectronWP(el, "tight"));
    SetBranch(t, "isPrompt",    isPrompt);
    SetBranch(t, "leptonType",  leptonType);
    SetBranch(t, "isData",      (bool)IsDATA);
    SetBranch(t, "pass_ele8",   passEle8);
    SetBranch(t, "pass_ele12",  passEle12);
    SetBranch(t, "pass_ele23",  passEle23);
    SetBranch(t, "w_ele8",      wEle8);
    SetBranch(t, "w_ele12",     wEle12);
    SetBranch(t, "w_ele23",     wEle23);
    SetBranch(t, "mcweight",    IsDATA ? 1.f : (float)MCweight());
    SetBranch(t, "run",         ev.run());
    SetBranch(t, "lumi",        ev.lumi());
    SetBranch(t, "evt",         ev.event());

    //---- electron kinematics / IP / isolation (전부 data 에서도 얻을 수 있는 양)
    //     FR 맵의 축은 (ptcorr, |scEta|) 이므로 el_conept / el_abssceta 가 그 둘이다
    SetBranch(t, "el_pt",          (float)el.Pt());
    SetBranch(t, "el_conept",      ptCorr);
    SetBranch(t, "el_eta",         (float)el.Eta());
    SetBranch(t, "el_abseta",      (float)fabs(el.Eta()));
    SetBranch(t, "el_sceta",       el.scEta());
    SetBranch(t, "el_abssceta",    (float)fabs(el.scEta()));
    SetBranch(t, "el_phi",         (float)el.Phi());
    SetBranch(t, "el_charge",      (int)el.Charge());
    SetBranch(t, "el_etaregion",   (int)el.etaRegion());
    SetBranch(t, "el_dxy",         el.dXY());
    SetBranch(t, "el_dxyerr",      el.dXYerr());
    SetBranch(t, "el_dxysig",      el.dXY() / max(el.dXYerr(), 1e-6f));
    SetBranch(t, "el_dz",          el.dZ());
    SetBranch(t, "el_dzerr",       el.dZerr());
    SetBranch(t, "el_dzsig",       el.dZ() / max(el.dZerr(), 1e-6f));
    SetBranch(t, "el_ip3d",        el.IP3D());
    SetBranch(t, "el_sip3d",       el.SIP3D());
    SetBranch(t, "el_miniiso",     el.MiniPFRelIso());
    SetBranch(t, "el_pfiso03",     el.PfRelIso03());
    SetBranch(t, "el_pfiso04",     el.PfRelIso04());
    SetBranch(t, "el_tkreliso",    el.TkRelIso());
    SetBranch(t, "el_sieie",       el.sieie());
    SetBranch(t, "el_hoe",         el.hoe());
    SetBranch(t, "el_einvminuspinv", el.eInvMinusPInv());
    SetBranch(t, "el_detainsc",    el.deltaEtaInSC());
    SetBranch(t, "el_dphiinsc",    el.deltaPhiInSC());
    SetBranch(t, "el_detainseed",  el.deltaEtaInSeed());
    SetBranch(t, "el_r9",          el.r9());
    SetBranch(t, "el_energyerr",   el.energyErr());
    SetBranch(t, "el_ecalpfclusteriso", el.ecalPFClusterIso());
    SetBranch(t, "el_hcalpfclusteriso", el.hcalPFClusterIso());
    SetBranch(t, "el_dr03tksumpt", el.dr03TkSumPt());
    SetBranch(t, "el_mvanoiso",    el.MvaNoIso());
    SetBranch(t, "el_mvaiso",      el.MvaIso());
    SetBranch(t, "el_mvaTTH",      el.MvaTTH());
    SetBranch(t, "el_mvanoiso_wp80", el.isMVANoIsoWP80());
    SetBranch(t, "el_mvanoiso_wp90", el.isMVANoIsoWP90());
    SetBranch(t, "el_mvaiso_wp90",   el.isMVAIsoWP90());
    SetBranch(t, "el_cutbased",    (int)el.CutBased());
    SetBranch(t, "el_convveto",    el.ConvVeto());
    SetBranch(t, "el_losthits",    (int)el.LostHits());
    SetBranch(t, "el_tightCharge", (int)el.TightCharge());
    SetBranch(t, "el_seedgain",    (int)el.SeedGain());
    SetBranch(t, "el_hasJet",      (bool)(mj != nullptr));

    //---- source jet 변수 (없으면 sentinel)
    SetBranch(t, "jet_pt",          jetPt);
    SetBranch(t, "jet_rawpt",       jetRawPt);
    SetBranch(t, "jet_eta",         jetEta);
    SetBranch(t, "jet_mass",        jetMass);
    SetBranch(t, "jet_area",        jetArea);
    SetBranch(t, "jet_dr_el",       jetDR);
    SetBranch(t, "jet_ptratio",     jetPtRatio);
    SetBranch(t, "jet_ptrel",       jetPtRel);
    SetBranch(t, "jet_chHEF",       jetChHEF);
    SetBranch(t, "jet_neHEF",       jetNeHEF);
    SetBranch(t, "jet_chEmEF",      jetChEmEF);
    SetBranch(t, "jet_neEmEF",      jetNeEmEF);
    SetBranch(t, "jet_muEF",        jetMuEF);
    SetBranch(t, "jet_muonSubtr",   jetMuonSubtr);
    SetBranch(t, "jet_nconst",      jetNConst);
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
    SetBranch(t, "sv_dr_el",     svDREl);

    //---- event / away jet context
    SetBranch(t, "ev_met",          (float)METv.Pt());
    SetBranch(t, "ev_mt",           MT);
    SetBranch(t, "ev_dphi_el_met",  (float)fabs(el.DeltaPhi(METv)));
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
