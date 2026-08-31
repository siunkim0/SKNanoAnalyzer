#include "closnff.h"

//==== 출력 tree 이름 (한 곳에서만 정의)
static const TString TREE = "clos";

//==== sentinel (fake.cc::fillTreeRow 와 같은 값 — skim.py 가 NaN 으로 바꾼다)
static const float FSENT = -999.;
static const int   ISENT = -99;

//==== Constructor and Destructor
closnff::closnff() {}
closnff::~closnff() {}

//==== Initialize variables
void closnff::initializeAnalyzer() {

    //==== era 별 DoubleMuon HLT 경로 (aa/TriLeptonBase.cc 와 글자 그대로 동일)
    if (DataEra == "2016preVFP") {
        DblMuTriggers = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
                         "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL"};
    } else if (DataEra == "2016postVFP") {
        DblMuTriggers = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                         "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
                         "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"};
    } else if (DataEra == "2017") {
        DblMuTriggers = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                         "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"};
    } else {
        DblMuTriggers = {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"};
    }

    //==== era 별 EMu HLT 경로 (aa/TriLeptonBase.cc 와 글자 그대로 동일).
    //==== 2016preVFP 만 _DZ 가 없다.
    if (DataEra == "2016preVFP") {
        EMuTriggers = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                       "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"};
    } else {
        EMuTriggers = {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                       "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"};
    }

    //==== loose electron WP 의 raw MVANoIso 컷 (elecfake.cc 와 같은 값)
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

    //==== era 별 jet |eta| 컷 (fake.cc 와 동일)
    if (DataEra == "2016preVFP" || DataEra == "2016postVFP") {
        cuts.jet_eta_max = 2.4;
    } else {
        cuts.jet_eta_max = 2.5;
    }

    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);

    NewTree(TREE);
}

//==============================================================
// Muon / electron / jet WP — fake.cc 의 구현을 글자 그대로 복사한 것.
// 여기서 시험하는 fake rate 가 그 WP 으로 측정됐으므로 한 글자도 달라선 안 된다.
// (framework 의 HcToWA ID 는 muon TkRelIso 컷이 no-op 이 되는 버그가 있어
//  fake.cc 가 자체 구현을 쓴다. 같은 이유로 여기서도 자체 구현을 쓴다.)
//==============================================================
bool closnff::PassMuonWP(const Muon &mu, const TString &wp) const {
    if (!mu.isPOGMediumId()) return false;
    if (!(fabs(mu.dZ()) < cuts.muon_dz_max)) return false;
    if (!(mu.TkRelIso() < cuts.muon_tkiso_max)) return false;

    if (wp == "tight") {
        if (!(fabs(mu.SIP3D()) < cuts.tight_sip3d_max)) return false;
        if (!(mu.MiniPFRelIso() < cuts.tight_miniiso_max)) return false;
    } else {
        if (!(fabs(mu.SIP3D()) < cuts.loose_sip3d_max)) return false;
        if (!(mu.MiniPFRelIso() < cuts.loose_miniiso_max)) return false;
    }
    return true;
}

//==== elecfake.cc::PassElectronWP 의 글자 그대로 복사.
//     1E2Mu 의 fake 후보가 electron 이고 그 FR 을 elecfake 가 이 WP 으로 쟀다.
//     (Electron::Pass_HcToWALooseRun2 는 raw-MVA 컷의 부호가 뒤집힌 버그가 있어
//      framework ID 를 쓰지 않는다.)
//     공통 컷은 Electron::Pass_HcToWABaseline 과 같은 네 줄이라, 이전 구현의
//     electron veto 와 동치다 -> 3Mu 결과는 이 교체로 달라지지 않는다.
bool closnff::PassElectronWP(const Electron &el, const TString &wp) const {
    if (!el.Pass_CaloIdL_TrackIdL_IsoVL()) return false;   // trigger emulation
    if (!el.ConvVeto()) return false;
    if (!(el.LostHits() <= cuts.electron_losthits_max)) return false;
    if (!(fabs(el.dZ()) < cuts.electron_dz_max)) return false;

    if (wp == "tight") {
        if (!el.isMVANoIsoWP90()) return false;
        if (!(fabs(el.SIP3D()) < cuts.tight_el_sip3d_max)) return false;
        if (!(el.MiniPFRelIso() < cuts.tight_el_miniiso_max)) return false;
    } else {
        //==== tight ⊂ loose 를 보장하기 위해 wp90 통과도 loose MVA 컷으로 인정
        float mvaCut = looseMvaEC;
        if (el.etaRegion() == Electron::ETAREGION::IB)      mvaCut = looseMvaIB;
        else if (el.etaRegion() == Electron::ETAREGION::OB) mvaCut = looseMvaOB;
        if (!(el.isMVANoIsoWP90() || el.MvaNoIso() > mvaCut)) return false;
        if (!(fabs(el.SIP3D()) < cuts.loose_el_sip3d_max)) return false;
        if (!(el.MiniPFRelIso() < cuts.loose_el_miniiso_max)) return false;
    }
    return true;
}

bool closnff::PassVetoMapJet(const Jet &jet, const RVec<Muon> &muons) const {
    if (jet.chEmEF() + jet.neEmEF() > 0.9) return true;
    for (const auto &mu : muons)
        if (jet.DeltaR(mu) < 0.2) return true;
    return !myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), "jetvetomap");
}

//==============================================================
// ClosFakeRate::configureChargeOf 와 동일: 같은 부호 두 개를 앞에 두고
// 반대 부호 하나를 뒤에 둔 (ss1, ss2, os) 를 돌려준다.
//==============================================================
std::tuple<Muon, Muon, Muon> closnff::configureChargeOf(const RVec<Muon> &muons) const {
    const Muon &mu1 = muons.at(0);
    const Muon &mu2 = muons.at(1);
    const Muon &mu3 = muons.at(2);
    if (mu1.Charge() == mu2.Charge()) return std::make_tuple(mu1, mu2, mu3);
    if (mu1.Charge() == mu3.Charge()) return std::make_tuple(mu1, mu3, mu2);
    return std::make_tuple(mu2, mu3, mu1);   // |q1+q2+q3|==1 이므로 여기 도달하면 mu2/mu3 가 ss
}

void closnff::executeEvent() {

    Event ev = GetEvent();
    rawJets = GetAllJets();

    if (!PassNoiseFilter(rawJets, ev)) return;

    RVec<Muon> rawMuons = GetAllMuons();

    if (!PassVetoMap(rawJets, rawMuons, "jetvetomap")) return;

    RVec<Electron> rawElectrons = GetAllElectrons();
    sort(rawMuons.begin(), rawMuons.end(), PtComparing);

    //==== Step 1: loose / tight muon (tight 는 loose 의 부분집합)
    //   ClosFakeRate 는 veto == loose 를 요구한다. 우리 WP 에서는 veto muon 이
    //   곧 loose muon 이므로 (fake.cc 와 같은 정의) 조건이 자동으로 만족된다.
    looseMuons.clear();
    tightMuons.clear();
    for (const auto &mu : rawMuons) {
        if (mu.Pt() < cuts.muon_pt_min || fabs(mu.Eta()) > cuts.muon_eta_max) continue;
        if (!PassMuonWP(mu, "loose")) continue;
        looseMuons.push_back(mu);
        if (PassMuonWP(mu, "tight")) tightMuons.push_back(mu);
    }

    //==== Step 2: electron 두 단계 (ClosFakeRate 와 같은 pT 문턱)
    //   veto  : loose WP, pT > 10  -> 채널 판정과 jet cleaning
    //   loose : loose WP, pT > 15  -> 1E2Mu 의 fake 후보 (tight ⊂ loose ⊂ veto)
    sort(rawElectrons.begin(), rawElectrons.end(), PtComparing);
    vetoElectrons.clear();
    looseElectrons.clear();
    tightElectrons.clear();
    for (const auto &el : rawElectrons) {
        if (!(fabs(el.Eta()) < cuts.electron_eta_max)) continue;
        if (!(el.Pt() > cuts.electron_veto_pt_min)) continue;
        if (!PassElectronWP(el, "loose")) continue;
        vetoElectrons.push_back(el);
        if (!(el.Pt() > cuts.electron_pt_min)) continue;
        looseElectrons.push_back(el);
        if (PassElectronWP(el, "tight")) tightElectrons.push_back(el);
    }

    //==== Step 3: 채널 판정. multiplicity 가 배타적이라 한 이벤트는 최대 한 채널.
    //   veto muon == loose muon 이므로 (fake.cc 와 같은 WP) muon 쪽은 크기만 본다.
    const bool is3Mu   = (looseMuons.size() == 3 && vetoElectrons.size() == 0);
    const bool is1E2Mu = (looseMuons.size() == 2 && vetoElectrons.size() == 1
                                                 && looseElectrons.size() == 1);
    if (!is3Mu && !is1E2Mu) return;
    const int channel = is3Mu ? 0 : 1;

    //==== Step 4: trigger (채널마다 다르다)
    bool passTrig = false;
    for (const auto &trig : (is3Mu ? DblMuTriggers : EMuTriggers))
        if (ev.PassTrigger(trig)) { passTrig = true; break; }
    if (!passTrig) return;

    //==== Step 5: lepton pT (raw pT 로 자른다 — ClosFakeRate 도 스케일 전에 자른다)
    if (is3Mu) {
        if (!(looseMuons.at(0).Pt() > cuts.lead_mu_pt)) return;
        if (!(looseMuons.at(1).Pt() > cuts.sub_mu_pt))  return;
        if (!(looseMuons.at(2).Pt() > cuts.sub_mu_pt))  return;
    } else {
        //==== 비대칭 leg 두 개의 OR (EMu trigger 의 두 경로에 대응). mu2 는
        //==== object 문턱(10) 말고 따로 컷이 없다 — ClosFakeRate 와 같다.
        const float mu1pt = looseMuons.at(0).Pt();
        const float elpt  = looseElectrons.at(0).Pt();
        const bool passLeadMu  = (mu1pt > cuts.emu_lead_mu_pt && elpt > cuts.emu_sub_el_pt);
        const bool passLeadEle = (mu1pt > cuts.emu_sub_mu_pt  && elpt > cuts.emu_lead_el_pt);
        if (!passLeadMu && !passLeadEle) return;
    }

    //==== Step 6: jet 선택 (fake.cc 와 같은 순서, ClosFakeRate 와 같은 20 GeV)
    //   cleaning 은 **veto** electron (pT > 10) 으로 한다 — ClosFakeRate 와 동일.
    jets = SelectJets(rawJets, "tight", cuts.jet_pt_min, cuts.jet_eta_max);
    jets = JetsVetoLeptonInside(jets, vetoElectrons, looseMuons, cuts.jet_lep_dr);
    RVec<Jet> jetsNoPuId;
    if (Run == 2) {
        RVec<Jet> vetoMapped;
        for (const auto &jet : jets)
            if (PassVetoMapJet(jet, rawMuons)) vetoMapped.push_back(jet);
        jetsNoPuId = vetoMapped;
        jets = SelectJets(vetoMapped, "loosePuId", cuts.jet_pt_min, cuts.jet_eta_max);
    }
    sort(jets.begin(), jets.end(), PtComparing);

    //==== b-jet: DeepJet Medium (ClosFakeRate 와 동일)
    const float bwp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet,
                                            JetTagging::JetFlavTaggerWP::Medium);
    bjets.clear();
    for (const auto &jet : jets)
        if (jet.GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet) > bwp)
            bjets.push_back(jet);

    if ((int)jets.size()  < cuts.n_jet_min)  return;
    if ((int)bjets.size() < cuts.n_bjet_min) return;

    //==== Step 7: MET
    METv = ApplyTypeICorrection(ev.GetMETVector(Event::MET_Type::PUPPI),
                                rawJets, rawElectrons, rawMuons);

    gens    = IsDATA ? RVec<Gen>() : GetAllGens();
    genJets = IsDATA ? RVec<GenJet>() : GetAllGenJets();

    //==== Step 8: cone-corrected 4-vector 로 스케일한 lepton 사본
    //   AN Fig. 50 의 M(mu mu) 와 pT(l^fake) 는 이 스케일된 4-vector 로 그린
    //   것이다 (ClosFakeRate::selectEvent 가 charge/mass 컷 직전에
    //   GetPTCorrScaled{Muons,Electrons} 로 컬렉션을 통째로 갈아끼운다).
    //   electron 쪽은 질량을 **0** 으로 두는데 (muon 은 원래 질량을 유지한다)
    //   이것도 GetPTCorrScaledElectrons 와 같게 맞춘 것이다.
    //   원래 pT 도 mu{i}_pt / el{i}_pt 로 남기므로 offline 에서 둘 다 쓸 수 있다.
    RVec<Muon> scaled;
    for (const auto &mu : looseMuons) {
        Muon s = mu;
        const float ptCorr = mu.Pt() * (1. + max(0.f, mu.MiniPFRelIso() - cuts.tight_miniiso_max));
        s.SetPtEtaPhiM(ptCorr, mu.Eta(), mu.Phi(), mu.M());
        scaled.emplace_back(s);
    }
    Electron scaledEl;
    if (is1E2Mu) {
        const Electron &el = looseElectrons.at(0);
        scaledEl = el;
        const float ptCorr = el.Pt() * (1. + max(0.f, el.MiniPFRelIso() - cuts.tight_el_miniiso_max));
        scaledEl.SetPtEtaPhiM(ptCorr, el.Eta(), el.Phi(), 0.);
    }

    //==== Step 9: charge + OS pair mass (스케일된 lepton 으로)
    //   Fig. 50 의 관측량. 안 쓰는 쪽 채널의 브랜치는 sentinel 로 남긴다.
    Particle ZCand, nZCand, pair_lowM, pair_highM, pair_emu, trilep;
    float nonpromptPt = FSENT, nonpromptEta = FSENT, nonpromptPhi = FSENT;
    int nonpromptSlot = ISENT;
    int nonpromptFlavor = ISENT;   // 0 = muon, 1 = electron

    if (is3Mu) {
        float qsum = 0.;
        for (const auto &mu : scaled) qsum += mu.Charge();
        if (!(fabs(qsum) == 1.)) return;

        auto [mu_ss1, mu_ss2, mu_os] = configureChargeOf(scaled);
        Particle pair1 = mu_ss1 + mu_os;
        Particle pair2 = mu_ss2 + mu_os;
        if (!(pair1.M() > cuts.pair_mass_min)) return;
        if (!(pair2.M() > cuts.pair_mass_min)) return;

        //   ZCand   = 91.2 에 가까운 OS pair, nZCand = 나머지
        //   nonprompt = nZCand 를 이루는 same-sign muon (= fake 후보)
        Muon nonprompt;
        if (fabs(pair1.M() - cuts.z_mass) < fabs(pair2.M() - cuts.z_mass)) {
            ZCand = pair1; nZCand = pair2; nonprompt = mu_ss2;
        } else {
            ZCand = pair2; nZCand = pair1; nonprompt = mu_ss1;
        }
        pair_lowM  = (pair1.M() < pair2.M()) ? pair1 : pair2;
        pair_highM = (pair1.M() > pair2.M()) ? pair1 : pair2;
        trilep = scaled.at(0) + scaled.at(1) + scaled.at(2);

        nonpromptPt  = nonprompt.Pt();
        nonpromptEta = nonprompt.Eta();
        nonpromptPhi = nonprompt.Phi();
        nonpromptFlavor = 0;
        //   nonprompt 가 어느 slot 인지 — offline 에서 그 muon 의 feature 를 찾을
        //   때 쓴다. scaled 와 looseMuons 는 같은 순서라 (eta, phi) 로 충분하다.
        for (size_t i = 0; i < scaled.size(); i++)
            if (fabs(scaled.at(i).Phi() - nonprompt.Phi()) < 1e-6 &&
                fabs(scaled.at(i).Eta() - nonprompt.Eta()) < 1e-6) nonpromptSlot = (int)i + 1;
    } else {
        //   1E2Mu: OS muon pair 하나, M > 12. Z veto 도 두 번째 pair 도 없다.
        if (!(scaled.at(0).Charge() + scaled.at(1).Charge() == 0)) return;
        pair_emu = scaled.at(0) + scaled.at(1);
        if (!(pair_emu.M() > cuts.pair_mass_min)) return;
        trilep = scaled.at(0) + scaled.at(1) + scaledEl;

        //   fake 후보는 **electron** 이다 — ClosFakeRate 는 1E2Mu 에서 어느
        //   lepton 이 실제로 fail 했는지 보지 않고 electrons[0] 을 nonprompt 로
        //   놓는다. eta 는 scEta 가 아니라 Eta() 를 쓴다 (FR 맵의 축은
        //   absScEta 인데 그림은 Eta 로 그린다 — AN 쪽의 불일치를 그대로 따른다).
        nonpromptPt  = scaledEl.Pt();      // 스케일된 = conePt
        nonpromptEta = scaledEl.Eta();
        nonpromptPhi = scaledEl.Phi();
        nonpromptFlavor = 1;
        nonpromptSlot = 1;
    }

    //==== Step 10: SR / SB
    //   SR = loose lepton 이 모두 tight. tight 는 loose 의 부분집합이므로 크기
    //   비교로 충분하다 (ClosFakeRate 와 같은 판정). electron 이 fail 해도
    //   muon 이 fail 한 것과 똑같이 세고 부호도 같다.
    const int nFail = ((int)looseMuons.size()     - (int)tightMuons.size())
                    + ((int)looseElectrons.size() - (int)tightElectrons.size());
    const int region = (nFail == 0) ? 0 : 1;

    //==== Step 11: weight — ClosFakeRate::getWeights 와 동일한 구성
    //   lepton ID SF / b-tag SF / trigger SF 는 일부러 넣지 않는다 (MC 대 MC
    //   closure 이므로 양쪽에 똑같이 곱해져 상쇄된다). top pT reweight 은
    //   ClosFakeRate 가 적용하지 않으므로 여기서도 weight 에 넣지 않고
    //   w_toppt 로 따로 남겨 offline 에서 선택할 수 있게 한다.
    float weight = 1., wTopPt = 1.;
    if (!IsDATA) {
        weight  = MCweight() * ev.GetTriggerLumi("Full");
        weight *= GetL1PrefireWeight();
        weight *= myCorr->GetPUWeight(ev.nTrueInt());
        if (MCSample.Contains("TTLL") || MCSample.Contains("TTLJ"))
            wTopPt = myCorr->GetTopPtReweight(gens);
    }

    //==== Step 12: 이벤트 요약량
    float HT = 0.;
    int nJet40 = 0;
    for (const auto &jet : jets) {
        HT += jet.Pt();
        if (jet.Pt() > 40.) nJet40++;   // 모델의 ev_njet40 과 같은 정의
    }

    //==== 이하 tree 기록. **모든 브랜치를 매 row 무조건 쓴다** —
    //==== SetBranch 가 deque 에 push 하고 FillTrees 가 비우는 구조라
    //==== 조건부로 건너뛰면 다음 row 의 주소가 어긋난다.
    SetBranch(TREE, "channel",  channel);          // 0 = 3Mu, 1 = 1E2Mu
    SetBranch(TREE, "region",   region);
    SetBranch(TREE, "n_fail",   nFail);
    SetBranch(TREE, "n_loose",  (int)(looseMuons.size() + looseElectrons.size()));
    SetBranch(TREE, "n_tight",  (int)(tightMuons.size() + tightElectrons.size()));
    SetBranch(TREE, "n_loose_mu", (int)looseMuons.size());
    SetBranch(TREE, "n_loose_el", (int)looseElectrons.size());
    SetBranch(TREE, "isData",   IsDATA);
    SetBranch(TREE, "run",      ev.run());
    SetBranch(TREE, "lumi",     ev.lumi());
    SetBranch(TREE, "evt",      ev.event());
    SetBranch(TREE, "weight",   weight);
    SetBranch(TREE, "w_toppt",  wTopPt);
    SetBranch(TREE, "mcweight", IsDATA ? 1.f : MCweight());

    SetBranch(TREE, "ev_met",     METv.Pt());
    SetBranch(TREE, "ev_met_phi", METv.Phi());
    SetBranch(TREE, "ev_ht",      HT);
    SetBranch(TREE, "ev_njet",    (int)jets.size());
    SetBranch(TREE, "ev_njet40",  nJet40);
    SetBranch(TREE, "ev_nbjet",   (int)bjets.size());
    SetBranch(TREE, "ev_nPV",     (int)PV_npvs);
    SetBranch(TREE, "ev_nTrueInt", IsDATA ? -1.f : ev.nTrueInt());

    //==== AN Fig. 50 의 관측량. 채널마다 쓰는 것이 다르고, 안 쓰는 쪽은
    //==== 기본생성된 Particle 이라 M() = pt() = 0 이 들어간다 (histogram 은
    //==== channel 로 걸러 그린다).
    //   3Mu   row : pair_lowM_mass, pair_highM_mass, nonprompt_pt
    //   1E2Mu row : pair_mass,      nonprompt_pt,    nonprompt_eta
    SetBranch(TREE, "ZCand_mass",     ZCand.M());
    SetBranch(TREE, "ZCand_pt",       ZCand.Pt());
    SetBranch(TREE, "nZCand_mass",    nZCand.M());
    SetBranch(TREE, "nZCand_pt",      nZCand.Pt());
    SetBranch(TREE, "pair_lowM_mass",  pair_lowM.M());
    SetBranch(TREE, "pair_highM_mass", pair_highM.M());
    SetBranch(TREE, "pair_mass",       pair_emu.M());     // 1E2Mu 의 단일 OS pair
    SetBranch(TREE, "pair_pt",         pair_emu.Pt());
    SetBranch(TREE, "m3l",             trilep.M());
    SetBranch(TREE, "nonprompt_pt",     nonpromptPt);     // 스케일된 = conePt
    SetBranch(TREE, "nonprompt_eta",    nonpromptEta);
    SetBranch(TREE, "nonprompt_phi",    nonpromptPhi);
    SetBranch(TREE, "nonprompt_slot",   nonpromptSlot);
    SetBranch(TREE, "nonprompt_flavor", nonpromptFlavor); // 0 = muon, 1 = electron

    //==== lepton slot 의 FR 입력 feature.
    //   브랜치 목록은 채널과 무관하게 항상 mu1/mu2/mu3 + el1 로 고정한다.
    //   3Mu 는 el1 이, 1E2Mu 는 mu3 이 sentinel 이다 (nullptr -> sentinel).
    for (int i = 0; i < 3; i++)
        fillMuonSlot(i + 1, i < (int)looseMuons.size() ? &looseMuons.at(i) : nullptr);
    fillElectronSlot(1, is1E2Mu ? &looseElectrons.at(0) : nullptr);

    FillTrees(TREE);
}

//==============================================================
// muon 한 개의 feature 를 mu{slot}_* 로 기록.
// source jet 탐색과 SV 블록은 fake.cc::fillTreeRow 에서 그대로 옮겼다 —
// 모델이 학습한 입력과 정의가 한 글자도 달라선 안 되기 때문이다
// (대표 SV = dlenSig 최대, jet 축 기준 dR < 0.4 매칭).
//==============================================================
// source jet flavour -> fake.cc 가 쓰는 것과 같은 5-class 정수 코드
//   0=b 1=c 2=uds 3=g 4=pileup, -1=source jet 없음(jet-less), -2=DATA
// fake::FlavorTag + fake::FlavorCode 를 그대로 옮겼다. 우선순위가 중요하다:
// pileup(genJetIdx<0) 을 hadronFlavour 보다 FIRST 로 본다 — 공식 코드
// MeasFakeRateV4::getMotherJetFlavour 가 그렇게 하고 fake.cc 도 그렇게 한다.
// 하나라도 어긋나면 per-event flavour 분해를 per-muon 쪽과 비교할 수 없다.
//==============================================================
int closnff::FlavorCode(const Jet *jet) const {
    if (IsDATA) return -2;
    if (!jet)   return -1;                     // jetIdx 가 없거나 컬렉션에서 못 찾음
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
void closnff::fillMuonSlot(int slot, const Muon *mu) {

    const TString p = Form("mu%d_", slot);

    //==== 이 slot 이 안 쓰이는 채널이면 (1E2Mu 의 mu3) 값만 sentinel 로 두고
    //==== 브랜치 목록은 아래에서 똑같이 다 쓴다.
    const float ptCorr = mu ? mu->Pt() * (1. + max(0.f, mu->MiniPFRelIso() - cuts.tight_miniiso_max)) : FSENT;
    const float MT = mu ? sqrt(2. * mu->Pt() * METv.Pt() * (1. - cos(mu->DeltaPhi(METv)))) : FSENT;

    //==== source jet: lepton cleaning 이전 컬렉션에서 NanoAOD jetIdx 로 찾는다
    const short jetIdx = mu ? mu->JetIdx() : -1;
    const Jet *mj = nullptr;
    if (mu) {
        for (const auto &jet : rawJets) {
            if (jet.OriginalIndex() != jetIdx) continue;
            mj = &jet;
            break;
        }
    }

    float jetPt = FSENT, jetEta = FSENT, jetAbsEta = FSENT, jetDR = FSENT;
    float jetPtRatio = FSENT, jetMuEF = FSENT, jetDeepB = FSENT;
    int   jetNConst = ISENT, jetNMuons = ISENT;
    int   jetChMult = ISENT, jetNeMult = ISENT;   // Run3 전용
    if (mj) {
        jetPt      = mj->Pt();
        jetEta     = mj->Eta();
        jetAbsEta  = fabs(mj->Eta());
        jetDR      = mj->DeltaR(*mu);
        jetPtRatio = (mj->Pt() > 0.) ? mu->Pt() / mj->Pt() : FSENT;
        jetMuEF    = mj->muEF();
        jetNConst  = mj->nConstituents();
        jetNMuons  = mj->nMuons();
        jetDeepB   = mj->GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        //==== charged / neutral multiplicity 는 Run3 (NanoAODv12+) 전용.
        //     Run2 v9 에는 branch 자체가 없고 Jet 생성자도 이 두 멤버를
        //     초기화하지 않아 Run2 에서 읽으면 쓰레기 값이 나온다
        //     (fake.cc::fillTreeRow 와 완전히 같은 처리).
        if (Run == 3) {
            jetChMult = int(mj->chMultiplicity());
            jetNeMult = int(mj->neMultiplicity());
        }
    }

    //==== source jet 안의 secondary vertex (IVF) — fake.cc 와 동일한 로직
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
            svMass    = SV_mass[bestIdx];
            svPt      = SV_pt[bestIdx];
            svPtRatio = (jetPt > 0.) ? SV_pt[bestIdx] / jetPt : FSENT;
            svDlenSig = SV_dlenSig[bestIdx];
            svDxySig  = SV_dxySig[bestIdx];
            svPAngle  = SV_pAngle[bestIdx];
            svChi2    = SV_chi2[bestIdx];
            svNTracks = (int)SV_ntracks[bestIdx];
            float dphiJ = SV_phi[bestIdx] - jetPhi;
            while (dphiJ >  M_PI) dphiJ -= 2. * M_PI;
            while (dphiJ < -M_PI) dphiJ += 2. * M_PI;
            const float detaJ = SV_eta[bestIdx] - jetEta;
            svDRJet = sqrt(detaJ * detaJ + dphiJ * dphiJ);
            float dphiM = SV_phi[bestIdx] - mu->Phi();
            while (dphiM >  M_PI) dphiM -= 2. * M_PI;
            while (dphiM < -M_PI) dphiM += 2. * M_PI;
            const float detaM = SV_eta[bestIdx] - mu->Eta();
            svDRMu = sqrt(detaM * detaM + dphiM * dphiM);
        }
    }

    const int leptonType = (mu && !IsDATA) ? GetLeptonType(*mu, gens) : ISENT;

    SetBranch(TREE, p + "pt",         mu ? (float)mu->Pt() : FSENT);
    SetBranch(TREE, p + "conept",     ptCorr);
    SetBranch(TREE, p + "eta",        mu ? (float)mu->Eta() : FSENT);
    SetBranch(TREE, p + "abseta",     mu ? (float)fabs(mu->Eta()) : FSENT);
    SetBranch(TREE, p + "phi",        mu ? (float)mu->Phi() : FSENT);
    SetBranch(TREE, p + "charge",     mu ? (int)mu->Charge() : ISENT);
    SetBranch(TREE, p + "isTight",    mu ? PassMuonWP(*mu, "tight") : false);
    SetBranch(TREE, p + "isPrompt",   mu && !IsDATA && leptonType > 0);
    SetBranch(TREE, p + "leptonType", leptonType);
    SetBranch(TREE, p + "mt",         MT);
    SetBranch(TREE, p + "miniiso",    mu ? mu->MiniPFRelIso() : FSENT);
    SetBranch(TREE, p + "sip3d",      mu ? mu->SIP3D() : FSENT);
    SetBranch(TREE, p + "hasJet",     (bool)(mj != nullptr));
    SetBranch(TREE, p + "flavor",     mu ? FlavorCode(mj) : ISENT);

    SetBranch(TREE, p + "jet_pt",      jetPt);
    SetBranch(TREE, p + "jet_eta",     jetEta);
    SetBranch(TREE, p + "jet_abseta",  jetAbsEta);
    SetBranch(TREE, p + "jet_dr_mu",   jetDR);
    SetBranch(TREE, p + "jet_ptratio", jetPtRatio);
    SetBranch(TREE, p + "jet_muEF",    jetMuEF);
    SetBranch(TREE, p + "jet_nconst",  jetNConst);
    SetBranch(TREE, p + "jet_chmult",  jetChMult);   // Run3 전용
    SetBranch(TREE, p + "jet_nemult",  jetNeMult);   // Run3 전용
    SetBranch(TREE, p + "jet_nmuons",  jetNMuons);
    SetBranch(TREE, p + "jet_deepjet_b", jetDeepB);

    SetBranch(TREE, p + "sv_n",       svN);
    SetBranch(TREE, p + "sv_mass",    svMass);
    SetBranch(TREE, p + "sv_masssum", svMassSum);
    SetBranch(TREE, p + "sv_pt",      svPt);
    SetBranch(TREE, p + "sv_ptratio", svPtRatio);
    SetBranch(TREE, p + "sv_dlensig", svDlenSig);
    SetBranch(TREE, p + "sv_dxysig",  svDxySig);
    SetBranch(TREE, p + "sv_pangle",  svPAngle);
    SetBranch(TREE, p + "sv_chi2",    svChi2);
    SetBranch(TREE, p + "sv_ntracks", svNTracks);
    SetBranch(TREE, p + "sv_dr_jet",  svDRJet);
    SetBranch(TREE, p + "sv_dr_mu",   svDRMu);
}

//==============================================================
// electron 한 개의 feature 를 el{slot}_* 로 기록 (1E2Mu 의 fake 후보).
// 이름과 정의는 elecfake.cc::fillTreeRow 와 맞춘다 — electron neural FR 이
// 그 컬럼으로 학습되므로 (el_conept, el_abssceta, sv_masssum, jet_nconst)
// 한 글자라도 다르면 모델에 다른 입력이 들어간다.
// 3Mu 이벤트에서는 el == nullptr 로 불려 전부 sentinel 이 된다.
//==============================================================
void closnff::fillElectronSlot(int slot, const Electron *el) {

    const TString p = Form("el%d_", slot);

    const float ptCorr = el ? el->Pt() * (1. + max(0.f, el->MiniPFRelIso() - cuts.tight_el_miniiso_max)) : FSENT;
    const float MT = el ? sqrt(2. * el->Pt() * METv.Pt() * (1. - cos(el->DeltaPhi(METv)))) : FSENT;

    //==== source jet: elecfake 와 같이 NanoAOD jetIdx 로 찾는다
    const short jetIdx = el ? el->JetIdx() : -1;
    const Jet *mj = nullptr;
    if (el && jetIdx >= 0) {
        for (const auto &jet : rawJets) {
            if (jet.OriginalIndex() != jetIdx) continue;
            mj = &jet;
            break;
        }
    }

    float jetPt = FSENT, jetEta = FSENT, jetAbsEta = FSENT, jetDR = FSENT;
    float jetPtRatio = FSENT, jetMuEF = FSENT, jetDeepB = FSENT;
    int   jetNConst = ISENT, jetNMuons = ISENT;
    int   jetChMult = ISENT, jetNeMult = ISENT;   // Run3 전용
    if (mj) {
        jetPt      = mj->Pt();
        jetEta     = mj->Eta();
        jetAbsEta  = fabs(mj->Eta());
        jetDR      = mj->DeltaR(*el);
        jetPtRatio = (mj->Pt() > 0.) ? el->Pt() / mj->Pt() : FSENT;
        jetMuEF    = mj->muEF();
        jetNConst  = mj->nConstituents();
        jetNMuons  = mj->nMuons();
        jetDeepB   = mj->GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        //==== charged / neutral multiplicity 는 Run3 (NanoAODv12+) 전용.
        //     Run2 v9 에는 branch 자체가 없고 Jet 생성자도 이 두 멤버를
        //     초기화하지 않아 Run2 에서 읽으면 쓰레기 값이 나온다
        //     (fake.cc::fillTreeRow 와 완전히 같은 처리).
        if (Run == 3) {
            jetChMult = int(mj->chMultiplicity());
            jetNeMult = int(mj->neMultiplicity());
        }
    }

    //==== source jet 안의 secondary vertex (IVF) — elecfake / fake 와 같은 로직
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
            svMass    = SV_mass[bestIdx];
            svPt      = SV_pt[bestIdx];
            svPtRatio = (jetPt > 0.) ? SV_pt[bestIdx] / jetPt : FSENT;
            svDlenSig = SV_dlenSig[bestIdx];
            svDxySig  = SV_dxySig[bestIdx];
            svPAngle  = SV_pAngle[bestIdx];
            svChi2    = SV_chi2[bestIdx];
            svNTracks = (int)SV_ntracks[bestIdx];
            float dphiJ = SV_phi[bestIdx] - jetPhi;
            while (dphiJ >  M_PI) dphiJ -= 2. * M_PI;
            while (dphiJ < -M_PI) dphiJ += 2. * M_PI;
            const float detaJ = SV_eta[bestIdx] - jetEta;
            svDRJet = sqrt(detaJ * detaJ + dphiJ * dphiJ);
            float dphiE = SV_phi[bestIdx] - el->Phi();
            while (dphiE >  M_PI) dphiE -= 2. * M_PI;
            while (dphiE < -M_PI) dphiE += 2. * M_PI;
            const float detaE = SV_eta[bestIdx] - el->Eta();
            svDREl = sqrt(detaE * detaE + dphiE * dphiE);
        }
    }

    const int leptonType = (el && !IsDATA) ? GetLeptonType(*el, gens) : ISENT;

    SetBranch(TREE, p + "pt",         el ? (float)el->Pt() : FSENT);
    SetBranch(TREE, p + "conept",     ptCorr);
    SetBranch(TREE, p + "eta",        el ? (float)el->Eta() : FSENT);
    SetBranch(TREE, p + "abseta",     el ? (float)fabs(el->Eta()) : FSENT);
    //   FR 맵과 모델의 축은 supercluster eta 다 (그림의 eta 와 다르다)
    SetBranch(TREE, p + "sceta",      el ? el->scEta() : FSENT);
    SetBranch(TREE, p + "abssceta",   el ? (float)fabs(el->scEta()) : FSENT);
    SetBranch(TREE, p + "phi",        el ? (float)el->Phi() : FSENT);
    SetBranch(TREE, p + "charge",     el ? (int)el->Charge() : ISENT);
    SetBranch(TREE, p + "isTight",    el ? PassElectronWP(*el, "tight") : false);
    SetBranch(TREE, p + "isPrompt",   el && !IsDATA && leptonType > 0);
    SetBranch(TREE, p + "leptonType", leptonType);
    SetBranch(TREE, p + "mt",         MT);
    SetBranch(TREE, p + "miniiso",    el ? el->MiniPFRelIso() : FSENT);
    SetBranch(TREE, p + "sip3d",      el ? el->SIP3D() : FSENT);
    SetBranch(TREE, p + "mvanoiso",   el ? el->MvaNoIso() : FSENT);
    SetBranch(TREE, p + "hasJet",     (bool)(mj != nullptr));
    SetBranch(TREE, p + "flavor",     el ? FlavorCode(mj) : ISENT);

    SetBranch(TREE, p + "jet_pt",      jetPt);
    SetBranch(TREE, p + "jet_eta",     jetEta);
    SetBranch(TREE, p + "jet_abseta",  jetAbsEta);
    SetBranch(TREE, p + "jet_dr_el",   jetDR);
    SetBranch(TREE, p + "jet_ptratio", jetPtRatio);
    SetBranch(TREE, p + "jet_muEF",    jetMuEF);
    SetBranch(TREE, p + "jet_nconst",  jetNConst);
    SetBranch(TREE, p + "jet_chmult",  jetChMult);   // Run3 전용
    SetBranch(TREE, p + "jet_nemult",  jetNeMult);   // Run3 전용
    SetBranch(TREE, p + "jet_nmuons",  jetNMuons);
    SetBranch(TREE, p + "jet_deepjet_b", jetDeepB);

    SetBranch(TREE, p + "sv_n",       svN);
    SetBranch(TREE, p + "sv_mass",    svMass);
    SetBranch(TREE, p + "sv_masssum", svMassSum);
    SetBranch(TREE, p + "sv_pt",      svPt);
    SetBranch(TREE, p + "sv_ptratio", svPtRatio);
    SetBranch(TREE, p + "sv_dlensig", svDlenSig);
    SetBranch(TREE, p + "sv_dxysig",  svDxySig);
    SetBranch(TREE, p + "sv_pangle",  svPAngle);
    SetBranch(TREE, p + "sv_chi2",    svChi2);
    SetBranch(TREE, p + "sv_ntracks", svNTracks);
    SetBranch(TREE, p + "sv_dr_jet",  svDRJet);
    SetBranch(TREE, p + "sv_dr_el",   svDREl);
}
