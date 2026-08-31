#include "chig.h"

//==== Constructor and Destructor
chig::chig() {}
chig::~chig() {}

//==== Initialize variables
void chig::initializeAnalyzer() {

    //==== era 별 unprescaled trigger + trigger-safe pT 컷
    if (DataEra == "2016preVFP" || DataEra == "2016postVFP") {
        muTriggers = {"HLT_IsoMu24", "HLT_IsoTkMu24"};
        elTriggers = {"HLT_Ele27_WPTight_Gsf"};
        cuts.mu_pt_min = 26.;
        cuts.el_pt_min = 30.;
    } else if (DataEra == "2017") {
        muTriggers = {"HLT_IsoMu27"};
        elTriggers = {"HLT_Ele32_WPTight_Gsf_L1DoubleEG"};
        cuts.mu_pt_min = 29.;
        cuts.el_pt_min = 35.;
    } else {   // 2018
        muTriggers = {"HLT_IsoMu24"};
        elTriggers = {"HLT_Ele32_WPTight_Gsf"};
        cuts.mu_pt_min = 26.;
        cuts.el_pt_min = 35.;
    }
    //==== Had channel: 일단 HT trigger 만 (high mass 용).
    //==== HLT_Path.json 에서 "active": true 인 것만 PassTrigger 가 동작한다
    //==== TODO: low mass 는 multijet + btag trigger study 필요
    hadTriggers = {"HLT_PFHT1050"};

    //==== lepton ID: POG cut-based (tight = analysis, loose = veto)
    muonTightID = "POGTight";            // + POGPfIsoTight
    muonVetoID = "POGLoose";             // + POGPfIsoLoose
    electronTightID = "POGTight";
    electronVetoID = "POGLoose";

    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);

    //==== private Delphes-mimic 신호 샘플: PU/prefire/jet-PUID 입력이 없다
    isDelphes = MCSample.BeginsWith("chig_HcToTB");

    //==== b-tagging: DeepJet WP cut 값 (Medium = tag, Loose = 1b1L dataset 용)
    bTagWPMedium = myCorr->GetBTaggingWP(bTagger, JetTagging::JetFlavTaggerWP::Medium);
    bTagWPLoose = myCorr->GetBTaggingWP(bTagger, JetTagging::JetFlavTaggerWP::Loose);

    //==== BDT reweighting 학습용 flat tree
    NewTree("tree");
}

//==============================================================
// b-candidate 분류 (HIG-20-012 방식):
// DeepJet score 상위 2개 jet 를 먼저 b-candidate 로 고정한 뒤 분류
//   2b   : 둘 다 Medium
//   1b1L : score 1위가 Medium, 2위는 Loose 통과 + Medium 실패
//   ""   : 그 외 (reject)
// score 정렬 덕에 "Medium 1개 + Loose-only 1개" 는 항상 (1위 M, 2위 L)
//==============================================================
TString chig::ClassifyBCands() {
    bCands.clear();
    otherJets.clear();
    if (jets.size() < 2) return "";

    //==== score 순 정렬 index (동점이면 pT 높은 쪽 —
    //==== Delphes 신호는 score 가 binary 라 동점이 흔하다: unstable sort 로
    //==== b-candidate 가 빌드마다 뒤바뀌는 것을 막는 결정적 tie-break)
    RVec<int> idx(jets.size());
    for (int i = 0; i < (int)jets.size(); i++) idx[i] = i;
    sort(idx.begin(), idx.end(), [&](int a, int b) {
        const float sa = jets.at(a).GetBTaggerResult(bTagger);
        const float sb = jets.at(b).GetBTaggerResult(bTagger);
        if (sa != sb) return sa > sb;
        return jets.at(a).Pt() > jets.at(b).Pt();
    });

    bCands.push_back(jets.at(idx[0]));
    bCands.push_back(jets.at(idx[1]));
    for (int k = 2; k < (int)idx.size(); k++) otherJets.push_back(jets.at(idx[k]));
    sort(otherJets.begin(), otherJets.end(), PtComparing);

    const float s1 = bCands.at(0).GetBTaggerResult(bTagger);
    const float s2 = bCands.at(1).GetBTaggerResult(bTagger);

    if (s1 > bTagWPMedium && s2 > bTagWPMedium) return "2b";
    if (s1 > bTagWPMedium && s2 > bTagWPLoose && s2 <= bTagWPMedium) return "1b1L";
    return "";
}

//==============================================================
// analysis region: |m_top - 172.5| < 25 / 45 / 90 -> SR / VR / CR
//==============================================================
TString chig::GetRegion(float mtop) const {
    const float dm = fabs(mtop - cuts.top_mass);
    if (dm < cuts.sr_dm) return "SR";
    if (dm < cuts.vr_dm) return "VR";
    if (dm < cuts.cr_dm) return "CR";
    return "";
}

//==============================================================
// W -> l nu 의 nu pz 복원: (l + nu)^2 = mW^2 constraint 의 2차방정식
// 판별식 < 0 (mT > mW) 이면 실수부 하나만 반환
//==============================================================
RVec<float> chig::SolveNeutrinoPz(const Particle &lep, const Particle &met) const {
    const float A = cuts.w_mass * cuts.w_mass / 2.
                    + lep.Px() * met.Px() + lep.Py() * met.Py();
    const float pt2 = lep.Pt() * lep.Pt();
    const float disc = A * A - met.Pt() * met.Pt() * (lep.E() * lep.E() - lep.Pz() * lep.Pz());

    RVec<float> sols;
    if (disc < 0.) {
        sols.push_back(A * lep.Pz() / pt2);   // 실수부
    } else {
        const float sq = sqrt(disc);
        sols.push_back((A * lep.Pz() + lep.E() * sq) / pt2);
        sols.push_back((A * lep.Pz() - lep.E() * sq) / pt2);
    }
    return sols;
}

void chig::executeEvent() {

    //==== 이벤트와 raw physics object 읽기
    Event ev = GetEvent();
    RVec<Jet> rawJets = GetAllJets();

    //==== noise (MET) filter 통과 요구
    if (!PassNoiseFilter(rawJets, ev)) return;

    RVec<Muon> rawMuons = GetAllMuons();
    RVec<Electron> rawElectrons = GetAllElectrons();

    //==== Step 1: lepton 선택 (tight = analysis, veto = loose ID + loose iso)
    vetoMuons = SelectMuons(rawMuons, muonVetoID, cuts.veto_lep_pt_min, cuts.mu_eta_max);
    vetoMuons = SelectMuons(vetoMuons, "POGPfIsoLoose", cuts.veto_lep_pt_min, cuts.mu_eta_max);
    tightMuons = SelectMuons(vetoMuons, muonTightID, cuts.mu_pt_min, cuts.mu_eta_max);
    tightMuons = SelectMuons(tightMuons, "POGPfIsoTight", cuts.mu_pt_min, cuts.mu_eta_max);

    vetoElectrons = SelectElectrons(rawElectrons, electronVetoID,
                                    cuts.veto_lep_pt_min, cuts.el_eta_max);
    tightElectrons = SelectElectrons(vetoElectrons, electronTightID,
                                     cuts.el_pt_min, cuts.el_eta_max);

    sort(tightMuons.begin(), tightMuons.end(), PtComparing);
    sort(tightElectrons.begin(), tightElectrons.end(), PtComparing);

    //==== Step 2: jet 선택 - tight ID (+ Run2 는 loose PU ID), lepton cleaning
    jets = SelectJets(rawJets, "tight", cuts.jet_pt_min, cuts.jet_eta_max);
    if (Run == 2 && !isDelphes)
        jets = SelectJets(jets, "loosePuId", cuts.jet_pt_min, cuts.jet_eta_max);
    jets = JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, cuts.jet_lep_dr);
    sort(jets.begin(), jets.end(), PtComparing);

    //==== Step 3: PUPPI MET + Type-I correction
    METv = ApplyTypeICorrection(ev.GetMETVector(Event::MET_Type::PUPPI),
                                rawJets, rawElectrons, rawMuons);

    //==== Step 4: weight (MC only)
    //==== TODO: lepton ID/trigger SF, b-tag SF 는 v3 에서 추가
    float weight = 1.;
    if (!IsDATA) {
        weight = MCweight() * ev.GetTriggerLumi("Full");
        if (!isDelphes) {
            weight *= GetL1PrefireWeight();
            weight *= myCorr->GetPUWeight(ev.nTrueInt());
        }
    }

    const int nTight = tightMuons.size() + tightElectrons.size();
    const int nVeto = vetoMuons.size() + vetoElectrons.size();

    //==============================================================
    // Lep channel: W -> l nu
    //   trigger + 딱 1개의 tight lepton (추가 veto lepton 없음) + MET
    //   + >=2 jets + b-candidate 분류 (2b / 1b1L)
    //==============================================================
    bool passMuTrig = false, passElTrig = false;
    for (const auto &t : muTriggers) if (ev.PassTrigger(t)) { passMuTrig = true; break; }
    for (const auto &t : elTriggers) if (ev.PassTrigger(t)) { passElTrig = true; break; }

    //==== channel 결정 (+ 데이터는 stream 으로 중복 카운트 방지)
    TString lepFlav = "";
    Particle lepton;
    if (nTight == 1 && nVeto == 1) {
        if (tightMuons.size() == 1 && passMuTrig
            && (!IsDATA || DataStream.Contains("SingleMuon"))) {
            lepFlav = "mu";
            lepton = tightMuons.at(0);
        } else if (tightElectrons.size() == 1 && passElTrig
                   && (!IsDATA || DataStream.Contains("EGamma"))) {
            lepFlav = "el";
            lepton = tightElectrons.at(0);
        }
    }

    if (lepFlav != "") {
        const TString pre = "Lep/" + lepFlav;
        FillHist(pre + "/cutflow", 0.5, weight, 10, 0., 10.);   // 0: 1 tight lep + trig

        if (METv.Pt() > cuts.met_min) {
            FillHist(pre + "/cutflow", 1.5, weight, 10, 0., 10.);   // 1: MET

            if ((int)jets.size() >= cuts.njet_lep_min) {
                FillHist(pre + "/cutflow", 2.5, weight, 10, 0., 10.);   // 2: njet

                const TString cat = ClassifyBCands();
                if (cat != "") {
                    FillHist(pre + "/cutflow", 3.5, weight, 10, 0., 10.);   // 3: 2b|1b1L

                    //==== W 복원: nu pz 후보 x b-candidate 조합에서 |m(Wb)-mt| 최소
                    Particle bestW, bestTop;
                    int bestB = -1;
                    float bestDiff = 1e9;
                    for (const float pz : SolveNeutrinoPz(lepton, METv)) {
                        Particle nu;
                        nu.SetPxPyPzE(METv.Px(), METv.Py(), pz,
                                      sqrt(METv.Pt() * METv.Pt() + pz * pz));
                        const Particle W = lepton + nu;
                        for (int ib = 0; ib < 2; ib++) {
                            const float diff = fabs((W + bCands.at(ib)).M() - cuts.top_mass);
                            if (diff < bestDiff) {
                                bestDiff = diff;
                                bestW = W;
                                bestTop = W + bCands.at(ib);
                                bestB = ib;
                            }
                        }
                    }
                    const Jet &bTop = bCands.at(bestB);
                    const Jet &bHc = bCands.at(1 - bestB);
                    const Particle Hc = bestTop + bHc;

                    //==== region 나누기 전 m_top (region 경계 튜닝용)
                    FillHist(pre + "/" + cat + "/m_top_inclusive", bestTop.M(),
                             weight, 100, 0., 500.);

                    const TString region = GetRegion(bestTop.M());
                    if (region != "") {
                        FillHist(pre + "/cutflow", 4.5, weight, 10, 0., 10.);   // 4: region

                        const float mT = sqrt(2. * lepton.Pt() * METv.Pt()
                                              * (1. - cos(lepton.DeltaPhi(METv))));
                        const int chId = (lepFlav == "mu") ? 0 : 1;
                        const int catId = (cat == "2b") ? 2 : 1;
                        const int regId = (region == "SR") ? 0 : (region == "VR") ? 1 : 2;
                        fillOutputs(pre + "/" + cat + "/" + region, chId, catId, regId,
                                    bestW, bestTop, Hc, bTop, bHc,
                                    lepton.Pt(), mT, weight);
                    }
                }
            }
        }
    }

    //==============================================================
    // Had channel: W -> qq'
    //   trigger + lepton veto + >=4 jets + b-candidate 분류
    //   W = b-candidate 제외 jet pair 중 m(jj) 가 80.4 에 가장 가까운 것
    //==============================================================
    bool passHadTrig = false;
    for (const auto &t : hadTriggers) if (ev.PassTrigger(t)) { passHadTrig = true; break; }

    //==== offline HT: PFHT1050 plateau 요구 (모든 샘플 공통)
    float htOffline = 0.;
    for (const auto &jet : jets) htOffline += jet.Pt();

    if (passHadTrig && nVeto == 0 && htOffline > cuts.had_ht_min
        && (!IsDATA || DataStream.Contains("JetHT"))) {
        FillHist("Had/cutflow", 0.5, weight, 10, 0., 10.);   // 0: trig + HT + lep veto

        if ((int)jets.size() >= cuts.njet_had_min) {
            FillHist("Had/cutflow", 1.5, weight, 10, 0., 10.);   // 1: njet

            const TString cat = ClassifyBCands();
            if (cat != "" && otherJets.size() >= 2) {
                FillHist("Had/cutflow", 2.5, weight, 10, 0., 10.);   // 2: 2b|1b1L

                //==== W candidate: b-candidate 제외 jet pair, m(jj) closest to mW
                int wj1 = -1, wj2 = -1;
                float bestDiff = 1e9;
                for (int i = 0; i < (int)otherJets.size(); i++) {
                    for (int j = i + 1; j < (int)otherJets.size(); j++) {
                        const float m = (otherJets.at(i) + otherJets.at(j)).M();
                        const float diff = fabs(m - cuts.w_mass);
                        if (diff < bestDiff) { bestDiff = diff; wj1 = i; wj2 = j; }
                    }
                }

                if (bestDiff < cuts.w_window) {
                    FillHist("Had/cutflow", 3.5, weight, 10, 0., 10.);   // 3: W window
                    const Particle W = otherJets.at(wj1) + otherJets.at(wj2);

                    //==== top pairing: 2개 b-candidate 중 |m(Wb)-mt| 최소
                    const int bestB =
                        (fabs((W + bCands.at(0)).M() - cuts.top_mass) <
                         fabs((W + bCands.at(1)).M() - cuts.top_mass)) ? 0 : 1;
                    const Jet &bTop = bCands.at(bestB);
                    const Jet &bHc = bCands.at(1 - bestB);
                    const Particle top = W + bTop;
                    const Particle Hc = top + bHc;

                    FillHist("Had/" + cat + "/m_top_inclusive", top.M(),
                             weight, 100, 0., 500.);

                    const TString region = GetRegion(top.M());
                    if (region != "") {
                        FillHist("Had/cutflow", 4.5, weight, 10, 0., 10.);   // 4: region
                        const int catId = (cat == "2b") ? 2 : 1;
                        const int regId = (region == "SR") ? 0 : (region == "VR") ? 1 : 2;
                        fillOutputs("Had/" + cat + "/" + region, 2, catId, regId,
                                    W, top, Hc, bTop, bHc, 0., 0., weight);
                    }
                }
            }
        }
    }
}

//==============================================================
// 히스토그램 (<prefix> = <ch>/<cat>/<region>) + BDT 학습용 tree
// tree 변수는 HIG-20-012 의 BDT input 을 tb 로 번역한 것:
//   질량 (m_Hc, m_top), b pT, 후보 pT/eta, dR, boundary 변수 + HT/njet/MET
//==============================================================
void chig::fillOutputs(const TString &prefix, int channel, int cat, int region,
                       const Particle &W, const Particle &top, const Particle &Hc,
                       const Jet &bTop, const Jet &bHc,
                       float lepPt, float mT, float weight) {
    float ht = 0.;
    for (const auto &jet : jets) ht += jet.Pt();
    const float dR_bb = bTop.DeltaR(bHc);
    const float dR_WbTop = W.DeltaR(bTop);

    //==== 히스토그램
    FillHist(prefix + "/count", 0.5, weight, 1, 0., 1.);
    FillHist(prefix + "/m_W", W.M(), weight, 60, 0., 300.);
    FillHist(prefix + "/m_top", top.M(), weight, 100, 0., 500.);
    FillHist(prefix + "/m_Hc", Hc.M(), weight, 175, 0., 3500.);
    FillHist(prefix + "/pt_W", W.Pt(), weight, 150, 0., 1500.);
    FillHist(prefix + "/pt_top", top.Pt(), weight, 150, 0., 1500.);
    FillHist(prefix + "/pt_Hc", Hc.Pt(), weight, 100, 0., 1000.);
    FillHist(prefix + "/pt_bTop", bTop.Pt(), weight, 150, 0., 1500.);
    FillHist(prefix + "/pt_bHc", bHc.Pt(), weight, 150, 0., 1500.);
    FillHist(prefix + "/dR_bb", dR_bb, weight, 60, 0., 6.);
    FillHist(prefix + "/dR_WbTop", dR_WbTop, weight, 60, 0., 6.);
    FillHist(prefix + "/nJets", jets.size(), weight, 12, 0., 12.);
    FillHist(prefix + "/MET", METv.Pt(), weight, 100, 0., 1000.);
    FillHist(prefix + "/HT", ht, weight, 150, 0., 3000.);
    if (channel != 2) {
        FillHist(prefix + "/lep_pt", lepPt, weight, 200, 0., 1000.);
        FillHist(prefix + "/mT", mT, weight, 100, 0., 500.);
    }

    //==== flat tree (GBReweighter 학습용)
    SetBranch("tree", "channel", channel);   // 0=mu, 1=el, 2=had
    SetBranch("tree", "cat", cat);           // 2 = 2b, 1 = 1b1L
    SetBranch("tree", "region", region);     // 0=SR, 1=VR, 2=CR
    SetBranch("tree", "m_top", (float)top.M());
    SetBranch("tree", "m_Hc", (float)Hc.M());
    SetBranch("tree", "m_W", (float)W.M());
    SetBranch("tree", "dM_topHc", (float)(Hc.M() - top.M()));   // boundary 변수
    SetBranch("tree", "pt_bTop", (float)bTop.Pt());
    SetBranch("tree", "eta_bTop", (float)bTop.Eta());
    SetBranch("tree", "pt_bHc", (float)bHc.Pt());
    SetBranch("tree", "eta_bHc", (float)bHc.Eta());
    SetBranch("tree", "pt_W", (float)W.Pt());
    SetBranch("tree", "eta_W", (float)W.Eta());
    SetBranch("tree", "pt_top", (float)top.Pt());
    SetBranch("tree", "eta_top", (float)top.Eta());
    SetBranch("tree", "pt_Hc", (float)Hc.Pt());
    SetBranch("tree", "eta_Hc", (float)Hc.Eta());
    SetBranch("tree", "dR_bb", dR_bb);
    SetBranch("tree", "dR_WbTop", dR_WbTop);
    SetBranch("tree", "HT", ht);
    SetBranch("tree", "nJets", (int)jets.size());
    SetBranch("tree", "MET", (float)METv.Pt());
    SetBranch("tree", "lep_pt", lepPt);
    SetBranch("tree", "mT", mT);
    SetBranch("tree", "weight", weight);
    FillTrees("tree");
}
