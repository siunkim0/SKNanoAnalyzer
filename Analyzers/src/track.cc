#include "track.h"

//==== Constructor and Destructor
track::track() {}
track::~track() {}

//==============================================================
// Initialize
//==============================================================
void track::initializeAnalyzer() {

    //==== era 별 single muon trigger
    if (DataEra == "2016preVFP" || DataEra == "2016postVFP" || DataEra == "2018" ||
        DataEra == "2022" || DataEra == "2022EE") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2017") {
        IsoMuTriggerName = "HLT_IsoMu27";
        TriggerSafePtCut = 29.;
    } else {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    }

    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);
}

//==============================================================
// Event loop
//==============================================================
void track::executeEvent() {

    const TString prefix = "Central/";

    //==== 이벤트와 raw jet/muon 읽기
    Event ev = GetEvent();
    RVec<Jet> rawJets = GetAllJets();
    RVec<Muon> rawMuons = GetAllMuons();

    //==== event weight (data = 1, MC = MCweight * trigger lumi)
    float weight = 1.;
    if (!IsDATA) {
        weight *= MCweight();
        weight *= ev.GetTriggerLumi("Full");
    }

    //==== (v6) offline muon impact parameter — selection 없음 (QCD MuEnriched 스터디용)
    //==== (v7) + muon 소속 jet flavour별 분해 (NoSelFlav/): IP + iso + MT + MET
    Particle METv = ev.GetMETVector(Event::MET_Type::PUPPI);

    //==== (v12) undef 트레이스용 gen 컬렉션 (MC 전용, 루프 밖에서 1회)
    RVec<Gen> gens;
    if (!IsDATA) gens = GetAllGens();

    for (const auto &mu : rawMuons) {
        FillHist("NoSel/Muon_Pt",     mu.Pt(),     weight, 200,  0.,  200.);
        FillHist("NoSel/Muon_Eta",    mu.Eta(),    weight,  50, -2.5,   2.5);
        FillHist("NoSel/Muon_dXY",    mu.dXY(),    weight, 200, -0.5,   0.5);
        FillHist("NoSel/Muon_dXYerr", mu.dXYerr(), weight, 100,  0.,    0.05);
        FillHist("NoSel/Muon_dZ",     mu.dZ(),     weight, 200, -1.,    1.);
        FillHist("NoSel/Muon_dZerr",  mu.dZerr(),  weight, 100,  0.,    0.1);
        FillHist("NoSel/Muon_IP3D",   mu.IP3D(),   weight, 200,  0.,    0.5);
        FillHist("NoSel/Muon_SIP3D",  mu.SIP3D(),  weight, 200,  0.,   50.);

        if (IsDATA) continue;   // flavour truth는 MC 전용

        //---- (v16) muon–jet 매칭 통계 (prompt veto '이전', 전 MC muon 대상):
        //     0 = JetIdx 없음 (NanoAOD Muon_jetIdx < 0),
        //     1 = JetIdx는 있으나 저장 jet과 OriginalIndex 미매칭 (소속 jet<15GeV → NanoAOD 미저장),
        //     2 = 매칭 성공. weighted + raw(weight=1) 둘 다 기록.
        {
            int mstat;
            if (mu.JetIdx() < 0) {
                mstat = 0;
            } else {
                int jj = -1;
                for (int ij = 0; ij < (int)rawJets.size(); ij++)
                    if (rawJets[ij].OriginalIndex() == mu.JetIdx()) { jj = ij; break; }
                mstat = (jj >= 0) ? 2 : 1;
            }
            FillHist("NoSelFlavIdx/MuonJetMatchStat",    mstat, weight, 3, -0.5, 2.5);
            FillHist("NoSelFlavIdx/MuonJetMatchStatRaw", mstat, 1.f,    3, -0.5, 2.5);
        }

        //---- (v15) prompt muon 판정: undef/gluon이 정말 prompt(+uds/g mistag)인지 검증용.
        //     GetLeptonType>0 = prompt(비-fake, 비-hadronic 기원: EW prompt/tau daughter/internal conv),
        //     <=0 = fake(hadron 기원/external conv/미매칭). fake.cc·elecfake.cc와 동일 컨벤션.
        //     (v17) v15처럼 skip하지 않고 두 디렉토리에 동시 저장:
        //       NoSelFlavIdx/     = 전 muon        (= v13)
        //       NoSelFlavIdxFake/ = prompt 제거    (= v15)
        //     → 한 번의 job으로 v13/v15 비교가 모두 가능.
        const bool isPrompt = (GetLeptonType(mu, gens) > 0);

        //---- flavour 카테고리 (유저 정의: b/c=hadronFlavour, u/d/s/g=partonFlavour)
        //     (v17) _pu   = partonFlavour 0 & hadronFlavour 0 → gen parton 미매칭 = pileup jet
        //           _undef = 나머지 (주로 |partonFlavour| 4/5 인데 hadronFlavour 0 = ghost 매칭 실패)
        auto flavCat = [&](int idx) -> TString {
            if (idx < 0) return "_nojet";
            const Jet &j = rawJets[idx];
            int hf  = j.hadronFlavour();
            int apf = std::abs(j.partonFlavour());
            if (hf == 5)   return "_b";
            if (hf == 4)   return "_c";
            if (apf == 2)  return "_u";
            if (apf == 1)  return "_d";
            if (apf == 3)  return "_s";
            if (apf == 21) return "_g";
            if (apf == 0)  return "_pu";
            return "_undef";
        };

        //---- (v12) 진짜 매칭 = NanoAOD PF 연관: Muon_jetIdx ↔ Jet::OriginalIndex
        //     (muon이 실제 그 jet의 constituent일 때만; 연관 없으면 _nojet =
        //      소속 jet이 15 GeV 저장 임계 미만이라 NanoAOD에 없는 경우)
        int jidx = -1;
        if (mu.JetIdx() >= 0) {
            for (int ij = 0; ij < (int)rawJets.size(); ij++) {
                if (rawJets[ij].OriginalIndex() == mu.JetIdx()) { jidx = ij; break; }
            }
        }

        //---- (v8) FR 파라미터화 변수: near-jet pT, cone-corrected pT (fake README 컨벤션)
        float conePt = mu.Pt() * (1.f + std::max(0.f, mu.MiniPFRelIso() - 0.1f));
        float mt = std::sqrt(2.f * mu.Pt() * METv.Pt() * (1.f - std::cos(mu.DeltaPhi(METv))));

        auto fillMuVars = [&](const TString &pf, const TString &cat) -> void {
            FillHist(pf + "Pt"       + cat, mu.Pt(),           weight, 200,  0.,  200.);
            FillHist(pf + "Eta"      + cat, mu.Eta(),          weight,  50, -2.5,   2.5);
            FillHist(pf + "dXY"      + cat, mu.dXY(),          weight, 200, -0.5,   0.5);
            FillHist(pf + "dXYerr"   + cat, mu.dXYerr(),       weight, 100,  0.,    0.05);
            FillHist(pf + "dZ"       + cat, mu.dZ(),           weight, 200, -1.,    1.);
            FillHist(pf + "dZerr"    + cat, mu.dZerr(),        weight, 100,  0.,    0.1);
            FillHist(pf + "IP3D"     + cat, mu.IP3D(),         weight, 200,  0.,    0.5);
            FillHist(pf + "SIP3D"    + cat, mu.SIP3D(),        weight, 200,  0.,   50.);
            FillHist(pf + "MiniIso"  + cat, mu.MiniPFRelIso(), weight, 100,  0.,    2.);
            FillHist(pf + "RelIso"   + cat, mu.PfRelIso04(),   weight, 100,  0.,    2.);
            FillHist(pf + "TrackIso" + cat, mu.TkRelIso(),     weight, 100,  0.,    2.);
            FillHist(pf + "MT"       + cat, mt,                weight, 100,  0.,  200.);
            FillHist(pf + "MET"      + cat, METv.Pt(),         weight, 100,  0.,  200.);
            FillHist(pf + "ConePt"   + cat, conePt,            weight, 200,  0.,  400.);
        };

        //---- PF 연관 jet flavour 기준으로 전 변수 fill
        TString cat = flavCat(jidx);

        //---- (v17) 한 muon을 dir 하나에 통째로 fill하는 람다.
        //     "NoSelFlavIdx/"(전 muon) + prompt가 아니면 "NoSelFlavIdxFake/"에도 fill.
        auto fillBlock = [&](const TString &dir) -> void {
            fillMuVars(dir + "Muon_", cat);
            if (jidx >= 0) {
                const Jet &mj = rawJets[jidx];
                FillHist(dir + "Muon_MatchedJetPt" + cat, mj.Pt(),       weight, 100, 0., 1000.);
                FillHist(dir + "Muon_MatchedJetDR" + cat, mu.DeltaR(mj), weight, 100, 0.,   10.);

                //---- (v13) 소속 jet의 track 다중도 (total/charged/neutral) — QCD vs TTJJ 비교용
                //     charged/neutral branch는 Run3(2022+) NanoAODv12+에서만 정상 fill (Run2 v9엔 없음)
                const int   nTot = mj.nConstituents();
                const int   nCh  = mj.chMultiplicity();
                const int   nNe  = mj.neMultiplicity();
                const float mjpt = mj.Pt();
                auto fillTrk = [&](const TString &suf) {
                    FillHist(dir + "Muon_MatchedJetNConst" + suf, nTot, weight, 100, 0., 100.);
                    FillHist(dir + "Muon_MatchedJetChMult" + suf, nCh,  weight, 100, 0., 100.);
                    FillHist(dir + "Muon_MatchedJetNeMult" + suf, nNe,  weight, 100, 0., 100.);
                    //---- (v13) ⟨N⟩ vs jet pT 프로파일용 2D (v6 profile_flavour 스타일)
                    FillHist(dir + "Muon_MatchedJetNConst_vs_Pt" + suf, mjpt, nTot, weight, 300, 0., 3000., 100, 0., 100.);
                    FillHist(dir + "Muon_MatchedJetChMult_vs_Pt" + suf, mjpt, nCh,  weight, 300, 0., 3000., 100, 0., 100.);
                    FillHist(dir + "Muon_MatchedJetNeMult_vs_Pt" + suf, mjpt, nNe,  weight, 300, 0., 3000., 100, 0., 100.);
                };
                fillTrk("");     // inclusive
                fillTrk(cat);    // flavour split
            }

            //---- (v12) undef 트레이스: 어떤 flavour 코드/pdgId가 fall-through 되는지
            //     (v17) pileup(_pu)이 분리됐으므로 여기 남는 건 주로 |parton|=4/5 & hadron=0
            if (cat == "_undef") {
                //-- 소속 jet의 raw flavour 코드 (flavCat 어느 case에도 안 걸린 값)
                if (jidx >= 0) {
                    FillHist(dir + "Undef_JetPartonFlav", rawJets[jidx].partonFlavour(), weight, 101, -50.5, 50.5);
                    FillHist(dir + "Undef_JetHadronFlav", rawJets[jidx].hadronFlavour(), weight,  11,  -0.5, 10.5);
                }
                //-- muon 자체의 gen 기원: NanoAOD genPartFlav 요약 (0=nomatch,1=prompt,3=light,4=c,5=b,15=tau)
                FillHist(dir + "Undef_MuonGenPartFlav", (int)mu.GenPartFlav(), weight, 16, -0.5, 15.5);
                //-- muon gen-link(Muon_genPartIdx)으로 실제 pdgId 해석 (|PID|, 하드론이면 큰 값)
                int gPID = 0;
                short gIdx = mu.GenPartIdx();
                if (gIdx >= 0) {
                    for (const auto &g : gens) {
                        if (g.Index() == gIdx) { gPID = g.PID(); break; }
                    }
                }
                FillHist(dir + "Undef_MuonGenPID", std::abs(gPID), weight, 1000, 0., 1000.);
            }
        };

        fillBlock("NoSelFlavIdx/");                    // = v13 (전 muon)
        if (!isPrompt) fillBlock("NoSelFlavIdxFake/"); // = v15 (prompt 제거)
    }

    //==== (v14) raw jet (무선택) track 다중도 vs pT — TTLJ 스터디 (MC 전용, flavour truth)
    //     muon-매칭과 무관하게 GetAllJets() 전체를 flavour별로 (b/c=hadronFlavour, u/d/s/g=partonFlavour)
    if (!IsDATA) {
        for (const auto &j : rawJets) {
            int hf  = j.hadronFlavour();
            int apf = std::abs(j.partonFlavour());
            TString fcat = "_undef";
            if      (hf == 5)   fcat = "_b";
            else if (hf == 4)   fcat = "_c";
            else if (apf == 2)  fcat = "_u";
            else if (apf == 1)  fcat = "_d";
            else if (apf == 3)  fcat = "_s";
            else if (apf == 21) fcat = "_g";
            else if (apf == 0)  fcat = "_pu";   // (v17) parton/hadron 둘 다 0 = gen 미매칭 pileup jet
            const float jpt = j.Pt();
            const int   nT  = j.nConstituents();
            const int   nC  = j.chMultiplicity();
            const int   nN  = j.neMultiplicity();
            auto fillRaw = [&](const TString &suf) {
                FillHist("RawJet/NConst_vs_Pt" + suf, jpt, nT, weight, 300, 0., 3000., 100, 0., 100.);
                FillHist("RawJet/ChMult_vs_Pt" + suf, jpt, nC, weight, 300, 0., 3000., 100, 0., 100.);
                FillHist("RawJet/NeMult_vs_Pt" + suf, jpt, nN, weight, 300, 0., 3000., 100, 0., 100.);
            };
            fillRaw("");      // inclusive
            fillRaw(fcat);    // flavour split
        }
    }
}
/*
    FillHist(prefix + "CutFlow", 0.5, weight, 10, 0., 10.);  // 0: all events

    //==== noise (MET) filter
    if (!PassNoiseFilter(rawJets, ev)) return;
    FillHist(prefix + "CutFlow", 1.5, weight, 10, 0., 10.);  // 1: noise filter

    //==== single muon trigger 요구
    if (!ev.PassTrigger(IsoMuTriggerName)) return;
    FillHist(prefix + "CutFlow", 2.5, weight, 10, 0., 10.);  // 2: trigger

    //==== jet veto map (event 단위)
    if (!PassVetoMap(rawJets, rawMuons, "jetvetomap")) return;
    FillHist(prefix + "CutFlow", 3.5, weight, 10, 0., 10.);  // 3: jet veto map

    //==== MET > 20 GeV (W->munu enriched; METv는 NoSel 블록에서 이미 읽음)
    if (METv.Pt() < cuts.met_pt_min) return;
    FillHist(prefix + "CutFlow", 4.5, weight, 10, 0., 10.);  // 4: MET > 20 GeV
    FillHist(prefix + "MET", METv.Pt(), weight, 100, 0., 200.);

    //==== jet 선택: Tight ID, pT > 30, |eta| < 2.4
    RVec<Jet> jets = SelectJets(rawJets, Jet::JetID::TIGHT, cuts.jet_pt_min, cuts.jet_eta_max);
    if (jets.size() == 0) return;
    FillHist(prefix + "CutFlow", 5.5, weight, 10, 0., 10.);  // 5: >= 1 good jet

    //==== event 단위 jet 개수
    FillHist(prefix + "NJets", jets.size(), weight, 20, 0., 20.);

    //==== jet 단위 분포 (pT축 TeV까지 확장)
    for (const auto &jet : jets) {
        const float pt   = jet.Pt();
        const int   nTot = jet.nConstituents();
        const int   nCh  = jet.chMultiplicity();   // Run3(v12+) 전용 branch — 2022부터 사용
        const int   nNe  = jet.neMultiplicity();

        //---- inclusive kinematics (data + MC)
        FillHist(prefix + "Jet_Pt",  pt,        weight, 300, 0., 3000.);
        FillHist(prefix + "Jet_Eta", jet.Eta(), weight,  50, -2.5, 2.5);
        FillHist(prefix + "Jet_qgl", jet.qgl(), weight,  50, 0.,   1.);   // Run2 전용 (Run3 = -999 → underflow)

        //---- total / charged / neutral 다중도 3종 (1D + vs pT 2D)
        auto fillMult = [&](const TString &suf) {
            FillHist(prefix + "Jet_nConstituents"  + suf, nTot, weight, 100, 0., 100.);
            FillHist(prefix + "Jet_chMultiplicity" + suf, nCh,  weight, 100, 0., 100.);
            FillHist(prefix + "Jet_neMultiplicity" + suf, nNe,  weight, 100, 0., 100.);
            FillHist(prefix + "Jet_nConstituents_vs_Pt"  + suf, pt, nTot, weight, 300, 0., 3000., 100, 0., 100.);
            FillHist(prefix + "Jet_chMultiplicity_vs_Pt" + suf, pt, nCh,  weight, 300, 0., 3000., 100, 0., 100.);
            FillHist(prefix + "Jet_neMultiplicity_vs_Pt" + suf, pt, nNe,  weight, 300, 0., 3000., 100, 0., 100.);
        };
        fillMult("");   // inclusive (data + MC)

        //---- flavour split (MC truth only): v6 = u/d/s 개별 분리 (pdg: d=1, u=2, s=3)
        if (!IsDATA) {
            int apf = std::abs(jet.partonFlavour());
            TString flav = "";
            if      (apf == 1)  flav = "_d";
            else if (apf == 2)  flav = "_u";
            else if (apf == 3)  flav = "_s";
            else if (apf == 4)  flav = "_c";
            else if (apf == 5)  flav = "_b";
            else if (apf == 21) flav = "_gluon";
            if (flav != "") fillMult(flav);
        }
    }
}

        */