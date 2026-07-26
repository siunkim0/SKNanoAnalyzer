#include "Wtag.h"
#include <set>

//==== Constructor and Destructor
Wtag::Wtag() {}
Wtag::~Wtag() {}

//==============================================================
// Initialize
//==============================================================
void Wtag::initializeAnalyzer() {
    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);
    NewTree("jets");   // v2: per-AK8-jet ntuple for tagger training
}

//==============================================================
// Event loop  (MC only — gen-truth W→qq' study)
//==============================================================
void Wtag::executeEvent() {

    if (IsDATA) return;   // needs gen truth

    RVec<Gen>    gens    = GetAllGens();
    RVec<Jet>    rawJets = GetAllJets();
    RVec<FatJet> fatJets = GetAllFatJets();

    //==== event weight (sign × genWeight, xsec-normalized). No trigger/selection:
    //     fraction curves are ratios within a sample so weight cancels; kept for
    //     correct NLO sign handling and cross-sample stacking.
    float weight = MCweight();

    //==== good AK4 jets (Tight ID, pT > 20, |eta| < 2.5)
    RVec<Jet> jets = SelectJets(rawJets, Jet::JetID::TIGHT, cuts.ak4_pt_min, cuts.ak4_eta_max);

    const int   NWPT    = Wtag::NWPT;
    const float WPT_MAX = Wtag::WPT_MAX;

    //---- per-jet observable filler (Stage 3): track mult / mass / nSVs vs jet pT
    auto fillJetVars = [&](const TString &cat, const Jet &j) {
        const float pt  = j.Pt();
        const int   nT  = j.nConstituents();
        const int   nC  = j.chMultiplicity();   // charged ≈ # tracks (Run3 v12+)
        const int   nN  = j.neMultiplicity();
        const int   nSV = j.nSVs();
        const float m   = j.M();
        //-- 1D
        FillHist("Jet/Pt"      + cat, pt,  weight, 300, 0., 3000.);
        FillHist("Jet/Mass"    + cat, m,   weight, 120, 0.,  300.);
        FillHist("Jet/NConst"  + cat, nT,  weight, 100, 0.,  100.);
        FillHist("Jet/ChMult"  + cat, nC,  weight, 100, 0.,  100.);
        FillHist("Jet/NeMult"  + cat, nN,  weight, 100, 0.,  100.);
        FillHist("Jet/NSV"     + cat, nSV, weight,  10, 0.,   10.);
        //-- vs jet pT (2D profiles — p_T-matched comparison)
        FillHist("Jet/NConst_vs_Pt" + cat, pt, nT,  weight, 300, 0., 3000., 100, 0., 100.);
        FillHist("Jet/ChMult_vs_Pt" + cat, pt, nC,  weight, 300, 0., 3000., 100, 0., 100.);
        FillHist("Jet/NeMult_vs_Pt" + cat, pt, nN,  weight, 300, 0., 3000., 100, 0., 100.);
        FillHist("Jet/NSV_vs_Pt"    + cat, pt, nSV, weight, 300, 0., 3000.,  10, 0.,  10.);
        FillHist("Jet/Mass_vs_Pt"   + cat, pt, m,   weight, 300, 0., 3000., 120, 0., 300.);
    };

    //==== find hadronic W: last-copy W with exactly two quark daughters
    std::set<int> usedJetIdx;     // AK4 jets assigned to a W (excluded from control)
    std::set<int> usedFatJetIdx;  // AK8 fatjets assigned to a W

    for (const auto &w : gens) {
        if (std::abs(w.PID()) != 24) continue;
        if (!w.isLastCopy())         continue;

        //---- direct quark daughters of this W
        RVec<const Gen*> qd;
        for (const auto &g : gens) {
            if (g.MotherIndex() != w.Index()) continue;
            int apid = std::abs(g.PID());
            if (apid >= 1 && apid <= 6) qd.push_back(&g);
        }
        if (qd.size() != 2) continue;   // require hadronic W→qq'

        const Gen &q1 = *qd[0];
        const Gen &q2 = *qd[1];
        const float wpt  = w.Pt();
        const float dRqq = q1.DeltaR(q2);

        //==== Stage 1 — gen truth: boost vs opening angle
        FillHist("Gen/W_Pt_all",   wpt,       weight, NWPT, 0., WPT_MAX);
        FillHist("Gen/dRqq_vs_WPt", wpt, dRqq, weight, NWPT, 0., WPT_MAX, 100, 0., 2.0);
        if (dRqq < cuts.ak4_dr) FillHist("Gen/W_Pt_dRqq_lt04", wpt, weight, NWPT, 0., WPT_MAX);
        if (dRqq < cuts.ak8_dr) FillHist("Gen/W_Pt_dRqq_lt08", wpt, weight, NWPT, 0., WPT_MAX);

        //==== Stage 2 — detector: match each quark to nearest AK4 jet (ΔR<0.4)
        auto matchJet = [&](const Gen &q) -> int {
            int   best = -1;
            float bestdr = cuts.ak4_dr;
            for (int ij = 0; ij < (int)jets.size(); ij++) {
                float dr = jets[ij].DeltaR(q);
                if (dr < bestdr) { bestdr = dr; best = ij; }
            }
            return best;
        };
        int j1 = matchJet(q1);
        int j2 = matchJet(q2);
        bool bothMatched = (j1 >= 0 && j2 >= 0);
        bool merged      = bothMatched && (j1 == j2);
        bool resolved    = bothMatched && (j1 != j2);

        if (bothMatched) FillHist("Reco/W_Pt_ak4bothmatched", wpt, weight, NWPT, 0., WPT_MAX);
        if (merged)      FillHist("Reco/W_Pt_ak4merged",      wpt, weight, NWPT, 0., WPT_MAX);
        if (resolved)    FillHist("Reco/W_Pt_ak4resolved",    wpt, weight, NWPT, 0., WPT_MAX);

        //---- overlap of the two resolved jets
        if (resolved) {
            float dRjj = jets[j1].DeltaR(jets[j2]);
            FillHist("Reco/dRjj_vs_WPt", wpt, dRjj, weight, NWPT, 0., WPT_MAX, 100, 0., 4.0);
            if (dRjj < cuts.ak8_dr) FillHist("Reco/W_Pt_ak4overlap", wpt, weight, NWPT, 0., WPT_MAX);
        }

        //==== AK8 capture: one good FatJet within ΔR<0.8 of both quarks
        int fjIdx = -1;
        for (int ifj = 0; ifj < (int)fatJets.size(); ifj++) {
            const FatJet &fj = fatJets[ifj];
            if (fj.Pt() < cuts.ak8_pt_min || std::abs(fj.Eta()) > cuts.ak8_eta_max) continue;
            if (fj.DeltaR(q1) < cuts.ak8_dr && fj.DeltaR(q2) < cuts.ak8_dr) { fjIdx = ifj; break; }
        }
        bool captured = (fjIdx >= 0);
        if (captured) FillHist("Reco/W_Pt_ak8captured", wpt, weight, NWPT, 0., WPT_MAX);

        //==== Stage 3 — jet kinematics of the W jet(s)
        if (merged) {
            usedJetIdx.insert(j1);
            fillJetVars("_mergedW", jets[j1]);
        } else if (resolved) {
            usedJetIdx.insert(j1);
            usedJetIdx.insert(j2);
            fillJetVars("_singleW", jets[j1]);
            fillJetVars("_singleW", jets[j2]);
        }

        //==== AK8 substructure of the captured W fatjet
        if (captured) {
            usedFatJetIdx.insert(fjIdx);
            const FatJet &fj = fatJets[fjIdx];
            const float tau21 = (fj.Tau1() > 0.f) ? fj.Tau2() / fj.Tau1() : -1.f;
            const float wvq   = fj.ParticleNetWithMass_WvsQCD();
            FillHist("FatJet/Pt_capturedW",     fj.Pt(),           weight, 300, 0., 3000.);
            FillHist("FatJet/SDMass_capturedW", fj.SDMass(),       weight, 120, 0.,  300.);
            FillHist("FatJet/NConst_capturedW", fj.NConstituents(), weight, 150, 0., 150.);
            FillHist("FatJet/Tau21_capturedW",  tau21,             weight, 100, 0.,   1.);
            FillHist("FatJet/WvsQCD_capturedW", wvq,               weight, 100, 0.,   1.);
            FillHist("FatJet/SDMass_vs_Pt_capturedW", fj.Pt(), fj.SDMass(),        weight, 300, 0., 3000., 120, 0., 300.);
            FillHist("FatJet/NConst_vs_Pt_capturedW", fj.Pt(), fj.NConstituents(), weight, 300, 0., 3000., 150, 0., 150.);
            FillHist("FatJet/Tau21_vs_Pt_capturedW",  fj.Pt(), tau21,               weight, 300, 0., 3000., 100, 0.,   1.);
            FillHist("FatJet/WvsQCD_vs_Pt_capturedW", fj.Pt(), wvq,                 weight, 300, 0., 3000., 100, 0.,   1.);
        }
    }

    //==== Stage 3 control: light/gluon (and c,b) jets NOT assigned to any W
    for (int ij = 0; ij < (int)jets.size(); ij++) {
        if (usedJetIdx.count(ij)) continue;
        int apf = std::abs(jets[ij].partonFlavour());
        TString fc;
        if      (apf >= 1 && apf <= 3) fc = "_lightuds";
        else if (apf == 21)           fc = "_gluon";
        else if (apf == 4)            fc = "_c";
        else if (apf == 5)            fc = "_b";
        else continue;   // undefined flavour
        fillJetVars(fc, jets[ij]);
    }

    //==== AK8 control: fatjets not matched to any W (QCD-like)
    for (int ifj = 0; ifj < (int)fatJets.size(); ifj++) {
        if (usedFatJetIdx.count(ifj)) continue;
        const FatJet &fj = fatJets[ifj];
        if (fj.Pt() < cuts.ak8_pt_min || std::abs(fj.Eta()) > cuts.ak8_eta_max) continue;
        const float tau21 = (fj.Tau1() > 0.f) ? fj.Tau2() / fj.Tau1() : -1.f;
        const float wvq   = fj.ParticleNetWithMass_WvsQCD();
        FillHist("FatJet/Pt_other",     fj.Pt(),            weight, 300, 0., 3000.);
        FillHist("FatJet/SDMass_other", fj.SDMass(),        weight, 120, 0.,  300.);
        FillHist("FatJet/NConst_other", fj.NConstituents(), weight, 150, 0., 150.);
        FillHist("FatJet/Tau21_other",  tau21,              weight, 100, 0.,   1.);
        FillHist("FatJet/WvsQCD_other", wvq,                weight, 100, 0.,   1.);
        FillHist("FatJet/SDMass_vs_Pt_other", fj.Pt(), fj.SDMass(),        weight, 300, 0., 3000., 120, 0., 300.);
        FillHist("FatJet/NConst_vs_Pt_other", fj.Pt(), fj.NConstituents(), weight, 300, 0., 3000., 150, 0., 150.);
        FillHist("FatJet/Tau21_vs_Pt_other",  fj.Pt(), tau21,               weight, 300, 0., 3000., 100, 0.,   1.);
        FillHist("FatJet/WvsQCD_vs_Pt_other", fj.Pt(), wvq,                 weight, 300, 0., 3000., 100, 0.,   1.);
    }

    //==============================================================
    // v2 — per-AK8-jet ntuple for tagger training
    //   label: 1 = clean hadronic-W jet (both quarks in cone, no extra
    //          hard b), 2 = W + hard b in cone (t→bqq̄' capture),
    //          0 = no gen-W match (background candidate; use QCD-sample
    //          jets for training, TT/VV unmatched jets are spectators)
    //==============================================================

    //---- hadronic gen Ws (last copy, two quark daughters)
    struct HadW { const Gen* w; const Gen* q1; const Gen* q2; };
    std::vector<HadW> hadWs;
    for (const auto &w : gens) {
        if (std::abs(w.PID()) != 24 || !w.isLastCopy()) continue;
        RVec<const Gen*> qd;
        for (const auto &g : gens) {
            if (g.MotherIndex() != w.Index()) continue;
            int apid = std::abs(g.PID());
            if (apid >= 1 && apid <= 6) qd.push_back(&g);
        }
        if (qd.size() == 2) hadWs.push_back({&w, qd[0], qd[1]});
    }

    //---- hard-process b quarks (top decay) for the contamination flag
    std::vector<const Gen*> hardBs;
    for (const auto &g : gens) {
        if (std::abs(g.PID()) != 5)                  continue;
        if (!g.fromHardProcess() || !g.isLastCopy()) continue;
        hardBs.push_back(&g);
    }

    for (int ifj = 0; ifj < (int)fatJets.size(); ifj++) {
        const FatJet &fj = fatJets[ifj];
        if (fj.Pt() < cuts.ak8_pt_min || std::abs(fj.Eta()) > cuts.ak8_eta_max) continue;
        // FatJet::PassTight() is a no-op (SetJetID never fills j_jetId) —
        // framework bug, worked around with the raw branch (index-aligned)
        if (!(FatJet_jetId[ifj] & 2)) continue;

        //---- gen-W match: both quark daughters inside the AK8 cone
        const HadW *mW = nullptr;
        float bestdr = 1e9;
        for (const auto &hw : hadWs) {
            if (fj.DeltaR(*hw.q1) > cuts.ak8_dr || fj.DeltaR(*hw.q2) > cuts.ak8_dr) continue;
            float dr = fj.DeltaR(*hw.w);
            if (dr < bestdr) { bestdr = dr; mW = &hw; }
        }
        //---- extra hard b in cone (not a daughter of the matched W — W→cb̄ exists)
        int hasB = 0;
        for (const auto *b : hardBs) {
            if (mW && (b == mW->q1 || b == mW->q2)) continue;
            if (fj.DeltaR(*b) < cuts.ak8_dr) { hasB = 1; break; }
        }
        int label = mW ? (hasB ? 2 : 1) : 0;

        //==== v2b: raw-branch shape comparison, merged W (label 1) vs QCD (label 0).
        //     label 2 (W+b/top) is a spectator and is not drawn.
        //     Class is by label here, but the QCD class must be taken from the QCD
        //     samples only at merge time — label 0 also exists in ttbar/diboson.
        //
        //     Two histogram sets, because the m_SD window changes the answer rather
        //     than just tightening it. "Branch/" keeps m_SD in [20,250] (the v2 ML
        //     dataset boundary, and the regime a deployed tagger actually runs in);
        //     "BranchNoMSD/" drops it. That window keeps 98% of merged W but only
        //     41% of QCD, so ranking the mass branches inside it is circular —
        //     msoftdrop is being scored on a sample msoftdrop already purified.
        //     Filling both makes every branch's conditionality measurable, including
        //     the ParticleNet scores, which are not written to the `jets` tree and so
        //     cannot be re-measured offline the way tau21/n2b1 can.
        const bool pass_pt  = fj.Pt() > 200. && fj.Pt() < 1200.;
        const bool pass_msd = fj.SDMass() > 20. && fj.SDMass() < 250.;
        for (int pass = 0; pass < 2 && label != 2 && pass_pt; ++pass) {
            if (pass == 0 && !pass_msd) continue;

            const TString c = TString(pass == 0 ? "Branch/" : "BranchNoMSD/")
                            + (label == 1 ? "W_" : "QCD_");
            // Scores live in [-1,1], but two branches escape it: the ParticleNet
            // vs-QCD ratios use **-10** (not -1) for "not evaluated", and lsf3 has a
            // tail out to ~33. On a plain [-1,1] axis both land in under/overflow,
            // which Integral() and the separation loop ignore — so the two classes
            // get normalised over different subsets (measured: 1.5% of W vs 3.3% of
            // QCD undefined for the PNet_x* heads). Pad by 10 bins each side and
            // bucket out-of-range values there: physical binning is unchanged at
            // 0.01/bin, nothing is dropped, and the undefined spike sits at -1.05
            // where it cannot be mistaken for a real score.
            auto S = [&](const TString &n, float v) {
                float x = (v < -1.f) ? -1.05f : ((v > 1.f) ? 1.05f : v);
                FillHist(c + n, x, weight, 220, -1.1, 1.1);
            };
            auto F = [&](const TString &n, float v) { FillHist(c + n, v, weight, 100,  0., 1.); };
            auto N = [&](const TString &n, float v, int hi) { FillHist(c + n, v, weight, hi, 0., (float)hi); };
            auto M = [&](const TString &n, float v) { FillHist(c + n, v, weight, 120,  0., 300.); };

            //-- kinematics / bulk
            FillHist(c + "pt",   fj.Pt(),  weight, 100, 200., 1200.);
            FillHist(c + "eta",  fj.Eta(), weight,  60,  -3.,    3.);
            FillHist(c + "phi",  fj.Phi(), weight,  64, -3.2,  3.2);
            FillHist(c + "area", fj.Area(), weight, 100, 0.,     4.);
            M("mass",      fj.M());
            M("msoftdrop", fj.SDMass());
            M("mreg",      fj.ParticleNet_MassCorr() * fj.M());
            FillHist(c + "masscorr", fj.ParticleNet_MassCorr(), weight, 100, 0., 2.);

            //-- substructure
            F("tau1", fj.Tau1());  F("tau2", fj.Tau2());
            F("tau3", fj.Tau3());  F("tau4", fj.Tau4());
            S("tau21", fj.Tau1() > 0. ? fj.Tau2() / fj.Tau1() : -1.f);
            S("tau32", fj.Tau2() > 0. ? fj.Tau3() / fj.Tau2() : -1.f);
            S("tau43", fj.Tau3() > 0. ? fj.Tau4() / fj.Tau3() : -1.f);
            S("n2b1", fj.N2b1());   // -1 sentinel: defined only for raw pT > 250
            // N3 is NOT bounded by 1 the way N2 is: measured range [0.04, 4.21],
            // ~94% of defined entries exceed 1. With the [-1,1] score axis they
            // all landed in overflow, which Integral()/separation() ignore, and
            // n3b1 wrongly read as S=0 rather than its true 0.024. Own axis.
            FillHist(c + "n3b1", fj.N3b1(), weight, 250, -1.25, 5.);
            S("lsf3", fj.LSF3());

            //-- multiplicities & energy fractions
            N("nconst", fj.NConstituents(), 150);
            N("chmult", fj.ChMultiplicity(), 150);
            N("nemult", fj.NeMultiplicity(), 150);
            N("nbhad",  fj.NBHadrons(), 10);   // MC truth: QCD composition only
            N("nchad",  fj.NCHadrons(), 10);
            F("chhef",  fj.ChHEF());   F("nehef",  fj.NeHEF());
            F("chemef", fj.ChEmEF());  F("neemef", fj.NeEmEF());
            F("muef",   fj.MuEF());

            //-- double-b / Hbb taggers (loader arrays: no getters in FatJet)
            S("ddbvl",      FatJet_btagDDBvLV2[ifj]);
            S("ddcvb",      FatJet_btagDDCvBV2[ifj]);
            S("ddcvl",      FatJet_btagDDCvLV2[ifj]);
            S("btag_deepb", FatJet_btagDeepB[ifj]);
            S("btag_hbb",   FatJet_btagHbb[ifj]);

            //-- ParticleNet. PNetM_ = with-mass, PNet_ = mass-decorrelated; the
            //   two "QCD" sums collide without the prefix. Run-3-only arrays.
            if (Run == 3) {
                S("PNetM_h4qvsqcd", FatJet_particleNetWithMass_H4qvsQCD[ifj]);
                S("PNetM_hbbvsqcd", FatJet_particleNetWithMass_HbbvsQCD[ifj]);
                S("PNetM_hccvsqcd", FatJet_particleNetWithMass_HccvsQCD[ifj]);
                S("PNetM_wvsqcd",   FatJet_particleNetWithMass_WvsQCD[ifj]);
                S("PNetM_zvsqcd",   FatJet_particleNetWithMass_ZvsQCD[ifj]);
                S("PNetM_tvsqcd",   FatJet_particleNetWithMass_TvsQCD[ifj]);
                S("PNetM_qcd",      FatJet_particleNetWithMass_QCD[ifj]);
                S("PNet_qcd",       FatJet_particleNet_QCD[ifj]);
                S("PNet_qcd0hf",    FatJet_particleNet_QCD0HF[ifj]);
                S("PNet_qcd1hf",    FatJet_particleNet_QCD1HF[ifj]);
                S("PNet_qcd2hf",    FatJet_particleNet_QCD2HF[ifj]);
                S("PNet_xbbvsqcd",  FatJet_particleNet_XbbVsQCD[ifj]);
                S("PNet_xccvsqcd",  FatJet_particleNet_XccVsQCD[ifj]);
                S("PNet_xqqvsqcd",  FatJet_particleNet_XqqVsQCD[ifj]);
                S("PNet_xggvsqcd",  FatJet_particleNet_XggVsQCD[ifj]);
                S("PNet_xtevsqcd",  FatJet_particleNet_XteVsQCD[ifj]);
                S("PNet_xtmvsqcd",  FatJet_particleNet_XtmVsQCD[ifj]);
                S("PNet_xttvsqcd",  FatJet_particleNet_XttVsQCD[ifj]);
            }
        }

        //---- softdrop subjets (loader arrays, indexed from the fatjet)
        float sj1_pt = -1., sj1_mass = -1., sj1_btag = -1.;
        float sj2_pt = -1., sj2_mass = -1., sj2_btag = -1.;
        float sj_dr = -1., sj_z = -1.;
        const int i1 = fj.SubJetIdx1(), i2 = fj.SubJetIdx2();
        if (i1 >= 0 && i1 < nSubJet) { sj1_pt = SubJet_pt[i1]; sj1_mass = SubJet_mass[i1]; sj1_btag = SubJet_btagDeepB[i1]; }
        if (i2 >= 0 && i2 < nSubJet) { sj2_pt = SubJet_pt[i2]; sj2_mass = SubJet_mass[i2]; sj2_btag = SubJet_btagDeepB[i2]; }
        if (sj1_pt > 0. && sj2_pt > 0.) {
            const float dEta = SubJet_eta[i1] - SubJet_eta[i2];
            const float dPhi = TVector2::Phi_mpi_pi(SubJet_phi[i1] - SubJet_phi[i2]);
            sj_dr = std::sqrt(dEta*dEta + dPhi*dPhi);
            sj_z  = std::min(sj1_pt, sj2_pt) / (sj1_pt + sj2_pt);
        }

        //---- secondary vertices inside the AK8 cone
        int   nsv = 0, sv_sumntracks = 0;
        float sv_summass = 0., sv_maxdlensig = -1.;
        RVec<int> svIdx;
        for (int is = 0; is < nSV; is++) {
            const float dEta = fj.Eta() - SV_eta[is];
            const float dPhi = TVector2::Phi_mpi_pi(fj.Phi() - SV_phi[is]);
            if (dEta*dEta + dPhi*dPhi > cuts.ak8_dr*cuts.ak8_dr) continue;
            nsv++;
            sv_summass    += SV_mass[is];
            sv_sumntracks += SV_ntracks[is];
            sv_maxdlensig  = std::max(sv_maxdlensig, SV_dlenSig[is]);
            svIdx.push_back(is);
        }
        //---- leading 3 in-cone SVs by decay-length significance (ParT-lite tokens)
        std::sort(svIdx.begin(), svIdx.end(),
                  [&](int a, int b) { return SV_dlenSig[a] > SV_dlenSig[b]; });
        float svtok[3][4];   // pt, mass, dlensig, ntracks
        for (int k = 0; k < 3; k++) {
            if (k < (int)svIdx.size()) {
                const int is = svIdx[k];
                svtok[k][0] = SV_pt[is];  svtok[k][1] = SV_mass[is];
                svtok[k][2] = SV_dlenSig[is]; svtok[k][3] = (float)SV_ntracks[is];
            } else {
                svtok[k][0] = svtok[k][1] = svtok[k][2] = svtok[k][3] = -1.f;
            }
        }

        SetBranch("jets", "label",  label);
        SetBranch("jets", "weight", weight);
        //-- kinematics
        SetBranch("jets", "pt",       (float)fj.Pt());
        SetBranch("jets", "eta",      (float)fj.Eta());
        SetBranch("jets", "phi",      (float)fj.Phi());
        SetBranch("jets", "mass",     (float)fj.M());
        SetBranch("jets", "sdmass",   fj.SDMass());
        SetBranch("jets", "masscorr", fj.ParticleNet_MassCorr());
        //-- substructure
        SetBranch("jets", "tau1", fj.Tau1());
        SetBranch("jets", "tau2", fj.Tau2());
        SetBranch("jets", "tau3", fj.Tau3());
        SetBranch("jets", "tau4", fj.Tau4());
        SetBranch("jets", "n2b1", fj.N2b1());
        SetBranch("jets", "n3b1", fj.N3b1());
        SetBranch("jets", "lsf3", fj.LSF3());
        //-- multiplicities & energy fractions
        SetBranch("jets", "nconst", (int)fj.NConstituents());
        SetBranch("jets", "chmult", (int)fj.ChMultiplicity());
        SetBranch("jets", "nemult", (int)fj.NeMultiplicity());
        SetBranch("jets", "chhef",  fj.ChHEF());
        SetBranch("jets", "nehef",  fj.NeHEF());
        SetBranch("jets", "chemef", fj.ChEmEF());
        SetBranch("jets", "neemef", fj.NeEmEF());
        SetBranch("jets", "muef",   fj.MuEF());
        //-- subjets
        SetBranch("jets", "sj1_pt",   sj1_pt);
        SetBranch("jets", "sj1_mass", sj1_mass);
        SetBranch("jets", "sj1_btag", sj1_btag);
        SetBranch("jets", "sj2_pt",   sj2_pt);
        SetBranch("jets", "sj2_mass", sj2_mass);
        SetBranch("jets", "sj2_btag", sj2_btag);
        SetBranch("jets", "sj_dr",    sj_dr);
        SetBranch("jets", "sj_z",     sj_z);
        //-- secondary vertices
        SetBranch("jets", "nsv",           nsv);
        SetBranch("jets", "sv_summass",    sv_summass);
        SetBranch("jets", "sv_maxdlensig", sv_maxdlensig);
        SetBranch("jets", "sv_sumntracks", sv_sumntracks);
        for (int k = 0; k < 3; k++) {
            const TString p = TString::Format("sv%d_", k + 1);
            SetBranch("jets", p + "pt",      svtok[k][0]);
            SetBranch("jets", p + "mass",    svtok[k][1]);
            SetBranch("jets", p + "dlensig", svtok[k][2]);
            SetBranch("jets", p + "ntracks", svtok[k][3]);
        }
        //-- stored-tagger benchmarks
        SetBranch("jets", "pnet_wvsqcd",   fj.ParticleNetWithMass_WvsQCD());
        SetBranch("jets", "pnet_tvsqcd",   fj.ParticleNetWithMass_TvsQCD());
        SetBranch("jets", "pnet_xqqvsqcd", fj.ParticleNet_XqqVsQCD());
        SetBranch("jets", "pnet_xggvsqcd", fj.ParticleNet_XggVsQCD());
        //-- gen truth of the match
        SetBranch("jets", "gen_wpt",  mW ? (float)mW->w->Pt() : -1.f);
        SetBranch("jets", "gen_drqq", mW ? (float)mW->q1->DeltaR(*mW->q2) : -1.f);
        SetBranch("jets", "gen_drwj", mW ? (float)fj.DeltaR(*mW->w) : -1.f);
        SetBranch("jets", "gen_hasb", hasB);
        FillTrees("jets");
    }
}
