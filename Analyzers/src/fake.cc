#include "fake.h"

fake::fake() : MeasFakeMu8(false), MeasFakeMu17(false), RunSyst(false) {}
fake::~fake() {}

void fake::initializeAnalyzer() {
    // Userflags
    MeasFakeMu8 = HasFlag("MeasFakeMu8");
    MeasFakeMu17 = HasFlag("MeasFakeMu17");
    RunSyst = HasFlag("RunSyst");

    // If no path is specified, measure both
    if (!MeasFakeMu8 && !MeasFakeMu17) {
        MeasFakeMu8 = true;
        MeasFakeMu17 = true;
    }

    // Prescaled single muon triggers with the pT cuts and
    // cone-corrected pT regions from the analysis note (Table 38)
    if (MeasFakeMu8)  paths.push_back({"Mu8",  "HLT_Mu8_TrkIsoVVL",  10., 10.});
    if (MeasFakeMu17) paths.push_back({"Mu17", "HLT_Mu17_TrkIsoVVL", 20., 30.});

    systs = {"Central"};
    if (RunSyst) {
        systs.push_back("AwayJetPt30");
        systs.push_back("AwayJetPt60");
    }

    // Fake rate binning: cone-corrected pT x |eta|
    ptCorrBins = {10., 15., 20., 30., 50., 100.};
    absEtaBins = {0., 0.9, 1.6, 2.4};

    // Loose electron ID for the lepton veto
    electronLooseID = (Run == 2) ? "HcToWALooseRun2" : "HcToWALooseRun3";

    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);
}

// WPs from the analysis note:
// common      : POG medium ID, |dz| < 0.1 cm, rel. tracker iso (R03) < 0.4 (trigger emulation)
// tight       : |SIP3D| < 3, rel. mini-iso < 0.1
// loose (Run2): |SIP3D| < 5, rel. mini-iso < 0.6
bool fake::PassMuonWP(const Muon &mu, const TString &wp) const {
    if (!mu.isPOGMediumId()) return false;
    if (!(fabs(mu.dZ()) < 0.1)) return false;
    if (!(mu.TkRelIso() < 0.4)) return false;
    if (wp == "tight") {
        if (!(fabs(mu.SIP3D()) < 3.)) return false;
        if (!(mu.MiniPFRelIso() < 0.1)) return false;
    } else {
        if (!(fabs(mu.SIP3D()) < 5.)) return false;
        if (!(mu.MiniPFRelIso() < 0.6)) return false;
    }
    return true;
}

void fake::executeEvent() {
    Event ev = GetEvent();
    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;

    RVec<Muon> rawMuons = GetAllMuons();
    RVec<Electron> rawElectrons = GetAllElectrons();
    sort(rawMuons.begin(), rawMuons.end(), PtComparing);

    RVec<Muon> looseMuons, tightMuons;
    for (const auto &mu : rawMuons) {
        if (mu.Pt() < 10. || fabs(mu.Eta()) > 2.4) continue;
        if (!PassMuonWP(mu, "loose")) continue;
        looseMuons.push_back(mu);
        if (PassMuonWP(mu, "tight")) tightMuons.push_back(mu);
    }
    RVec<Electron> looseElectrons = SelectElectrons(rawElectrons, electronLooseID, 15., 2.5);

    const float maxJetEta = DataEra.Contains("2016") ? 2.4 : 2.5;
    RVec<Jet> jets = SelectJets(rawJets, "tight", 25., maxJetEta);
    if (Run == 2) jets = SelectJets(jets, "loosePuId", 25., maxJetEta);
    jets = JetsVetoLeptonInside(jets, looseElectrons, looseMuons, 0.4);
    sort(jets.begin(), jets.end(), PtComparing);

    Particle METv = ApplyTypeICorrection(ev.GetMETVector(Event::MET_Type::PUPPI),
                                         rawJets, rawElectrons, rawMuons);

    RVec<Gen> gens = IsDATA ? RVec<Gen>() : GetAllGens();

    for (const auto &path : paths) {
        if (!ev.PassTrigger(path.trigger)) continue;

        float weight = 1.;
        if (!IsDATA) {
            // MC normalisation is only provisional here (prescales are not in
            // the trigger lumi); the final normalisation per trigger path is
            // extracted offline from the ZEnriched region
            weight = MCweight() * ev.GetTriggerLumi(path.trigger);
            weight *= GetL1PrefireWeight();
            weight *= myCorr->GetPUWeight(ev.nTrueInt());
        }

        measureFakeRate(path, looseMuons, looseElectrons, jets, METv, weight, gens);
        fillZEnriched(path, ev, looseMuons, tightMuons, looseElectrons, jets, weight, gens);
    }
}

void fake::measureFakeRate(const TriggerPath &path,
                           const RVec<Muon> &looseMuons,
                           const RVec<Electron> &looseElectrons,
                           const RVec<Jet> &jets,
                           const Particle &METv,
                           float weight,
                           const RVec<Gen> &gens) {
    // exactly one loose ID lepton
    if (looseMuons.size() != 1) return;
    if (looseElectrons.size() != 0) return;

    const Muon &mu = looseMuons.at(0);
    if (mu.Pt() < path.ptCut) return;

    // cone-corrected pT: pt * (1 + max(0, I_mini - I_mini^tight))
    const float ptCorr = mu.Pt() * (1. + max(0.f, mu.MiniPFRelIso() - 0.1f));
    if (ptCorr < path.ptCorrCut) return;

    const bool isTight = PassMuonWP(mu, "tight");
    const TString ltype = IsDATA ? "" : LeptonTypeToString(GetLeptonType(mu, gens));

    const float MT = sqrt(2. * mu.Pt() * METv.Pt() * (1. - cos(mu.DeltaPhi(METv))));

    for (const auto &syst : systs) {
        float awayJetPtCut = 40.;
        if (syst == "AwayJetPt30") awayJetPtCut = 30.;
        else if (syst == "AwayJetPt60") awayJetPtCut = 60.;

        // require at least one away jet: pT > cut && dR(mu, j) > 0.7
        int nJets = 0;
        bool hasAwayJet = false;
        for (const auto &jet : jets) {
            if (jet.Pt() < awayJetPtCut) continue;
            nJets++;
            if (jet.DeltaR(mu) > 0.7) hasAwayJet = true;
        }
        if (!hasAwayJet) continue;

        // QCD-enriched validation region: MET / MT cuts removed
        if (syst == "Central") {
            fillMuonKinematics(path.name + "/Inclusive/loose", mu, ptCorr, METv, MT, nJets, weight);
            if (!IsDATA)
                fillMuonKinematics(path.name + "/Inclusive/loose/" + ltype, mu, ptCorr, METv, MT, nJets, weight);
            if (isTight) {
                fillMuonKinematics(path.name + "/Inclusive/tight", mu, ptCorr, METv, MT, nJets, weight);
                if (!IsDATA)
                    fillMuonKinematics(path.name + "/Inclusive/tight/" + ltype, mu, ptCorr, METv, MT, nJets, weight);
            }
        }

        // measurement region
        if (!(METv.Pt() < 25.)) continue;
        if (!(MT < 25.)) continue;

        // keep overflow in the last bin
        const float xval = min(ptCorr, ptCorrBins.back() - 0.1f);
        const float yval = fabs(mu.Eta());

        const TString base = path.name + "/MeasReg/" + syst;
        FillHist(base + "/loose/fake_ptcorr_abseta", xval, yval, weight, ptCorrBins, absEtaBins);
        if (!IsDATA)
            FillHist(base + "/loose/" + ltype + "/fake_ptcorr_abseta", xval, yval, weight, ptCorrBins, absEtaBins);
        if (isTight) {
            FillHist(base + "/tight/fake_ptcorr_abseta", xval, yval, weight, ptCorrBins, absEtaBins);
            if (!IsDATA)
                FillHist(base + "/tight/" + ltype + "/fake_ptcorr_abseta", xval, yval, weight, ptCorrBins, absEtaBins);
        }

        if (syst == "Central") {
            fillMuonKinematics(base + "/loose", mu, ptCorr, METv, MT, nJets, weight);
            if (!IsDATA)
                fillMuonKinematics(base + "/loose/" + ltype, mu, ptCorr, METv, MT, nJets, weight);
            if (isTight) {
                fillMuonKinematics(base + "/tight", mu, ptCorr, METv, MT, nJets, weight);
                if (!IsDATA)
                    fillMuonKinematics(base + "/tight/" + ltype, mu, ptCorr, METv, MT, nJets, weight);
            }
        }
    }
}

// Z-enriched region to extract the MC normalisation for each
// prescaled trigger path: OS tight pair, |M(ll) - 91.2| < 15,
// at least one jet with pT > 40
void fake::fillZEnriched(const TriggerPath &path,
                         const Event &ev,
                         const RVec<Muon> &looseMuons,
                         const RVec<Muon> &tightMuons,
                         const RVec<Electron> &looseElectrons,
                         const RVec<Jet> &jets,
                         float weight,
                         const RVec<Gen> &gens) {
    if (!(looseMuons.size() == 2 && tightMuons.size() == 2)) return;
    if (!looseElectrons.empty()) return;

    const Muon &mu1 = tightMuons.at(0);
    const Muon &mu2 = tightMuons.at(1);
    if (mu1.Charge() + mu2.Charge() != 0) return;
    if (mu1.Pt() < path.ptCut || mu2.Pt() < 10.) return;

    bool hasJet = false;
    for (const auto &jet : jets) {
        if (jet.Pt() > 40.) { hasJet = true; break; }
    }
    if (!hasJet) return;

    const float mll = (mu1 + mu2).M();
    if (fabs(mll - 91.2) > 15.) return;

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

void fake::fillMuonKinematics(const TString &prefix,
                              const Muon &mu, float ptCorr,
                              const Particle &METv, float MT,
                              int nJets, float weight) {
    FillHist(prefix + "/pt", mu.Pt(), weight, 200, 0., 200.);
    FillHist(prefix + "/ptcorr", ptCorr, weight, 200, 0., 200.);
    FillHist(prefix + "/eta", mu.Eta(), weight, 48, -2.4, 2.4);
    FillHist(prefix + "/phi", mu.Phi(), weight, 64, -3.2, 3.2);
    FillHist(prefix + "/MET", METv.Pt(), weight, 100, 0., 100.);
    FillHist(prefix + "/MT", MT, weight, 100, 0., 100.);
    FillHist(prefix + "/nJets", nJets, weight, 10, 0., 10.);
}

// LeptonType > 0: prompt (EW/BSM prompt, tau daughter, internal conversion)
// LeptonType <= 0: hadronic origin / external conversion / unmatched
TString fake::LeptonTypeToString(int leptonType) const {
    return (leptonType > 0) ? "prompt" : "fake";
}
