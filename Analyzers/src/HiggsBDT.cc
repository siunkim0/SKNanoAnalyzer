#include "HiggsBDT.h"

#include <TLorentzVector.h>
#include <TVector3.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

// ---------------------------------------------------------------------------
// Helicity-angle helpers — mirror src/features.py::add_helicity_angles in the
// BDT project. Sign conventions follow arXiv:1208.4018. A flipped sign here
// will silently shift the BDT score; check the parity CSV from the notebook.
// ---------------------------------------------------------------------------

namespace {

// Return the spatial 3-vector of `p` after boosting into `parent`'s rest frame.
TVector3 boostInto(const TLorentzVector& p, const TLorentzVector& parent) {
    TLorentzVector q = p;
    q.Boost(-parent.BoostVector());
    return q.Vect();
}

// Signed angle from n1 to n2 around `axis`. All inputs assumed unit vectors.
double signedAngle(const TVector3& n1, const TVector3& n2, const TVector3& axis) {
    double cosA = std::clamp(n1.Dot(n2), -1.0, 1.0);
    double s    = axis.Dot(n1.Cross(n2));
    if (s == 0.0) s = 1.0;     // matches np.where(sign==0, 1, sign) in Python
    return std::copysign(std::acos(cosA), s);
}

}  // namespace

// ---------------------------------------------------------------------------
HiggsBDT::HiggsBDT()  {}
HiggsBDT::~HiggsBDT() {}

void HiggsBDT::initializeAnalyzer() {

    RunSyst = HasFlag("RunSyst");

    // Trigger OR — must match config/selection.yaml::triggers in the BDT project.
    BDTTriggers = {
        "HLT_IsoMu24"
    };

    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);

    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", DataStream, DataEra);
    } else {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/ExampleSystematic.yaml", MCSample, DataEra);
    }

    // Locate the ONNX model. By convention we drop it under
    //   $SKNANO_DATA/<DataEra>/HZZ4mu/bdt_v1.onnx
    // (the export instructions in the BDT-project README §7.1).
    const char* sknano_data = getenv("SKNANO_DATA");
    if (sknano_data == nullptr) {
        std::cerr << "[HiggsBDT] SKNANO_DATA env var not set; BDT inference disabled\n";
        return;
    }
    bdtModelPath = TString(sknano_data) + "/" + DataEra + "/BDT/HZZ4mu/bdt_v1_v1.onnx";

    std::ifstream probe(bdtModelPath.Data());
    if (!probe.good()) {
        std::cerr << "[HiggsBDT] BDT model not found at " << bdtModelPath
                  << " — BDT inference disabled, selection-only mode\n";
        return;
    }
    // Read the feature order exported alongside the model. This is the single
    // source of truth shared with src/features.py::FEATURES — never hardcode it.
    TString featListPath = bdtModelPath;
    featListPath.ReplaceAll(".onnx", "_features.txt");
    std::ifstream flist(featListPath.Data());
    if (!flist.good()) {
        std::cerr << "[HiggsBDT] feature list not found at " << featListPath
                  << " — BDT inference disabled\n";
        return;
    }
    featureNames.clear();
    std::string line;
    while (std::getline(flist, line)) {
        // Trim trailing CR/space (handles CRLF and stray whitespace).
        while (!line.empty() &&
               (line.back() == '\r' || line.back() == ' ' || line.back() == '\t')) {
            line.pop_back();
        }
        if (!line.empty()) featureNames.push_back(line);
    }
    if (featureNames.empty()) {
        std::cerr << "[HiggsBDT] feature list " << featListPath
                  << " is empty — BDT inference disabled\n";
        return;
    }

    bdt = std::make_unique<MLHelper>(bdtModelPath.Data(), MLHelper::ModelType::ONNX);
    std::cout << "[HiggsBDT] loaded BDT model from " << bdtModelPath
              << " with " << featureNames.size() << " features from "
              << featListPath << std::endl;
}

void HiggsBDT::executeEvent() {
    AllMuons = GetAllMuons();
    AllJets  = GetAllJets();

    for (const auto& syst_dummy : *systHelper) {
        executeEventFromParameter();
    }
}

void HiggsBDT::executeEventFromParameter() {
    float weight = 1.0;
    Event ev = GetEvent();
    const TString this_syst = systHelper->getCurrentSysName();

    if (!IsDATA) {
        weight *= MCweight();
        weight *= ev.GetTriggerLumi("Full");
        weight *= myCorr->GetPUWeight(ev.nTrueInt());
    }

    // 1. cutflow bin 0 = nocuts
    FillHist(this_syst + "/CutFlow", 0, weight, 10, 0, 10);

    // 2. trigger OR — same list as the training skim
    bool pass_trig = false;
    for (const auto& t : BDTTriggers) {
        if (ev.PassTrigger(t)) { pass_trig = true; break; }
    }
    if (!pass_trig) return;
    FillHist(this_syst + "/CutFlow", 1, weight, 10, 0, 10);

    // 3. muon ID + iso, exactly 4
    auto muons_tight = SelectMuons(AllMuons, "POGTight", cuts.muon_pt_other, cuts.muon_eta);
    RVec<Muon> selectedMuons;
    // src/skim.py uses pfRelIso04_all < 0.35 (HZZ tight).
    for (const auto& mu : muons_tight) {
        if (mu.PfRelIso04() < cuts.muon_iso_max) selectedMuons.push_back(mu);
    }
    if (selectedMuons.size() != 4) return;
    sort(selectedMuons.begin(), selectedMuons.end(), PtComparing);
    FillHist(this_syst + "/CutFlow", 2, weight, 10, 0, 10);

    // 4. pT thresholds (20/10/5/5) and total charge 0
    if (selectedMuons[0].Pt() < cuts.muon_pt_lead)    return;
    if (selectedMuons[1].Pt() < cuts.muon_pt_sublead) return;
    if (selectedMuons[2].Pt() < cuts.muon_pt_other)   return;
    if (selectedMuons[3].Pt() < cuts.muon_pt_other)   return;
    if (selectedMuons[0].Charge() + selectedMuons[1].Charge() +
        selectedMuons[2].Charge() + selectedMuons[3].Charge() != 0) return;
    FillHist(this_syst + "/CutFlow", 3, weight, 10, 0, 10);

    // 5. Z1 = OS pair closest to mZ
    constexpr double mZ_PDG = 91.1876;
    int i1 = -1, j1 = -1;
    double best = 1e9;
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            if (selectedMuons[i].Charge() + selectedMuons[j].Charge() != 0) continue;
            double m = (selectedMuons[i] + selectedMuons[j]).M();
            double diff = std::fabs(m - mZ_PDG);
            if (diff < best) { best = diff; i1 = i; j1 = j; }
        }
    }
    if (i1 < 0) return;

    int i2 = -1, j2 = -1;
    for (int k = 0; k < 4; ++k) {
        if (k == i1 || k == j1) continue;
        if (i2 < 0) i2 = k; else j2 = k;
    }

    TLorentzVector p4_Z1 = selectedMuons[i1] + selectedMuons[j1];
    TLorentzVector p4_Z2 = selectedMuons[i2] + selectedMuons[j2];
    TLorentzVector p4_H  = p4_Z1 + p4_Z2;

    double mZ1   = p4_Z1.M();
    double mZ2   = p4_Z2.M();
    double m4l   = p4_H.M();

    if (mZ1 < cuts.mZ1_min || mZ1 > cuts.mZ1_max) return;
    if (mZ2 < cuts.mZ2_min || mZ2 > cuts.mZ2_max) return;
    if (m4l < cuts.m4l_min || m4l > cuts.m4l_max) return;
    FillHist(this_syst + "/CutFlow", 4, weight, 10, 0, 10);

    // 6. Compute every known variable, then assemble the input tensor in the
    //    order given by featureNames (read from <model>_features.txt, which
    //    mirrors src/features.py::FEATURES). No index bookkeeping here — adding
    //    or reordering features in features.py needs no change to this code.

    // Helicity angles (5). Convention: ℓ⁻ in each Z is the reference.
    const Muon& mZ1a = selectedMuons[i1];
    const Muon& mZ1b = selectedMuons[j1];
    const Muon& mZ2a = selectedMuons[i2];
    const Muon& mZ2b = selectedMuons[j2];
    const TLorentzVector& z1m = (mZ1a.Charge() == -1) ? (TLorentzVector)mZ1a : (TLorentzVector)mZ1b;
    const TLorentzVector& z2m = (mZ2a.Charge() == -1) ? (TLorentzVector)mZ2a : (TLorentzVector)mZ2b;

    // +z lab-frame beam reference (massless 4-vector along beam axis).
    TLorentzVector beam_lab(0.0, 0.0, 1.0, 1.0);

    TVector3 n_Z1_H   = boostInto(p4_Z1,   p4_H ).Unit();
    TVector3 n_Z2_H   = boostInto(p4_Z2,   p4_H ).Unit();
    TVector3 n_z1m_H  = boostInto(z1m,     p4_H ).Unit();
    TVector3 n_z2m_H  = boostInto(z2m,     p4_H ).Unit();
    TVector3 n_beam_H = boostInto(beam_lab,p4_H ).Unit();

    TVector3 n_l1m_Z1   = boostInto(z1m,   p4_Z1).Unit();
    TVector3 n_negZ2_Z1 = (-boostInto(p4_Z2, p4_Z1)).Unit();
    TVector3 n_l2m_Z2   = boostInto(z2m,   p4_Z2).Unit();
    TVector3 n_negZ1_Z2 = (-boostInto(p4_Z1, p4_Z2)).Unit();

    double cos_theta_star = n_Z1_H.Dot(n_beam_H);
    double cos_theta1     = n_l1m_Z1.Dot(n_negZ2_Z1);
    double cos_theta2     = n_l2m_Z2.Dot(n_negZ1_Z2);

    TVector3 n_plane_Z1 = (n_Z1_H.Cross(n_z1m_H)).Unit();
    TVector3 n_plane_Z2 = (n_Z2_H.Cross(n_z2m_H)).Unit();
    TVector3 n_plane_sc = (n_Z1_H.Cross(n_beam_H)).Unit();

    double Phi  = signedAngle(n_plane_Z1, n_plane_Z2, n_Z1_H);
    double Phi1 = signedAngle(n_plane_Z1, n_plane_sc, n_Z1_H);

    // Named feature registry. Keys MUST match the column names in
    // src/features.py (i.e. the strings in features.py::FEATURES). To add a new
    // input variable: register it here once and add it to features.py — the run
    // -time order then comes from the exported _features.txt automatically.
    const std::map<std::string, float> featMap = {
        {"m4l",            static_cast<float>(m4l)},
        {"mZ1",            static_cast<float>(mZ1)},
        {"mZ2",            static_cast<float>(mZ2)},
        {"pt4l",           static_cast<float>(p4_H.Pt())},
        {"eta4l",          static_cast<float>(p4_H.Eta())},
        {"pt_Z1",          static_cast<float>(p4_Z1.Pt())},
        {"pt_Z2",          static_cast<float>(p4_Z2.Pt())},
        {"dR_Z1Z2",        static_cast<float>(p4_Z1.DeltaR(p4_Z2))},
        {"pt_mu1",         static_cast<float>(selectedMuons[0].Pt())},
        {"pt_mu2",         static_cast<float>(selectedMuons[1].Pt())},
        {"pt_mu3",         static_cast<float>(selectedMuons[2].Pt())},
        {"pt_mu4",         static_cast<float>(selectedMuons[3].Pt())},
        {"eta_mu1",        static_cast<float>(selectedMuons[0].Eta())},
        {"eta_mu2",        static_cast<float>(selectedMuons[1].Eta())},
        {"eta_mu3",        static_cast<float>(selectedMuons[2].Eta())},
        {"eta_mu4",        static_cast<float>(selectedMuons[3].Eta())},
        {"pt_Z1_over_m4l", static_cast<float>(p4_Z1.Pt() / m4l)},
        {"pt_Z2_over_m4l", static_cast<float>(p4_Z2.Pt() / m4l)},
        {"cos_theta_star", static_cast<float>(cos_theta_star)},
        {"cos_theta1",     static_cast<float>(cos_theta1)},
        {"cos_theta2",     static_cast<float>(cos_theta2)},
        {"Phi",            static_cast<float>(Phi)},
        {"Phi1",           static_cast<float>(Phi1)},
    };

    // Assemble the input vector in the exact training order.
    FloatArray feats(featureNames.size());
    for (size_t k = 0; k < featureNames.size(); ++k) {
        auto it = featMap.find(featureNames[k]);
        if (it == featMap.end()) {
            std::cerr << "[HiggsBDT] feature '" << featureNames[k]
                      << "' from the feature list is not registered in "
                      << "HiggsBDT.cc — add it to featMap\n";
            return;
        }
        feats[k] = it->second;
    }

    // 7. Selection-stage histograms (always filled, regardless of BDT availability).
    FillHist(this_syst + "/H_Mass",  m4l, weight, 125, 0., 250.);
    FillHist(this_syst + "/Z1_Mass", mZ1, weight, 100, 0., 200.);
    FillHist(this_syst + "/Z2_Mass", mZ2, weight, 100, 0., 200.);
    FillHist(this_syst + "/H_Pt",    p4_H.Pt(), weight, 50, 0., 200.);

    // Path-1: only score events inside the training SR. v5 never saw
    // m4l outside [105,140]; scoring there is extrapolation noise.
    if (m4l < cuts.m4l_sr_min || m4l > cuts.m4l_sr_max) return;

    // 8. BDT inference. If the model wasn't loaded, we stop here.
    if (!bdt) return;

    std::unordered_map<std::string, VariousArray> in {
        { bdtInputName.Data(), feats }
    };
    std::unordered_map<std::string, IntArray> in_shape {
        { bdtInputName.Data(), IntArray{1, static_cast<int>(featureNames.size())} }
    };
    auto out = bdt->Run_ONNX_Model(in, in_shape);

    auto it = out.find(bdtOutputName.Data());
    if (it == out.end()) {
        // First event: print what the model actually exposes so the user can
        // fix `bdtOutputName` once and rebuild.
        std::cerr << "[HiggsBDT] output node '" << bdtOutputName << "' not found. Got: ";
        for (const auto& kv : out) std::cerr << kv.first << " ";
        std::cerr << std::endl;
        return;
    }
    // zipmap=False export → tensor of shape (1, 2); column 1 is P(signal).
    if (it->second.size() < 2) {
        std::cerr << "[HiggsBDT] expected probability vector of size >= 2, got "
                  << it->second.size() << std::endl;
        return;
    }
    float p_sig = it->second[1];

    FillHist(this_syst + "/CutFlow",  5, weight, 10, 0, 10);
    FillHist(this_syst + "/BDT_Score", p_sig, weight, 50, 0., 1.);
    FillHist(this_syst + "/BDT_vs_m4l", p_sig, m4l, weight, 50, 0., 1., 60, 70., 250.);

    // 9. Working-point scan — one mass histogram per BDT threshold, single pass.
    static const double wps[] = {0.50, 0.55, 0.60, 0.65, 0.70,
                                 0.75, 0.80, 0.85, 0.90, 0.95};
    for (double wp : wps) {
        if (p_sig > wp) {
            int tag = (int)std::lround(wp * 100);
            FillHist(this_syst + Form("/H_Mass_BDT%03d", tag),
                     m4l, weight, 125, 0., 250.);
        }
    }
}
