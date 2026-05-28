#include "DiPhoton.h"

DiPhoton::DiPhoton() {}
DiPhoton::~DiPhoton() {}

void DiPhoton::initializeAnalyzer() {
    // HLT_DoublePhoton70 OR HLT_DoublePhoton85
    DiPhoTriggers = {"HLT_DoublePhoton70"};

    // cut-based Loose, pT > 20 GeV, |eta|<2.5 (EB-EE gap vetoed inside PassID)
    PhotonID     = "POGLoose";
    PhotonPtCut  = 20.;
    PhotonEtaCut = 2.5;

    // MyCorrection needed for MCweight() / trigger lumi bookkeeping
    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);

    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", DataStream, DataEra);
    } else {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/ExampleSystematic.yaml", MCSample, DataEra);
    }

    cout << "[DiPhoton::initializeAnalyzer] DataEra = " << DataEra << endl;
    cout << "[DiPhoton::initializeAnalyzer] Triggers = HLT_DoublePhoton70 || HLT_DoublePhoton85" << endl;
    cout << "[DiPhoton::initializeAnalyzer] Photon ID = " << PhotonID
         << ", pT > " << PhotonPtCut << ", |eta| < " << PhotonEtaCut << endl;
}

void DiPhoton::executeEvent() {
    for (const auto &syst_dummy : *systHelper) {
        executeEventFromParameter();
    }
}

void DiPhoton::executeEventFromParameter() {
    const TString this_syst = systHelper->getCurrentSysName();
    Event ev = GetEvent();
    float weight = 1.0;

    // Trigger selection
    if (!ev.PassTrigger(DiPhoTriggers)) return;

    // Build photon collection with ID + kinematic selection
    RVec<Photon> photons = GetPhotons(PhotonID, PhotonPtCut, PhotonEtaCut);

    // Need at least 2 photons
    if (photons.size() < 2) return;


    // Sort in pT-descending order
    sort(photons.begin(), photons.end(),
         [](const Photon &a, const Photon &b) { return a.Pt() > b.Pt(); });

    const Photon &g1 = photons.at(0);
    const Photon &g2 = photons.at(1);

    if (g1.Pt() < 90. || g2.Pt() < 75.) return;

    // Event weight
    if (!IsDATA) {
        weight = MCweight() * ev.GetTriggerLumi("Full");
    }

    // Number of selected photons
    FillHist(this_syst + "/CutFlow", 2, 1 , 10, 0 , 10 );

    FillHist(this_syst + "/nPhotons", photons.size(), weight, 10, 0., 10.);

    // Per-photon kinematics (leading = 1, sub-leading = 2)
    FillHist(this_syst + "/photon1/pt",     g1.Pt(),       weight, 200, 0., 1000.);
    FillHist(this_syst + "/photon1/eta",    g1.Eta(),      weight, 50, -2.5, 2.5);
    FillHist(this_syst + "/photon1/scEta",  g1.scEta(),    weight, 50, -2.5, 2.5);
    FillHist(this_syst + "/photon1/phi",    g1.Phi(),      weight, 64, -3.2, 3.2);
    FillHist(this_syst + "/photon1/energy", g1.energy(),   weight, 200, 0., 2000.);

    FillHist(this_syst + "/photon2/pt",     g2.Pt(),       weight, 200, 0., 1000.);
    FillHist(this_syst + "/photon2/eta",    g2.Eta(),      weight, 50, -2.5, 2.5);
    FillHist(this_syst + "/photon2/scEta",  g2.scEta(),    weight, 50, -2.5, 2.5);
    FillHist(this_syst + "/photon2/phi",    g2.Phi(),      weight, 64, -3.2, 3.2);
    FillHist(this_syst + "/photon2/energy", g2.energy(),   weight, 200, 0., 2000.);

    // Diphoton system
    Particle pair = g1 + g2;
    FillHist(this_syst + "/pair/mass", pair.M(),   weight, 300, 0., 3000.);
    FillHist(this_syst + "/pair/pt",   pair.Pt(),  weight, 200, 0., 1000.);
    FillHist(this_syst + "/pair/eta",  pair.Eta(), weight, 100, -5., 5.);
    FillHist(this_syst + "/pair/phi",  pair.Phi(), weight, 64, -3.2, 3.2);

    // Angular separations between the two photons
    const float dR   = g1.DeltaR(g2);
    const float dPhi = fabs(g1.DeltaPhi(g2));
    const float dEta = fabs(g1.Eta() - g2.Eta());
    FillHist(this_syst + "/pair/deltaR",   dR,   weight, 100, 0., 10.);
    FillHist(this_syst + "/pair/deltaPhi", dPhi, weight, 64, 0., 3.2);
    FillHist(this_syst + "/pair/deltaEta", dEta, weight, 100, 0., 5.);
}
