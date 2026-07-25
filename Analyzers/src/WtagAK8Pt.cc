#include "WtagAK8Pt.h"
#include <cmath>

//==== Constructor and Destructor
WtagAK8Pt::WtagAK8Pt() {}
WtagAK8Pt::~WtagAK8Pt() {}

//==============================================================
// Initialize
//==============================================================
void WtagAK8Pt::initializeAnalyzer() {}

//==============================================================
// Event loop  (MC only — inclusive AK8 pT spectrum)
//==============================================================
void WtagAK8Pt::executeEvent() {

    if (IsDATA) return;

    RVec<FatJet> fatJets = GetAllFatJets();

    //==== event weight (sign x genWeight, xsec-normalized) — same as Wtag,
    //     so spectra stack consistently across pT-hat bins.
    float weight = MCweight();

    for (const auto &fj : fatJets) {
        //---- inclusive: every stored AK8 jet (NanoAOD floor pT >= 170)
        FillHist("FatJet/Pt_all",   fj.Pt(), weight, 300, 0., 3000.);
        //---- within tracker/tagger acceptance (matches the v2 jets-tree cut)
        if (std::abs(fj.Eta()) < ak8_eta_max)
            FillHist("FatJet/Pt_eta24", fj.Pt(), weight, 300, 0., 3000.);
    }
}
