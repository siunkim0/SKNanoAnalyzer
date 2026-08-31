#ifndef WtagAK8Pt_h
#define WtagAK8Pt_h

#include "AnalyzerCore.h"

//==============================================================
// Inclusive AK8 (FatJet) pT spectrum — light companion to Wtag.
//
//   Purpose: the "all QCD jets" denominator for the boosted-W
//   fraction study. NanoAOD stores FatJets only down to pT = 170
//   GeV (hard storage floor), so this IS the full reconstructable
//   AK8 spectrum. To populate the 170-200 region without the
//   pT-hat>170 bias, run over the *inclusive* QCD_Pt suite
//   including the low-pT-hat bins (80to120, 120to170) that feed
//   AK8 jets up across 170 GeV.
//
//   No gen truth, no ntuple — just the xsec-weighted pT spectrum,
//   so it is fast and light enough to run the huge low-pT-hat bins.
//==============================================================
class WtagAK8Pt : public AnalyzerCore {
public:
    WtagAK8Pt();
    ~WtagAK8Pt();

    void initializeAnalyzer();
    void executeEvent();

    float ak8_eta_max = 2.4;   // tracker/tagger acceptance
};

#endif
