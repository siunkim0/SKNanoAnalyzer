#ifndef HiggsBDT_h
#define HiggsBDT_h

// Example analyzer: H -> ZZ -> 4mu selection (mirrors Higgs.cc) + xgboost BDT
// inference via MLHelper. The BDT is trained in /data6/Users/snuintern2/BDT
// and exported to ONNX (notebooks/bdt_tutorial.ipynb section 11). Run-time
// model lookup expects the file at:
//   $SKNANO_DATA/<DataEra>/HZZ4mu/bdt_v1.onnx
//
// To enable in the build:
//   1) Analyzers/CMakeLists.txt:44  -> uncomment
//        target_link_libraries(Analyzers PUBLIC MLHelper)
//   2) Analyzers/include/AnalyzersLinkDef.hpp -> add
//        #pragma link C++ class HiggsBDT+;
//   3) ./scripts/rebuild.sh
//
// Run:
//   SKNano.py -a HiggsBDT -i 'GluGluHto2Zto4L*' -e 2022EE -n 10

#include "AnalyzerCore.h"
#include "SystematicHelper.h"
#include "MLHelper.h"

#include <memory>

class HiggsBDT : public AnalyzerCore {
public:
    HiggsBDT();
    ~HiggsBDT();

    void initializeAnalyzer();
    void executeEvent();
    void executeEventFromParameter();

    bool RunSyst;

    // BDT trigger OR — same list as the training skim (config/selection.yaml).
    RVec<TString> BDTTriggers;

    // Same kinematic cuts as the BDT skim, NOT Higgs.cc's IsoMu24 turn-on.
    struct AnalysisCuts {
        float muon_pt_lead     = 26.0;
        float muon_pt_sublead  = 10.0;
        float muon_pt_other    =  5.0;
        float muon_eta         =  2.4;
        float muon_iso_max     =  0.35;
        float mZ1_min = 40.0,   mZ1_max = 120.0;
        float mZ2_min = 12.0,   mZ2_max = 120.0;
        float m4l_min = 0.0,   m4l_max = 1000.0;
        // Higgs signal region used for the high-purity BDT histograms.
        float m4l_sr_min = 105.0;
        float m4l_sr_max = 140.0;
    } cuts;

    // The BDT (v5, Path 1 / SR-restricted) expects 23 features in the order
    // defined in /data6/Users/snuintern2/BDT/src/features.py::FEATURES
    // (also see bdt_v5_features.txt next to the .onnx).
    static constexpr int N_FEAT = 23;

    // Per-event work areas.
    RVec<Muon> AllMuons;
    RVec<Jet>  AllJets;

    // Owned by initializeAnalyzer; nullptr if the model file isn't present
    // (executeEventFromParameter then falls back to selection-only mode).
    std::unique_ptr<MLHelper> bdt;
    TString bdtModelPath;
    TString bdtInputName  = "features";       // matches the export in section 11.2
    TString bdtOutputName = "probabilities";  // verify with sess.get_outputs() once

    unique_ptr<SystematicHelper> systHelper;
};

#endif
