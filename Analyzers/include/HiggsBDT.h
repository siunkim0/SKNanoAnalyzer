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
#include <string>
#include <vector>

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
        float muon_iso_max     =  0.35;  // pfRelIso04
        float muon_dxy_max     =  0.5;   // cm   } mirror src/skim.py::apply_muon_id
        float muon_dz_max      =  1.0;   // cm   } (HZZ loose muon ID + IP cuts)
        float muon_sip3d_max   =  4.0;   //      }
        float mZ1_min = 40.0,   mZ1_max = 120.0;
        float mZ2_min = 12.0,   mZ2_max = 120.0;
        float m4l_min = 70.0,   m4l_max = 1000.0;  // training window (selection.yaml)
        // Higgs signal region used for the high-purity BDT histograms.
        float m4l_sr_min = 105.0;
        float m4l_sr_max = 140.0;
    } cuts;

    // Feature order is NOT hardcoded. It is read at run time from the
    // <model_stem>_features.txt file written next to the .onnx by the BDT
    // project's scripts/export_onnx.py. That file is the single source of
    // truth (it mirrors /data6/Users/snuintern2/BDT/src/features.py::FEATURES).
    // executeEventFromParameter computes every known variable into a name->value
    // map and assembles the input tensor in this order, so changing the feature
    // set in features.py needs only a re-export + redeploy, no C++ edit.
    std::vector<std::string> featureNames;

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
