#ifndef Higgs_h
#define Higgs_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"

class Higgs : public AnalyzerCore {
public:
    Higgs();
    ~Higgs();

    void initializeAnalyzer();
    void executeEvent();
    void executeEventFromParameter();

    // Analysis flags
    bool RunSyst;
    
    // Trigger settings
    TString IsoMuTriggerName;
    float TriggerSafePtCut;
    float TriggerSafePtCut_Leading;
    float TriggerSafePtCut_Subleading;
    
    // Object ID settings
    RVec<Muon::MuonID> MuonIDs;
    RVec<TString> MuonIDSFKeys;
    
    
    // Physics objects
    RVec<Muon> AllMuons;
    RVec<Jet> AllJets;
    
    RVec<Muon> selectedMuons;
    float dilepton_mass;
    
    // Analysis cuts
    // Mirror HiggsBDT.h::AnalysisCuts so the cut-based preselection starts
    // from the SAME event population as the BDT path (fair comparison).
    struct AnalysisCuts {
        float muon_pt_lead    = 26.0;
        float muon_pt_sublead = 10.0;
        float muon_pt_other   =  5.0;
        float muon_eta        =  2.4;
        float muon_iso_max    =  0.35;   // pfRelIso04
        float muon_dxy_max    =  0.5;    // cm   } mirror src/skim.py::apply_muon_id
        float muon_dz_max     =  1.0;    // cm   } (HZZ loose muon ID + IP cuts)
        float muon_sip3d_max  =  4.0;    //      }
        float mZ1_min = 40.0,  mZ1_max = 120.0;
        float mZ2_min = 12.0,  mZ2_max = 120.0;
        float m4l_min = 70.0,  m4l_max = 1000.0;
    } cuts;
    
    // Systematic helper
    unique_ptr<SystematicHelper> systHelper;
    
    
    // Helper functions
    // RVec<Muon> SelectMuons(const RVec<Muon>& muons);

};

#endif