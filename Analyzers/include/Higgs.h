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
    struct AnalysisCuts {
        float muon_pt_lead = 26.0;
        float muon_pt_sublead = 10.0;
        float muon_pt_sub2lead = 7.0;
        float muon_pt_sub3lead = 4.0;
        float muon_eta = 2.4;
    } cuts;
    
    // Systematic helper
    unique_ptr<SystematicHelper> systHelper;
    
    
    // Helper functions
    // RVec<Muon> SelectMuons(const RVec<Muon>& muons);

};

#endif