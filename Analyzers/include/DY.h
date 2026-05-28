#ifndef DY_h
#define DY_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"

class DY : public AnalyzerCore {
public:
    DY();
    ~DY();

    void initializeAnalyzer();
    void executeEvent();
    void executeEventFromParameter();

    // Analysis flags
    bool RunSyst;
    
    // Trigger settings
    TString IsoMuTriggerName;
    float TriggerSafePtCut;

    TString IsoEleTriggerName;
    float EleTriggerSafePtCut;
    
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
        float muon_pt_lead = 28.0;
        float muon_pt_sublead = 15.0;
        float muon_eta = 2.4;
    } cuts;
    
    // Systematic helper
    unique_ptr<SystematicHelper> systHelper;
    
    
    // Helper functions
    // RVec<Muon> SelectMuons(const RVec<Muon>& muons);

};

#endif