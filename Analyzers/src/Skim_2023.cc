#include "Skim_2023.h"    
#include <utility>

Skim_2023::Skim_2023() {
    newtree = NULL; 
}

Skim_2023::~Skim_2023() {}

void Skim_2023::initializeAnalyzer() {
    GetOutfile()->cd();
    newtree = fChain->CloneTree(0);

    if (DataEra == "2018")
    {
        mu_set.Muon_Trigger = {"HLT_IsoMu24"}; 
        mu_set.Muon_Trigger_Safe_Pt_Cut = 26.;
        // el_set.Ele_Trigger = {"HLT_Ele32_WPTight_Gsf","HLT_Photon200","HLT_Ele115_CaloIdVT_GsfTrkIdT"};
        // el_set.Ele_Trigger_Safe_Pt_Cut = 35.;  
    }
    if (DataEra == "2022")
    {
        mu_set.Muon_Trigger = {"HLT_IsoMu24"}; 
        mu_set.Muon_Trigger_Safe_Pt_Cut = 26.;
        // el_set.Ele_Trigger = {"HLT_Ele32_WPTight_Gsf","HLT_Photon200","HLT_Ele115_CaloIdVT_GsfTrkIdT"};
        // el_set.Ele_Trigger_Safe_Pt_Cut = 35.;  
    }
    if (DataEra == "2022EE")
    {
        mu_set.Muon_Trigger = {"HLT_IsoMu24"}; 
        mu_set.Muon_Trigger_Safe_Pt_Cut = 26.;
        // el_set.Ele_Trigger = {"HLT_Ele32_WPTight_Gsf","HLT_Photon200","HLT_Ele115_CaloIdVT_GsfTrkIdT"};
        // el_set.Ele_Trigger_Safe_Pt_Cut = 35.; 
    }
    if (DataEra == "2023")
    {
        mu_set.Muon_Trigger = {"HLT_IsoMu24"}; 
        mu_set.Muon_Trigger_Safe_Pt_Cut = 26.;
        // el_set.Ele_Trigger = {"HLT_Ele32_WPTight_Gsf","HLT_Photon200","HLT_Ele115_CaloIdVT_GsfTrkIdT"};
        // el_set.Ele_Trigger_Safe_Pt_Cut = 35.; 
    }
    if (DataEra == "2023BPix")
    {
        mu_set.Muon_Trigger = {"HLT_IsoMu24"}; 
        mu_set.Muon_Trigger_Safe_Pt_Cut = 26.;
        // el_set.Ele_Trigger = {"HLT_Ele32_WPTight_Gsf","HLT_Photon200","HLT_Ele115_CaloIdVT_GsfTrkIdT"};
        // el_set.Ele_Trigger_Safe_Pt_Cut = 35.; 
    }

    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);

    
    // Initialize systematic helper
    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", DataStream, DataEra);
    } else {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/ExampleSystematic.yaml", MCSample, DataEra);
    }

}

void Skim_2023::executeEvent() {
    executeEventFromParameter();
}

void Skim_2023::executeEventFromParameter() {
    Event ev = GetEvent();

    bool pass_trig_muon = ev.PassTrigger(mu_set.Muon_Trigger);
    bool pass_trig_elec = ev.PassTrigger(el_set.Ele_Trigger);

    /////////Trugger pass ///////////
    if ((!pass_trig_muon) and (!pass_trig_elec)) return;

    
    newtree->Fill();

}

void Skim_2023::WriteHist() {
    //GetOutfile()->mkdir("Events");
    GetOutfile()->cd();
    newtree->Write();
    //GetOutfile()->cd();   
}


