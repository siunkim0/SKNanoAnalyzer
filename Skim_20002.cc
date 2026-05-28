#include "Skim_20002.h"    

Skim_20002::Skim_20002() {
    newtree = NULL; 
}

Skim_20002::~Skim_20002() {}

void Skim_20002::initializeAnalyzer() {
    GetOutfile()->cd();
    newtree = fChain->CloneTree(0);

    el_set.AllElectrons.clear();
    mu_set.AllMuons.clear();
    jet_set.AllJets.clear();
    fatjet_set.AllFatJets.clear();
    gen_set.gens.clear();

    mu_set.Muon_Trigger.clear();
    mu_set.Muon_Trigger_Safe_Pt_Cut = 0.;
    el_set.Ele_Trigger.clear();
    el_set.Ele_Trigger_Safe_Pt_Cut = 0.;
    
    if (DataEra == "2022")
    {
        mu_set.Muon_Trigger = {"HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100"}; 
        mu_set.Muon_Trigger_Safe_Pt_Cut = 52.;
        el_set.Ele_Trigger = {"HLT_Ele32_WPTight_Gsf","HLT_Photon200","HLT_Ele115_CaloIdVT_GsfTrkIdT"};
        el_set.Ele_Trigger_Safe_Pt_Cut = 35.;  
    }
    if (DataEra == "2022EE")
    {
        mu_set.Muon_Trigger = {"HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100"}; 
        mu_set.Muon_Trigger_Safe_Pt_Cut = 52.;
        el_set.Ele_Trigger = {"HLT_Ele32_WPTight_Gsf","HLT_Photon200","HLT_Ele115_CaloIdVT_GsfTrkIdT"};
        el_set.Ele_Trigger_Safe_Pt_Cut = 35.; 
    }
    if (DataEra == "2023")
    {
        mu_set.Muon_Trigger = {"HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100"}; 
        mu_set.Muon_Trigger_Safe_Pt_Cut = 52.;
        el_set.Ele_Trigger = {"HLT_Ele32_WPTight_Gsf","HLT_Photon200","HLT_Ele115_CaloIdVT_GsfTrkIdT"};
        el_set.Ele_Trigger_Safe_Pt_Cut = 35.; 
    }
    if (DataEra == "2023BPix")
    {
        mu_set.Muon_Trigger = {"HLT_Mu50", "HLT_CascadeMu100", "HLT_HighPtTkMu100"}; 
        mu_set.Muon_Trigger_Safe_Pt_Cut = 52.;
        el_set.Ele_Trigger = {"HLT_Ele32_WPTight_Gsf","HLT_Photon200","HLT_Ele115_CaloIdVT_GsfTrkIdT"};
        el_set.Ele_Trigger_Safe_Pt_Cut = 35.; 
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

void Skim_20002::executeEvent() {

    el_set.AllElectrons =  GetAllElectrons();
    mu_set.AllMuons = GetAllMuons();
    jet_set.AllJets = GetAllJets();
    fatjet_set.AllFatJets = GetAllFatJets();
    gen_set.gens = GetAllGens();
    // no tune p muon 
    

    // DY pt reweight
    if(MCSample.Contains("DY")){
    // wtf ?
    }
    for (const auto &syst_dummy : *systHelper) {
        executeEventFromParameter();
    }
    // Nvtx ( PU )


    //Parameter -no function of parameter -> use struct 
    

    //Run syst  # 457




}

void Skim_20002::executeEventFromParameter() {
    const TString this_syst = systHelper->getCurrentSysName();

    Event ev = GetEvent();

    RVec<Electron> electrons = el_set.AllElectrons;
    RVec<Muon> muons = mu_set.AllMuons;
    RVec<Jet> jets = jet_set.AllJets;
    RVec<FatJet> fatjets = fatjet_set.AllFatJets;
    if (!PassNoiseFilter(jets,ev,Event::MET_Type::PUPPI)) return;


    bool pass_trig_muon = ev.PassTrigger(mu_set.Muon_Trigger);
    bool pass_trig_elec = ev.PassTrigger(el_set.Ele_Trigger);

    /////////Trugger pass ///////////
    if ((!pass_trig_muon) and (!pass_trig_elec)) return;

    RVec<Electron> my_electrons = SelectElectrons(electrons, "NOCUT" , el_set.Electron_MinPt, 2.4); //ID in 429
    RVec<Muon> my_muons = SelectMuons(muons, "NOCUT" , mu_set.Muon_MinPt, 2.4); //ID in 446

    RVec<Electron *> Loose_electrons , Tight_electrons;
    RVec<Muon *> Loose_muons , Tight_muons;
    RVec<Lepton *> Tight_leps_el , Tight_leps_mu , Tight_leps;
    RVec<Lepton *> Loose_leps_el , Loose_leps_mu , Loose_leps;

    for (unsigned int i=0 ; i< my_electrons.size(); i ++) {
        Electron & el = my_electrons.at(i);
        if (el.PassID(el_set.Electron_Tight_ID[0])) {
            Tight_electrons.push_back(&el);
            Tight_leps_el.push_back( &el);
            Tight_leps.push_back(&el);
        }
        // Loose ID
        if (el_set.isPassCustomLooseID(el)){
            Loose_electrons.push_back(&el);
            Loose_leps_el.push_back(&el);
            Loose_leps.push_back(&el);
        }
    }
    
    for (unsigned int i=0 ; i< my_muons.size(); i ++) {
        Muon & mu = my_muons.at(i);
        float tkRelIso = mu.TkRelIso();
        if ((mu.PassID(mu_set.Muon_Tight_ID[0]))&&( tkRelIso < 0.1) ){ //global high pt id 
            Tight_muons.push_back(&mu);
            Tight_leps_mu.push_back( &mu);
            Tight_leps.push_back(&mu);
        }
        if (mu.PassID(mu_set.Muon_Loose_ID[0])) {
            Loose_muons.push_back(&mu);
            Loose_leps_mu.push_back(&mu);
            Loose_leps.push_back(&mu);
        }
    }
    int n_Loose_leptons  = Loose_electrons.size() + Loose_muons.size();
    int n_Tight_leptons  = Tight_electrons.size() + Tight_muons.size();
    jet_set.cleanedjet_with_tight_leptons = Clean_jet_with_tight_leptons(jet_set.AllJets, Tight_leps);
    RVec<Jet> selected_jets = SelectJets(jet_set.cleanedjet_with_tight_leptons, jet_set.Jet_ID[0] , jet_set.Jet_MinPt, jet_set.Jet_MaxEta);
    
    if (Tight_leps.size() < 1) return;
    newtree->Fill();

}

void Skim_20002::WriteHist() {
    //GetOutfile()->mkdir("Events");
    GetOutfile()->cd();
    newtree->Write();
    //GetOutfile()->cd();   
}

bool Skim_20002::Electrons::isPassCustomLooseID(const Electron& el) const {
    if (!(el.hoe() < 0.5)) return false;
    
    if (fabs(el.scEta()) <= 1.479){
        if (!(el.sieie() < 0.0112)) return false;
        if (!(fabs(el.deltaEtaInSC()) < 0.00377)) return false;
        if (!(fabs(el.deltaPhiInSeed()) < 0.0884)) return false;
        if (!(fabs(el.eInvMinusPInv()) < 0.193)) return false;
        if (!(el.LostHits() <= 1)) return false;
        if (!(el.ConvVeto())) return false;
        return true;
    }
    else {
        if (!(el.sieie() < 0.0425)) return false;
        if (!(fabs(el.deltaEtaInSC()) < 0.00674)) return false;
        if (!(fabs(el.deltaPhiInSeed()) < 0.169)) return false;
        if (!(fabs(el.eInvMinusPInv()) < 0.111)) return false;
        if (!(el.LostHits() <= 1)) return false;
        if (!(el.ConvVeto())) return false;
        return true;
    }
    return true ;
}

RVec<Jet> Skim_20002::Clean_jet_with_tight_leptons(const RVec<Jet> & jets, const RVec<Lepton *> & tight_leps) {
    RVec<Jet> cleanedjets;
    for (unsigned int i=0 ; i< jets.size(); i ++) {
        Jet jet = jets.at(i);
        bool isDRtoTightLepton(false);
        for (unsigned int j=0 ; j< tight_leps.size(); j ++) {
            Lepton * lep = tight_leps.at(j);
            if ( jet.DeltaR(*lep) < dR_Separation ) {
                isDRtoTightLepton = true;
                break;
            }
        }
        if (!isDRtoTightLepton) {
            cleanedjets.push_back(jet);
        }
    }
    return cleanedjets;
}


