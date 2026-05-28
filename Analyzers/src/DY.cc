#include "DY.h"
#include <utility>

//==== Constructor and Destructor
DY::DY() {}
DY::~DY() {}

void DY::initializeAnalyzer() {

    RunSyst = HasFlag("RunSyst");
    MuonIDs.clear();
    MuonIDs.push_back(Muon::MuonID::POG_TIGHT);
    MuonIDs.push_back(Muon::MuonID::POG_PFISO_TIGHT);
    

    //==== Trigger settings
    if (DataEra == "2022") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2022EE") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2023") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2023BPix") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2024") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;  
    } else if (DataEra == "2018") {  // <--- 2018년 조건 추가!
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;     
    } else {
        cerr << "[DY::initializeAnalyzer] DataEra is not set properly: " << DataEra << endl;
        exit(EXIT_FAILURE);
    }
    
    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);
    

    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", DataStream, DataEra);
    } else {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/ExampleSystematic.yaml", MCSample, DataEra);
    }
    
    
}

void DY::executeEvent() {

    AllMuons = GetAllMuons();
    for (const auto &syst_dummy : *systHelper) {
        executeEventFromParameter();
    }
}

void DY::executeEventFromParameter() {
    // === Getting systematic  
    const TString this_syst = systHelper->getCurrentSysName();
    Event ev = GetEvent();
    float weight = 1.0;

    FillHist(this_syst + "/CutFlow", 1, 1 , 10, 0 , 10 );

    //==== Trigger
        //==== see Dataformat/src/Event.cc
    if (!ev.PassTrigger(IsoMuTriggerName)) return;
 
    //Cutflow after trigger
    FillHist(this_syst + "/CutFlow", 2, 1 , 10, 0 , 10 );

    auto muons_tight = SelectMuons(AllMuons,"POGTight", 15., 2.4);
    auto muons_iso = SelectMuons(muons_tight,"POGPfIsoTight", 15., 2.4); 

    FillHist(this_syst + "/CutFlow", 3, 1, 10, 0, 10);

    // 분석에서 사용할 메인 뮤온 벡터 정의
    RVec<Muon> selectedMuons = muons_iso; 
    if (selectedMuons.size() != 2) return; // 2개 이상일 때만 통과
    sort(selectedMuons.begin(), selectedMuons.end(), PtComparing);

    //Cutflow after delete single muon events
    FillHist(this_syst + "/CutFlow", 4, 1 , 10, 0 , 10 );

    if ((selectedMuons[0].Pt() < cuts.muon_pt_lead)) return;
    if ((selectedMuons[1].Pt() < cuts.muon_pt_sublead)) return;
    
    //Cutflow after leading and subleading pt cuts
    FillHist(this_syst + "/CutFlow", 5, 1 , 10, 0 , 10 );

    if (selectedMuons[0].Charge() * selectedMuons[1].Charge() > 0) return; //opposite sign
    
    //Cutflow after OS selection
    FillHist(this_syst + "/CutFlow", 6, 1 , 10, 0 , 10 );

    // --- [Correction 추가 영역] ---
    
    weight = 1.0;  // <--- 'float'를 지워서 재선언이 아닌 '값 대입'으로 변경!
    
    if (!IsDATA) {
        // MC
        weight *= MCweight();
        weight *= ev.GetTriggerLumi("Full");
        
        // 1. Muon RECO SF (전체 선택된 뮤온 벡터에 대해 적용)
        weight *= myCorr->GetMuonRECOSF(selectedMuons);

        // 2. Muon ID SF 
        // "NUM_TightID_DEN_TrackerMuons"는 예시입니다. 사용하시는 ID에 맞는 Key를 넣으세요.
        weight *= myCorr->GetMuonIDSF("NUM_TightID_DEN_TrackerMuons", selectedMuons);
        weight *= myCorr->GetMuonIDSF("NUM_TightPFIso_DEN_TightID", selectedMuons);
        weight *= myCorr->GetMuonTriggerSF("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight", selectedMuons);
        
        // 3. Pileup Weight
        // nTrueInt는 시뮬레이션된 실제 Pileup 수입니다.
        weight *= myCorr->GetPUWeight(ev.nTrueInt());
        
        // ----------------------------
    }
    float dilepton_mass = (selectedMuons[0] + selectedMuons[1]).M();

    FillHist(this_syst + "/DileptonMass", dilepton_mass, weight, 3000, 0., 3000.);
}

