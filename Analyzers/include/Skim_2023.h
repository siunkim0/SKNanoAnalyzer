#ifndef Skim_2023_h
#define Skim_2023_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"
#include "Muon.h"
#include "Electron.h"

class Skim_2023 : public AnalyzerCore {
public:
    Skim_2023();
    ~Skim_2023();

    void initializeAnalyzer();
    void executeEvent();
    void executeEventFromParameter();
    TTree *newtree;
    void WriteHist();

    bool RunSyst;
    bool RunWRCut;
    unique_ptr<SystematicHelper> systHelper;
    
    
    struct flags{
        bool RunFake;
        bool RunCF ;
        bool RunSyst ;
        bool PromptLeptononly ;
        bool ApplyDYPtReweight ;
        bool ApplyDYReshape ;
        bool RunXsecSyst ;
        bool Signal ;
        bool CalculateAverageKFactor ;
        bool SignalElectronOnly ;
        bool SignalMuonOnly ;
        bool UseJetPtRwg ;
        bool UseDYCR1Reshape ;
    } flag;

    struct parameters{
        std::map<std::string ,float> Zpt;

        parameters(){
            Zpt["ZPTReweight"] = 1.0;
            Zpt["ZPTReweight_Up"] = 1.0;
            Zpt["ZPTReweight_Down"] = 1.0;
            Zpt["ZPTReweight_QCDScaleUp"] = 1.0;
            Zpt["ZPTReweight_QCDScaleDown"] = 1.0;
            Zpt["ZPTReweight_QCDPDFErrorUp"] = 1.0;
            Zpt["ZPTReweight_QCDPDFErrorDown"] = 1.0;
            Zpt["ZPTReweight_QCDPDFAlphaSUp"] = 1.0;
            Zpt["ZPTReweight_QCDPDFAlphaSDown"] = 1.0;
            Zpt["ZPtEWCorr"] = 1.0;
            Zpt["ZPtEWCorr_E1Up"] = 1.0;
            Zpt["ZPtEWCorr_E1Down"] = 1.0;
            Zpt["ZPtEWCorr_E2Up"] = 1.0;
            Zpt["ZPtEWCorr_E2Down"] = 1.0;
            Zpt["ZPtEWCorr_E3Up"] = 1.0;
            Zpt["ZPtEWCorr_E3Down"] = 1.0;
        }
    } param;
    
    struct Electrons {
        RVec<Electron> AllElectrons;

        std::string TriggerNameForSF_Electron ="WREGammaTrigger";
        RVec<Electron::ElectronID> Electron_Tight_ID = {Electron::ElectronID::POG_HEEP};
        RVec<Electron::ElectronID> Electron_Loose_ID = {Electron::ElectronID::POG_LOOSE};
        //set to custom loose id
        float Electron_MinPt = 53.;

        RVec<std::string> Ele_Trigger;
        float Ele_Trigger_Safe_Pt_Cut;
        
        std::string Electron_IS_SF_Key = "HEEP";
        std::string Electron_Trigger_SF_Key = "HEEP";
        std::string Electron_FR_ID = "HNWR";
        std::string Electron_FR_Key = "AwayJetPt40";
        std::string Electron_CF_ID = "HNWRTight";
        std::string Electron_CF_Key = "ZToLL";
        bool Electron_UseMini = false;
        bool Electron_UsePtCone = false;
        bool isPassCustomLooseID(const Electron& el) const;
    }el_set;

    struct Muons {
        RVec<Muon> AllMuons;

        std::string TriggerNameForSF_Muon ="Mu50";
        RVec<Muon::MuonID> Muon_Tight_ID = {Muon::MuonID::POG_GLOBAL_HIGH_PT};
        RVec<Muon::MuonID> Muon_Loose_ID = {Muon::MuonID::POG_GLOBAL_HIGH_PT};
        float Muon_MinPt = 53.;

        RVec<std::string> Muon_Trigger;
        float Muon_Trigger_Safe_Pt_Cut;
        
        float Muon_TkrelISO ;
        float Muon_Iso_Cut = 0.1;
        std::string Muon_IS_SF_Key = "HighPtMuonIDSF";
        std::string Muon_Trigger_SF_Key = "POGHighPtLooseTrkIso";
        std::string Muon_FR_ID = "HNWR";
        std::string Muon_FR_Key = "AwayJetPt40";
        std::string Muon_CF_ID = "HNWRTight";
        std::string Muon_CF_Key = "ZToLL";
        bool Muon_UseMini = false;
        bool Muon_UsePtCone = false;
        bool Muon_UseTuneP = true;
    }mu_set;

    struct Jets{
        RVec<Jet> AllJets;
        //TString Jet_ID ="TIGHTLEPVETO";
        RVec<Jet::JetID> Jet_ID = {Jet::JetID::TIGHTLEPVETO};
        float Jet_MinPt = 40.;
        float Jet_MaxEta = 2.4;
        RVec<Jet::JetID> JetIds;
        RVec<Jet>cleanedjet_with_tight_leptons;
        RVec<Jet>cleanedjet_with_loose_leptons;
    }jet_set;

    struct FatJets{
        RVec<FatJet> AllFatJets;
        TString FatJet_ID ="Tight";
        float Fatjet_LSF = 0.75;
        float FatJet_MinPt = 200.;
        float FatJet_MaxEta = 2.4;
        float FatJet_SDM = 40;
    }fatjet_set;

    struct Gens{
        RVec<Gen> gens;
    }gen_set;

    RVec<Jet> Clean_jet_with_tight_leptons(const RVec<Jet> & jets, const RVec<Lepton *> & tight_leps) ;
    float dR_Separation = 0.4;
    float FatJet_dR_Separation = 0.8;




};

#endif 
