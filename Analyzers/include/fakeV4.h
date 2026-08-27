#ifndef fakeV4_h
#define fakeV4_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"
#include "MyCorrection.h"

class fakeV4: public AnalyzerCore {
public:
    fakeV4();
    virtual ~fakeV4();

    void initializeAnalyzer();
    void executeEvent();

    enum class Channel {
        NONE,
        INCLUSIVE,
        QCDENRICHED,
        WENRICHED,
        ZENRICHED
    };

    inline TString channelToString(Channel ch) {
        if (ch == Channel::INCLUSIVE) return "Inclusive";
        if (ch == Channel::QCDENRICHED) return "QCDEnriched";
        if (ch == Channel::WENRICHED) return "WEnriched";
        if (ch == Channel::ZENRICHED) return "ZEnriched";
        return "NONE";
    }

    enum class LeptonType {
        NONE,
        MUON,
        ELECTRON
    };

    // TriggerInfo struct: stores trigger info including event-level fired status
    struct TriggerInfo {
        TString name;           // HLT path: "HLT_Mu8_TrkIsoVVL"
        TString prefix;         // Histogram prefix: "Mu8"
        float trigSafePtCut;    // pT threshold (10, 15, 20, 25 GeV)
        bool fired;             // Event-level: did this trigger fire?
    };

    struct RecoObjects {
        RVec<Muon> vetoMuons;                // 10 GeV, loose ID (for event veto)
        RVec<Muon> looseMuons;               // 10 GeV, loose ID (for measurement)
        RVec<Muon> tightMuons;               // 10 GeV, tight ID (for measurement)
        RVec<Electron> vetoElectrons;        // 10 GeV, loose ID (for event veto)
        RVec<Electron> looseElectrons;       // 15 GeV, loose ID (for measurement)
        RVec<Electron> tightElectrons;       // 15 GeV, tight ID (for measurement)
        RVec<int> looseMuonJetFlavours;      // mother jet hadron flavour for loose muons
        RVec<int> tightMuonJetFlavours;      // mother jet hadron flavour for tight muons
        RVec<int> looseElectronJetFlavours;  // mother jet hadron flavour for loose electrons
        RVec<int> tightElectronJetFlavours;  // mother jet hadron flavour for tight electrons
        RVec<Jet> tightJets;
        RVec<Jet> tightJets_noPUID;  // Jets before PUID for SF calculation (Run 2)
        RVec<Jet> bjets;
        RVec<GenJet> genJets;
        Particle METv;
    };

    struct WeightInfo {
        float genWeight;
        float prefireWeight;
        float pileupWeight;
        float topPtWeight;
        float muonRecoSF;
        float eleRecoSF;
        float btagSF;
        float pileupIDSF;
    };

private:
    // Configuration flags
    bool MeasFakeMu, MeasFakeEl;
    bool RunSyst;
    bool RunNoHEMVeto;

    // Analysis configuration
    LeptonType leptonType;

    // Triggers for the chosen lepton type
    RVec<TriggerInfo> triggers;

    // Loosest pT cut across all triggers (for initial event selection)
    float lowestPtCut;

    // Binning
    RVec<float> ptcorr_bins;
    RVec<float> abseta_bins;

    // IDs
    IDContainer *MuonIDs, *ElectronIDs;

    // SystematicHelper
    std::unique_ptr<SystematicHelper> systHelper;

    // Core analysis methods
    Channel selectEvent(Event& ev, const RecoObjects& recoObjects, const TString& ID, const TString& syst);
    RecoObjects defineObjects(Event& ev, const RVec<Muon>& rawMuons,
                             const RVec<Electron>& rawElectrons,
                             const RVec<Jet>& rawJets,
                             const RVec<GenJet>& genJets,
                             const TString& ID,
                             const TString& syst = "Central");
    WeightInfo getWeights(const Channel& channel,
                          const TString& ID,
                          Event& event,
                          const RecoObjects& recoObjects,
                          const RVec<Gen>& genParts,
                          const TString& syst = "Central");

    void fillObjects(const Channel& channel,
                     const TString& ID,
                     const RecoObjects& recoObjects,
                     const WeightInfo& weights,
                     const TString& syst = "Central");

    // Helper methods
    TString getBinPrefix(const double ptcorr, const double abseta);
    float getJetPtCut(const TString& selection);
    template<typename T>
    int getMotherJetFlavour(const T& lep, const RVec<Jet>& allJets);
    TString getFlavourSubdir(int flavour);

    // Cutflow functionality
    enum class CutStage {
        Initial = 0,
        NoiseFilter = 1,
        VetoMap = 2,
        AnyTrigger = 3,
        LeptonSelection = 4,
        JetRequirements = 5,
        AwayJetRequirements = 6,
        ZMassWindow = 7,
        Final = 8
    };

    void fillCutflow(CutStage stage, const Channel& channel, const TString& ID, float weight, const TString& syst);

    //==== 아래는 fake.cc 에서 그대로 옮겨온 flat ntuple 부분 (userflag MakeTree).
    //     측정 로직(defineObjects/selectEvent/fillObjects)은 공식 MeasFakeRateV4
    //     코드 그대로다 — 그래야 tree 의 population 이 공식 FR map 을 재현하는
    //     선택과 정확히 같다. 우리 fake.cc 는 jet 개수에서 공식과 어긋난다.
    bool MakeTree;
    struct { float awayjet_dr = 0.7; float tight_miniiso_max = 0.1; } treecuts;
    RVec<Muon> looseMuons;
    RVec<Electron> looseElectrons;
    RVec<Jet> rawJets, jets;
    RVec<Gen> gens;
    RVec<GenJet> genJets;
    Particle METv;

    //==== loose WP 만 직접 구현한다: 우리 fork 의 Muon 에는
    //     Pass_HcToWALooseRun2/Run3 자체가 없고, Electron 은 있지만 Run3 값이
    //     낡았다 (SIP3D<8, MVA 0.5/-0.8/-0.5). tight 은 두 repo 가 완전히 같아서
    //     PassID("HcToWATight") 를 그대로 쓴다. 값은 공식 DataFormats 와 대조 완료.
    bool PassLooseMuon(const Muon &mu) const;
    bool PassLooseElectron(const Electron &el) const;
    RVec<Muon> PickMuons(const RVec<Muon> &muons, const TString &wp,
                         float ptmin, float fetamax) const;
    RVec<Electron> PickElectrons(const RVec<Electron> &els, const TString &wp,
                                 float ptmin, float fetamax, bool vetoHEM) const;

    TString LeptonTypeToString(int leptonType) const;
    TString FlavorTag(int partonFlavour, int hadronFlavour, bool isPileup) const;
    TString RecoJetFlavor(const Muon &mu) const;
    TString GenJetFlavor(const Muon &mu) const;
    int FlavorCode(const TString &flav) const;
    void fillTreeRow(Event &ev, float evtSF);
};

#endif
