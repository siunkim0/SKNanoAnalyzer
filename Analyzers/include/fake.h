#ifndef fake_h
#define fake_h

#include "AnalyzerCore.h"

// Muon fake rate measurement with prescaled single muon triggers
// (HLT_Mu8_TrkIsoVVL / HLT_Mu17_TrkIsoVVL), Run 2 2016.
// Userflags:
//   MeasFakeMu8  - run only the Mu8 path
//   MeasFakeMu17 - run only the Mu17 path
//   (no flag)    - run both paths
//   RunSyst      - fill away-jet pT variations (30 / 60 GeV)
class fake : public AnalyzerCore {
public:
    fake();
    ~fake();

    void initializeAnalyzer();
    void executeEvent();

    // public: rootcling generates a dictionary for vector<fake::TriggerPath>
    struct TriggerPath {
        TString name;      // histogram prefix: Mu8, Mu17
        TString trigger;   // HLT path name
        float ptCut;       // offline muon pT cut
        float ptCorrCut;   // cone-corrected pT region cut
    };

private:
    // Userflags
    bool MeasFakeMu8, MeasFakeMu17, RunSyst;

    vector<TriggerPath> paths;
    RVec<TString> systs;   // Central (+ AwayJetPt30, AwayJetPt60 if RunSyst)
    RVec<float> ptCorrBins, absEtaBins;
    TString electronLooseID;

    // WP cuts from the analysis note (Run 2 loose / tight)
    bool PassMuonWP(const Muon &mu, const TString &wp) const;

    void measureFakeRate(const TriggerPath &path,
                         const RVec<Muon> &looseMuons,
                         const RVec<Electron> &looseElectrons,
                         const RVec<Jet> &jets,
                         const Particle &METv,
                         float weight,
                         const RVec<Gen> &gens);
    void fillZEnriched(const TriggerPath &path,
                       const Event &ev,
                       const RVec<Muon> &looseMuons,
                       const RVec<Muon> &tightMuons,
                       const RVec<Electron> &looseElectrons,
                       const RVec<Jet> &jets,
                       float weight,
                       const RVec<Gen> &gens);
    void fillMuonKinematics(const TString &prefix,
                            const Muon &mu, float ptCorr,
                            const Particle &METv, float MT,
                            int nJets, float weight);
    TString LeptonTypeToString(int leptonType) const;
};

#endif
