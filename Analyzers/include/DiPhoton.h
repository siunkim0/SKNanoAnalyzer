#ifndef DiPhoton_h
#define DiPhoton_h

#include "AnalyzerCore.h"
#include "SystematicHelper.h"

class DiPhoton : public AnalyzerCore {
public:
    DiPhoton();
    virtual ~DiPhoton();

    void initializeAnalyzer();
    void executeEvent();
    void executeEventFromParameter();

private:
    // HLT diphoton triggers (HLT_DoublePhoton70 OR HLT_DoublePhoton85)
    RVec<TString> DiPhoTriggers;

    // Photon selection thresholds
    TString PhotonID;
    float   PhotonPtCut;
    float   PhotonEtaCut;

    // Systematic helper
    unique_ptr<SystematicHelper> systHelper;
};

#endif
