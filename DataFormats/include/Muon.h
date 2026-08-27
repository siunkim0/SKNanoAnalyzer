#ifndef Muon_h
#define Muon_h

#include "TString.h"
#include "Lepton.h"

// Need update
// - TuneP object
// - momentum scale
// - rochester correction
// - chi2
// - truth matching info
// - Analysis dependent IDs

class Muon: public Lepton {
public:
    Muon();
    ~Muon();

    // Boolean IDs
    enum class BooleanID {NONE, LOOSE, MEDIUM, MEDIUMPROMPT, TIGHT, SOFT, SOFTMVA, TRIGGERLOOSE};

    void SetBIDBit(BooleanID id, bool idbit);
    inline bool isPOGTightId() const {return j_tightId;}
    inline bool isPOGMediumId() const {return j_mediumId;}
    inline bool isPOGMediumPromptId() const {return j_mediumPromptId;}
    inline bool isPOGLooseId() const {return j_looseId;}
    inline bool isPOGSoftId() const {return j_softId;}
    inline bool isPOGSoftMvaId() const {return j_softMvaId;}
    inline bool isPOGTriggerIdLoose() const {return j_triggerIdLoose;}

    // Muon type methods
    void SetIsTracker(bool isTracker) { j_isTracker = isTracker; }
    void SetIsStandalone(bool isStandalone) { j_isStandalone = isStandalone; }
    void SetIsGlobal(bool isGlobal) { j_isGlobal = isGlobal; }
    inline bool isTracker() const { return j_isTracker; }
    inline bool isStandalone() const { return j_isStandalone; }
    inline bool isGlobal() const { return j_isGlobal; }

    // Unsigned char IDs
    enum class WorkingPointID {NONE, HIGHPT, MINIISO, MULTIISO, MVAMU, PFISO, PUPPIISO, TKISO};
    enum class WorkingPoint {NONE, VLOOSE, LOOSE, MEDIUM, TIGHT, VTIGHT, VVTIGHT};

    enum class MuonID
    {
        NOCUT,
        POG_TIGHT,
        POG_MEDIUM,
        POG_MEDIUM_PROMPT,
        POG_LOOSE,
        POG_SOFT,
        POG_SOFT_MVA,
        POG_TRIGGER_LOOSE,
        POG_TRACKER_HIGH_PT,
        POG_GLOBAL_HIGH_PT,
        POG_MINISO_LOOSE,
        POG_MINISO_MEDIUM,
        POG_MINISO_TIGHT,
        POG_MINISO_VTIGHT,
        POG_MULTISO_LOOSE,
        POG_MULTISO_MEDIUM,
        POG_MVA_MU_MEDIUM,
        POG_MVA_MU_TIGHT,
        POG_PFISO_VLOOSE,
        POG_PFISO_LOOSE,
        POG_PFISO_MEDIUM,
        POG_PFISO_TIGHT,
        POG_PFISO_VTIGHT,
        POG_PFISO_VVTIGHT,
        POG_PUPPIISO_LOOSE,
        POG_PUPPIISO_MEDIUM,
        POG_PUPPIISO_TIGHT,
        POG_TKISO_LOOSE,
        POG_TKISO_TIGHT
    };

    void SetWIDBit(WorkingPointID id, unsigned char value);
    inline WorkingPoint HighPtId() const {return (WorkingPoint)j_highPtId;}
    inline WorkingPoint MiniIsoId() const {return (WorkingPoint)j_miniIsoId;}
    inline WorkingPoint MultiIsoId() const {return (WorkingPoint)j_multiIsoId;}
    inline WorkingPoint MvaMuId() const {return (WorkingPoint)j_mvaMuId;}
    //inline WorkingPoint MvaLowPtId() const {return (WorkingPoint)j_mvaLowPtId;}
    inline WorkingPoint PfIsoId() const {return (WorkingPoint)j_pfIsoId;}
    inline WorkingPoint PuppiIsoId() const {return (WorkingPoint)j_puppiIsoId;}
    inline WorkingPoint TkIsoId() const {return (WorkingPoint)j_tkIsoId;}

    void SetNTrackerLayers(int n) {j_nTrackerLayers = n;}
    inline int nTrackerLayers() const {return j_nTrackerLayers;}
    void SetOriginalPt(float pt) {j_miniAODPt = pt;}
    inline float OriginalPt() const {return j_miniAODPt;}
    void SetMomentumScaleUpDown(float up, float down) {j_momentumScaleUp = up; j_momentumScaleDown = down;}
    inline float MomentumScaleUp() const {return j_momentumScaleUp;}
    inline float MomentumScaleDown() const {return j_momentumScaleDown;}
    // MVA ID scores
    enum class MVAID {NONE, SOFTMVA, MVALOWPT, MVATTH};

    void SetMVAID(MVAID id, float score);
    inline float SoftMva() const {return j_softMva;}
    inline float MvaLowPt() const {return j_mvaLowPt;}
    inline float MvaTTH() const {return j_mvaTTH;}

    void SetGenPartIdx(short genPartIdx) { j_genPartIdx = genPartIdx; }
    inline short GenPartIdx() const { return j_genPartIdx; }

    void SetGenPartFlav(unsigned char genPartFlav) { j_genPartFlav = genPartFlav; }
    inline unsigned char GenPartFlav() const { return j_genPartFlav; }

    void SetJetIdx(short jetIdx) { j_jetIdx = jetIdx; }
    inline short JetIdx() const { return j_jetIdx; }

    // NanoAOD 에서 그대로 가져오는 추가 변수들 (fake flavor BDT 입력용)
    //   svIdx 는 NanoAODv12+ 에만 있으므로 Run2 에서는 -1 로 남는다
    void SetSVIdx(short svIdx) { j_svIdx = svIdx; }
    inline short SVIdx() const { return j_svIdx; }

    void SetTrackErrors(float ptErr, float tunepRelPt) { j_ptErr = ptErr; j_tunepRelPt = tunepRelPt; }
    inline float PtErr() const { return j_ptErr; }
    inline float TunepRelPt() const { return j_tunepRelPt; }

    void SetMuonQuality(float segmentComp, int nStations) { j_segmentComp = segmentComp; j_nStations = nStations; }
    inline float SegmentComp() const { return j_segmentComp; }
    inline int NStations() const { return j_nStations; }

    void SetChargedIso(float miniChg, float pf03Chg) { j_miniPFRelIso_chg = miniChg; j_pfRelIso03_chg = pf03Chg; }
    inline float MiniPFRelIsoChg() const { return j_miniPFRelIso_chg; }
    inline float PfRelIso03Chg() const { return j_pfRelIso03_chg; }

    void SetJetVars(float jetPtRelv2, float jetRelIso, int jetNDauCharged) {
        j_jetPtRelv2 = jetPtRelv2; j_jetRelIso = jetRelIso; j_jetNDauCharged = jetNDauCharged;
    }
    inline float JetPtRelv2() const { return j_jetPtRelv2; }
    inline float JetRelIso() const { return j_jetRelIso; }
    inline int JetNDauCharged() const { return j_jetNDauCharged; }

    void SetTrackFlags(bool isPFcand, bool highPurity, int tightCharge) {
        j_isPFcand = isPFcand; j_highPurity = highPurity; j_tightCharge = tightCharge;
    }
    inline bool isPFcand() const { return j_isPFcand; }
    inline bool highPurity() const { return j_highPurity; }
    inline int TightCharge() const { return j_tightCharge; }

    // ID helper functions
    bool PassID(const MuonID ID) const;
    bool PassID(const TString ID) const;

    // Private IDs
    bool Pass_HcToWATight() const;
    bool Pass_HcToWALoose() const;

private:
    bool j_isTracker, j_isStandalone, j_isGlobal;
    bool j_looseId, j_mediumId, j_mediumPromptId, j_tightId, j_softId, j_softMvaId, j_triggerIdLoose;
    unsigned char j_highPtId, j_miniIsoId, j_multiIsoId, j_mvaMuId, j_pfIsoId, j_puppiIsoId, j_tkIsoId;
    float j_softMva, j_mvaLowPt, j_mvaTTH;
    int j_nTrackerLayers;
    float j_miniAODPt, j_momentumScaleUp, j_momentumScaleDown;
    short j_genPartIdx;
    unsigned char j_genPartFlav;
    short j_jetIdx;
    short j_svIdx;
    float j_ptErr, j_tunepRelPt;
    float j_segmentComp;
    int j_nStations;
    float j_miniPFRelIso_chg, j_pfRelIso03_chg;
    float j_jetPtRelv2, j_jetRelIso;
    int j_jetNDauCharged;
    bool j_isPFcand, j_highPurity;
    int j_tightCharge;
    ClassDef(Muon, 2);
};

#endif
