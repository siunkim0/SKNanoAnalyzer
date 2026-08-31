#ifndef Wtag_h
#define Wtag_h

#include "AnalyzerCore.h"

//==============================================================
// Boosted W kinematics → boosted-W tagger study (MC only)
//
//   Hadronic W (W→qq') at gen level: find last-copy W with exactly
//   two quark daughters (|PID|≤6). Per W record:
//     - Stage 1 (gen truth): W pT, ΔR(q,q') → when do the two quarks
//       fall within one AK4 (0.4) / AK8 (0.8) cone.
//     - Stage 2 (detector): match the two quarks to reco AK4 jets.
//       merged = both quarks → same jet; resolved = two distinct jets;
//       overlap = ΔR(jet,jet) < 0.8. AK8 captured = one FatJet within
//       0.8 of both quarks. → fraction-vs-WpT curves.
//     - Stage 3 (jet kinematics): track multiplicity (total/charged/
//       neutral), mass, nSVs vs jet pT for
//         mergedW jet  vs  single-quark W jet  vs  light/gluon control.
//       (charged/neutral branches valid on Run3 2022+ NanoAODv12+)
//
//   v2: per-AK8-jet flat TTree ("jets") for ML tagger training —
//       label (1=W, 2=W+b/top, 0=unmatched), kinematics, substructure
//       (τ1–4, N2/N3, LSF3), multiplicities & energy fractions, softdrop
//       subjets, in-cone SVs, stored ParticleNet benchmarks, gen truth.
//==============================================================
class Wtag : public AnalyzerCore {
public:
    Wtag();
    ~Wtag();

    void initializeAnalyzer();
    void executeEvent();

    // Analysis cuts
    struct AnalysisCuts {
        float ak4_pt_min  =  20.0;
        float ak4_eta_max =   2.5;
        float ak8_pt_min  = 170.0;
        float ak8_eta_max =   2.4;
        float ak4_dr      =   0.4;   // quark→AK4 match / merge cone
        float ak8_dr      =   0.8;   // quark→AK8 capture cone
    } cuts;

    // W pT axis for fraction curves (150 bins × 10 GeV)
    static constexpr int   NWPT = 150;
    static constexpr float WPT_MAX = 1500.0;
};

#endif
