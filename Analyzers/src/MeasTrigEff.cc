#include "MeasTrigEff.h"

MeasTrigEff::MeasTrigEff() {}
MeasTrigEff::~MeasTrigEff() {}

void MeasTrigEff::initializeAnalyzer() {
    // ID settings
    MuonIDs = new IDContainer("HcToWATight", "HcToWALoose");
    
    // IsoMu24 trigger for all eras
    if (DataEra == "2016preVFP" || DataEra == "2016postVFP") {
        SglMuTrigs = {"HLT_IsoMu24", "HLT_IsoTkMu24"};
    } else if (DataEra == "2017") {
        SglMuTrigs = {"HLT_IsoMu27"};
    } else if (DataEra == "2018") {
        SglMuTrigs = {"HLT_IsoMu24"};
    } else {  // Run3
        SglMuTrigs = {"HLT_IsoMu24"};
    }
    
    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA?DataStream:MCSample ,IsDATA);
}

void MeasTrigEff::executeEvent() {
    Event ev = GetEvent();
    RVec<Jet> rawJets = GetAllJets();
    
    if (!PassNoiseFilter(rawJets, ev)) return;
    
    FillHist("CutFlow", 1, 1, 10, 0, 10);
    
    if (!ev.PassTrigger(SglMuTrigs)) return;
    FillHist("CutFlow", 2, 1, 10, 0, 10);
    
    // Object Definition
    RVec<Muon> rawMuons = GetMuons("NOCUT", 5., 2.4);
    sort(rawMuons.begin(), rawMuons.end(), PtComparing);
    
    RVec<Muon> looseMuons = SelectMuons(rawMuons, MuonIDs->GetID("loose"), 10., 2.4);
    RVec<Muon> tightMuons = SelectMuons(rawMuons, MuonIDs->GetID("tight"), 10., 2.4);
    RVec<TrigObj> trigObjs = GetAllTrigObjs();
    
    // Weights
    float weight = 1.;
    if (!IsDATA) {
        weight = MCweight()*ev.GetTriggerLumi("Full");
        weight *= GetL1PrefireWeight();
        weight *= myCorr->GetPUWeight(ev.nTrueInt());
    }
    
    measIsoMu24TrigEff(tightMuons, looseMuons, trigObjs, weight);
}

void MeasTrigEff::measIsoMu24TrigEff(RVec<Muon> &tightMuons, RVec<Muon> &looseMuons, RVec<TrigObj> &trigObjs, float &weight) {
    // Event Selection: Tag-and-Probe with Z->mumu
    if (tightMuons.size() < 1 || looseMuons.size() < 2) return;
    
    const Muon &tag = tightMuons.at(0);
    
    // Tag muon must pass trigger matching
    if (!PassIsoMuTrigger(tag, trigObjs)) return;
    
    // Tag muon kinematic cuts
    const float tag_pt_cut = (DataYear == 2017) ? 29. : 26.;
    if (tag.Pt() < tag_pt_cut) return;
    
    // Find probe muon (opposite charge, forms Z peak)
    for (const auto &probe : looseMuons) {
        if (&probe == &tag) continue;
        if (tag.Charge   () + probe.Charge() != 0) continue;
        
        double Mll = (tag + probe).M();
        if (fabs(Mll - 91.2) > 10.) continue;  // Still use wide Z window for selection
        if (tag.DeltaR(probe) < 0.3) continue;
        
        // Check if probe passes trigger
        const bool probePassTrig = PassIsoMuTrigger(probe, trigObjs);
        
        // ★ Overall (integrated, 86-96 GeV)
        if (Mll >= 86. && Mll < 96.) {
            FillHist("/TrigEff_IsoMu24_Overall_DENOM", fabs(probe.Eta()), probe.Pt(), weight, 
                     24, 0., 2.4, 100, 0., 200.);
            if (probePassTrig) {
                FillHist("/TrigEff_IsoMu24_Overall_NUM", fabs(probe.Eta()), probe.Pt(), weight, 
                         24, 0., 2.4, 100, 0., 200.);
            }
        }
        
        // ★ Mass bin별 (1 GeV bins: 86-96 GeV)
        int mass_bin = (int)Mll;  // 86, 87, 88, ..., 95
        
        if (mass_bin >= 86 && mass_bin < 96) {
            TString histname_denom = Form("/TrigEff_IsoMu24_MassBin%d_%d_DENOM", mass_bin, mass_bin+1);
            TString histname_num = Form("/TrigEff_IsoMu24_MassBin%d_%d_NUM", mass_bin, mass_bin+1);
            
            FillHist(histname_denom, fabs(probe.Eta()), probe.Pt(), weight, 
                     24, 0., 2.4, 100, 0., 200.);
            
            if (probePassTrig) {
                FillHist(histname_num, fabs(probe.Eta()), probe.Pt(), weight, 
                         24, 0., 2.4, 100, 0., 200.);
            }
        }
    }
}

bool MeasTrigEff::PassIsoMuTrigger(const Muon &mu, const RVec<TrigObj> &trigObjs) {
    const float trig_pt_cut = (DataYear == 2017) ? 27. : 24.;
    
    for (const auto &trigObj : trigObjs) {
        if (!trigObj.isMuon()) continue;
        if (trigObj.DeltaR(mu) > 0.3) continue;
        if (!trigObj.hasBit(3)) continue;
        if (trigObj.Pt() < trig_pt_cut) continue;
        return true;
    }
    return false;
}