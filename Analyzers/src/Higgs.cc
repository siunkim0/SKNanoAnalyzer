#include "Higgs.h"
#include <utility>

//==== Constructor and Destructor
Higgs::Higgs() {}
Higgs::~Higgs() {}

//==== Initialize variables
void Higgs::initializeAnalyzer() {

    RunSyst = HasFlag("RunSyst");

    // MuonIDs.clear();
    // MuonIDs.push_back(Muon::MuonID::POG_TIGHT);

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
    } else if (DataEra == "2018") {  // <--- DoubleMuon Trigger
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

void Higgs::executeEvent() {

    AllMuons = GetAllMuons();
    AllJets = GetAllJets();
    
    for (const auto &syst_dummy : *systHelper) {
        executeEventFromParameter();
    }
}

void Higgs::executeEventFromParameter() {
    float weight = 1.0; 
    Event ev = GetEvent();

    const TString this_syst = systHelper->getCurrentSysName();

    // MC인 경우만 weight 계산 및 Gen 분석
    if (!IsDATA) {
        weight *= MCweight();
        weight *= ev.GetTriggerLumi("Full");
        weight *= myCorr->GetPUWeight(ev.nTrueInt());
    
        // === Gen level 정보 확인 (MC only) ===
        auto gens = GetAllGens();
        RVec<Gen> higgs_muons;

        for (const auto& g : gens) {
            if (abs(g.PID()) == 13) {
                int mom_idx = g.MotherIndex();
                if (mom_idx >= 0 && abs(gens.at(mom_idx).PID()) == 23) {
                    int grand_mom_idx = gens.at(mom_idx).MotherIndex();
                    if (grand_mom_idx >= 0 && abs(gens.at(grand_mom_idx).PID()) == 25) {
                        higgs_muons.push_back(g);
                    }
                }
            }
        }

        // Gen 레벨 히스토그램 (selection 전)
        if (higgs_muons.size() == 4) {
            float gen_H_mass = (higgs_muons[0] + higgs_muons[1] + 
                                higgs_muons[2] + higgs_muons[3]).M();
            FillHist(this_syst + "/Gen_H_Mass", gen_H_mass, weight, 125, 0, 250);
            FillHist(this_syst + "/Gen_H4Mu_Event", 1, weight, 2, 0, 2);
        }
    }  // <-- if (!IsDATA) 블록 끝

    // 여기부터 Data/MC 공통 코드
    //nocuts
    FillHist(this_syst + "/CutFlow", 1, weight, 10, 0, 10 );

    //==== muons and jets for this event are one the device memory
		//==== select muons  for current ID POG Tight ID
    RVec<Muon> muons = AllMuons;
    RVec<Jet> jets = AllJets;

    //==== Trigger
		//==== see Dataformat/src/Event.cc
    if (!ev.PassTrigger(IsoMuTriggerName)) return;

    //Cutflow after trigger
    FillHist(this_syst + "/CutFlow", 2, weight, 10, 0, 10 );

    auto muons_tight = SelectMuons(AllMuons,"POGTight", 4., 2.4);
    auto muons_iso = SelectMuons(muons_tight,"POGPfIsoTight", 4., 2.4); 
    FillHist(this_syst + "/CutFlow", 3, weight, 10, 0, 10);
    RVec<Muon> selectedMuons = muons_iso;

    if ( selectedMuons.size() != 4 ) return;
    sort(selectedMuons.begin(), selectedMuons.end(), PtComparing);
    
    //Cutflow after delete single muon events
    FillHist(this_syst + "/CutFlow", 4, weight , 10, 0 , 10 );
    // mass cut 주기 전 후 MET, Jet & Z delta R, delta Phi, pt
    if ((selectedMuons[0].Pt() < cuts.muon_pt_lead)) return;
    if ((selectedMuons[1].Pt() < cuts.muon_pt_sublead)) return;
    if ((selectedMuons[2].Pt() < cuts.muon_pt_sub2lead)) return;
    if ((selectedMuons[3].Pt() < cuts.muon_pt_sub3lead)) return;
    
    //Cutflow after leading and subleading pt cuts
    FillHist(this_syst + "/CutFlow", 5, weight , 10, 0 , 10 );

    if (selectedMuons[0].Charge() + selectedMuons[1].Charge() +
        selectedMuons[2].Charge() + selectedMuons[3].Charge() != 0) return; //total charge zero
    
    //Cutflow after OS selection
    FillHist(this_syst + "/CutFlow", 6, weight , 10, 0 , 10 );

    double min_diff_Z1 = 999999.;
    int i1 = -1, j1 = -1;

    // 뮤온 짝 짓기
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            // 1. 반대 전하(OS)여야 함
            if (selectedMuons[i].Charge() + selectedMuons[j].Charge() != 0) continue;

            // 2. 두 뮤온의 불변 질량 계산
            float m_ll = (selectedMuons[i] + selectedMuons[j]).M();
            float diff = fabs(m_ll - 91.18);

            // 3. mZ(91.18)에 가장 가까운 쌍을 업데이트
            if (diff < min_diff_Z1) {
                min_diff_Z1 = diff;
                i1 = i; j1 = j;
            }
        }
    }

    // 만약 적절한 OS 쌍을 못 찾았다면 (매우 드물음) return
    if (i1 == -1) return;

    int i2 = -1, j2 = -1;
    for (int k = 0; k < 4; ++k) {
        if (k != i1 && k != j1) {
            if (i2 == -1) i2 = k;
            else j2 = k;
        }
    }

    // Z1, Z2 사차원 벡터 구성
    TLorentzVector p4_Z1 = selectedMuons[i1] + selectedMuons[j1];
    TLorentzVector p4_Z2 = selectedMuons[i2] + selectedMuons[j2];

    // 최종 4-muon(힉스 후보) 질량
    float mass_4mu = (p4_Z1 + p4_Z2).M();
    float Hpt = (p4_Z1 + p4_Z2).Pt();
    float Z1_pt = p4_Z1.Pt();
    float Z2_pt = p4_Z2.Pt();
    float Z1_mass = p4_Z1.M();
    float Z2_mass = p4_Z2.M();
    double dR_Z1Z2 = p4_Z1.DeltaR(p4_Z2);
    double dR_muons_Z1 = selectedMuons[i1].DeltaR(selectedMuons[j1]);
    double dR_muons_Z2 = selectedMuons[i2].DeltaR(selectedMuons[j2]);
    float Z1Muon1_pt = selectedMuons[i1].Pt();
    float Z1Muon2_pt = selectedMuons[j1].Pt();
    float Z2Muon1_pt = selectedMuons[i2].Pt();
    float Z2Muon2_pt = selectedMuons[j2].Pt();
    float Z1Muon1_eta = selectedMuons[i1].Eta();
    float Z1Muon2_eta = selectedMuons[j1].Eta();
    float Z2Muon1_eta = selectedMuons[i2].Eta();
    float Z2Muon2_eta = selectedMuons[j2].Eta();

    if (!IsDATA) {
        // 1. Muon ID SF, RECO SF
        // "NUM_TightID_DEN_TrackerMuons"는 예시입니다. 사용하시는 ID에 맞는 Key를 넣으세요.
        for (const auto& mu : selectedMuons) {
            weight *= myCorr->GetMuonRECOSF(mu); // RECO 효율 보정
            if (mu.Pt() > 15.0) {
                // for 2023, 2022EE
                //weight *= myCorr->GetMuonIDSF("NUM_MediumID_DEN_TrackerMuons", mu);
                //weight *= myCorr->GetMuonIDSF("NUM_LoosePFIso_DEN_MediumID", mu);

                // for 2023, 2022EE
                weight *= myCorr->GetMuonIDSF("NUM_TightID_DEN_TrackerMuons", mu);
                weight *= myCorr->GetMuonIDSF("NUM_TightPFIso_DEN_TightID", mu);

                // for 2017, 2018 Tight
                //weight *= myCorr->GetMuonIDSF("NUM_TightID_DEN_TrackerMuons", mu);
                //weight *= myCorr->GetMuonIDSF("NUM_TightRelIso_DEN_TightIDandIPCut", mu);

                // for 2018 Loose
                //weight *= myCorr->GetMuonIDSF("NUM_MediumID_DEN_TrackerMuons", mu);
                //weight *= myCorr->GetMuonIDSF("NUM_LooseRelIso_DEN_MediumID", mu);

                // weight *= myCorr->GetMuonIDSF("NUM_TightID_DEN_TrackerMuons", mu);
                // weight *= myCorr->GetMuonIDSF("NUM_TightPFIso_DEN_TightID", mu);
            } else {
                weight *= 1.0;
            }
        }
    }

    // Gen-level mother/grandmother PID study for Z1 in 50-80 GeV region (MC only)
    // Gen-level ancestor study (MC only)
    /*
    if (!IsDATA) {
        if (Z1_mass > 50 && Z1_mass < 80) {
            auto gens = GetAllGens();
        
            for (int muon_num = 0; muon_num < 2; ++muon_num) {
                int idx = (muon_num == 0) ? i1 : j1;
                const auto& mu = selectedMuons[idx];
                int gen_idx = mu.GenPartIdx();
            
                if (gen_idx >= 0) {
                    // 1대: Mother
                    int mom_idx = gens.at(gen_idx).MotherIndex();
                    if (mom_idx >= 0) {
                        FillHist(this_syst + "/MotherPID_Z1Region", abs(gens.at(mom_idx).PID()), weight, 30, 0, 30);
                    
                        // 2대: Grand-mother
                        int g_mom_idx = gens.at(mom_idx).MotherIndex();
                        if (g_mom_idx >= 0) {
                            FillHist(this_syst + "/GrandMotherPID_Z1Region", abs(gens.at(g_mom_idx).PID()), weight, 30, 0, 30);
                        
                            // 3대: Great-Grand-mother (할머니의 엄마)
                            int gg_mom_idx = gens.at(g_mom_idx).MotherIndex();
                            if (gg_mom_idx >= 0) {
                                int gg_mom_pid = abs(gens.at(gg_mom_idx).PID());
                                FillHist(this_syst + "/GreatGrandMotherPID_Z1Region", gg_mom_pid, weight, 30, 0, 30);
                            
                                // 4대: Great-Great-Grand-mother (증조할머니의 엄마)
                                int ggg_mom_idx = gens.at(gg_mom_idx).MotherIndex();
                                if (ggg_mom_idx >= 0) {
                                    int ggg_mom_pid = abs(gens.at(ggg_mom_idx).PID());
                                    FillHist(this_syst + "/GreatGreatGrandMotherPID_Z1Region", ggg_mom_pid, weight, 30, 0, 30);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    */
/*
    // 1. 논리 연산자 수정 (&& 사용)
    if (Z1_mass > 40 && Z1_mass < 80) {
    
        // Muon 1 Pt & Eta
        FillHist(this_syst + "/Z1Muon1_Pt_40-80",  Z1Muon1_pt,  weight, 100, 0., 200.);
        FillHist(this_syst + "/Z1Muon1_Eta_40-80", Z1Muon1_eta, weight, 100, -2.5, 2.5); // eta 범위는 보통 -2.5 ~ 2.5
    
        // Muon 2 Pt & Eta
        FillHist(this_syst + "/Z1Muon2_Pt_40-80",  Z1Muon2_pt,  weight, 100, 0., 200.);
        FillHist(this_syst + "/Z1Muon2_Eta_40-80", Z1Muon2_eta, weight, 100, -2.5, 2.5);
    }

    if (Z2_mass > 10 && Z2_mass < 30) {
        // Muon 1 Pt & Eta
        FillHist(this_syst + "/Z2Muon1_Pt_10-30",  Z2Muon1_pt,  weight, 100, 0., 200.);
        FillHist(this_syst + "/Z2Muon1_Eta_10-30", Z2Muon1_eta, weight, 100, -2.5, 2.5); // eta 범위는 보통 -2.5 ~ 2.5
    
        // Muon 2 Pt & Eta
        FillHist(this_syst + "/Z2Muon2_Pt_10-30",  Z2Muon2_pt,  weight, 100, 0., 200.);
        FillHist(this_syst + "/Z2Muon2_Eta_10-30", Z2Muon2_eta, weight, 100, -2.5, 2.5);
    }
*/
    // 히스토그램 채우기
    FillHist(this_syst + "/H_Mass", mass_4mu, weight, 125, 0., 250.);
    FillHist(this_syst + "/H_Pt", Hpt, weight, 50, 0., 200.);
    FillHist(this_syst + "/Z1_Mass", p4_Z1.M(), weight, 100, 0., 200.);
    FillHist(this_syst + "/Z2_Mass", p4_Z2.M(), weight, 100, 0., 200.);
    FillHist(this_syst + "/Z1_Pt", Z1_pt, weight, 50, 0., 200.);
    FillHist(this_syst + "/Z2_Pt", Z2_pt, weight, 50, 0., 200.);
    FillHist(this_syst + "/DeltaR_Z1Z2", dR_Z1Z2, weight, 50, 0., 5.);
    FillHist(this_syst + "/Z1_Muon_dR", dR_muons_Z1, weight, 50, 0., 5.0);
    FillHist(this_syst + "/Z2_Muon_dR", dR_muons_Z2, weight, 50, 0., 5.0);

    FillHist(this_syst + "/Z1Muon1_Pt", Z1Muon1_pt, weight, 100, 0., 200.);
    FillHist(this_syst + "/Z1Muon2_Pt", Z1Muon2_pt, weight, 100, 0., 200.);
    FillHist(this_syst + "/Z2Muon1_Pt", Z2Muon1_pt, weight, 100, 0., 200.);
    FillHist(this_syst + "/Z2Muon2_Pt", Z2Muon2_pt, weight, 100, 0., 200.);

    if (Z1_mass < 40.0 || Z1_mass > 120.0) return;
    if (Z2_mass < 12.0 || Z2_mass > 100.0) return;

    FillHist(this_syst + "/H_Mass_cut", mass_4mu, weight, 125, 0., 250.);
    FillHist(this_syst + "/Z1_Muon_dR_cut", dR_muons_Z1, weight, 50, 0., 5.0);
    FillHist(this_syst + "/Z2_Muon_dR_cut", dR_muons_Z2, weight, 50, 0., 5.0);

    if (mass_4mu < 70.0) return;
    FillHist(this_syst + "/H_Mass_final", mass_4mu, weight, 125, 0., 250.);
} 
    




