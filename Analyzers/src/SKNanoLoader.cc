#define SKNanoLoader_cxx
#include "SKNanoLoader.h"
#include <typeinfo>
using json = nlohmann::json;

SKNanoLoader::SKNanoLoader() {
    MaxEvent = -1;
    NSkipEvent = 0;
    LogEvery = 1000;
    IsDATA = false;
    SkimmingMode = false;
    DataStream = "";
    MCSample = "";
    SetEra("2018");
    xsec = 1.;
    sumW = 1.;
    sumSign = 1.;
    Userflags.clear();
}

SKNanoLoader::~SKNanoLoader() {
    for (auto& [key, value] : TriggerMap) {
        delete value.first;
    }
    // The reader proxies (owned by the filler lambdas) must be destroyed
    // before the reader, and the reader before the chain.
    fScalarFillers.clear();
    fArrayFillers.clear();
    fReader.reset();
    if (!fChain) return;
    if (fChain->GetCurrentFile()) fChain->GetCurrentFile()->Close();
    delete fChain;
}

void SKNanoLoader::Loop() {
    long nentries = fChain->GetEntries();
    if (MaxEvent > 0) nentries = std::min(nentries, MaxEvent);
    auto startTime = std::chrono::steady_clock::now();
    cout << "[SKNanoLoader::Loop] Event Loop Started" << endl;
    
    for (long jentry = 0; jentry < nentries; jentry++) {
        if (jentry < NSkipEvent) continue;

        // Log progress for every LogEvery events
        if (jentry % LogEvery == 0) {
            auto currentTime = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsedTime = currentTime - startTime;
            double timePerEvent = elapsedTime.count() / (jentry + 1);
            double estimatedRemaining = (nentries - jentry) * timePerEvent;

            cout << "[SKNanoLoader::Loop] Processing " << jentry << " / " << nentries
                 << " | Elapsed: " << std::fixed << std::setprecision(2) << elapsedTime.count() << "s, Remaining: " << estimatedRemaining << "s" << endl;
        }

        if (fReader) {
            if (fReader->SetEntry(jentry) != TTreeReader::kEntryValid) {
                cerr << "[SKNanoLoader::Loop] Error reading event " << jentry << endl;
                exit(1);
            }
            for (const auto &fill : fScalarFillers) fill();
        } else {
            if (fChain->GetEntry(jentry) < 0) {
                cerr << "[SKNanoLoader::Loop] Error reading event " << jentry << endl;
                exit(1);
            }
            // Legacy (SkimmingMode) path: NanoAODv13 stores these branches with
            // narrower integer types than the analyzer-facing members;
            // SetBranchAddress byte-copies without conversion, so they were read
            // into raw buffers at Init and are widened here. The TTreeReader path
            // reads the on-file type directly and does not need this.
            if (Run == 3) {
                for (int i = 0; i < nMuon && i < (int)Buf_Muon_nTrackerLayers.size(); i++)
                    Muon_nTrackerLayers[i] = Buf_Muon_nTrackerLayers[i];
                for (int i = 0; i < nJet && i < (int)Buf_Jet_chMultiplicity.size(); i++) {
                    Jet_chMultiplicity[i] = Buf_Jet_chMultiplicity[i];
                    Jet_neMultiplicity[i] = Buf_Jet_neMultiplicity[i];
                }
                for (int i = 0; i < nGenVisTau && i < (int)Buf_GenVisTau_charge.size(); i++) {
                    GenVisTau_charge[i] = Buf_GenVisTau_charge[i];
                    GenVisTau_genPartIdxMother[i] = Buf_GenVisTau_genPartIdxMother[i];
                    GenVisTau_status[i] = Buf_GenVisTau_status[i];
                }
            }
        }

        // make sure Run2 and Run3 variables are in sync
        if (Run == 2) {
            nLHEPart = static_cast<Int_t>(nLHEPart_RunII);
            nGenPart = static_cast<Int_t>(nGenPart_RunII);
            nGenJet = static_cast<Int_t>(nGenJet_RunII);
            nGenJetAK8 = static_cast<Int_t>(nGenJetAK8_RunII);
            nGenIsolatedPhoton = static_cast<Int_t>(nGenIsolatedPhoton_RunII);
            nGenDressedLepton = static_cast<Int_t>(nGenDressedLepton_RunII);
            nGenVisTau = static_cast<Int_t>(nGenVisTau_RunII);
            nMuon = static_cast<Int_t>(nMuon_RunII);
            nElectron = static_cast<Int_t>(nElectron_RunII);
            nTau = static_cast<Int_t>(nTau_RunII);
            nPhoton = static_cast<Int_t>(nPhoton_RunII);
            nJet = static_cast<Int_t>(nJet_RunII);
            nFatJet = static_cast<Int_t>(nFatJet_RunII);
            nTrigObj = static_cast<Int_t>(nTrigObj_RunII);
            nSubJet = static_cast<Int_t>(nSubJet_RunII);
            nSV = static_cast<Int_t>(nSV_RunII);
        }

        if (fReader) {
            for (const auto &fill : fArrayFillers) fill();
        }

        executeEvent();
    }
    cout << "[SKNanoLoader::Loop] Event Loop Finished"<< endl;
}

void SKNanoLoader::SetMaxLeafSize(){
    //check how much time it takes to read the tree
    //and set the maximum leaf size accordingly
    auto start = std::chrono::high_resolution_clock::now();

    // Get Maximum length of arrays from the leaf metadata (TLeaf::GetMaximum)
    // of every file in the chain -- no event loop needed.
    const std::vector<TString> counterBranches = {
        "nLHEPdfWeight", "nLHEScaleWeight", "nPSWeight", "nLHEPart",
        "nGenPart", "nGenJet", "nGenJetAK8", "nGenIsolatedPhoton",
        "nGenDressedLepton", "nGenVisTau", "nPhoton", "nJet", "nMuon",
        "nElectron", "nTau", "nFatJet", "nTrigObj", "nSubJet", "nSV"};
    std::map<TString, int> maxValues;
    for (const auto &branchName : counterBranches) maxValues[branchName] = 0;

    TObjArray *fileElements = fChain->GetListOfFiles();
    for (int i = 0; i < fileElements->GetEntries(); i++) {
        TChainElement *element = (TChainElement *)fileElements->At(i);
        TFile *file = TFile::Open(element->GetTitle());
        TTree *tree = (TTree *)file->Get(fChain->GetName());
        for (const auto &branchName : counterBranches) {
            TLeaf *leaf = tree ? tree->GetLeaf(branchName) : nullptr;
            if (leaf) maxValues[branchName] = std::max(maxValues[branchName], static_cast<int>(leaf->GetMaximum()));
        }
        file->Close();
        delete file;
    }

    auto getMaxBranchValue = [this, &maxValues](const TString &branchName) {
        if (!fChain->GetBranch(branchName)) {
            cout << "[SKNanoGenLoader::SetMaxLeafSize] Warning: Branch " << branchName << " not found" << endl;
            return 0;
        }
        int maxValue = maxValues[branchName];
        cout << "[SKNanoLoader::SetMaxLeafSize] Branch: " << branchName << ", Max Value: " << maxValue << endl;
        return maxValue;
    };

    const UInt_t kMaxLHEPdfWeight = getMaxBranchValue("nLHEPdfWeight");
    const UInt_t kMaxLHEScaleWeight = getMaxBranchValue("nLHEScaleWeight");
    const UInt_t kMaxPSWeight = getMaxBranchValue("nPSWeight");
    const UInt_t kMaxLHEPart = getMaxBranchValue("nLHEPart");
    const UInt_t kMaxGenPart = getMaxBranchValue("nGenPart");
    const UInt_t kMaxGenJet = getMaxBranchValue("nGenJet");
    const UInt_t kMaxGenJetAK8 = getMaxBranchValue("nGenJetAK8");
    const UInt_t kMaxGenIsolatedPhoton = getMaxBranchValue("nGenIsolatedPhoton");
    const UInt_t kMaxGenDressedLepton = getMaxBranchValue("nGenDressedLepton");
    const UInt_t kMaxGenVisTau = getMaxBranchValue("nGenVisTau");
    const UInt_t kMaxPhoton = getMaxBranchValue("nPhoton");
    const UInt_t kMaxJet = getMaxBranchValue("nJet");
    const UInt_t kMaxMuon = getMaxBranchValue("nMuon");
    const UInt_t kMaxElectron = getMaxBranchValue("nElectron");
    const UInt_t kMaxTau = getMaxBranchValue("nTau");
    const UInt_t kMaxFatJet = getMaxBranchValue("nFatJet");
    const UInt_t kMaxTrigObj = getMaxBranchValue("nTrigObj");
    const UInt_t kMaxSubJet = getMaxBranchValue("nSubJet");
    const UInt_t kMaxSV = getMaxBranchValue("nSV");
    cout << "[SKNanoLoader::SetMaxLeafSize] Maximum Leaf Size Set" << endl;
    auto RDataFrameFinishTime = std::chrono::high_resolution_clock::now();
    
    // Now missing parts will automatically shrink
    LHEPdfWeight.resize(kMaxLHEPdfWeight);
    LHEScaleWeight.resize(kMaxLHEScaleWeight);
    PSWeight.resize(kMaxPSWeight);

    // LHEPart
    LHEPart_pt.resize(kMaxLHEPart);
    LHEPart_eta.resize(kMaxLHEPart);
    LHEPart_phi.resize(kMaxLHEPart);
    LHEPart_mass.resize(kMaxLHEPart);
    LHEPart_incomingpz.resize(kMaxLHEPart);
    LHEPart_pdgId.resize(kMaxLHEPart);
    LHEPart_status.resize(kMaxLHEPart);
    LHEPart_spin.resize(kMaxLHEPart);

    // GenPart
    GenPart_eta.resize(kMaxGenPart);
    GenPart_mass.resize(kMaxGenPart);
    GenPart_pdgId.resize(kMaxGenPart);
    GenPart_phi.resize(kMaxGenPart);
    GenPart_pt.resize(kMaxGenPart);
    GenPart_status.resize(kMaxGenPart);
    GenPart_genPartIdxMother.resize(kMaxGenPart);
    GenPart_statusFlags.resize(kMaxGenPart);
    GenPart_genPartIdxMother_RunII.resize(kMaxGenPart);
    GenPart_statusFlags_RunII.resize(kMaxGenPart);

    // GenJet
    GenJet_pt.resize(kMaxGenJet);
    GenJet_eta.resize(kMaxGenJet);
    GenJet_phi.resize(kMaxGenJet);
    GenJet_mass.resize(kMaxGenJet);
    GenJet_partonFlavour.resize(kMaxGenJet);
    GenJet_hadronFlavour.resize(kMaxGenJet);
    GenJet_partonFlavour_RunII.resize(kMaxGenJet);
    
    // GenJetAK8
    GenJetAK8_pt.resize(kMaxGenJetAK8);
    GenJetAK8_eta.resize(kMaxGenJetAK8);
    GenJetAK8_phi.resize(kMaxGenJetAK8);
    GenJetAK8_mass.resize(kMaxGenJetAK8);
    GenJetAK8_partonFlavour.resize(kMaxGenJetAK8);
    GenJetAK8_hadronFlavour.resize(kMaxGenJetAK8);
    GenJetAK8_partonFlavour_RunII.resize(kMaxGenJetAK8);

    // GenDressedLepton
    GenDressedLepton_pt.resize(kMaxGenDressedLepton);
    GenDressedLepton_eta.resize(kMaxGenDressedLepton);
    GenDressedLepton_phi.resize(kMaxGenDressedLepton);
    GenDressedLepton_mass.resize(kMaxGenDressedLepton);
    GenDressedLepton_pdgId.resize(kMaxGenDressedLepton);
    GenDressedLepton_hasTauAnc.resize(kMaxGenDressedLepton);

    // GenIsolatedPhoton
    GenIsolatedPhoton_pt.resize(kMaxGenIsolatedPhoton);
    GenIsolatedPhoton_eta.resize(kMaxGenIsolatedPhoton);
    GenIsolatedPhoton_phi.resize(kMaxGenIsolatedPhoton);
    GenIsolatedPhoton_mass.resize(kMaxGenIsolatedPhoton);

    // GenVisTau
    GenVisTau_pt.resize(kMaxGenVisTau);
    GenVisTau_eta.resize(kMaxGenVisTau);
    GenVisTau_phi.resize(kMaxGenVisTau);
    GenVisTau_mass.resize(kMaxGenVisTau);
    GenVisTau_charge.resize(kMaxGenVisTau);
    GenVisTau_genPartIdxMother.resize(kMaxGenVisTau);
    GenVisTau_status.resize(kMaxGenVisTau);
    Buf_GenVisTau_charge.resize(Run == 3 ? kMaxGenVisTau : 0);
    Buf_GenVisTau_genPartIdxMother.resize(Run == 3 ? kMaxGenVisTau : 0);
    Buf_GenVisTau_status.resize(Run == 3 ? kMaxGenVisTau : 0);

    // Muon----------------------------
    Muon_charge.resize(kMaxMuon);
    Muon_dxy.resize(kMaxMuon);
    Muon_dxyErr.resize(kMaxMuon);
    Muon_dxybs.resize(kMaxMuon);
    Muon_dz.resize(kMaxMuon);
    Muon_dzErr.resize(kMaxMuon);
    Muon_eta.resize(kMaxMuon);
    Muon_highPtId.resize(kMaxMuon);
    Muon_ip3d.resize(kMaxMuon);
    Muon_nTrackerLayers.resize(kMaxMuon);
    Buf_Muon_nTrackerLayers.resize(Run == 3 ? kMaxMuon : 0);
    Muon_isGlobal.resize(kMaxMuon);
    Muon_isStandalone.resize(kMaxMuon);
    Muon_isTracker.resize(kMaxMuon);
    Muon_looseId.resize(kMaxMuon);
    Muon_mass.resize(kMaxMuon);
    Muon_mediumId.resize(kMaxMuon);
    Muon_mediumPromptId.resize(kMaxMuon);
    Muon_miniIsoId.resize(kMaxMuon);
    Muon_miniPFRelIso_all.resize(kMaxMuon);
    Muon_multiIsoId.resize(kMaxMuon);
    Muon_mvaLowPt.resize(kMaxMuon);
    Muon_mvaTTH.resize(kMaxMuon);
    Muon_pfIsoId.resize(kMaxMuon);
    Muon_pfRelIso03_all.resize(kMaxMuon);
    Muon_pfRelIso04_all.resize(kMaxMuon);
    Muon_phi.resize(kMaxMuon);
    Muon_pt.resize(kMaxMuon);
    Muon_puppiIsoId.resize(kMaxMuon);
    Muon_sip3d.resize(kMaxMuon);
    Muon_softId.resize(kMaxMuon);
    Muon_softMva.resize(kMaxMuon);
    Muon_softMvaId.resize(kMaxMuon);
    Muon_tightId.resize(kMaxMuon);
    Muon_tkIsoId.resize(kMaxMuon);
    Muon_tkRelIso.resize(kMaxMuon);
    Muon_triggerIdLoose.resize(kMaxMuon);
    Muon_genPartFlav.resize(kMaxMuon);
    if(Run == 3){ 
        Muon_mvaMuID_WP.resize(kMaxMuon);
        Muon_genPartIdx.resize(kMaxMuon);
        Muon_jetIdx.resize(kMaxMuon);
        Muon_genPartIdx_RunII.resize(0);
        Muon_jetIdx_RunII.resize(0);
        Muon_mvaId.resize(0);
    }
    else if(Run == 2){
        Muon_mvaMuID_WP.resize(0);
        Muon_jetIdx.resize(0);
        Muon_genPartIdx.resize(0);
        Muon_genPartIdx_RunII.resize(kMaxMuon);
        Muon_jetIdx_RunII.resize(kMaxMuon);
        Muon_mvaId.resize(kMaxMuon);
    }
    // Electron----------------------------
    Electron_charge.resize(kMaxElectron);
    Electron_convVeto.resize(kMaxElectron);
    Electron_cutBased_HEEP.resize(kMaxElectron);
    Electron_scEta.resize(kMaxElectron);
    Electron_deltaEtaInSC.resize(kMaxElectron);
    Electron_deltaEtaInSeed.resize(kMaxElectron);
    Electron_deltaPhiInSC.resize(kMaxElectron);
    Electron_deltaPhiInSeed.resize(kMaxElectron);
    Electron_ecalPFClusterIso.resize(kMaxElectron);
    Electron_hcalPFClusterIso.resize(kMaxElectron);
    Electron_dr03EcalRecHitSumEt.resize(kMaxElectron);
    Electron_dr03HcalDepth1TowerSumEt.resize(kMaxElectron);
    Electron_dr03TkSumPt.resize(kMaxElectron);
    Electron_dr03TkSumPtHEEP.resize(kMaxElectron);
    Electron_dxy.resize(kMaxElectron);
    Electron_dxyErr.resize(kMaxElectron);
    Electron_dz.resize(kMaxElectron);
    Electron_dzErr.resize(kMaxElectron);
    Electron_eInvMinusPInv.resize(kMaxElectron);
    Electron_energyErr.resize(kMaxElectron);
    Electron_eta.resize(kMaxElectron);
    Electron_hoe.resize(kMaxElectron);
    Electron_ip3d.resize(kMaxElectron);
    Electron_isPFcand.resize(kMaxElectron);
    Electron_jetNDauCharged.resize(kMaxElectron);
    Electron_jetPtRelv2.resize(kMaxElectron);
    Electron_jetRelIso.resize(kMaxElectron);
    Electron_lostHits.resize(kMaxElectron);
    Electron_mass.resize(kMaxElectron);
    Electron_miniPFRelIso_all.resize(kMaxElectron);
    Electron_miniPFRelIso_chg.resize(kMaxElectron);
    Electron_mvaTTH.resize(kMaxElectron);
    Electron_pdgId.resize(kMaxElectron);
    Electron_pfRelIso03_all.resize(kMaxElectron);
    Electron_pfRelIso03_chg.resize(kMaxElectron);
    Electron_phi.resize(kMaxElectron);
    Electron_pt.resize(kMaxElectron);
    Electron_r9.resize(kMaxElectron);
    Electron_scEtOverPt.resize(kMaxElectron);
    Electron_seedGain.resize(kMaxElectron);
    Electron_sieie.resize(kMaxElectron);
    Electron_sip3d.resize(kMaxElectron);
    Electron_genPartFlav.resize(kMaxElectron);
    if(Run == 3){
        Electron_cutBased.resize(kMaxElectron);
        Electron_genPartIdx.resize(kMaxElectron);
        Electron_jetIdx.resize(kMaxElectron);
        Electron_mvaIso.resize(kMaxElectron);
        Electron_mvaIso_WP80.resize(kMaxElectron);
        Electron_mvaIso_WP90.resize(kMaxElectron);
        Electron_mvaIso_WPL.resize(kMaxElectron);
        Electron_mvaNoIso.resize(kMaxElectron);
        Electron_mvaNoIso_WP80.resize(kMaxElectron);
        Electron_mvaNoIso_WP90.resize(kMaxElectron);
        Electron_mvaNoIso_WPL.resize(kMaxElectron);
        Electron_cutBased_RunII.resize(0);
        Electron_genPartIdx_RunII.resize(0);
        Electron_jetIdx_RunII.resize(0);
        Electron_mvaFall17V2Iso.resize(0);
        Electron_mvaFall17V2Iso_WP80.resize(0);
        Electron_mvaFall17V2Iso_WP90.resize(0);
        Electron_mvaFall17V2Iso_WPL.resize(0);
        Electron_mvaFall17V2noIso.resize(0);
        Electron_mvaFall17V2noIso_WP80.resize(0);
        Electron_mvaFall17V2noIso_WP90.resize(0);
        Electron_mvaFall17V2noIso_WPL.resize(0);
    }
    else if(Run == 2){
        Electron_cutBased.resize(0);
        Electron_genPartIdx.resize(0);
        Electron_jetIdx.resize(0);
        Electron_mvaIso.resize(0);
        Electron_mvaIso_WP80.resize(0);
        Electron_mvaIso_WP90.resize(0);
        Electron_mvaIso_WPL.resize(0);
        Electron_mvaNoIso.resize(0);
        Electron_mvaNoIso_WP80.resize(0);
        Electron_mvaNoIso_WP90.resize(0);
        Electron_mvaNoIso_WPL.resize(0);
        Electron_cutBased_RunII.resize(kMaxElectron);
        Electron_genPartIdx_RunII.resize(kMaxElectron);
        Electron_jetIdx_RunII.resize(kMaxElectron);
        Electron_mvaFall17V2Iso.resize(kMaxElectron);
        Electron_mvaFall17V2Iso_WP80.resize(kMaxElectron);
        Electron_mvaFall17V2Iso_WP90.resize(kMaxElectron);
        Electron_mvaFall17V2Iso_WPL.resize(kMaxElectron);
        Electron_mvaFall17V2noIso.resize(kMaxElectron);
        Electron_mvaFall17V2noIso_WP80.resize(kMaxElectron);
        Electron_mvaFall17V2noIso_WP90.resize(kMaxElectron);
        Electron_mvaFall17V2noIso_WPL.resize(kMaxElectron);
        Electron_dEsigmaUp.resize(kMaxElectron);
        Electron_dEsigmaDown.resize(kMaxElectron);
    } 
    
    //Photon----------------------------
    Photon_energyErr.resize(kMaxPhoton);
    Photon_eta.resize(kMaxPhoton);
    Photon_hoe.resize(kMaxPhoton);
    Photon_isScEtaEB.resize(kMaxPhoton);
    Photon_isScEtaEE.resize(kMaxPhoton);
    Photon_mvaID.resize(kMaxPhoton);
    Photon_mvaID_WP80.resize(kMaxPhoton);
    Photon_mvaID_WP90.resize(kMaxPhoton);
    Photon_phi.resize(kMaxPhoton);
    Photon_pt.resize(kMaxPhoton);
    Photon_sieie.resize(kMaxPhoton);
    if(Run == 3){
        Photon_cutBased.resize(kMaxPhoton);
        Photon_energyRaw.resize(kMaxPhoton);
        Photon_cutBased_RunII.resize(kMaxPhoton);
    }
    else if(Run == 2){
        Photon_cutBased.resize(0);
        Photon_energyRaw.resize(0);
        Photon_cutBased_RunII.resize(kMaxPhoton);
    }

    //Jet----------------------------
    Jet_area.resize(kMaxJet);
    Jet_btagDeepFlavB.resize(kMaxJet);
    Jet_btagDeepFlavCvB.resize(kMaxJet);
    Jet_btagDeepFlavCvL.resize(kMaxJet);
    Jet_btagDeepFlavQG.resize(kMaxJet);
    Jet_chEmEF.resize(kMaxJet);
    Jet_chHEF.resize(kMaxJet);
    Jet_eta.resize(kMaxJet);
    Jet_hfadjacentEtaStripsSize.resize(kMaxJet);
    Jet_hfcentralEtaStripSize.resize(kMaxJet);
    Jet_hfsigmaEtaEta.resize(kMaxJet);
    Jet_hfsigmaPhiPhi.resize(kMaxJet);
    Jet_mass.resize(kMaxJet);
    Jet_muEF.resize(kMaxJet);
    Jet_muonSubtrFactor.resize(kMaxJet);
    Jet_nConstituents.resize(kMaxJet);
    Jet_neEmEF.resize(kMaxJet);
    Jet_neHEF.resize(kMaxJet);
    Jet_phi.resize(kMaxJet);
    Jet_pt.resize(kMaxJet);
    Jet_rawFactor.resize(kMaxJet);
    if(Run == 3){
        Jet_PNetRegPtRawCorr.resize(kMaxJet);
        Jet_PNetRegPtRawCorrNeutrino.resize(kMaxJet);
        Jet_PNetRegPtRawRes.resize(kMaxJet);
        Jet_btagPNetB.resize(kMaxJet);
        Jet_btagPNetCvB.resize(kMaxJet);
        Jet_btagPNetCvL.resize(kMaxJet);
        Jet_btagPNetQvG.resize(kMaxJet);
        Jet_btagPNetTauVJet.resize(kMaxJet);
        Jet_btagRobustParTAK4B.resize(kMaxJet);
        Jet_btagRobustParTAK4CvB.resize(kMaxJet);
        Jet_btagRobustParTAK4CvL.resize(kMaxJet);
        Jet_btagRobustParTAK4QG.resize(kMaxJet);
        Jet_electronIdx1.resize(kMaxJet);
        Jet_electronIdx2.resize(kMaxJet);
        Jet_genJetIdx.resize(kMaxJet);
        Jet_hadronFlavour.resize(kMaxJet);
        Jet_jetId.resize(kMaxJet);
        Jet_muonIdx1.resize(kMaxJet);
        Jet_muonIdx2.resize(kMaxJet);
        Jet_nElectrons.resize(kMaxJet);
        Jet_nMuons.resize(kMaxJet);
        Jet_nSVs.resize(kMaxJet);
        Jet_partonFlavour.resize(kMaxJet);
        Jet_svIdx1.resize(kMaxJet);
        Jet_svIdx2.resize(kMaxJet);
        Jet_chMultiplicity.resize(kMaxJet);
        Jet_neMultiplicity.resize(kMaxJet);
        Buf_Jet_chMultiplicity.resize(kMaxJet);
        Buf_Jet_neMultiplicity.resize(kMaxJet);
        Jet_bRegCorr.resize(0);
        Jet_bRegRes.resize(0);
        Jet_btagCSVV2.resize(0);
        //Jet_btagDeepB.resize(0);
        //Jet_btagDeepCvB.resize(0);
        //Jet_btagDeepCvL.resize(0);
        Jet_cRegCorr.resize(0);
        Jet_cRegRes.resize(0);
        Jet_chFPV0EF.resize(0);
        Jet_cleanmask.resize(0);
        Jet_electronIdx1_RunII.resize(0);
        Jet_electronIdx2_RunII.resize(0);
        Jet_genJetIdx_RunII.resize(0);
        Jet_hadronFlavour_RunII.resize(0);
        Jet_jetId_RunII.resize(0);
        Jet_muonIdx1_RunII.resize(0);
        Jet_muonIdx2_RunII.resize(0);
        Jet_nElectrons_RunII.resize(0);
        Jet_nMuons_RunII.resize(0);
        Jet_partonFlavour_RunII.resize(0);
        Jet_puId.resize(0);
        Jet_puIdDisc.resize(0);
        Jet_qgl.resize(0);
    }
    else if(Run == 2){
        Jet_PNetRegPtRawCorr.resize(0);
        Jet_PNetRegPtRawCorrNeutrino.resize(0);
        Jet_PNetRegPtRawRes.resize(0);
        Jet_btagPNetB.resize(0);
        Jet_btagPNetCvB.resize(0);
        Jet_btagPNetCvL.resize(0);
        Jet_btagPNetQvG.resize(0);
        Jet_btagPNetTauVJet.resize(0);
        Jet_btagRobustParTAK4B.resize(0);
        Jet_btagRobustParTAK4CvB.resize(0);
        Jet_btagRobustParTAK4CvL.resize(0);
        Jet_btagRobustParTAK4QG.resize(0);
        Jet_electronIdx1.resize(0);
        Jet_electronIdx2.resize(0);
        Jet_genJetIdx.resize(0);
        Jet_hadronFlavour.resize(0);
        Jet_jetId.resize(0);
        Jet_muonIdx1.resize(0);
        Jet_muonIdx2.resize(0);
        Jet_nElectrons.resize(0);
        Jet_nMuons.resize(0);
        Jet_nSVs.resize(0);
        Jet_partonFlavour.resize(0);
        Jet_svIdx1.resize(0);
        Jet_svIdx2.resize(0);
        Jet_bRegCorr.resize(kMaxJet);
        Jet_bRegRes.resize(kMaxJet);
        Jet_btagCSVV2.resize(kMaxJet);
        //Jet_btagDeepB.resize(kMaxJet);
        //Jet_btagDeepCvB.resize(kMaxJet);
        //Jet_btagDeepCvL.resize(kMaxJet);
        Jet_cRegCorr.resize(kMaxJet);
        Jet_cRegRes.resize(kMaxJet);
        Jet_chFPV0EF.resize(kMaxJet);
        Jet_cleanmask.resize(kMaxJet);
        Jet_electronIdx1_RunII.resize(kMaxJet);
        Jet_electronIdx2_RunII.resize(kMaxJet);
        Jet_genJetIdx_RunII.resize(kMaxJet);
        Jet_hadronFlavour_RunII.resize(kMaxJet);
        Jet_jetId_RunII.resize(kMaxJet);
        Jet_muonIdx1_RunII.resize(kMaxJet);
        Jet_muonIdx2_RunII.resize(kMaxJet);
        Jet_nElectrons_RunII.resize(kMaxJet);
        Jet_nMuons_RunII.resize(kMaxJet);
        Jet_partonFlavour_RunII.resize(kMaxJet);
        Jet_puId.resize(kMaxJet);
        Jet_puIdDisc.resize(kMaxJet);
        Jet_qgl.resize(kMaxJet);
    }

    //Tau----------------------------
    Tau_dxy.resize(kMaxTau);
    Tau_dz.resize(kMaxTau);
    Tau_eta.resize(kMaxTau);
    Tau_genPartFlav.resize(kMaxTau);
    Tau_idDeepTau2017v2p1VSe.resize(kMaxTau);
    Tau_idDeepTau2017v2p1VSjet.resize(kMaxTau);
    Tau_idDeepTau2017v2p1VSmu.resize(kMaxTau);
    Tau_mass.resize(kMaxTau);
    Tau_phi.resize(kMaxTau);
    Tau_pt.resize(kMaxTau);
    // Run3
    if(Run == 3){
        Tau_charge.resize(kMaxTau);
        Tau_decayMode.resize(kMaxTau);
        Tau_genPartIdx.resize(kMaxTau);
        Tau_idDecayModeNewDMs.resize(kMaxTau);
        Tau_idDeepTau2018v2p5VSe.resize(kMaxTau);
        Tau_idDeepTau2018v2p5VSjet.resize(kMaxTau);
        Tau_idDeepTau2018v2p5VSmu.resize(kMaxTau);
        Tau_charge_RunII.resize(0);
        Tau_decayMode_RunII.resize(0);
        Tau_genPartIdx_RunII.resize(0);
    }
    else if(Run == 2){
        Tau_charge.resize(0);
        Tau_decayMode.resize(0);
        Tau_genPartIdx.resize(0);
        Tau_idDecayModeNewDMs.resize(0);
        Tau_idDeepTau2018v2p5VSe.resize(0);
        Tau_idDeepTau2018v2p5VSjet.resize(0);
        Tau_idDeepTau2018v2p5VSmu.resize(0);
        Tau_charge_RunII.resize(kMaxTau);
        Tau_decayMode_RunII.resize(kMaxTau);
        Tau_genPartIdx_RunII.resize(kMaxTau);
    }

    // FatJet----------------------------
    FatJet_area.resize(kMaxFatJet);
    FatJet_btagDDBvLV2.resize(kMaxFatJet);
    FatJet_btagDDCvBV2.resize(kMaxFatJet);
    FatJet_btagDDCvLV2.resize(kMaxFatJet);
    FatJet_btagDeepB.resize(kMaxFatJet);
    FatJet_btagHbb.resize(kMaxFatJet);
    FatJet_eta.resize(kMaxFatJet);
    FatJet_lsf3.resize(kMaxFatJet);
    FatJet_mass.resize(kMaxFatJet);
    FatJet_msoftdrop.resize(kMaxFatJet);
    FatJet_nBHadrons.resize(kMaxFatJet);
    FatJet_nCHadrons.resize(kMaxFatJet);
    FatJet_nConstituents.resize(kMaxFatJet);
    FatJet_particleNet_QCD.resize(kMaxFatJet);
    FatJet_phi.resize(kMaxFatJet);
    FatJet_pt.resize(kMaxFatJet);
    FatJet_tau1.resize(kMaxFatJet);
    FatJet_tau2.resize(kMaxFatJet);
    FatJet_tau3.resize(kMaxFatJet);
    FatJet_tau4.resize(kMaxFatJet);
    if(Run == 3){
        FatJet_genJetAK8Idx.resize(kMaxFatJet);
        FatJet_genJetAK8Idx_RunII.resize(0);
        FatJet_jetId.resize(kMaxFatJet);
        FatJet_particleNetWithMass_H4qvsQCD.resize(kMaxFatJet);
        FatJet_particleNetWithMass_HbbvsQCD.resize(kMaxFatJet);
        FatJet_particleNetWithMass_HccvsQCD.resize(kMaxFatJet);
        FatJet_particleNetWithMass_QCD.resize(kMaxFatJet);
        FatJet_particleNetWithMass_TvsQCD.resize(kMaxFatJet);
        FatJet_particleNetWithMass_WvsQCD.resize(kMaxFatJet);
        FatJet_particleNetWithMass_ZvsQCD.resize(kMaxFatJet);
        FatJet_particleNet_QCD0HF.resize(kMaxFatJet);
        FatJet_particleNet_QCD1HF.resize(kMaxFatJet);
        FatJet_particleNet_QCD2HF.resize(kMaxFatJet);
        FatJet_particleNet_XbbVsQCD.resize(kMaxFatJet);
        FatJet_particleNet_XccVsQCD.resize(kMaxFatJet);
        FatJet_particleNet_XggVsQCD.resize(kMaxFatJet);
        FatJet_particleNet_XqqVsQCD.resize(kMaxFatJet);
        FatJet_particleNet_XteVsQCD.resize(kMaxFatJet);
        FatJet_particleNet_XtmVsQCD.resize(kMaxFatJet);
        FatJet_particleNet_XttVsQCD.resize(kMaxFatJet);
        FatJet_particleNet_massCorr.resize(kMaxFatJet);
        FatJet_subJetIdx1.resize(kMaxFatJet);
        FatJet_subJetIdx2.resize(kMaxFatJet);
        FatJet_n2b1.resize(kMaxFatJet);
        FatJet_n3b1.resize(kMaxFatJet);
        FatJet_chMultiplicity.resize(kMaxFatJet);
        FatJet_neMultiplicity.resize(kMaxFatJet);
        FatJet_chHEF.resize(kMaxFatJet);
        FatJet_neHEF.resize(kMaxFatJet);
        FatJet_chEmEF.resize(kMaxFatJet);
        FatJet_neEmEF.resize(kMaxFatJet);
        FatJet_muEF.resize(kMaxFatJet);
        SubJet_pt.resize(kMaxSubJet);
        SubJet_eta.resize(kMaxSubJet);
        SubJet_phi.resize(kMaxSubJet);
        SubJet_mass.resize(kMaxSubJet);
        SubJet_btagDeepB.resize(kMaxSubJet);
        SV_pt.resize(kMaxSV);
        SV_eta.resize(kMaxSV);
        SV_phi.resize(kMaxSV);
        SV_mass.resize(kMaxSV);
        SV_dlenSig.resize(kMaxSV);
        SV_dxySig.resize(kMaxSV);
        SV_chi2.resize(kMaxSV);
        SV_pAngle.resize(kMaxSV);
        SV_ntracks.resize(kMaxSV);
        FatJet_jetId_RunII.resize(0);
        FatJet_particleNetMD_QCD.resize(0);
        FatJet_particleNetMD_Xbb.resize(0);
        FatJet_particleNetMD_Xcc.resize(0);
        FatJet_particleNetMD_Xqq.resize(0);
        FatJet_particleNet_H4qvsQCD.resize(0);
        FatJet_particleNet_HbbvsQCD.resize(0);
        FatJet_particleNet_HccvsQCD.resize(0);
        FatJet_particleNet_TvsQCD.resize(0);
        FatJet_particleNet_WvsQCD.resize(0);
        FatJet_particleNet_ZvsQCD.resize(0);
        FatJet_particleNet_mass.resize(0);
        FatJet_subJetIdx1_RunII.resize(0);
        FatJet_subJetIdx2_RunII.resize(0);
    }
    else if(Run == 2){
        FatJet_genJetAK8Idx.resize(0);
        FatJet_genJetAK8Idx_RunII.resize(kMaxFatJet);
        FatJet_jetId.resize(0);
        FatJet_particleNetWithMass_H4qvsQCD.resize(0);
        FatJet_particleNetWithMass_HbbvsQCD.resize(0);
        FatJet_particleNetWithMass_HccvsQCD.resize(0);
        FatJet_particleNetWithMass_QCD.resize(0);
        FatJet_particleNetWithMass_TvsQCD.resize(0);
        FatJet_particleNetWithMass_WvsQCD.resize(0);
        FatJet_particleNetWithMass_ZvsQCD.resize(0);
        FatJet_particleNet_QCD0HF.resize(0);
        FatJet_particleNet_QCD1HF.resize(0);
        FatJet_particleNet_QCD2HF.resize(0);
        FatJet_particleNet_XbbVsQCD.resize(0);
        FatJet_particleNet_XccVsQCD.resize(0);
        FatJet_particleNet_XggVsQCD.resize(0);
        FatJet_particleNet_XqqVsQCD.resize(0);
        FatJet_particleNet_XteVsQCD.resize(0);
        FatJet_particleNet_XtmVsQCD.resize(0);
        FatJet_particleNet_XttVsQCD.resize(0);
        FatJet_particleNet_massCorr.resize(0);
        FatJet_subJetIdx1.resize(0);
        FatJet_subJetIdx2.resize(0);
        FatJet_n2b1.resize(0);
        FatJet_n3b1.resize(0);
        FatJet_chMultiplicity.resize(0);
        FatJet_neMultiplicity.resize(0);
        FatJet_chHEF.resize(0);
        FatJet_neHEF.resize(0);
        FatJet_chEmEF.resize(0);
        FatJet_neEmEF.resize(0);
        FatJet_muEF.resize(0);
        SubJet_pt.resize(0);
        SubJet_eta.resize(0);
        SubJet_phi.resize(0);
        SubJet_mass.resize(0);
        SubJet_btagDeepB.resize(0);
        SV_pt.resize(0);
        SV_eta.resize(0);
        SV_phi.resize(0);
        SV_mass.resize(0);
        SV_dlenSig.resize(0);
        SV_dxySig.resize(0);
        SV_chi2.resize(0);
        SV_pAngle.resize(0);
        SV_ntracks.resize(0);
        FatJet_jetId_RunII.resize(kMaxFatJet);
        FatJet_particleNetMD_QCD.resize(kMaxFatJet);
        FatJet_particleNetMD_Xbb.resize(kMaxFatJet);
        FatJet_particleNetMD_Xcc.resize(kMaxFatJet);
        FatJet_particleNetMD_Xqq.resize(kMaxFatJet);
        FatJet_particleNet_H4qvsQCD.resize(kMaxFatJet);
        FatJet_particleNet_HbbvsQCD.resize(kMaxFatJet);
        FatJet_particleNet_HccvsQCD.resize(kMaxFatJet);
        FatJet_particleNet_TvsQCD.resize(kMaxFatJet);
        FatJet_particleNet_WvsQCD.resize(kMaxFatJet);
        FatJet_particleNet_ZvsQCD.resize(kMaxFatJet);
        FatJet_particleNet_mass.resize(kMaxFatJet);
        FatJet_subJetIdx1_RunII.resize(kMaxFatJet);
        FatJet_subJetIdx2_RunII.resize(kMaxFatJet);
    }

    // TrigObj----------------------------
    TrigObj_pt.resize(kMaxTrigObj);
    TrigObj_eta.resize(kMaxTrigObj);
    TrigObj_phi.resize(kMaxTrigObj);
    TrigObj_id.resize(kMaxTrigObj);
    TrigObj_id_RunII.resize(kMaxTrigObj);
    TrigObj_filterBits.resize(kMaxTrigObj);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    auto RDataFrameDuration = std::chrono::duration_cast<std::chrono::milliseconds>(RDataFrameFinishTime - start);
    auto resizingDuration = std::chrono::duration_cast<std::chrono::milliseconds>(end - RDataFrameFinishTime);
    cout << "[SKNanoLoader::SetMaxLeafSize] Resizing Time: " << resizingDuration.count() << " ms" << endl;
    cout << "[SKNanoLoader::SetMaxLeafSize] Max-value metadata scan Time: " << RDataFrameDuration.count() << " ms" << endl;
    cout << "[SKNanoLoader::SetMaxLeafSize] Time taken: " << duration.count() << " ms" << endl;
}

void SKNanoLoader::Init() {
    cout << "[SKNanoLoader::Init] Initializing. Era = " << DataEra << " Run =  " << Run << endl;
    if(fChain->GetEntries() == 0) {
        cout << "[SKNanoLoader::Init] No Entries in the Tree" << endl;
        cout << "[SKNanoLoader::Init] Exiting without make output..." << endl;
        exit(0);
    }

    // Skimmers copy raw events with fChain->CloneTree(0) + Fill(), which only
    // works when fChain->GetEntry() fills the SetBranchAddress buffers, so
    // they must use the legacy reading mode.
    if (!SkimmingMode && TString(typeid(*this).name()).Contains("Skim")) {
        cout << "[SKNanoLoader::Init] Skim analyzer detected: enabling SkimmingMode (legacy branch-address reading)" << endl;
        SkimmingMode = true;
    }

    if (SkimmingMode) InitLegacy();
    else InitTTreeReader();
}

namespace {
// TTreeReader requires the proxy template type to match the on-file leaf type
// exactly, unlike SetBranchAddress(void*) which copies raw bytes with no type
// check. The on-file type is therefore looked up at Init time and each value
// is static_cast into the member, keeping the member types (the interface the
// analyzers use) unchanged across NanoAOD versions.
template <typename FileT, typename MemT>
void MakeScalarFiller(TTreeReader &reader, const TString &branchName, MemT &dest,
                      std::vector<std::function<void()>> &fillers) {
    auto rdr = std::make_shared<TTreeReaderValue<FileT>>(reader, branchName.Data());
    fillers.push_back([rdr, &dest]() { dest = static_cast<MemT>(**rdr); });
}

template <typename FileT, typename MemT>
void MakeArrayFiller(TTreeReader &reader, const TString &branchName, RVec<MemT> &dest,
                     std::vector<std::function<void()>> &fillers) {
    auto rdr = std::make_shared<TTreeReaderArray<FileT>>(reader, branchName.Data());
    fillers.push_back([rdr, &dest]() {
        const size_t n = rdr->GetSize();
        dest.resize(n);
        for (size_t i = 0; i < n; i++) dest[i] = static_cast<MemT>((*rdr)[i]);
    });
}

template <typename MemT>
bool AddScalarFiller(TTreeReader &reader, const TString &branchName, const TString &leafType,
                     MemT &dest, std::vector<std::function<void()>> &fillers) {
    if      (leafType == "Float_t")   MakeScalarFiller<Float_t>(reader, branchName, dest, fillers);
    else if (leafType == "Double_t")  MakeScalarFiller<Double_t>(reader, branchName, dest, fillers);
    else if (leafType == "Int_t")     MakeScalarFiller<Int_t>(reader, branchName, dest, fillers);
    else if (leafType == "UInt_t")    MakeScalarFiller<UInt_t>(reader, branchName, dest, fillers);
    else if (leafType == "Bool_t")    MakeScalarFiller<Bool_t>(reader, branchName, dest, fillers);
    else if (leafType == "Char_t")    MakeScalarFiller<Char_t>(reader, branchName, dest, fillers);
    else if (leafType == "UChar_t")   MakeScalarFiller<UChar_t>(reader, branchName, dest, fillers);
    else if (leafType == "Short_t")   MakeScalarFiller<Short_t>(reader, branchName, dest, fillers);
    else if (leafType == "UShort_t")  MakeScalarFiller<UShort_t>(reader, branchName, dest, fillers);
    else if (leafType == "Long64_t")  MakeScalarFiller<Long64_t>(reader, branchName, dest, fillers);
    else if (leafType == "ULong64_t") MakeScalarFiller<ULong64_t>(reader, branchName, dest, fillers);
    else return false;
    return true;
}

template <typename MemT>
bool AddArrayFiller(TTreeReader &reader, const TString &branchName, const TString &leafType,
                    RVec<MemT> &dest, std::vector<std::function<void()>> &fillers) {
    if      (leafType == "Float_t")   MakeArrayFiller<Float_t>(reader, branchName, dest, fillers);
    else if (leafType == "Double_t")  MakeArrayFiller<Double_t>(reader, branchName, dest, fillers);
    else if (leafType == "Int_t")     MakeArrayFiller<Int_t>(reader, branchName, dest, fillers);
    else if (leafType == "UInt_t")    MakeArrayFiller<UInt_t>(reader, branchName, dest, fillers);
    else if (leafType == "Bool_t")    MakeArrayFiller<Bool_t>(reader, branchName, dest, fillers);
    else if (leafType == "Char_t")    MakeArrayFiller<Char_t>(reader, branchName, dest, fillers);
    else if (leafType == "UChar_t")   MakeArrayFiller<UChar_t>(reader, branchName, dest, fillers);
    else if (leafType == "Short_t")   MakeArrayFiller<Short_t>(reader, branchName, dest, fillers);
    else if (leafType == "UShort_t")  MakeArrayFiller<UShort_t>(reader, branchName, dest, fillers);
    else if (leafType == "Long64_t")  MakeArrayFiller<Long64_t>(reader, branchName, dest, fillers);
    else if (leafType == "ULong64_t") MakeArrayFiller<ULong64_t>(reader, branchName, dest, fillers);
    else return false;
    return true;
}
} // namespace

void SKNanoLoader::InitTTreeReader() {
    fReader = std::make_unique<TTreeReader>(fChain);
    fScalarFillers.clear();
    fArrayFillers.clear();

    // TTreeReader only reads the branches that have a proxy, so the legacy
    // SetBranchStatus("*", 0) optimization is automatic here.

    auto BindScalar = [this](const TString &branchName, auto &dest) {
        TBranch *branch = fChain->GetBranch(branchName);
        if (!branch) {
            cout << "[SKNanoGenLoader::Init] Warning:Branch " << branchName << " not found" << endl;
            return;
        }
        TLeaf *leaf = (TLeaf *)branch->GetListOfLeaves()->At(0);
        if (!AddScalarFiller(*fReader, branchName, leaf->GetTypeName(), dest, fScalarFillers))
            cout << "[SKNanoLoader::Init] Warning: Branch " << branchName << " has unsupported leaf type " << leaf->GetTypeName() << endl;
    };

    // In the legacy mode an array member whose branch is absent (or that was
    // never bound at all) stayed zero-filled at the collection's max size, so
    // in-range reads returned 0; reproduce that with a per-event zero-fill of
    // the current collection size.
    auto ZeroFillArray = [this](const TString & /*branchName*/, auto &dest, const Int_t *count, const TString &countBranchName) {
        if (!fChain->GetBranch(countBranchName)) return;
        fArrayFillers.push_back([&dest, count]() {
            dest.clear();
            dest.resize(std::max(*count, 0)); // value-initialized -> zeros
        });
    };

    auto BindArray = [this, &ZeroFillArray](const TString &branchName, auto &dest,
                                            const Int_t *count, const TString &countBranchName) {
        TBranch *branch = fChain->GetBranch(branchName);
        if (!branch) {
            cout << "[SKNanoGenLoader::Init] Warning:Branch " << branchName << " not found" << endl;
            if (count) ZeroFillArray(branchName, dest, count, countBranchName);
            return;
        }
        TLeaf *leaf = (TLeaf *)branch->GetListOfLeaves()->At(0);
        if (!AddArrayFiller(*fReader, branchName, leaf->GetTypeName(), dest, fArrayFillers))
            cout << "[SKNanoLoader::Init] Warning: Branch " << branchName << " has unsupported leaf type " << leaf->GetTypeName() << endl;
    };

    // For type conversion between Run2 and Run3
    auto BindScalarWithRunCheck = [this, &BindScalar](const TString &branchName, Int_t &run3Var, UInt_t &runIIVar) {
        if (Run == 3) {
            BindScalar(branchName, run3Var);
        } else {
            BindScalar(branchName, runIIVar);
        }
    };

    // Weights
    BindScalar("genWeight", genWeight);
    BindScalar("LHEWeight_originalXWGTUP", LHEWeight_originalXWGTUP);
    BindScalar("Generator_weight", Generator_weight);
    BindScalar("nLHEPdfWeight", nLHEPdfWeight);
    BindScalar("nLHEScaleWeight", nLHEScaleWeight);
    BindScalar("nPSWeight", nPSWeight);
    BindArray("LHEPdfWeight", LHEPdfWeight, &nLHEPdfWeight, "nLHEPdfWeight");
    BindArray("LHEScaleWeight", LHEScaleWeight, &nLHEScaleWeight, "nLHEScaleWeight");
    BindArray("PSWeight", PSWeight, &nPSWeight, "nPSWeight");

    // PDFs
    BindScalar("Generator_id1", Generator_id1);
    BindScalar("Generator_id2", Generator_id2);
    BindScalar("Generator_x1", Generator_x1);
    BindScalar("Generator_x2", Generator_x2);
    BindScalar("Generator_xpdf1", Generator_xpdf1);
    BindScalar("Generator_xpdf2", Generator_xpdf2);
    BindScalar("Generator_scalePDF", Generator_scalePDF);

    // LHE
    BindScalar("LHE_HT", LHE_HT);
    BindScalar("LHE_HTIncoming", LHE_HTIncoming);
    BindScalar("LHE_Vpt", LHE_Vpt);
    BindScalar("LHE_AlphaS", LHE_AlphaS);
    BindScalar("LHE_Njets", LHE_Njets);
    BindScalar("LHE_Nb", LHE_Nb);
    BindScalar("LHE_Nc", LHE_Nc);
    BindScalar("LHE_Nuds", LHE_Nuds);
    BindScalar("LHE_Nglu", LHE_Nglu);
    BindScalar("LHE_NpLO", LHE_NpLO);
    BindScalar("LHE_NpNLO", LHE_NpNLO);

    // LHEPart
    BindScalarWithRunCheck("nLHEPart", nLHEPart, nLHEPart_RunII);
    BindArray("LHEPart_pt", LHEPart_pt, &nLHEPart, "nLHEPart");
    BindArray("LHEPart_eta", LHEPart_eta, &nLHEPart, "nLHEPart");
    BindArray("LHEPart_phi", LHEPart_phi, &nLHEPart, "nLHEPart");
    BindArray("LHEPart_mass", LHEPart_mass, &nLHEPart, "nLHEPart");
    BindArray("LHEPart_pdgId", LHEPart_pdgId, &nLHEPart, "nLHEPart");
    BindArray("LHEPart_status", LHEPart_status, &nLHEPart, "nLHEPart");
    BindArray("LHEPart_spin", LHEPart_spin, &nLHEPart, "nLHEPart");
    BindArray("LHEPart_incomingpz", LHEPart_incomingpz, &nLHEPart, "nLHEPart");

    // GenPart
    BindScalarWithRunCheck("nGenPart", nGenPart, nGenPart_RunII);
    BindArray("GenPart_eta", GenPart_eta, &nGenPart, "nGenPart");
    BindArray("GenPart_mass", GenPart_mass, &nGenPart, "nGenPart");
    BindArray("GenPart_pdgId", GenPart_pdgId, &nGenPart, "nGenPart");
    BindArray("GenPart_phi", GenPart_phi, &nGenPart, "nGenPart");
    BindArray("GenPart_pt", GenPart_pt, &nGenPart, "nGenPart");
    BindArray("GenPart_status", GenPart_status, &nGenPart, "nGenPart");
    if(Run == 3) {
        BindArray("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &nGenPart, "nGenPart");
        BindArray("GenPart_statusFlags", GenPart_statusFlags, &nGenPart, "nGenPart");
        ZeroFillArray("GenPart_genPartIdxMother_RunII", GenPart_genPartIdxMother_RunII, &nGenPart, "nGenPart");
        ZeroFillArray("GenPart_statusFlags_RunII", GenPart_statusFlags_RunII, &nGenPart, "nGenPart");
    } else if(Run == 2) {
        BindArray("GenPart_genPartIdxMother", GenPart_genPartIdxMother_RunII, &nGenPart, "nGenPart");
        BindArray("GenPart_statusFlags", GenPart_statusFlags_RunII, &nGenPart, "nGenPart");
        ZeroFillArray("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &nGenPart, "nGenPart");
        ZeroFillArray("GenPart_statusFlags", GenPart_statusFlags, &nGenPart, "nGenPart");
    }

    // GenJet
    BindScalarWithRunCheck("nGenJet", nGenJet, nGenJet_RunII);
    BindArray("GenJet_eta", GenJet_eta, &nGenJet, "nGenJet");
    BindArray("GenJet_hadronFlavour", GenJet_hadronFlavour, &nGenJet, "nGenJet");
    BindArray("GenJet_mass", GenJet_mass, &nGenJet, "nGenJet");
    BindArray("GenJet_phi", GenJet_phi, &nGenJet, "nGenJet");
    BindArray("GenJet_pt", GenJet_pt, &nGenJet, "nGenJet");
    if(Run == 3){
        BindArray("GenJet_partonFlavour", GenJet_partonFlavour, &nGenJet, "nGenJet");
        ZeroFillArray("GenJet_partonFlavour_RunII", GenJet_partonFlavour_RunII, &nGenJet, "nGenJet");
    } else if(Run == 2) {
        BindArray("GenJet_partonFlavour", GenJet_partonFlavour_RunII, &nGenJet, "nGenJet");
        ZeroFillArray("GenJet_partonFlavour", GenJet_partonFlavour, &nGenJet, "nGenJet");
    }

    // GenJetAK8 (partonFlavour was never bound in the legacy mode -> zeros)
    BindScalarWithRunCheck("nGenJetAK8", nGenJetAK8, nGenJetAK8_RunII);
    BindArray("GenJetAK8_eta", GenJetAK8_eta, &nGenJetAK8, "nGenJetAK8");
    BindArray("GenJetAK8_hadronFlavour", GenJetAK8_hadronFlavour, &nGenJetAK8, "nGenJetAK8");
    BindArray("GenJetAK8_mass", GenJetAK8_mass, &nGenJetAK8, "nGenJetAK8");
    BindArray("GenJetAK8_phi", GenJetAK8_phi, &nGenJetAK8, "nGenJetAK8");
    BindArray("GenJetAK8_pt", GenJetAK8_pt, &nGenJetAK8, "nGenJetAK8");
    ZeroFillArray("GenJetAK8_partonFlavour", GenJetAK8_partonFlavour, &nGenJetAK8, "nGenJetAK8");
    ZeroFillArray("GenJetAK8_partonFlavour_RunII", GenJetAK8_partonFlavour_RunII, &nGenJetAK8, "nGenJetAK8");

    // GenMET
    BindScalar("GenMET_pt", GenMET_pt);
    BindScalar("GenMET_phi", GenMET_phi);

    // GenDressedLepton
    BindScalarWithRunCheck("nGenDressedLepton", nGenDressedLepton, nGenDressedLepton_RunII);
    BindArray("GenDressedLepton_pt", GenDressedLepton_pt, &nGenDressedLepton, "nGenDressedLepton");
    BindArray("GenDressedLepton_eta", GenDressedLepton_eta, &nGenDressedLepton, "nGenDressedLepton");
    BindArray("GenDressedLepton_phi", GenDressedLepton_phi, &nGenDressedLepton, "nGenDressedLepton");
    BindArray("GenDressedLepton_mass", GenDressedLepton_mass, &nGenDressedLepton, "nGenDressedLepton");
    BindArray("GenDressedLepton_pdgId", GenDressedLepton_pdgId, &nGenDressedLepton, "nGenDressedLepton");
    BindArray("GenDressedLepton_hasTauAnc", GenDressedLepton_hasTauAnc, &nGenDressedLepton, "nGenDressedLepton");

    // GenIsolatedPhoton
    BindScalarWithRunCheck("nGenIsolatedPhoton", nGenIsolatedPhoton, nGenIsolatedPhoton_RunII);
    BindArray("GenIsolatedPhoton_pt", GenIsolatedPhoton_pt, &nGenIsolatedPhoton, "nGenIsolatedPhoton");
    BindArray("GenIsolatedPhoton_eta", GenIsolatedPhoton_eta, &nGenIsolatedPhoton, "nGenIsolatedPhoton");
    BindArray("GenIsolatedPhoton_phi", GenIsolatedPhoton_phi, &nGenIsolatedPhoton, "nGenIsolatedPhoton");
    BindArray("GenIsolatedPhoton_mass", GenIsolatedPhoton_mass, &nGenIsolatedPhoton, "nGenIsolatedPhoton");

    // GenVisTau
    BindScalarWithRunCheck("nGenVisTau", nGenVisTau, nGenVisTau_RunII);
    BindArray("GenVisTau_pt", GenVisTau_pt, &nGenVisTau, "nGenVisTau");
    BindArray("GenVisTau_eta", GenVisTau_eta, &nGenVisTau, "nGenVisTau");
    BindArray("GenVisTau_phi", GenVisTau_phi, &nGenVisTau, "nGenVisTau");
    BindArray("GenVisTau_mass", GenVisTau_mass, &nGenVisTau, "nGenVisTau");
    BindArray("GenVisTau_charge", GenVisTau_charge, &nGenVisTau, "nGenVisTau");
    BindArray("GenVisTau_genPartIdxMother", GenVisTau_genPartIdxMother, &nGenVisTau, "nGenVisTau");
    BindArray("GenVisTau_status", GenVisTau_status, &nGenVisTau, "nGenVisTau");

    // GenVtx

    // PileUp & others
    BindScalar("Pileup_nPU", Pileup_nPU);
    BindScalar("Pileup_nTrueInt", Pileup_nTrueInt);
    BindScalar("genTtbarId", genTtbarId);

    // Muon----------------------------
    BindScalarWithRunCheck("nMuon", nMuon, nMuon_RunII);
    BindArray("Muon_charge", Muon_charge, &nMuon, "nMuon");
    BindArray("Muon_dxy", Muon_dxy, &nMuon, "nMuon");
    BindArray("Muon_dxyErr", Muon_dxyErr, &nMuon, "nMuon");
    BindArray("Muon_dxybs", Muon_dxybs, &nMuon, "nMuon");
    BindArray("Muon_dz", Muon_dz, &nMuon, "nMuon");
    BindArray("Muon_dzErr", Muon_dzErr, &nMuon, "nMuon");
    BindArray("Muon_eta", Muon_eta, &nMuon, "nMuon");
    BindArray("Muon_ip3d", Muon_ip3d, &nMuon, "nMuon");
    BindArray("Muon_nTrackerLayers", Muon_nTrackerLayers, &nMuon, "nMuon");
    BindArray("Muon_isGlobal", Muon_isGlobal, &nMuon, "nMuon");
    BindArray("Muon_highPtId", Muon_highPtId, &nMuon, "nMuon");
    BindArray("Muon_isStandalone", Muon_isStandalone, &nMuon, "nMuon");
    BindArray("Muon_isTracker", Muon_isTracker, &nMuon, "nMuon");
    BindArray("Muon_looseId", Muon_looseId, &nMuon, "nMuon");
    BindArray("Muon_mass", Muon_mass, &nMuon, "nMuon");
    BindArray("Muon_mediumId", Muon_mediumId, &nMuon, "nMuon");
    BindArray("Muon_mediumPromptId", Muon_mediumPromptId, &nMuon, "nMuon");
    BindArray("Muon_miniIsoId", Muon_miniIsoId, &nMuon, "nMuon");
    BindArray("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, &nMuon, "nMuon");
    BindArray("Muon_multiIsoId", Muon_multiIsoId, &nMuon, "nMuon");
    BindArray("Muon_mvaLowPt", Muon_mvaLowPt, &nMuon, "nMuon");
    BindArray("Muon_mvaTTH", Muon_mvaTTH, &nMuon, "nMuon");
    BindArray("Muon_pfIsoId", Muon_pfIsoId, &nMuon, "nMuon");
    BindArray("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &nMuon, "nMuon");
    BindArray("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &nMuon, "nMuon");
    BindArray("Muon_phi", Muon_phi, &nMuon, "nMuon");
    BindArray("Muon_pt", Muon_pt, &nMuon, "nMuon");
    BindArray("Muon_sip3d", Muon_sip3d, &nMuon, "nMuon");
    BindArray("Muon_softId", Muon_softId, &nMuon, "nMuon");
    BindArray("Muon_softMva", Muon_softMva, &nMuon, "nMuon");
    BindArray("Muon_softMvaId", Muon_softMvaId, &nMuon, "nMuon");
    BindArray("Muon_tightId", Muon_tightId, &nMuon, "nMuon");
    BindArray("Muon_tkIsoId", Muon_tkIsoId, &nMuon, "nMuon");
    BindArray("Muon_tkRelIso", Muon_tkRelIso, &nMuon, "nMuon");
    BindArray("Muon_triggerIdLoose", Muon_triggerIdLoose, &nMuon, "nMuon");
    BindArray("Muon_genPartFlav", Muon_genPartFlav, &nMuon, "nMuon");
    ZeroFillArray("Muon_puppiIsoId", Muon_puppiIsoId, &nMuon, "nMuon"); // never bound in legacy mode
    if (Run == 3) {
        BindArray("Muon_mvaMuID_WP", Muon_mvaMuID_WP, &nMuon, "nMuon");
        BindArray("Muon_jetIdx", Muon_jetIdx, &nMuon, "nMuon");
        BindArray("Muon_genPartIdx", Muon_genPartIdx, &nMuon, "nMuon");
    } else if(Run == 2) {
        BindArray("Muon_mvaId", Muon_mvaId, &nMuon, "nMuon");
        BindArray("Muon_jetIdx", Muon_jetIdx_RunII, &nMuon, "nMuon");
        BindArray("Muon_genPartIdx", Muon_genPartIdx_RunII, &nMuon, "nMuon");
    }

    //Electron----------------------------
    BindScalarWithRunCheck("nElectron", nElectron, nElectron_RunII);
    BindArray("Electron_charge", Electron_charge, &nElectron, "nElectron");
    BindArray("Electron_convVeto", Electron_convVeto, &nElectron, "nElectron");
    BindArray("Electron_cutBased_HEEP", Electron_cutBased_HEEP, &nElectron, "nElectron");
    BindArray("Electron_scEta", Electron_scEta, &nElectron, "nElectron");
    BindArray("Electron_deltaEtaInSC", Electron_deltaEtaInSC, &nElectron, "nElectron");
    BindArray("Electron_deltaEtaInSeed", Electron_deltaEtaInSeed, &nElectron, "nElectron");
    BindArray("Electron_deltaPhiInSC", Electron_deltaPhiInSC, &nElectron, "nElectron");
    BindArray("Electron_deltaPhiInSeed", Electron_deltaPhiInSeed, &nElectron, "nElectron");
    BindArray("Electron_ecalPFClusterIso", Electron_ecalPFClusterIso, &nElectron, "nElectron");
    BindArray("Electron_hcalPFClusterIso", Electron_hcalPFClusterIso, &nElectron, "nElectron");
    BindArray("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt, &nElectron, "nElectron");
    BindArray("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt, &nElectron, "nElectron");
    BindArray("Electron_dr03TkSumPt", Electron_dr03TkSumPt, &nElectron, "nElectron");
    BindArray("Electron_dr03TkSumPtHEEP", Electron_dr03TkSumPtHEEP, &nElectron, "nElectron");
    BindArray("Electron_dxy", Electron_dxy, &nElectron, "nElectron");
    BindArray("Electron_dxyErr", Electron_dxyErr, &nElectron, "nElectron");
    BindArray("Electron_dz", Electron_dz, &nElectron, "nElectron");
    BindArray("Electron_dzErr", Electron_dzErr, &nElectron, "nElectron");
    BindArray("Electron_eInvMinusPInv", Electron_eInvMinusPInv, &nElectron, "nElectron");
    BindArray("Electron_energyErr", Electron_energyErr, &nElectron, "nElectron");
    BindArray("Electron_eta", Electron_eta, &nElectron, "nElectron");
    BindArray("Electron_hoe", Electron_hoe, &nElectron, "nElectron");
    BindArray("Electron_ip3d", Electron_ip3d, &nElectron, "nElectron");
    BindArray("Electron_isPFcand", Electron_isPFcand, &nElectron, "nElectron");
    BindArray("Electron_jetNDauCharged", Electron_jetNDauCharged, &nElectron, "nElectron");
    BindArray("Electron_jetPtRelv2", Electron_jetPtRelv2, &nElectron, "nElectron");
    BindArray("Electron_jetRelIso", Electron_jetRelIso, &nElectron, "nElectron");
    BindArray("Electron_lostHits", Electron_lostHits, &nElectron, "nElectron");
    BindArray("Electron_mass", Electron_mass, &nElectron, "nElectron");
    BindArray("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, &nElectron, "nElectron");
    BindArray("Electron_miniPFRelIso_chg", Electron_miniPFRelIso_chg, &nElectron, "nElectron");
    BindArray("Electron_mvaTTH", Electron_mvaTTH, &nElectron, "nElectron");
    BindArray("Electron_pdgId", Electron_pdgId, &nElectron, "nElectron");
    BindArray("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &nElectron, "nElectron");
    BindArray("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, &nElectron, "nElectron");
    BindArray("Electron_phi", Electron_phi, &nElectron, "nElectron");
    BindArray("Electron_pt", Electron_pt, &nElectron, "nElectron");
    BindArray("Electron_r9", Electron_r9, &nElectron, "nElectron");
    BindArray("Electron_scEtOverPt", Electron_scEtOverPt, &nElectron, "nElectron");
    BindArray("Electron_seedGain", Electron_seedGain, &nElectron, "nElectron");
    BindArray("Electron_sieie", Electron_sieie, &nElectron, "nElectron");
    BindArray("Electron_sip3d", Electron_sip3d, &nElectron, "nElectron");
    BindArray("Electron_genPartFlav", Electron_genPartFlav, &nElectron, "nElectron");
    if (Run == 3) {
        BindArray("Electron_cutBased", Electron_cutBased, &nElectron, "nElectron");
        BindArray("Electron_genPartIdx", Electron_genPartIdx, &nElectron, "nElectron");
        BindArray("Electron_jetIdx", Electron_jetIdx, &nElectron, "nElectron");
        BindArray("Electron_mvaIso", Electron_mvaIso, &nElectron, "nElectron");
        BindArray("Electron_mvaIso_WP80", Electron_mvaIso_WP80, &nElectron, "nElectron");
        BindArray("Electron_mvaIso_WP90", Electron_mvaIso_WP90, &nElectron, "nElectron");
        BindArray("Electron_mvaNoIso", Electron_mvaNoIso, &nElectron, "nElectron");
        BindArray("Electron_mvaNoIso_WP80", Electron_mvaNoIso_WP80, &nElectron, "nElectron");
        BindArray("Electron_mvaNoIso_WP90", Electron_mvaNoIso_WP90, &nElectron, "nElectron");
        // WPL members were never bound in the legacy mode -> zeros
        ZeroFillArray("Electron_mvaIso_WPL", Electron_mvaIso_WPL, &nElectron, "nElectron");
        ZeroFillArray("Electron_mvaNoIso_WPL", Electron_mvaNoIso_WPL, &nElectron, "nElectron");
    } else if(Run == 2) {
        BindArray("Electron_cutBased", Electron_cutBased_RunII, &nElectron, "nElectron");
        BindArray("Electron_genPartIdx", Electron_genPartIdx_RunII, &nElectron, "nElectron");
        BindArray("Electron_jetIdx", Electron_jetIdx_RunII, &nElectron, "nElectron");
        BindArray("Electron_mvaFall17V2Iso", Electron_mvaFall17V2Iso, &nElectron, "nElectron");
        BindArray("Electron_mvaFall17V2Iso_WP80", Electron_mvaFall17V2Iso_WP80, &nElectron, "nElectron");
        BindArray("Electron_mvaFall17V2Iso_WP90", Electron_mvaFall17V2Iso_WP90, &nElectron, "nElectron");
        BindArray("Electron_mvaFall17V2Iso_WPL", Electron_mvaFall17V2Iso_WPL, &nElectron, "nElectron");
        BindArray("Electron_mvaFall17V2noIso", Electron_mvaFall17V2noIso, &nElectron, "nElectron");
        BindArray("Electron_mvaFall17V2noIso_WP80", Electron_mvaFall17V2noIso_WP80, &nElectron, "nElectron");
        BindArray("Electron_mvaFall17V2noIso_WP90", Electron_mvaFall17V2noIso_WP90, &nElectron, "nElectron");
        BindArray("Electron_mvaFall17V2noIso_WPLoose", Electron_mvaFall17V2noIso_WPL, &nElectron, "nElectron");
        BindArray("Electron_dEsigmaUp", Electron_dEsigmaUp, &nElectron, "nElectron");
        BindArray("Electron_dEsigmaDown", Electron_dEsigmaDown, &nElectron, "nElectron");
    }

    // Photon----------------------------
    BindScalarWithRunCheck("nPhoton", nPhoton, nPhoton_RunII);
    BindArray("Photon_eta", Photon_eta, &nPhoton, "nPhoton");
    BindArray("Photon_hoe", Photon_hoe, &nPhoton, "nPhoton");
    BindArray("Photon_isScEtaEB", Photon_isScEtaEB, &nPhoton, "nPhoton");
    BindArray("Photon_isScEtaEE", Photon_isScEtaEE, &nPhoton, "nPhoton");
    BindArray("Photon_mvaID", Photon_mvaID, &nPhoton, "nPhoton");
    BindArray("Photon_mvaID_WP80", Photon_mvaID_WP80, &nPhoton, "nPhoton");
    BindArray("Photon_mvaID_WP90", Photon_mvaID_WP90, &nPhoton, "nPhoton");
    BindArray("Photon_phi", Photon_phi, &nPhoton, "nPhoton");
    BindArray("Photon_pt", Photon_pt, &nPhoton, "nPhoton");
    BindArray("Photon_sieie", Photon_sieie, &nPhoton, "nPhoton");
    ZeroFillArray("Photon_energyErr", Photon_energyErr, &nPhoton, "nPhoton"); // never bound in legacy mode
    if (Run == 3) {
        BindArray("Photon_energyRaw", Photon_energyRaw, &nPhoton, "nPhoton");
        BindArray("Photon_cutBased", Photon_cutBased, &nPhoton, "nPhoton");
        ZeroFillArray("Photon_cutBased_RunII", Photon_cutBased_RunII, &nPhoton, "nPhoton");
    } else if(Run == 2) {
        BindArray("Photon_cutBased", Photon_cutBased_RunII, &nPhoton, "nPhoton");
    }

    //Jet----------------------------
    BindScalarWithRunCheck("nJet", nJet, nJet_RunII);
    BindArray("Jet_area", Jet_area, &nJet, "nJet");
    BindArray("Jet_btagDeepFlavB", Jet_btagDeepFlavB, &nJet, "nJet");
    BindArray("Jet_btagDeepFlavCvB", Jet_btagDeepFlavCvB, &nJet, "nJet");
    BindArray("Jet_btagDeepFlavCvL", Jet_btagDeepFlavCvL, &nJet, "nJet");
    BindArray("Jet_btagDeepFlavQG", Jet_btagDeepFlavQG, &nJet, "nJet");
    BindArray("Jet_chEmEF", Jet_chEmEF, &nJet, "nJet");
    BindArray("Jet_chHEF", Jet_chHEF, &nJet, "nJet");
    BindArray("Jet_eta", Jet_eta, &nJet, "nJet");
    BindArray("Jet_hfadjacentEtaStripsSize", Jet_hfadjacentEtaStripsSize, &nJet, "nJet");
    BindArray("Jet_hfcentralEtaStripSize", Jet_hfcentralEtaStripSize, &nJet, "nJet");
    BindArray("Jet_hfsigmaEtaEta", Jet_hfsigmaEtaEta, &nJet, "nJet");
    BindArray("Jet_hfsigmaPhiPhi", Jet_hfsigmaPhiPhi, &nJet, "nJet");
    BindArray("Jet_mass", Jet_mass, &nJet, "nJet");
    BindArray("Jet_muEF", Jet_muEF, &nJet, "nJet");
    BindArray("Jet_muonSubtrFactor", Jet_muonSubtrFactor, &nJet, "nJet");
    BindArray("Jet_nConstituents", Jet_nConstituents, &nJet, "nJet");
    BindArray("Jet_neEmEF", Jet_neEmEF, &nJet, "nJet");
    BindArray("Jet_neHEF", Jet_neHEF, &nJet, "nJet");
    BindArray("Jet_phi", Jet_phi, &nJet, "nJet");
    BindArray("Jet_pt", Jet_pt, &nJet, "nJet");
    BindArray("Jet_rawFactor", Jet_rawFactor, &nJet, "nJet");
    if (Run == 3) {
        BindArray("Jet_PNetRegPtRawCorr", Jet_PNetRegPtRawCorr, &nJet, "nJet");
        BindArray("Jet_PNetRegPtRawCorrNeutrino", Jet_PNetRegPtRawCorrNeutrino, &nJet, "nJet");
        BindArray("Jet_PNetRegPtRawRes", Jet_PNetRegPtRawRes, &nJet, "nJet");
        BindArray("Jet_btagPNetB", Jet_btagPNetB, &nJet, "nJet");
        BindArray("Jet_btagPNetCvB", Jet_btagPNetCvB, &nJet, "nJet");
        BindArray("Jet_btagPNetCvL", Jet_btagPNetCvL, &nJet, "nJet");
        BindArray("Jet_btagPNetQvG", Jet_btagPNetQvG, &nJet, "nJet");
        BindArray("Jet_btagPNetTauVJet", Jet_btagPNetTauVJet, &nJet, "nJet");
        BindArray("Jet_btagRobustParTAK4B", Jet_btagRobustParTAK4B, &nJet, "nJet");
        BindArray("Jet_btagRobustParTAK4CvB", Jet_btagRobustParTAK4CvB, &nJet, "nJet");
        BindArray("Jet_btagRobustParTAK4CvL", Jet_btagRobustParTAK4CvL, &nJet, "nJet");
        BindArray("Jet_btagRobustParTAK4QG", Jet_btagRobustParTAK4QG, &nJet, "nJet");
        BindArray("Jet_electronIdx1", Jet_electronIdx1, &nJet, "nJet");
        BindArray("Jet_electronIdx2", Jet_electronIdx2, &nJet, "nJet");
        BindArray("Jet_genJetIdx", Jet_genJetIdx, &nJet, "nJet");
        BindArray("Jet_hadronFlavour", Jet_hadronFlavour, &nJet, "nJet");
        BindArray("Jet_jetId", Jet_jetId, &nJet, "nJet");
        BindArray("Jet_muonIdx1", Jet_muonIdx1, &nJet, "nJet");
        BindArray("Jet_muonIdx2", Jet_muonIdx2, &nJet, "nJet");
        BindArray("Jet_nElectrons", Jet_nElectrons, &nJet, "nJet");
        BindArray("Jet_nMuons", Jet_nMuons, &nJet, "nJet");
        BindArray("Jet_nSVs", Jet_nSVs, &nJet, "nJet");
        BindArray("Jet_partonFlavour", Jet_partonFlavour, &nJet, "nJet");
        BindArray("Jet_svIdx1", Jet_svIdx1, &nJet, "nJet");
        BindArray("Jet_svIdx2", Jet_svIdx2, &nJet, "nJet");
        BindArray("Jet_chMultiplicity", Jet_chMultiplicity, &nJet, "nJet");
        BindArray("Jet_neMultiplicity", Jet_neMultiplicity, &nJet, "nJet");
    } else if (Run == 2) {
        BindArray("Jet_bRegCorr", Jet_bRegCorr, &nJet, "nJet");
        BindArray("Jet_bRegRes", Jet_bRegRes, &nJet, "nJet");
        BindArray("Jet_btagCSVV2", Jet_btagCSVV2, &nJet, "nJet");
        //BindArray("Jet_btagDeepB", Jet_btagDeepB, &nJet, "nJet");
        //BindArray("Jet_btagDeepCvB", Jet_btagDeepCvB, &nJet, "nJet");
        //BindArray("Jet_btagDeepCvL", Jet_btagDeepCvL, &nJet, "nJet");
        BindArray("Jet_cRegCorr", Jet_cRegCorr, &nJet, "nJet");
        BindArray("Jet_cRegRes", Jet_cRegRes, &nJet, "nJet");
        BindArray("Jet_chFPV0EF", Jet_chFPV0EF, &nJet, "nJet");
        BindArray("Jet_cleanmask", Jet_cleanmask, &nJet, "nJet");
        BindArray("Jet_electronIdx1", Jet_electronIdx1_RunII, &nJet, "nJet");
        BindArray("Jet_electronIdx2", Jet_electronIdx2_RunII, &nJet, "nJet");
        BindArray("Jet_genJetIdx", Jet_genJetIdx_RunII, &nJet, "nJet");
        BindArray("Jet_hadronFlavour", Jet_hadronFlavour_RunII, &nJet, "nJet");
        BindArray("Jet_jetId", Jet_jetId_RunII, &nJet, "nJet");
        BindArray("Jet_muonIdx1", Jet_muonIdx1_RunII, &nJet, "nJet");
        BindArray("Jet_muonIdx2", Jet_muonIdx2_RunII, &nJet, "nJet");
        BindArray("Jet_nElectrons", Jet_nElectrons_RunII, &nJet, "nJet");
        BindArray("Jet_nMuons", Jet_nMuons_RunII, &nJet, "nJet");
        BindArray("Jet_partonFlavour", Jet_partonFlavour_RunII, &nJet, "nJet");
        BindArray("Jet_puId", Jet_puId, &nJet, "nJet");
        BindArray("Jet_puIdDisc", Jet_puIdDisc, &nJet, "nJet");
        BindArray("Jet_qgl", Jet_qgl, &nJet, "nJet");
    }

    //Tau----------------------------
    BindScalarWithRunCheck("nTau", nTau, nTau_RunII);
    BindArray("Tau_dxy", Tau_dxy, &nTau, "nTau");
    BindArray("Tau_dz", Tau_dz, &nTau, "nTau");
    BindArray("Tau_eta", Tau_eta, &nTau, "nTau");
    BindArray("Tau_genPartFlav", Tau_genPartFlav, &nTau, "nTau");
    // These branch names do not exist (same as in the legacy mode); the bind
    // just reproduces the legacy warning. No zero-fill fallback (count=nullptr)
    // because the target members were kept at size 0 in Run2.
    BindArray("Tau_genPartidDeepTau2017v2p1VSe", Tau_idDeepTau2018v2p5VSe, nullptr, "");
    BindArray("Tau_genPartidDeepTau2017v2p1VSjet", Tau_idDeepTau2018v2p5VSjet, nullptr, "");
    BindArray("Tau_genPartidDeepTau2017v2p1VSmu", Tau_idDeepTau2018v2p5VSmu, nullptr, "");
    BindArray("Tau_mass", Tau_mass, &nTau, "nTau");
    BindArray("Tau_phi", Tau_phi, &nTau, "nTau");
    BindArray("Tau_pt", Tau_pt, &nTau, "nTau");
    // The 2017v2p1 members were never bound in the legacy mode -> zeros
    ZeroFillArray("Tau_idDeepTau2017v2p1VSe", Tau_idDeepTau2017v2p1VSe, &nTau, "nTau");
    ZeroFillArray("Tau_idDeepTau2017v2p1VSjet", Tau_idDeepTau2017v2p1VSjet, &nTau, "nTau");
    ZeroFillArray("Tau_idDeepTau2017v2p1VSmu", Tau_idDeepTau2017v2p1VSmu, &nTau, "nTau");
    if (Run == 3) {
        BindArray("Tau_charge", Tau_charge, &nTau, "nTau");
        BindArray("Tau_decayMode", Tau_decayMode, &nTau, "nTau");
        BindArray("Tau_genPartIdx", Tau_genPartIdx, &nTau, "nTau");
        BindArray("Tau_idDecayModeNewDMs", Tau_idDecayModeNewDMs, &nTau, "nTau");
        BindArray("Tau_idDeepTau2018v2p5VSe", Tau_idDeepTau2018v2p5VSe, &nTau, "nTau");
        BindArray("Tau_idDeepTau2018v2p5VSjet", Tau_idDeepTau2018v2p5VSjet, &nTau, "nTau");
        BindArray("Tau_idDeepTau2018v2p5VSmu", Tau_idDeepTau2018v2p5VSmu, &nTau, "nTau");
    } else if(Run == 2) {
        BindArray("Tau_charge", Tau_charge_RunII, &nTau, "nTau");
        BindArray("Tau_decayMode", Tau_decayMode_RunII, &nTau, "nTau");
        BindArray("Tau_genPartIdx", Tau_genPartIdx_RunII, &nTau, "nTau");
    }

    //FatJet----------------------------
    BindScalarWithRunCheck("nFatJet", nFatJet, nFatJet_RunII);
    BindArray("FatJet_area", FatJet_area, &nFatJet, "nFatJet");
    BindArray("FatJet_btagDDBvLV2", FatJet_btagDDBvLV2, &nFatJet, "nFatJet");
    BindArray("FatJet_btagDDCvBV2", FatJet_btagDDCvBV2, &nFatJet, "nFatJet");
    BindArray("FatJet_btagDDCvLV2", FatJet_btagDDCvLV2, &nFatJet, "nFatJet");
    BindArray("FatJet_btagDeepB", FatJet_btagDeepB, &nFatJet, "nFatJet");
    BindArray("FatJet_btagHbb", FatJet_btagHbb, &nFatJet, "nFatJet");
    BindArray("FatJet_eta", FatJet_eta, &nFatJet, "nFatJet");
    BindArray("FatJet_lsf3", FatJet_lsf3, &nFatJet, "nFatJet");
    BindArray("FatJet_mass", FatJet_mass, &nFatJet, "nFatJet");
    BindArray("FatJet_msoftdrop", FatJet_msoftdrop, &nFatJet, "nFatJet");
    BindArray("FatJet_nBHadrons", FatJet_nBHadrons, &nFatJet, "nFatJet");
    BindArray("FatJet_nCHadrons", FatJet_nCHadrons, &nFatJet, "nFatJet");
    BindArray("FatJet_nConstituents", FatJet_nConstituents, &nFatJet, "nFatJet");
    BindArray("FatJet_particleNet_QCD", FatJet_particleNet_QCD, &nFatJet, "nFatJet");
    BindArray("FatJet_phi", FatJet_phi, &nFatJet, "nFatJet");
    BindArray("FatJet_pt", FatJet_pt, &nFatJet, "nFatJet");
    BindArray("FatJet_tau1", FatJet_tau1, &nFatJet, "nFatJet");
    BindArray("FatJet_tau2", FatJet_tau2, &nFatJet, "nFatJet");
    BindArray("FatJet_tau3", FatJet_tau3, &nFatJet, "nFatJet");
    BindArray("FatJet_tau4", FatJet_tau4, &nFatJet, "nFatJet");
    if (Run == 3) {
        BindArray("FatJet_genJetAK8Idx", FatJet_genJetAK8Idx, &nFatJet, "nFatJet");
        BindArray("FatJet_jetId", FatJet_jetId, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNetWithMass_H4qvsQCD", FatJet_particleNetWithMass_H4qvsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNetWithMass_HbbvsQCD", FatJet_particleNetWithMass_HbbvsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNetWithMass_HccvsQCD", FatJet_particleNetWithMass_HccvsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNetWithMass_QCD", FatJet_particleNetWithMass_QCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNetWithMass_TvsQCD", FatJet_particleNetWithMass_TvsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNetWithMass_WvsQCD", FatJet_particleNetWithMass_WvsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNetWithMass_ZvsQCD", FatJet_particleNetWithMass_ZvsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_QCD0HF", FatJet_particleNet_QCD0HF, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_QCD1HF", FatJet_particleNet_QCD1HF, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_QCD2HF", FatJet_particleNet_QCD2HF, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_XbbVsQCD", FatJet_particleNet_XbbVsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_XccVsQCD", FatJet_particleNet_XccVsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_XggVsQCD", FatJet_particleNet_XggVsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_XqqVsQCD", FatJet_particleNet_XqqVsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_XteVsQCD", FatJet_particleNet_XteVsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_XtmVsQCD", FatJet_particleNet_XtmVsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_XttVsQCD", FatJet_particleNet_XttVsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_massCorr", FatJet_particleNet_massCorr, &nFatJet, "nFatJet");
        BindArray("FatJet_subJetIdx1", FatJet_subJetIdx1, &nFatJet, "nFatJet");
        BindArray("FatJet_subJetIdx2", FatJet_subJetIdx2, &nFatJet, "nFatJet");
        BindArray("FatJet_n2b1", FatJet_n2b1, &nFatJet, "nFatJet");
        BindArray("FatJet_n3b1", FatJet_n3b1, &nFatJet, "nFatJet");
        BindArray("FatJet_chMultiplicity", FatJet_chMultiplicity, &nFatJet, "nFatJet");
        BindArray("FatJet_neMultiplicity", FatJet_neMultiplicity, &nFatJet, "nFatJet");
        BindArray("FatJet_chHEF", FatJet_chHEF, &nFatJet, "nFatJet");
        BindArray("FatJet_neHEF", FatJet_neHEF, &nFatJet, "nFatJet");
        BindArray("FatJet_chEmEF", FatJet_chEmEF, &nFatJet, "nFatJet");
        BindArray("FatJet_neEmEF", FatJet_neEmEF, &nFatJet, "nFatJet");
        BindArray("FatJet_muEF", FatJet_muEF, &nFatJet, "nFatJet");
        BindScalarWithRunCheck("nSubJet", nSubJet, nSubJet_RunII);
        BindArray("SubJet_pt", SubJet_pt, &nSubJet, "nSubJet");
        BindArray("SubJet_eta", SubJet_eta, &nSubJet, "nSubJet");
        BindArray("SubJet_phi", SubJet_phi, &nSubJet, "nSubJet");
        BindArray("SubJet_mass", SubJet_mass, &nSubJet, "nSubJet");
        BindArray("SubJet_btagDeepB", SubJet_btagDeepB, &nSubJet, "nSubJet");
        BindScalarWithRunCheck("nSV", nSV, nSV_RunII);
        BindArray("SV_pt", SV_pt, &nSV, "nSV");
        BindArray("SV_eta", SV_eta, &nSV, "nSV");
        BindArray("SV_phi", SV_phi, &nSV, "nSV");
        BindArray("SV_mass", SV_mass, &nSV, "nSV");
        BindArray("SV_dlenSig", SV_dlenSig, &nSV, "nSV");
        BindArray("SV_dxySig", SV_dxySig, &nSV, "nSV");
        BindArray("SV_chi2", SV_chi2, &nSV, "nSV");
        BindArray("SV_pAngle", SV_pAngle, &nSV, "nSV");
        BindArray("SV_ntracks", SV_ntracks, &nSV, "nSV");
    } else if(Run == 2) {
        BindArray("FatJet_genJetAK8Idx", FatJet_genJetAK8Idx_RunII, &nFatJet, "nFatJet");
        BindArray("FatJet_jetId", FatJet_jetId_RunII, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNetMD_QCD", FatJet_particleNetMD_QCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNetMD_Xbb", FatJet_particleNetMD_Xbb, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNetMD_Xcc", FatJet_particleNetMD_Xcc, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNetMD_Xqq", FatJet_particleNetMD_Xqq, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_H4qvsQCD", FatJet_particleNet_H4qvsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_HbbvsQCD", FatJet_particleNet_HbbvsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_HccvsQCD", FatJet_particleNet_HccvsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_TvsQCD", FatJet_particleNet_TvsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_WvsQCD", FatJet_particleNet_WvsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_ZvsQCD", FatJet_particleNet_ZvsQCD, &nFatJet, "nFatJet");
        BindArray("FatJet_particleNet_mass", FatJet_particleNet_mass, &nFatJet, "nFatJet");
        BindArray("FatJet_subJetIdx1", FatJet_subJetIdx1_RunII, &nFatJet, "nFatJet");
        BindArray("FatJet_subJetIdx2", FatJet_subJetIdx2_RunII, &nFatJet, "nFatJet");
    }

    // MET----------------------------
    BindScalar("MET_pt", MET_pt);
    BindScalar("MET_phi", MET_phi);
    BindScalar("PuppiMET_pt", PuppiMET_pt);
    BindScalar("PuppiMET_phi", PuppiMET_phi);
    BindScalar("PuppiMET_ptUnclusteredUp", PuppiMET_ptUnclusteredUp);
    BindScalar("PuppiMET_phiUnclusteredUp", PuppiMET_phiUnclusteredUp);
    BindScalar("PuppiMET_ptUnclusteredDown", PuppiMET_ptUnclusteredDown);
    BindScalar("PuppiMET_phiUnclusteredDown", PuppiMET_phiUnclusteredDown);

    //Rho----------------------------
    if(Run == 3) {
        BindScalar("Rho_fixedGridRhoFastjetAll", fixedGridRhoFastjetAll);
    } else if(Run == 2) {
        BindScalar("fixedGridRhoFastjetAll", fixedGridRhoFastjetAll);
    }

    // PV----------------------------
    BindScalar("PV_chi2", PV_chi2);
    BindScalar("PV_ndof", PV_ndof);
    BindScalar("PV_score", PV_score);
    BindScalar("PV_x", PV_x);
    BindScalar("PV_y", PV_y);
    BindScalar("PV_z", PV_z);
    if (Run==3) {
        BindScalar("PV_npvs", PV_npvs);
        BindScalar("PV_npvsGood", PV_npvsGood);
    } else if(Run==2) {
        BindScalar("PV_npvs", PV_npvs_RunII);
        BindScalar("PV_npvsGood", PV_npvsGood_RunII);
    }

    //L1PreFireweight----------------------------
    BindScalar("L1PreFiringWeight_Nom", L1PreFiringWeight_Nom);
    BindScalar("L1PreFiringWeight_Dn", L1PreFiringWeight_Dn);
    BindScalar("L1PreFiringWeight_Up", L1PreFiringWeight_Up);

    //Flags----------------------------
    BindScalar("Flag_METFilters", Flag_METFilters);
    BindScalar("Flag_goodVertices", Flag_goodVertices);
    BindScalar("Flag_globalSuperTightHalo2016Filter", Flag_globalSuperTightHalo2016Filter);
    BindScalar("Flag_HBHENoiseFilter", Flag_HBHENoiseFilter);
    BindScalar("Flag_HBHENoiseIsoFilter", Flag_HBHENoiseIsoFilter);
    BindScalar("Flag_EcalDeadCellTriggerPrimitiveFilter", Flag_EcalDeadCellTriggerPrimitiveFilter);
    BindScalar("Flag_BadPFMuonFilter", Flag_BadPFMuonFilter);
    BindScalar("Flag_BadPFMuonDzFilter", Flag_BadPFMuonDzFilter);
    BindScalar("Flag_hfNoisyHitsFilter", Flag_hfNoisyHitsFilter);
    BindScalar("Flag_ecalBadCalibFilter", Flag_ecalBadCalibFilter);
    BindScalar("Flag_eeBadScFilter", Flag_eeBadScFilter);
    BindScalar("run", RunNumber);
    BindScalar("luminosityBlock", LumiBlock);
    BindScalar("event", EventNumber);

    // TrigObj----------------------------
    BindScalarWithRunCheck("nTrigObj", nTrigObj, nTrigObj_RunII);
    BindArray("TrigObj_pt", TrigObj_pt, &nTrigObj, "nTrigObj");
    BindArray("TrigObj_eta", TrigObj_eta, &nTrigObj, "nTrigObj");
    BindArray("TrigObj_phi", TrigObj_phi, &nTrigObj, "nTrigObj");
    if (Run == 3) {
        BindArray("TrigObj_id", TrigObj_id, &nTrigObj, "nTrigObj");
        ZeroFillArray("TrigObj_id_RunII", TrigObj_id_RunII, &nTrigObj, "nTrigObj");
    } else {
        BindArray("TrigObj_id", TrigObj_id_RunII, &nTrigObj, "nTrigObj");
        ZeroFillArray("TrigObj_id", TrigObj_id, &nTrigObj, "nTrigObj");
    }
    BindArray("TrigObj_filterBits", TrigObj_filterBits, &nTrigObj, "nTrigObj");

    // For some data files, the branch is not in all files, especially for
    // triggers; a trigger read through TTreeReader must exist in every file
    // of the chain (same check as the legacy SuperSafeSetBranchAddress).
    auto BranchInAllFiles = [this](const TString &branchName) {
        TObjArray* fileElements = fChain->GetListOfFiles();
        for (int i = 0; i < fileElements->GetEntries(); i++) {
            TChainElement* element = (TChainElement*)fileElements->At(i);
            TString fileName = element->GetTitle();
            TFile* file = TFile::Open(fileName);
            TTree* tree = (TTree*)file->Get(fChain->GetName());
            if (!tree->GetBranch(branchName)) {
                cout << "[SKNanoGenLoader::Init] Warning: Branch " << branchName << " not found in file " << fileName << endl;
                file->Close();
                return false;
            }
            file->Close();
        }
        return true;
    };

    string json_path = string(getenv("SKNANO_DATA")) + "/" + DataEra.Data() + "/Trigger/HLT_Path.json";
    ifstream json_file(json_path);
    if (json_file.is_open()) {
        cout << "[SKNanoLoader::Init] Loading HLT Paths in " << json_path << endl;
        json j;
        json_file >> j;
        RVec<TString> not_in_tree;
        for (auto& [key, value] : j.items()) {
            if (!value.contains("active")) continue;
            if (!value["active"]) continue;
            Bool_t* passHLT = new Bool_t();
            TString key_str = key;
            TriggerMap[key_str].first = passHLT;
            TriggerMap[key_str].second = value["lumi"];
            //if key_str is in tree, set up a reader for it
            if (fChain->GetBranch(key_str)) {
                // In some data file, part of the trigger set is missing (changed during the run?)
                if (!IsDATA || BranchInAllFiles(key_str)) {
                    auto rdr = std::make_shared<TTreeReaderValue<Bool_t>>(*fReader, key_str.Data());
                    fScalarFillers.push_back([rdr, passHLT]() { *passHLT = **rdr; });
                }
                // if missing in some data file: stays false, as in the legacy mode
            } else if(key_str=="Full") {
                *TriggerMap[key_str].first = true;
            } else{
                not_in_tree.push_back(key_str);
                TriggerMap.erase(key_str);
            }
        }
        if (not_in_tree.size() > 0) {
            cout << "\033[1;33m[SKNanoLoader::Init] Following HLT Paths are not in the tree\033[0m" << endl;
            for (auto &path : not_in_tree) {
                cout << "\033[1;33m" << path << "\033[0m" << endl;
            }
        }
    }
    else cerr << "[SKNanoLoader::Init] Cannot open " << json_path << endl;
}

void SKNanoLoader::InitLegacy() {
    // Helper function to safely set branch address
    auto SafeSetBranchAddress = [this](const TString &branchName, void* address) {
        TBranch* branch = fChain->GetBranch(branchName);
        if (!branch) {
            cout << "[SKNanoGenLoader::Init] Warning:Branch " << branchName << " not found" << endl;
            return;
        }
        fChain->SetBranchStatus(branchName, 1);
        fChain->SetBranchAddress(branchName, address);
    };
    // For some data files, the branch is not in all files, especially for triggers
    auto SuperSafeSetBranchAddress = [this](const TString &branchName, void* address) {
        TObjArray* fileElements = fChain->GetListOfFiles();
        for (int i = 0; i < fileElements->GetEntries(); i++) {
            TChainElement* element = (TChainElement*)fileElements->At(i);
            TString fileName = element->GetTitle();
            TFile* file = TFile::Open(fileName);
            TTree* tree = (TTree*)file->Get("Events");
            if (!tree->GetBranch(branchName)) {
                cout << "[SKNanoGenLoader::Init] Warning: Branch " << branchName << " not found in file " << fileName << endl;
                file->Close();
                return;
            }
            file->Close();
        }
        fChain->SetBranchStatus(branchName, 1);
        fChain->SetBranchAddress(branchName, address);
    };
    // For type conversion between Run2 and Run3
    auto SetBranchWithRunCheck = [this, &SafeSetBranchAddress](const TString &branchName, Int_t &run3Var, UInt_t &runIIVar) {
        if (Run == 3) {
            SafeSetBranchAddress(branchName, &run3Var);
        } else {
            SafeSetBranchAddress(branchName, &runIIVar);
        }
    };

    cout << "[SKNanoLoader::Init] Using legacy branch-address reading (SkimmingMode)" << endl;

    SetMaxLeafSize();
    fChain->SetBranchStatus("*", 0);

    // Weights
    SafeSetBranchAddress("genWeight", &genWeight);
    SafeSetBranchAddress("LHEWeight_originalXWGTUP", &LHEWeight_originalXWGTUP);
    SafeSetBranchAddress("Generator_weight", &Generator_weight);
    SafeSetBranchAddress("nLHEPdfWeight", &nLHEPdfWeight);
    SafeSetBranchAddress("nLHEScaleWeight", &nLHEScaleWeight);
    SafeSetBranchAddress("nPSWeight", &nPSWeight);
    SafeSetBranchAddress("LHEPdfWeight", LHEPdfWeight.data());
    SafeSetBranchAddress("LHEScaleWeight", LHEScaleWeight.data());
    SafeSetBranchAddress("PSWeight", PSWeight.data());

    // PDFs
    SafeSetBranchAddress("Generator_id1", &Generator_id1);
    SafeSetBranchAddress("Generator_id2", &Generator_id2);
    SafeSetBranchAddress("Generator_x1", &Generator_x1);
    SafeSetBranchAddress("Generator_x2", &Generator_x2);
    SafeSetBranchAddress("Generator_xpdf1", &Generator_xpdf1);
    SafeSetBranchAddress("Generator_xpdf2", &Generator_xpdf2);
    SafeSetBranchAddress("Generator_scalePDF", &Generator_scalePDF);

    // LHE
    SafeSetBranchAddress("LHE_HT", &LHE_HT);
    SafeSetBranchAddress("LHE_HTIncoming", &LHE_HTIncoming);
    SafeSetBranchAddress("LHE_Vpt", &LHE_Vpt);
    SafeSetBranchAddress("LHE_AlphaS", &LHE_AlphaS);
    SafeSetBranchAddress("LHE_Njets", &LHE_Njets);
    SafeSetBranchAddress("LHE_Nb", &LHE_Nb);
    SafeSetBranchAddress("LHE_Nc", &LHE_Nc);
    SafeSetBranchAddress("LHE_Nuds", &LHE_Nuds);
    SafeSetBranchAddress("LHE_Nglu", &LHE_Nglu);
    SafeSetBranchAddress("LHE_NpLO", &LHE_NpLO);
    SafeSetBranchAddress("LHE_NpNLO", &LHE_NpNLO);

    // LHEPart
    SetBranchWithRunCheck("nLHEPart", nLHEPart, nLHEPart_RunII); 
    SafeSetBranchAddress("LHEPart_pt", LHEPart_pt.data());
    SafeSetBranchAddress("LHEPart_eta", LHEPart_eta.data());
    SafeSetBranchAddress("LHEPart_phi", LHEPart_phi.data());
    SafeSetBranchAddress("LHEPart_mass", LHEPart_mass.data());
    SafeSetBranchAddress("LHEPart_pdgId", LHEPart_pdgId.data());
    SafeSetBranchAddress("LHEPart_status", LHEPart_status.data());
    SafeSetBranchAddress("LHEPart_spin", LHEPart_spin.data());
    SafeSetBranchAddress("LHEPart_incomingpz", LHEPart_incomingpz.data());

    // GenPart
    SetBranchWithRunCheck("nGenPart", nGenPart, nGenPart_RunII);
    SafeSetBranchAddress("GenPart_eta", GenPart_eta.data());
    SafeSetBranchAddress("GenPart_mass", GenPart_mass.data());
    SafeSetBranchAddress("GenPart_pdgId", GenPart_pdgId.data());
    SafeSetBranchAddress("GenPart_phi", GenPart_phi.data());
    SafeSetBranchAddress("GenPart_pt", GenPart_pt.data());
    SafeSetBranchAddress("GenPart_status", GenPart_status.data());
    if(Run == 3) {
        SafeSetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother.data());
        SafeSetBranchAddress("GenPart_statusFlags", GenPart_statusFlags.data());
    } else if(Run == 2) {
        SafeSetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother_RunII.data());
        SafeSetBranchAddress("GenPart_statusFlags", GenPart_statusFlags_RunII.data());
    }

    // GenJet
    SetBranchWithRunCheck("nGenJet", nGenJet, nGenJet_RunII);
    SafeSetBranchAddress("GenJet_eta", GenJet_eta.data());
    SafeSetBranchAddress("GenJet_hadronFlavour", GenJet_hadronFlavour.data());
    SafeSetBranchAddress("GenJet_mass", GenJet_mass.data());
    SafeSetBranchAddress("GenJet_phi", GenJet_phi.data());
    SafeSetBranchAddress("GenJet_pt", GenJet_pt.data());
    if(Run == 3){
        SafeSetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour.data());
    } else if(Run == 2) {
        SafeSetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour_RunII.data());
    }

    // GenJetAK8
    SetBranchWithRunCheck("nGenJetAK8", nGenJetAK8, nGenJetAK8_RunII);
    SafeSetBranchAddress("GenJetAK8_eta", GenJetAK8_eta.data());
    SafeSetBranchAddress("GenJetAK8_hadronFlavour", GenJetAK8_hadronFlavour.data());
    SafeSetBranchAddress("GenJetAK8_mass", GenJetAK8_mass.data());
    SafeSetBranchAddress("GenJetAK8_phi", GenJetAK8_phi.data());
    SafeSetBranchAddress("GenJetAK8_pt", GenJetAK8_pt.data());

    // GenMET
    SafeSetBranchAddress("GenMET_pt", &GenMET_pt);
    SafeSetBranchAddress("GenMET_phi", &GenMET_phi);

    // GenDressedLepton
    SetBranchWithRunCheck("nGenDressedLepton", nGenDressedLepton, nGenDressedLepton_RunII);
    SafeSetBranchAddress("GenDressedLepton_pt", GenDressedLepton_pt.data());
    SafeSetBranchAddress("GenDressedLepton_eta", GenDressedLepton_eta.data());
    SafeSetBranchAddress("GenDressedLepton_phi", GenDressedLepton_phi.data());
    SafeSetBranchAddress("GenDressedLepton_mass", GenDressedLepton_mass.data());
    SafeSetBranchAddress("GenDressedLepton_pdgId", GenDressedLepton_pdgId.data());
    SafeSetBranchAddress("GenDressedLepton_hasTauAnc", GenDressedLepton_hasTauAnc.data());
    
    // GenIsolatedPhoton
    SetBranchWithRunCheck("nGenIsolatedPhoton", nGenIsolatedPhoton, nGenIsolatedPhoton_RunII);
    SafeSetBranchAddress("GenIsolatedPhoton_pt", GenIsolatedPhoton_pt.data());
    SafeSetBranchAddress("GenIsolatedPhoton_eta", GenIsolatedPhoton_eta.data());
    SafeSetBranchAddress("GenIsolatedPhoton_phi", GenIsolatedPhoton_phi.data());
    SafeSetBranchAddress("GenIsolatedPhoton_mass", GenIsolatedPhoton_mass.data());

    // GenVisTau
    SetBranchWithRunCheck("nGenVisTau", nGenVisTau, nGenVisTau_RunII);
    SafeSetBranchAddress("GenVisTau_pt", GenVisTau_pt.data());
    SafeSetBranchAddress("GenVisTau_eta", GenVisTau_eta.data());
    SafeSetBranchAddress("GenVisTau_phi", GenVisTau_phi.data());
    SafeSetBranchAddress("GenVisTau_mass", GenVisTau_mass.data());
    if (Run == 3) {
        // v13 leaf types: charge/genPartIdxMother are Short_t, status is UChar_t
        SafeSetBranchAddress("GenVisTau_charge", Buf_GenVisTau_charge.data());
        SafeSetBranchAddress("GenVisTau_genPartIdxMother", Buf_GenVisTau_genPartIdxMother.data());
        SafeSetBranchAddress("GenVisTau_status", Buf_GenVisTau_status.data());
    } else {
        SafeSetBranchAddress("GenVisTau_charge", GenVisTau_charge.data());
        SafeSetBranchAddress("GenVisTau_genPartIdxMother", GenVisTau_genPartIdxMother.data());
        SafeSetBranchAddress("GenVisTau_status", GenVisTau_status.data());
    }

    // GenVtx

    // PileUp & others
    SafeSetBranchAddress("Pileup_nPU", &Pileup_nPU);
    SafeSetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt); 
    SafeSetBranchAddress("FatJet_nBHadrons",FatJet_nBHadrons.data());
    SafeSetBranchAddress("FatJet_nCHadrons",FatJet_nCHadrons.data());
    SafeSetBranchAddress("genTtbarId", &genTtbarId);

    // Muon----------------------------
    SetBranchWithRunCheck("nMuon", nMuon, nMuon_RunII);
    SafeSetBranchAddress("Muon_charge", Muon_charge.data());
    SafeSetBranchAddress("Muon_dxy", Muon_dxy.data());
    SafeSetBranchAddress("Muon_dxyErr", Muon_dxyErr.data());
    SafeSetBranchAddress("Muon_dxybs", Muon_dxybs.data());
    SafeSetBranchAddress("Muon_dz", Muon_dz.data());
    SafeSetBranchAddress("Muon_dzErr", Muon_dzErr.data());
    SafeSetBranchAddress("Muon_eta", Muon_eta.data());
    SafeSetBranchAddress("Muon_ip3d", Muon_ip3d.data());
    if (Run == 3) // v13 leaf type is UChar_t
        SafeSetBranchAddress("Muon_nTrackerLayers", Buf_Muon_nTrackerLayers.data());
    else
        SafeSetBranchAddress("Muon_nTrackerLayers", Muon_nTrackerLayers.data());
    SafeSetBranchAddress("Muon_isGlobal", Muon_isGlobal.data());
    SafeSetBranchAddress("Muon_highPtId", Muon_highPtId.data());
    SafeSetBranchAddress("Muon_isStandalone", Muon_isStandalone.data());
    SafeSetBranchAddress("Muon_isTracker", Muon_isTracker.data());
    SafeSetBranchAddress("Muon_looseId", Muon_looseId.data());
    SafeSetBranchAddress("Muon_mass", Muon_mass.data());
    SafeSetBranchAddress("Muon_mediumId", Muon_mediumId.data());
    SafeSetBranchAddress("Muon_mediumPromptId", Muon_mediumPromptId.data());
    SafeSetBranchAddress("Muon_miniIsoId", Muon_miniIsoId.data());
    SafeSetBranchAddress("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all.data());
    SafeSetBranchAddress("Muon_multiIsoId", Muon_multiIsoId.data());
    SafeSetBranchAddress("Muon_mvaLowPt", Muon_mvaLowPt.data());
    SafeSetBranchAddress("Muon_mvaTTH", Muon_mvaTTH.data());
    SafeSetBranchAddress("Muon_pfIsoId", Muon_pfIsoId.data());
    SafeSetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all.data());
    SafeSetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all.data());
    SafeSetBranchAddress("Muon_phi", Muon_phi.data());
    SafeSetBranchAddress("Muon_pt", Muon_pt.data());
    SafeSetBranchAddress("Muon_sip3d", Muon_sip3d.data());
    SafeSetBranchAddress("Muon_softId", Muon_softId.data());
    SafeSetBranchAddress("Muon_softMva", Muon_softMva.data());
    SafeSetBranchAddress("Muon_softMvaId", Muon_softMvaId.data());
    SafeSetBranchAddress("Muon_tightId", Muon_tightId.data());
    SafeSetBranchAddress("Muon_tkIsoId", Muon_tkIsoId.data());
    SafeSetBranchAddress("Muon_tkRelIso", Muon_tkRelIso.data());
    SafeSetBranchAddress("Muon_triggerIdLoose", Muon_triggerIdLoose.data());
    SafeSetBranchAddress("Muon_genPartFlav", Muon_genPartFlav.data());
    if (Run == 3) {
        SafeSetBranchAddress("Muon_mvaMuID_WP", Muon_mvaMuID_WP.data());
        SafeSetBranchAddress("Muon_jetIdx", Muon_jetIdx.data());
        SafeSetBranchAddress("Muon_genPartIdx", Muon_genPartIdx.data());
    } else if(Run == 2) {
        SafeSetBranchAddress("Muon_mvaId", Muon_mvaId.data());
        SafeSetBranchAddress("Muon_jetIdx", Muon_jetIdx_RunII.data());
        SafeSetBranchAddress("Muon_genPartIdx", Muon_genPartIdx_RunII.data());
    }

    //Electron----------------------------
    SetBranchWithRunCheck("nElectron", nElectron, nElectron_RunII);
    SafeSetBranchAddress("Electron_charge", Electron_charge.data());
    SafeSetBranchAddress("Electron_convVeto", Electron_convVeto.data());
    SafeSetBranchAddress("Electron_cutBased_HEEP", Electron_cutBased_HEEP.data());
    SafeSetBranchAddress("Electron_scEta", Electron_scEta.data());
    SafeSetBranchAddress("Electron_deltaEtaInSC", Electron_deltaEtaInSC.data());
    SafeSetBranchAddress("Electron_deltaEtaInSeed", Electron_deltaEtaInSeed.data());
    SafeSetBranchAddress("Electron_deltaPhiInSC", Electron_deltaPhiInSC.data());
    SafeSetBranchAddress("Electron_deltaPhiInSeed", Electron_deltaPhiInSeed.data());
    SafeSetBranchAddress("Electron_ecalPFClusterIso", Electron_ecalPFClusterIso.data());
    SafeSetBranchAddress("Electron_hcalPFClusterIso", Electron_hcalPFClusterIso.data());
    SafeSetBranchAddress("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt.data());
    SafeSetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt.data());
    SafeSetBranchAddress("Electron_dr03TkSumPt", Electron_dr03TkSumPt.data());
    SafeSetBranchAddress("Electron_dr03TkSumPtHEEP", Electron_dr03TkSumPtHEEP.data());
    SafeSetBranchAddress("Electron_dxy", Electron_dxy.data());
    SafeSetBranchAddress("Electron_dxyErr", Electron_dxyErr.data());
    SafeSetBranchAddress("Electron_dz", Electron_dz.data());
    SafeSetBranchAddress("Electron_dzErr", Electron_dzErr.data());
    SafeSetBranchAddress("Electron_eInvMinusPInv", Electron_eInvMinusPInv.data());
    SafeSetBranchAddress("Electron_energyErr", Electron_energyErr.data());
    SafeSetBranchAddress("Electron_eta", Electron_eta.data());
    SafeSetBranchAddress("Electron_hoe", Electron_hoe.data());
    SafeSetBranchAddress("Electron_ip3d", Electron_ip3d.data());
    SafeSetBranchAddress("Electron_isPFcand", Electron_isPFcand.data());
    SafeSetBranchAddress("Electron_jetNDauCharged", Electron_jetNDauCharged.data());
    SafeSetBranchAddress("Electron_jetPtRelv2", Electron_jetPtRelv2.data());
    SafeSetBranchAddress("Electron_jetRelIso", Electron_jetRelIso.data());
    SafeSetBranchAddress("Electron_lostHits", Electron_lostHits.data());
    SafeSetBranchAddress("Electron_mass", Electron_mass.data());
    SafeSetBranchAddress("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all.data());
    SafeSetBranchAddress("Electron_miniPFRelIso_chg", Electron_miniPFRelIso_chg.data());
    SafeSetBranchAddress("Electron_mvaTTH", Electron_mvaTTH.data());
    SafeSetBranchAddress("Electron_pdgId", Electron_pdgId.data());
    SafeSetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all.data());
    SafeSetBranchAddress("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg.data());
    SafeSetBranchAddress("Electron_phi", Electron_phi.data());
    SafeSetBranchAddress("Electron_pt", Electron_pt.data());
    SafeSetBranchAddress("Electron_r9", Electron_r9.data());
    SafeSetBranchAddress("Electron_scEtOverPt", Electron_scEtOverPt.data());
    SafeSetBranchAddress("Electron_seedGain", Electron_seedGain.data());
    SafeSetBranchAddress("Electron_sieie", Electron_sieie.data());
    SafeSetBranchAddress("Electron_sip3d", Electron_sip3d.data());
    SafeSetBranchAddress("Electron_genPartFlav", Electron_genPartFlav.data());
    if (Run == 3) {
        SafeSetBranchAddress("Electron_cutBased", Electron_cutBased.data());
        SafeSetBranchAddress("Electron_genPartIdx", Electron_genPartIdx.data());
        SafeSetBranchAddress("Electron_jetIdx", Electron_jetIdx.data());
        SafeSetBranchAddress("Electron_mvaIso", Electron_mvaIso.data());
        SafeSetBranchAddress("Electron_mvaIso_WP80", Electron_mvaIso_WP80.data());
        SafeSetBranchAddress("Electron_mvaIso_WP90", Electron_mvaIso_WP90.data());
        SafeSetBranchAddress("Electron_mvaNoIso", Electron_mvaNoIso.data());
        SafeSetBranchAddress("Electron_mvaNoIso_WP80", Electron_mvaNoIso_WP80.data());
        SafeSetBranchAddress("Electron_mvaNoIso_WP90", Electron_mvaNoIso_WP90.data());
    } else if(Run == 2) {
        SafeSetBranchAddress("Electron_cutBased", Electron_cutBased_RunII.data());
        SafeSetBranchAddress("Electron_genPartIdx", Electron_genPartIdx_RunII.data());
        SafeSetBranchAddress("Electron_jetIdx", Electron_jetIdx_RunII.data());
        SafeSetBranchAddress("Electron_mvaFall17V2Iso", Electron_mvaFall17V2Iso.data());
        SafeSetBranchAddress("Electron_mvaFall17V2Iso_WP80", Electron_mvaFall17V2Iso_WP80.data());
        SafeSetBranchAddress("Electron_mvaFall17V2Iso_WP90", Electron_mvaFall17V2Iso_WP90.data());
        SafeSetBranchAddress("Electron_mvaFall17V2Iso_WPL", Electron_mvaFall17V2Iso_WPL.data());
        SafeSetBranchAddress("Electron_mvaFall17V2noIso", Electron_mvaFall17V2noIso.data());
        SafeSetBranchAddress("Electron_mvaFall17V2noIso_WP80", Electron_mvaFall17V2noIso_WP80.data());
        SafeSetBranchAddress("Electron_mvaFall17V2noIso_WP90", Electron_mvaFall17V2noIso_WP90.data());
        SafeSetBranchAddress("Electron_mvaFall17V2noIso_WPLoose", Electron_mvaFall17V2noIso_WPL.data());
        SafeSetBranchAddress("Electron_dEsigmaUp", Electron_dEsigmaUp.data());
        SafeSetBranchAddress("Electron_dEsigmaDown", Electron_dEsigmaDown.data());
    }

    // Photon----------------------------
    SetBranchWithRunCheck("nPhoton", nPhoton, nPhoton_RunII);
    SafeSetBranchAddress("Photon_eta", Photon_eta.data());
    SafeSetBranchAddress("Photon_hoe", Photon_hoe.data());
    SafeSetBranchAddress("Photon_isScEtaEB", Photon_isScEtaEB.data());
    SafeSetBranchAddress("Photon_isScEtaEE", Photon_isScEtaEE.data());
    SafeSetBranchAddress("Photon_mvaID", Photon_mvaID.data());
    SafeSetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80.data());
    SafeSetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90.data());
    SafeSetBranchAddress("Photon_phi", Photon_phi.data());
    SafeSetBranchAddress("Photon_pt", Photon_pt.data());
    SafeSetBranchAddress("Photon_sieie", Photon_sieie.data());
    if (Run == 3) {
        SafeSetBranchAddress("Photon_energyRaw", Photon_energyRaw.data());
        SafeSetBranchAddress("Photon_cutBased", Photon_cutBased.data());
    } else if(Run == 2) {
        SafeSetBranchAddress("Photon_cutBased", Photon_cutBased_RunII.data());
    }

    //Jet----------------------------
    SetBranchWithRunCheck("nJet", nJet, nJet_RunII);
    SafeSetBranchAddress("Jet_area", Jet_area.data());
    SafeSetBranchAddress("Jet_btagDeepFlavB", Jet_btagDeepFlavB.data());
    SafeSetBranchAddress("Jet_btagDeepFlavCvB", Jet_btagDeepFlavCvB.data());
    SafeSetBranchAddress("Jet_btagDeepFlavCvL", Jet_btagDeepFlavCvL.data());
    SafeSetBranchAddress("Jet_btagDeepFlavQG", Jet_btagDeepFlavQG.data());
    SafeSetBranchAddress("Jet_chEmEF", Jet_chEmEF.data());
    SafeSetBranchAddress("Jet_chHEF", Jet_chHEF.data());
    SafeSetBranchAddress("Jet_eta", Jet_eta.data());
    SafeSetBranchAddress("Jet_hfadjacentEtaStripsSize", Jet_hfadjacentEtaStripsSize.data());
    SafeSetBranchAddress("Jet_hfcentralEtaStripSize", Jet_hfcentralEtaStripSize.data());
    SafeSetBranchAddress("Jet_hfsigmaEtaEta", Jet_hfsigmaEtaEta.data());
    SafeSetBranchAddress("Jet_hfsigmaPhiPhi", Jet_hfsigmaPhiPhi.data());
    SafeSetBranchAddress("Jet_mass", Jet_mass.data());
    SafeSetBranchAddress("Jet_muEF", Jet_muEF.data());
    SafeSetBranchAddress("Jet_muonSubtrFactor", Jet_muonSubtrFactor.data());
    SafeSetBranchAddress("Jet_nConstituents", Jet_nConstituents.data());
    SafeSetBranchAddress("Jet_neEmEF", Jet_neEmEF.data());
    SafeSetBranchAddress("Jet_neHEF", Jet_neHEF.data());
    SafeSetBranchAddress("Jet_phi", Jet_phi.data());
    SafeSetBranchAddress("Jet_pt", Jet_pt.data());
    SafeSetBranchAddress("Jet_rawFactor", Jet_rawFactor.data());
    if (Run == 3) {
        SafeSetBranchAddress("Jet_PNetRegPtRawCorr", Jet_PNetRegPtRawCorr.data());
        SafeSetBranchAddress("Jet_PNetRegPtRawCorrNeutrino", Jet_PNetRegPtRawCorrNeutrino.data());
        SafeSetBranchAddress("Jet_PNetRegPtRawRes", Jet_PNetRegPtRawRes.data());
        SafeSetBranchAddress("Jet_btagPNetB", Jet_btagPNetB.data());
        SafeSetBranchAddress("Jet_btagPNetCvB", Jet_btagPNetCvB.data());
        SafeSetBranchAddress("Jet_btagPNetCvL", Jet_btagPNetCvL.data());
        SafeSetBranchAddress("Jet_btagPNetQvG", Jet_btagPNetQvG.data());
        SafeSetBranchAddress("Jet_btagPNetTauVJet", Jet_btagPNetTauVJet.data());
        SafeSetBranchAddress("Jet_btagRobustParTAK4B", Jet_btagRobustParTAK4B.data());
        SafeSetBranchAddress("Jet_btagRobustParTAK4CvB", Jet_btagRobustParTAK4CvB.data());
        SafeSetBranchAddress("Jet_btagRobustParTAK4CvL", Jet_btagRobustParTAK4CvL.data());
        SafeSetBranchAddress("Jet_btagRobustParTAK4QG", Jet_btagRobustParTAK4QG.data());
        SafeSetBranchAddress("Jet_electronIdx1", Jet_electronIdx1.data());
        SafeSetBranchAddress("Jet_electronIdx2", Jet_electronIdx2.data());
        SafeSetBranchAddress("Jet_genJetIdx", Jet_genJetIdx.data());
        SafeSetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour.data());
        SafeSetBranchAddress("Jet_jetId", Jet_jetId.data());
        SafeSetBranchAddress("Jet_muonIdx1", Jet_muonIdx1.data());
        SafeSetBranchAddress("Jet_muonIdx2", Jet_muonIdx2.data());
        SafeSetBranchAddress("Jet_nElectrons", Jet_nElectrons.data());
        SafeSetBranchAddress("Jet_nMuons", Jet_nMuons.data());
        SafeSetBranchAddress("Jet_nSVs", Jet_nSVs.data());
        SafeSetBranchAddress("Jet_partonFlavour", Jet_partonFlavour.data());
        SafeSetBranchAddress("Jet_svIdx1", Jet_svIdx1.data());
        SafeSetBranchAddress("Jet_svIdx2", Jet_svIdx2.data());
        // v13 leaf type is UChar_t; widened into Jet_chMultiplicity/Jet_neMultiplicity in Loop()
        SafeSetBranchAddress("Jet_chMultiplicity", Buf_Jet_chMultiplicity.data());
        SafeSetBranchAddress("Jet_neMultiplicity", Buf_Jet_neMultiplicity.data());
    } else if (Run == 2) {
        SafeSetBranchAddress("Jet_bRegCorr", Jet_bRegCorr.data());
        SafeSetBranchAddress("Jet_bRegRes", Jet_bRegRes.data());
        SafeSetBranchAddress("Jet_btagCSVV2", Jet_btagCSVV2.data());
        //SafeSetBranchAddress("Jet_btagDeepB", Jet_btagDeepB.data());
        //SafeSetBranchAddress("Jet_btagDeepCvB", Jet_btagDeepCvB.data());
        //SafeSetBranchAddress("Jet_btagDeepCvL", Jet_btagDeepCvL.data());
        SafeSetBranchAddress("Jet_cRegCorr", Jet_cRegCorr.data());
        SafeSetBranchAddress("Jet_cRegRes", Jet_cRegRes.data());
        SafeSetBranchAddress("Jet_chFPV0EF", Jet_chFPV0EF.data());
        SafeSetBranchAddress("Jet_cleanmask", Jet_cleanmask.data());
        SafeSetBranchAddress("Jet_electronIdx1", Jet_electronIdx1_RunII.data());
        SafeSetBranchAddress("Jet_electronIdx2", Jet_electronIdx2_RunII.data());
        SafeSetBranchAddress("Jet_genJetIdx", Jet_genJetIdx_RunII.data());
        SafeSetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour_RunII.data());
        SafeSetBranchAddress("Jet_jetId", Jet_jetId_RunII.data());
        SafeSetBranchAddress("Jet_muonIdx1", Jet_muonIdx1_RunII.data());
        SafeSetBranchAddress("Jet_muonIdx2", Jet_muonIdx2_RunII.data());
        SafeSetBranchAddress("Jet_nElectrons", Jet_nElectrons_RunII.data());
        SafeSetBranchAddress("Jet_nMuons", Jet_nMuons_RunII.data());
        SafeSetBranchAddress("Jet_partonFlavour", Jet_partonFlavour_RunII.data());
        SafeSetBranchAddress("Jet_puId", Jet_puId.data());
        SafeSetBranchAddress("Jet_puIdDisc", Jet_puIdDisc.data());
        SafeSetBranchAddress("Jet_qgl", Jet_qgl.data());
    }

    //Tau---------------------------- 
    SetBranchWithRunCheck("nTau", nTau, nTau_RunII);
    SafeSetBranchAddress("Tau_dxy", Tau_dxy.data());
    SafeSetBranchAddress("Tau_dz", Tau_dz.data());
    SafeSetBranchAddress("Tau_eta", Tau_eta.data());
    SafeSetBranchAddress("Tau_genPartFlav", Tau_genPartFlav.data());
    SafeSetBranchAddress("Tau_genPartidDeepTau2017v2p1VSe", Tau_idDeepTau2018v2p5VSe.data());
    SafeSetBranchAddress("Tau_genPartidDeepTau2017v2p1VSjet", Tau_idDeepTau2018v2p5VSjet.data());
    SafeSetBranchAddress("Tau_genPartidDeepTau2017v2p1VSmu", Tau_idDeepTau2018v2p5VSmu.data());
    SafeSetBranchAddress("Tau_mass", Tau_mass.data());
    SafeSetBranchAddress("Tau_phi", Tau_phi.data());
    SafeSetBranchAddress("Tau_pt", Tau_pt.data());
    if (Run == 3) {
        SafeSetBranchAddress("Tau_charge", Tau_charge.data());
        SafeSetBranchAddress("Tau_decayMode", Tau_decayMode.data());
        SafeSetBranchAddress("Tau_genPartIdx", Tau_genPartIdx.data());
        SafeSetBranchAddress("Tau_idDecayModeNewDMs", Tau_idDecayModeNewDMs.data());
        SafeSetBranchAddress("Tau_idDeepTau2018v2p5VSe", Tau_idDeepTau2018v2p5VSe.data());
        SafeSetBranchAddress("Tau_idDeepTau2018v2p5VSjet", Tau_idDeepTau2018v2p5VSjet.data());
        SafeSetBranchAddress("Tau_idDeepTau2018v2p5VSmu", Tau_idDeepTau2018v2p5VSmu.data());
    } else if(Run == 2) {
        SafeSetBranchAddress("Tau_charge", Tau_charge_RunII.data());
        SafeSetBranchAddress("Tau_decayMode", Tau_decayMode_RunII.data());
        SafeSetBranchAddress("Tau_genPartIdx", Tau_genPartIdx_RunII.data());
    }

    //FatJet----------------------------
    SetBranchWithRunCheck("nFatJet", nFatJet, nFatJet_RunII);
    SafeSetBranchAddress("FatJet_area", FatJet_area.data());
    SafeSetBranchAddress("FatJet_btagDDBvLV2", FatJet_btagDDBvLV2.data());
    SafeSetBranchAddress("FatJet_btagDDCvBV2", FatJet_btagDDCvBV2.data());
    SafeSetBranchAddress("FatJet_btagDDCvLV2", FatJet_btagDDCvLV2.data());
    SafeSetBranchAddress("FatJet_btagDeepB", FatJet_btagDeepB.data());
    SafeSetBranchAddress("FatJet_btagHbb", FatJet_btagHbb.data());
    SafeSetBranchAddress("FatJet_eta", FatJet_eta.data());
    SafeSetBranchAddress("FatJet_lsf3", FatJet_lsf3.data());
    SafeSetBranchAddress("FatJet_mass", FatJet_mass.data());
    SafeSetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop.data());
    SafeSetBranchAddress("FatJet_nBHadrons", FatJet_nBHadrons.data());
    SafeSetBranchAddress("FatJet_nCHadrons", FatJet_nCHadrons.data());
    SafeSetBranchAddress("FatJet_nConstituents", FatJet_nConstituents.data());
    SafeSetBranchAddress("FatJet_particleNet_QCD", FatJet_particleNet_QCD.data());
    SafeSetBranchAddress("FatJet_phi", FatJet_phi.data());
    SafeSetBranchAddress("FatJet_pt", FatJet_pt.data());
    SafeSetBranchAddress("FatJet_tau1", FatJet_tau1.data());
    SafeSetBranchAddress("FatJet_tau2", FatJet_tau2.data());
    SafeSetBranchAddress("FatJet_tau3", FatJet_tau3.data());
    SafeSetBranchAddress("FatJet_tau4", FatJet_tau4.data());
    if (Run == 3) {
        SafeSetBranchAddress("FatJet_genJetAK8Idx", FatJet_genJetAK8Idx.data());
        SafeSetBranchAddress("FatJet_jetId", FatJet_jetId.data());
        SafeSetBranchAddress("FatJet_particleNetWithMass_H4qvsQCD", FatJet_particleNetWithMass_H4qvsQCD.data());
        SafeSetBranchAddress("FatJet_particleNetWithMass_HbbvsQCD", FatJet_particleNetWithMass_HbbvsQCD.data());
        SafeSetBranchAddress("FatJet_particleNetWithMass_HccvsQCD", FatJet_particleNetWithMass_HccvsQCD.data());
        SafeSetBranchAddress("FatJet_particleNetWithMass_QCD", FatJet_particleNetWithMass_QCD.data());
        SafeSetBranchAddress("FatJet_particleNetWithMass_TvsQCD", FatJet_particleNetWithMass_TvsQCD.data());
        SafeSetBranchAddress("FatJet_particleNetWithMass_WvsQCD", FatJet_particleNetWithMass_WvsQCD.data());
        SafeSetBranchAddress("FatJet_particleNetWithMass_ZvsQCD", FatJet_particleNetWithMass_ZvsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_QCD0HF", FatJet_particleNet_QCD0HF.data());
        SafeSetBranchAddress("FatJet_particleNet_QCD1HF", FatJet_particleNet_QCD1HF.data());
        SafeSetBranchAddress("FatJet_particleNet_QCD2HF", FatJet_particleNet_QCD2HF.data());
        SafeSetBranchAddress("FatJet_particleNet_XbbVsQCD", FatJet_particleNet_XbbVsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_XccVsQCD", FatJet_particleNet_XccVsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_XggVsQCD", FatJet_particleNet_XggVsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_XqqVsQCD", FatJet_particleNet_XqqVsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_XteVsQCD", FatJet_particleNet_XteVsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_XtmVsQCD", FatJet_particleNet_XtmVsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_XttVsQCD", FatJet_particleNet_XttVsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_massCorr", FatJet_particleNet_massCorr.data());
        SafeSetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1.data());
        SafeSetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2.data());
        SafeSetBranchAddress("FatJet_n2b1", FatJet_n2b1.data());
        SafeSetBranchAddress("FatJet_n3b1", FatJet_n3b1.data());
        SafeSetBranchAddress("FatJet_chMultiplicity", FatJet_chMultiplicity.data());
        SafeSetBranchAddress("FatJet_neMultiplicity", FatJet_neMultiplicity.data());
        SafeSetBranchAddress("FatJet_chHEF", FatJet_chHEF.data());
        SafeSetBranchAddress("FatJet_neHEF", FatJet_neHEF.data());
        SafeSetBranchAddress("FatJet_chEmEF", FatJet_chEmEF.data());
        SafeSetBranchAddress("FatJet_neEmEF", FatJet_neEmEF.data());
        SafeSetBranchAddress("FatJet_muEF", FatJet_muEF.data());
        SetBranchWithRunCheck("nSubJet", nSubJet, nSubJet_RunII);
        SafeSetBranchAddress("SubJet_pt", SubJet_pt.data());
        SafeSetBranchAddress("SubJet_eta", SubJet_eta.data());
        SafeSetBranchAddress("SubJet_phi", SubJet_phi.data());
        SafeSetBranchAddress("SubJet_mass", SubJet_mass.data());
        SafeSetBranchAddress("SubJet_btagDeepB", SubJet_btagDeepB.data());
        SetBranchWithRunCheck("nSV", nSV, nSV_RunII);
        SafeSetBranchAddress("SV_pt", SV_pt.data());
        SafeSetBranchAddress("SV_eta", SV_eta.data());
        SafeSetBranchAddress("SV_phi", SV_phi.data());
        SafeSetBranchAddress("SV_mass", SV_mass.data());
        SafeSetBranchAddress("SV_dlenSig", SV_dlenSig.data());
        SafeSetBranchAddress("SV_dxySig", SV_dxySig.data());
        SafeSetBranchAddress("SV_chi2", SV_chi2.data());
        SafeSetBranchAddress("SV_pAngle", SV_pAngle.data());
        SafeSetBranchAddress("SV_ntracks", SV_ntracks.data());
    } else if(Run == 2) {
        SafeSetBranchAddress("FatJet_genJetAK8Idx", FatJet_genJetAK8Idx_RunII.data());
        SafeSetBranchAddress("FatJet_jetId", FatJet_jetId_RunII.data());
        SafeSetBranchAddress("FatJet_particleNetMD_QCD", FatJet_particleNetMD_QCD.data());
        SafeSetBranchAddress("FatJet_particleNetMD_Xbb", FatJet_particleNetMD_Xbb.data());
        SafeSetBranchAddress("FatJet_particleNetMD_Xcc", FatJet_particleNetMD_Xcc.data());
        SafeSetBranchAddress("FatJet_particleNetMD_Xqq", FatJet_particleNetMD_Xqq.data());
        SafeSetBranchAddress("FatJet_particleNet_H4qvsQCD", FatJet_particleNet_H4qvsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_HbbvsQCD", FatJet_particleNet_HbbvsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_HccvsQCD", FatJet_particleNet_HccvsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_TvsQCD", FatJet_particleNet_TvsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_WvsQCD", FatJet_particleNet_WvsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_ZvsQCD", FatJet_particleNet_ZvsQCD.data());
        SafeSetBranchAddress("FatJet_particleNet_mass", FatJet_particleNet_mass.data());
        SafeSetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1_RunII.data());
        SafeSetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2_RunII.data());
    }
    // MET----------------------------
    SafeSetBranchAddress("MET_pt", &MET_pt);
    SafeSetBranchAddress("MET_phi", &MET_phi);
    SafeSetBranchAddress("PuppiMET_pt", &PuppiMET_pt);
    SafeSetBranchAddress("PuppiMET_phi", &PuppiMET_phi);
    SafeSetBranchAddress("PuppiMET_ptUnclusteredUp", &PuppiMET_ptUnclusteredUp);
    SafeSetBranchAddress("PuppiMET_phiUnclusteredUp", &PuppiMET_phiUnclusteredUp);
    SafeSetBranchAddress("PuppiMET_ptUnclusteredDown", &PuppiMET_ptUnclusteredDown);
    SafeSetBranchAddress("PuppiMET_phiUnclusteredDown", &PuppiMET_phiUnclusteredDown); 

    //Rho----------------------------
    if(Run == 3) {
        SafeSetBranchAddress("Rho_fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll);
    } else if(Run == 2) {
        SafeSetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll);
    }

    // PV----------------------------
    SafeSetBranchAddress("PV_chi2", &PV_chi2);
    SafeSetBranchAddress("PV_ndof", &PV_ndof);
    SafeSetBranchAddress("PV_score", &PV_score);
    SafeSetBranchAddress("PV_x", &PV_x);
    SafeSetBranchAddress("PV_y", &PV_y);
    SafeSetBranchAddress("PV_z", &PV_z);
    if (Run==3) { 
        SafeSetBranchAddress("PV_npvs", &PV_npvs);
        SafeSetBranchAddress("PV_npvsGood", &PV_npvsGood);
    } else if(Run==2) {
        SafeSetBranchAddress("PV_npvs", &PV_npvs_RunII);
        SafeSetBranchAddress("PV_npvsGood", &PV_npvsGood_RunII);
    }

    //L1PreFireweight----------------------------
    SafeSetBranchAddress("L1PreFiringWeight_Nom", &L1PreFiringWeight_Nom);
    SafeSetBranchAddress("L1PreFiringWeight_Dn", &L1PreFiringWeight_Dn);
    SafeSetBranchAddress("L1PreFiringWeight_Up", &L1PreFiringWeight_Up);

    //Flags----------------------------
    SafeSetBranchAddress("Flag_METFilters", &Flag_METFilters);
    SafeSetBranchAddress("Flag_goodVertices", &Flag_goodVertices);
    SafeSetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter);
    SafeSetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter);
    SafeSetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter);
    SafeSetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter);
    SafeSetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter);
    SafeSetBranchAddress("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter);
    SafeSetBranchAddress("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter);
    SafeSetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter);
    SafeSetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter);
    SafeSetBranchAddress("run", &RunNumber);
    SafeSetBranchAddress("luminosityBlock", &LumiBlock);
    SafeSetBranchAddress("event", &EventNumber);

    // TrigObj----------------------------
    SetBranchWithRunCheck("nTrigObj", nTrigObj, nTrigObj_RunII);
    SafeSetBranchAddress("TrigObj_pt", TrigObj_pt.data());
    SafeSetBranchAddress("TrigObj_eta", TrigObj_eta.data());
    SafeSetBranchAddress("TrigObj_phi", TrigObj_phi.data());
    if (Run == 3) {
        SafeSetBranchAddress("TrigObj_id", TrigObj_id.data());
    } else {
        SafeSetBranchAddress("TrigObj_id", TrigObj_id_RunII.data());
    }
    SafeSetBranchAddress("TrigObj_filterBits", TrigObj_filterBits.data());

    string json_path = string(getenv("SKNANO_DATA")) + "/" + DataEra.Data() + "/Trigger/HLT_Path.json";
    ifstream json_file(json_path);
    if (json_file.is_open()) {
        cout << "[SKNanoLoader::Init] Loading HLT Paths in " << json_path << endl;
        json j;
        json_file >> j;
        RVec<TString> not_in_tree;
        for (auto& [key, value] : j.items()) {
            if (!value.contains("active")) continue;
            if (!value["active"]) continue;
            Bool_t* passHLT = new Bool_t();
            TString key_str = key;
            TriggerMap[key_str].first = passHLT;
            TriggerMap[key_str].second = value["lumi"];
            //if key_str is in tree, set branch address
            if (fChain->GetBranch(key_str)) {
                // In some data file, part of the trigger set is missing (changed during the run?)
                if (IsDATA) {
                    SuperSafeSetBranchAddress(key_str, TriggerMap[key_str].first);
                } else {
                    SafeSetBranchAddress(key_str, TriggerMap[key_str].first);
                }
            } else if(key_str=="Full") {
                *TriggerMap[key_str].first = true;
            } else{
                not_in_tree.push_back(key_str);
                TriggerMap.erase(key_str);
            }   
        }
        if (not_in_tree.size() > 0) {
            cout << "\033[1;33m[SKNanoLoader::Init] Following HLT Paths are not in the tree\033[0m" << endl;
            for (auto &path : not_in_tree) {
                cout << "\033[1;33m" << path << "\033[0m" << endl;
            }
        }
    }
    else cerr << "[SKNanoLoader::Init] Cannot open " << json_path << endl;
    
}
