#include "fakeV4.h"

fakeV4::fakeV4() : leptonType(LeptonType::NONE), lowestPtCut(10.) {}

fakeV4::~fakeV4() {}

void fakeV4::initializeAnalyzer() {
    // Set flags
    MeasFakeMu = HasFlag("MeasFakeMu");
    MeasFakeEl = HasFlag("MeasFakeEl");
    RunSyst = HasFlag("RunSyst");
    RunNoHEMVeto = HasFlag("RunNoHEMVeto");
    MakeTree = HasFlag("MakeTree");

    // Determine lepton type and configure triggers
    if (MeasFakeMu) {
        leptonType = LeptonType::MUON;
        ptcorr_bins = {10., 12., 14., 17., 20., 30., 50., 100., 200.};
        abseta_bins = {0., 0.9, 1.6, 2.4};
        triggers = {
            {"HLT_Mu8_TrkIsoVVL",  "Mu8",  10., false},
            {"HLT_Mu17_TrkIsoVVL", "Mu17", 20., false}
        };
        lowestPtCut = 10.;
    } else if (MeasFakeEl) {
        leptonType = LeptonType::ELECTRON;
        ptcorr_bins = {15., 17., 20., 25., 35., 50., 100., 200.};
        abseta_bins = {0., 0.8, 1.479, 2.5};
        triggers = {
            {"HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",  "El8",  10., false},
            {"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", "El12", 15., false},
            {"HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", "El23", 25., false}
        };
        lowestPtCut = 10.;
    } else {
        throw std::runtime_error("[fakeV4::initializeAnalyzer] No lepton type specified by flags. Use MeasFakeMu or MeasFakeEl.");
    }

    // Set IDs
    MuonIDs = new IDContainer("HcToWATight", ((Run == 2) ? "HcToWALooseRun2" : "HcToWALooseRun3"));
    ElectronIDs = new IDContainer("HcToWATight", ((Run == 2) ? "HcToWALooseRun2" : "HcToWALooseRun3"));

    // Initialize correction
    myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA ? DataStream : MCSample, IsDATA);

    if (MakeTree) NewTree("mu");

    // Initialize SystematicHelper
    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/AnalyzerTools/FakeSystematics.data.yaml", DataStream, DataEra);
    } else {
        if (RunSyst) {
            systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/AnalyzerTools/FakeSystematics.mc.yaml", MCSample, DataEra);
        } else {
            systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/AnalyzerTools/noSyst.yaml", DataStream, DataEra);
        }
    }
}

void fakeV4::executeEvent() {
    Event ev = GetEvent();

    // Initial cutflow entry
    float initialWeight = IsDATA ? 1.0 : MCweight() * ev.GetTriggerLumi("Full");
    fillCutflow(CutStage::Initial, Channel::NONE, "event", initialWeight, "Central");

    RVec<Jet> rawJets = GetAllJets();
    if (!PassNoiseFilter(rawJets, ev)) return;
    fillCutflow(CutStage::NoiseFilter, Channel::NONE, "event", initialWeight, "Central");

    RVec<Muon> rawMuons = GetAllMuons();
    if (!PassVetoMap(rawJets, rawMuons, "jetvetomap")) return;
    fillCutflow(CutStage::VetoMap, Channel::NONE, "event", initialWeight, "Central");

    RVec<Electron> rawElectrons = GetAllElectrons();

    // Check which triggers fired (store in triggers[i].fired)
    bool anyFired = false;
    for (auto& trig : triggers) {
        trig.fired = ev.PassTrigger(trig.name);
        anyFired |= trig.fired;
    }
    if (!anyFired) return;
    fillCutflow(CutStage::AnyTrigger, Channel::NONE, "event", initialWeight, "Central");

    RVec<Gen> genParts = !IsDATA ? GetAllGens() : RVec<Gen>();
    RVec<GenJet> genJets = !IsDATA ? GetAllGenJets() : RVec<GenJet>();

    // Loop over IDs (loose and tight)
    RVec<TString> IDs = {"loose", "tight"};

    for (const auto& ID : IDs) {
        if (RunSyst && systHelper) {
            // Process Central objects and weight-only systematics
            RecoObjects centralObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, ID, "Central");
            Channel selectedChannel = selectEvent(ev, centralObjects, ID, "Central");

            if (selectedChannel != Channel::NONE) {
                fillCutflow(CutStage::Final, selectedChannel, ID, initialWeight, "Central");
                WeightInfo centralWeights = getWeights(selectedChannel, ID, ev, centralObjects, genParts, "Central");
                fillObjects(selectedChannel, ID, centralObjects, centralWeights, "Central");

                // process weight-only systematics with central objects
                vector<string> weightOnlySysts = systHelper->getWeightOnlySystematics();
                for (const auto &systName : weightOnlySysts) {
                    TString systNameUp = systName + "_Up";
                    WeightInfo weightsUp = getWeights(selectedChannel, ID, ev, centralObjects, genParts, systNameUp);
                    fillObjects(selectedChannel, ID, centralObjects, weightsUp, systNameUp);

                    TString systNameDown = systName + "_Down";
                    WeightInfo weightsDown = getWeights(selectedChannel, ID, ev, centralObjects, genParts, systNameDown);
                    fillObjects(selectedChannel, ID, centralObjects, weightsDown, systNameDown);
                }
            }

            // Process systematics requiring evtLoopAgain
            for (const auto& syst : *systHelper) {
                TString systName = syst.iter_name;

                // Skip Central and weight-only systematics
                if (systName == "Central" || (!systHelper->findSystematic(syst.syst_name)->evtLoopAgain)) continue;

                RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, ID, systName);
                Channel systChannel = selectEvent(ev, recoObjects, ID, systName);

                if (systChannel != Channel::NONE) {
                    WeightInfo weights = getWeights(systChannel, ID, ev, recoObjects, genParts, systName);
                    fillObjects(systChannel, ID, recoObjects, weights, systName);
                }
            }
        } else {
            // systematics are off
            RecoObjects recoObjects = defineObjects(ev, rawMuons, rawElectrons, rawJets, genJets, ID, "Central");
            Channel selectedChannel = selectEvent(ev, recoObjects, ID, "Central");

            if (selectedChannel != Channel::NONE) {
                fillCutflow(CutStage::Final, selectedChannel, ID, initialWeight, "Central");
                WeightInfo weights = getWeights(selectedChannel, ID, ev, recoObjects, genParts, "Central");
                fillObjects(selectedChannel, ID, recoObjects, weights, "Central");

                //==== flat ntuple: loose pass 의 Inclusive 영역에서만 한 줄.
                //     ID loop 가 loose/tight 두 번 도므로 tight 에서도 부르면
                //     같은 muon 이 두 줄로 중복된다.
                if (MakeTree && ID == "loose" && selectedChannel == Channel::INCLUSIVE) {
                    looseMuons     = recoObjects.looseMuons;
                    looseElectrons = recoObjects.looseElectrons;
                    jets           = recoObjects.tightJets;
                    this->rawJets  = rawJets;
                    this->gens     = genParts;
                    this->genJets  = genJets;
                    METv           = recoObjects.METv;
                    fillTreeRow(ev, 1.0);   // topPt=1 (QCD), pileup-jet-ID SF is Run2 only
                }
            }
        }
    }
}

fakeV4::RecoObjects fakeV4::defineObjects(Event& ev,
                                                     const RVec<Muon>& rawMuons,
                                                     const RVec<Electron>& rawElectrons,
                                                     const RVec<Jet>& rawJets,
                                                     const RVec<GenJet>& genJets,
                                                     const TString& ID,
                                                     const TString& syst) {
    // Create copies for systematic variations
    RVec<Muon> allMuons = rawMuons;
    RVec<Electron> allElectrons = rawElectrons;
    RVec<Jet> allJets = rawJets;

    // Apply systematic variations
    if (syst.Contains("ElectronEn")) {
        TString variation = syst.Contains("Up") ? "up" : "down";
        allElectrons = ScaleElectrons(ev, allElectrons, variation);
    } else if (syst.Contains("ElectronRes")) {
        TString variation = syst.Contains("Up") ? "up" : "down";
        allElectrons = SmearElectrons(allElectrons, variation);
    } else if (syst.Contains("MuonEn")) {
        TString variation = syst.Contains("Up") ? "up" : "down";
        allMuons = ScaleMuons(allMuons, variation);
    } else if (syst.Contains("JetEn")) {
        TString variation = syst.Contains("Up") ? "up" : "down";
        TString systSource = "total";
        allJets = ScaleJets(allJets, variation, systSource);
    } else if (syst.Contains("JetRes")) {
        TString variation = syst.Contains("Up") ? "up" : "down";
        allJets = SmearJets(allJets, genJets, variation);
    } else {
        // No scale variation
    }

    // Get MET with unclustered energy variation if applicable
    Particle METv_default;
    if (syst.Contains("UnclusteredEn")) {
        Event::MET_Syst variation = syst.Contains("Up") ? Event::MET_Syst::UE_UP : Event::MET_Syst::UE_DOWN;
        METv_default = ev.GetMETVector(Event::MET_Type::PUPPI, variation);
    } else {
        METv_default = ev.GetMETVector(Event::MET_Type::PUPPI);
    }
    Particle METv = ApplyTypeICorrection(METv_default, allJets, allElectrons, allMuons);

    // Sort objects by pT
    sort(allMuons.begin(), allMuons.end(), [](const Muon& a, const Muon& b) { return a.Pt() > b.Pt(); });
    sort(allElectrons.begin(), allElectrons.end(), [](const Electron& a, const Electron& b) { return a.Pt() > b.Pt(); });
    sort(allJets.begin(), allJets.end(), [](const Jet& a, const Jet& b) { return a.Pt() > b.Pt(); });

    // Select objects based on ID
    // HEM veto: apply to measurement electrons in 2018 unless RunNoHEMVeto flag is set
    bool applyHEMVeto = DataEra.Contains("2018") && !RunNoHEMVeto;

    // Veto leptons: 10 GeV, loose ID (for event selection veto, no HEM veto)
    RVec<Muon> vetoMuons = PickMuons(allMuons, "loose", 10., 2.4);
    RVec<Electron> vetoElectrons = PickElectrons(allElectrons, "loose", 10., 2.5, false);

    // Measurement leptons: muons at 10 GeV, electrons at 15 GeV (with HEM veto for 2018)
    RVec<Muon> looseMuons = PickMuons(allMuons, "loose", 10., 2.4);
    RVec<Electron> looseElectrons = PickElectrons(vetoElectrons, "loose", 15., 2.5, applyHEMVeto);

    RVec<Muon> tightMuons;
    RVec<Electron> tightElectrons;

    if (ID == "loose") {
        tightMuons = looseMuons;
        tightElectrons = looseElectrons;
    } else if (ID == "tight") {
        tightMuons = PickMuons(allMuons, "tight", 10., 2.4);
        tightElectrons = PickElectrons(looseElectrons, "tight", 15., 2.5, applyHEMVeto);
    }

    // Jet selection with potential selection variations
    // Order: tight -> vetoLep -> vetoMap -> PUID
    const float jetPtCut = getJetPtCut(syst);
    const float max_jeteta = DataEra.Contains("2016") ? 2.4 : 2.5;
    RVec<Jet> tightJets = SelectJets(allJets, "tight", jetPtCut, max_jeteta);
    tightJets = JetsVetoLeptonInside(tightJets, vetoElectrons, vetoMuons, 0.4);
    RVec<Jet> tightJets_noPUID;  // Jets before PUID for SF calculation
    if (Run == 2) {
        RVec<Jet> tightJets_vetoMap;
        for (const auto &jet : tightJets) {
            if (PassVetoMap(jet, allMuons, "jetvetomap")) tightJets_vetoMap.push_back(jet);
        }
        // Store jets before PUID for SF calculation, then apply PUID
        tightJets_noPUID = tightJets_vetoMap;
        tightJets = SelectJets(tightJets_vetoMap, "loosePuId", jetPtCut, max_jeteta);
    }

    // B-jet selection
    RVec<Jet> bjets;
    float wp = myCorr->GetBTaggingWP(JetTagging::JetFlavTagger::DeepJet, JetTagging::JetFlavTaggerWP::Medium);
    for (const auto& jet : tightJets) {
        float btagScore = jet.GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        if (btagScore > wp) bjets.emplace_back(jet);
    }

    // Compute mother jet flavours for leptons (MC only)
    RVec<int> looseMuonJetFlavours, tightMuonJetFlavours;
    RVec<int> looseElectronJetFlavours, tightElectronJetFlavours;

    if (!IsDATA) {
        for (const auto& mu : looseMuons) {
            looseMuonJetFlavours.push_back(getMotherJetFlavour(mu, allJets));
        }
        for (const auto& mu : tightMuons) {
            tightMuonJetFlavours.push_back(getMotherJetFlavour(mu, allJets));
        }
        for (const auto& el : looseElectrons) {
            looseElectronJetFlavours.push_back(getMotherJetFlavour(el, allJets));
        }
        for (const auto& el : tightElectrons) {
            tightElectronJetFlavours.push_back(getMotherJetFlavour(el, allJets));
        }
    }

    RecoObjects objects;
    objects.vetoMuons = vetoMuons;
    objects.looseMuons = looseMuons;
    objects.tightMuons = tightMuons;
    objects.vetoElectrons = vetoElectrons;
    objects.looseElectrons = looseElectrons;
    objects.tightElectrons = tightElectrons;
    objects.looseMuonJetFlavours = looseMuonJetFlavours;
    objects.tightMuonJetFlavours = tightMuonJetFlavours;
    objects.looseElectronJetFlavours = looseElectronJetFlavours;
    objects.tightElectronJetFlavours = tightElectronJetFlavours;
    objects.tightJets = tightJets;
    objects.tightJets_noPUID = tightJets_noPUID;
    objects.bjets = bjets;
    objects.genJets = genJets;
    objects.METv = METv;

    return objects;
}

fakeV4::Channel fakeV4::selectEvent(Event& ev, const RecoObjects& recoObjects, const TString& ID, const TString& syst) {
    const RVec<Muon>& muons = (ID == "loose") ? recoObjects.looseMuons : recoObjects.tightMuons;
    const RVec<Electron>& electrons = (ID == "loose") ? recoObjects.looseElectrons : recoObjects.tightElectrons;
    const RVec<Muon>& vetoMuons = recoObjects.vetoMuons;
    const RVec<Electron>& vetoElectrons = recoObjects.vetoElectrons;
    const RVec<Jet>& jets = recoObjects.tightJets;
    const RVec<Jet>& bjets = recoObjects.bjets;

    float weight = IsDATA ? 1.0 : MCweight() * ev.GetTriggerLumi("Full");

    if (leptonType == LeptonType::MUON) {
        const bool sglMu = (muons.size() == 1 && vetoMuons.size() == 1 &&
                            electrons.size() == 0 && vetoElectrons.size() == 0);
        const bool dblMu = (muons.size() == 2 && vetoMuons.size() == 2 &&
                            electrons.size() == 0 && vetoElectrons.size() == 0);

        if (! (sglMu || dblMu)) return Channel::NONE;
        fillCutflow(CutStage::LeptonSelection, Channel::INCLUSIVE, ID, weight, syst);

        // Use loosest pT cut for initial selection (trigger-specific cut applied in fillObjects)
        if (! (muons[0].Pt() > lowestPtCut)) return Channel::NONE;
        if (! (jets.size() > 0)) return Channel::NONE;
        if (syst.Contains("RequireHeavyTag") && bjets.size() == 0) return Channel::NONE;
        fillCutflow(CutStage::JetRequirements, Channel::INCLUSIVE, ID, weight, syst);

        if (sglMu) {
            // Check for away jet
            bool existAwayJet = false;
            for (const auto& jet : jets) {
                if (jet.DeltaR(muons[0]) > 0.7) {
                    existAwayJet = true;
                    break;
                }
            }
            if (!existAwayJet) return Channel::NONE;
            fillCutflow(CutStage::AwayJetRequirements, Channel::INCLUSIVE, ID, weight, syst);
            return Channel::INCLUSIVE;
        } else { // doubleMu
            Particle ZCand = muons[0] + muons[1];
            bool isOnZ = (fabs(ZCand.M() - 91.2) < 15.);
            if (!isOnZ) return Channel::NONE;
            fillCutflow(CutStage::ZMassWindow, Channel::ZENRICHED, ID, weight, syst);
            return Channel::ZENRICHED;
        }
    } else { // ELECTRON
        const bool sglEl = (electrons.size() == 1 && vetoElectrons.size() == 1 &&
                        muons.size() == 0 && vetoMuons.size() == 0);
        const bool dblEl = (electrons.size() == 2 && vetoElectrons.size() == 2 &&
                        muons.size() == 0 && vetoMuons.size() == 0);

        if (! (sglEl || dblEl)) return Channel::NONE;
        fillCutflow(CutStage::LeptonSelection, Channel::INCLUSIVE, ID, weight, syst);

        // Use loosest pT cut for initial selection (trigger-specific cut applied in fillObjects)
        if (! (electrons[0].Pt() > lowestPtCut)) return Channel::NONE;
        if (! (jets.size() > 0)) return Channel::NONE;
        if (syst.Contains("RequireHeavyTag") && bjets.size() == 0) return Channel::NONE;
        fillCutflow(CutStage::JetRequirements, Channel::INCLUSIVE, ID, weight, syst);

        if (sglEl) {
            // Check for away jet
            bool existAwayJet = false;
            for (const auto& jet : jets) {
                if (jet.DeltaR(electrons[0]) > 0.7) {
                    existAwayJet = true;
                    break;
                }
            }
            if (!existAwayJet) return Channel::NONE;
            fillCutflow(CutStage::AwayJetRequirements, Channel::INCLUSIVE, ID, weight, syst);
            return Channel::INCLUSIVE;
        } else { // dblEl
            Particle ZCand = electrons[0] + electrons[1];
            bool isOnZ = (fabs(ZCand.M() - 91.2) < 15.);
            if (!isOnZ) return Channel::NONE;
            fillCutflow(CutStage::ZMassWindow, Channel::ZENRICHED, ID, weight, syst);
            return Channel::ZENRICHED;
        }
    }

    return Channel::NONE;
}

fakeV4::WeightInfo fakeV4::getWeights(const Channel& channel,
                                                  const TString& ID,
                                                  Event& event,
                                                  const RecoObjects& recoObjects,
                                                  const RVec<Gen>& genParts,
                                                  const TString& syst) {
    WeightInfo weights;
    weights.genWeight = 1.0;
    weights.prefireWeight = 1.0;
    weights.pileupWeight = 1.0;
    weights.topPtWeight = 1.0;
    weights.muonRecoSF = 1.0;
    weights.eleRecoSF = 1.0;
    weights.btagSF = 1.0;
    weights.pileupIDSF = 1.0;

    if (!IsDATA) {
        weights.genWeight = MCweight() * event.GetTriggerLumi("Full");

        // Determine systematic variation
        MyCorrection::variation var = MyCorrection::variation::nom;
        if (syst.Contains("_Up")) var = MyCorrection::variation::up;
        else if (syst.Contains("_Down")) var = MyCorrection::variation::down;

        // L1 Prefire
        if (syst.Contains("L1Prefire")) {
            weights.prefireWeight = GetL1PrefireWeight(var);
        } else {
            weights.prefireWeight = GetL1PrefireWeight(MyCorrection::variation::nom);
        }

        // Top pT reweight
        if (MCSample.Contains("TTLL") || MCSample.Contains("TTLJ")) {
            weights.topPtWeight = myCorr->GetTopPtReweight(genParts);
        }

        // Lepton reconstruction SF
        const RVec<Muon>& muons = (ID == "loose") ? recoObjects.looseMuons : recoObjects.tightMuons;
        const RVec<Electron>& electrons = (ID == "loose") ? recoObjects.looseElectrons : recoObjects.tightElectrons;

        if (syst.Contains("MuonRecoSF")) {
            weights.muonRecoSF = myCorr->GetMuonRECOSF(muons, var);
        } else {
            weights.muonRecoSF = myCorr->GetMuonRECOSF(muons, MyCorrection::variation::nom);
        }

        if (syst.Contains("ElectronRecoSF")) {
            weights.eleRecoSF = myCorr->GetElectronRECOSF(electrons, var);
        } else {
            weights.eleRecoSF = myCorr->GetElectronRECOSF(electrons, MyCorrection::variation::nom);
        }

        // B-tagging SF for RequireHeavyTag
        if (syst.Contains("RequireHeavyTag")) {
            weights.btagSF = 1.0;   // our fork's MyCorrection lacks
            // GetBTaggingReweightMethod1a. Only the RequireHeavyTag systematic
            // uses it and Central never does, so 1.0 is exact here.
        }

        // Jet PUID SF for Run2
        // Use jets before PUID selection for SF calculation
        if (Run == 2) {
            const RVec<Jet>& jets = recoObjects.tightJets_noPUID;
            const RVec<GenJet>& genJets = recoObjects.genJets;

            unordered_map<int, int> matched_idx = GenJetMatching(jets, genJets, fixedGridRhoFastjetAll, 0.4, 10.);
            if (syst.Contains("PileupJetIDSF")) {
                weights.pileupIDSF = myCorr->GetPileupJetIDSF(jets, matched_idx, "loose", var);
            } else {
                weights.pileupIDSF = myCorr->GetPileupJetIDSF(jets, matched_idx, "loose", MyCorrection::variation::nom);
            }
        }

        // Pileup reweighting
        if (syst.Contains("PileupReweight")) {
            weights.pileupWeight = myCorr->GetPUWeight(event.nTrueInt(), var);
        } else {
            weights.pileupWeight = myCorr->GetPUWeight(event.nTrueInt(), MyCorrection::variation::nom);
        }
    }

    return weights;
}

void fakeV4::fillObjects(const Channel& channel,
                               const TString& ID,
                               const RecoObjects& recoObjects,
                               const WeightInfo& weights,
                               const TString& syst) {
    float totalWeight = 1.;
    if (!IsDATA) {
        totalWeight = weights.genWeight;
        totalWeight *= weights.prefireWeight;
        totalWeight *= weights.pileupWeight;
        totalWeight *= weights.topPtWeight;
        totalWeight *= weights.muonRecoSF;
        totalWeight *= weights.eleRecoSF;
        if (syst == "RequireHevayTag") totalWeight *= weights.btagSF;
        if (Run == 2) totalWeight *= weights.pileupIDSF;
    }

    // Memory optimization: For promptNormOnly systematics, only fill ZEnriched Z mass
    // This reduces histogram count for weight-only systematic variations
    // Full measurement systematics (MotherJetPt, RequireHeavyTag, etc.) need all histograms
    if (syst != "Central") {
        // Check if this systematic is promptNormOnly
        bool isPromptNormOnly = false;
        if (systHelper) {
            // Extract base syst name (remove _Up/_Down suffix)
            TString baseSyst = syst;
            baseSyst.ReplaceAll("_Up", "").ReplaceAll("_Down", "");
            SystematicHelper::SYST* systInfo = systHelper->findSystematic(baseSyst.Data());
            if (systInfo) {
                isPromptNormOnly = false;   // our fork's SYST has no promptNormOnly field;
                                            // irrelevant for Central-only running
            }
        }

        if (isPromptNormOnly) {
            if (channel != Channel::ZENRICHED) return;  // Skip Inclusive channel

            // Only fill Z mass histogram for normalization
            for (const auto& trig : triggers) {
                if (!trig.fired) continue;

                TString basePrefix = channelToString(channel) + "/" + trig.prefix + "/" + ID + "/" + syst;

                if (leptonType == LeptonType::MUON) {
                    const RVec<Muon>& muons = (ID == "loose") ? recoObjects.looseMuons : recoObjects.tightMuons;
                    float ptcorr_lead = muons[0].Pt()*(1.+max(0., muons[0].MiniPFRelIso()-0.1));
                    if (muons[0].Pt() < trig.trigSafePtCut) continue;

                    const Particle ZCand = muons[0] + muons[1];
                    FillHist(basePrefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
                } else { // ELECTRON
                    const RVec<Electron>& electrons = (ID == "loose") ? recoObjects.looseElectrons : recoObjects.tightElectrons;
                    float ptcorr_lead = electrons[0].Pt()*(1.+max(0., electrons[0].MiniPFRelIso()-0.1));
                    if (electrons[0].Pt() < trig.trigSafePtCut) continue;

                    const Particle ZCand = electrons[0] + electrons[1];
                    FillHist(basePrefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
                }
            }
            return;  // Skip all other histograms for promptNormOnly systematics
        }
        // For non-promptNormOnly systematics, fall through to full histogram filling
    }

    // Central systematic: fill all histograms
    const RVec<Jet> &jets = recoObjects.tightJets;
    const RVec<Jet> &bjets = recoObjects.bjets;
    const Particle METv = recoObjects.METv;

    // Loop over triggers that fired and fill histograms for each
    for (const auto& trig : triggers) {
        if (!trig.fired) continue;

        // Build base prefix: Channel/TrigPrefix/ID/Syst
        TString basePrefix = channelToString(channel) + "/" + trig.prefix + "/" + ID + "/" + syst;

        // Fill histograms based on channel and lepton type
        if (channel == Channel::INCLUSIVE) {
            if (leptonType == LeptonType::MUON) {
                const Muon& mu = (ID == "loose") ? recoObjects.looseMuons[0] : recoObjects.tightMuons[0];
                float ptcorr = mu.Pt()*(1.+max(0., mu.MiniPFRelIso()-0.1));

                // Check if lepton passes this trigger's pT threshold
                if (mu.Pt() < trig.trigSafePtCut) continue;

                float abseta = fabs(mu.Eta());
                float mT = TMath::Sqrt(2.*mu.Pt()*METv.Pt()*(1.-TMath::Cos(mu.DeltaPhi(METv))));
                float mTfix = TMath::Sqrt(2.*35.*METv.Pt()*(1.-TMath::Cos(mu.DeltaPhi(METv))));
                TString binName = getBinPrefix(ptcorr, abseta);

                // Fill inclusive channel histograms
                FillHist(basePrefix + "/muon/pt", mu.Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/muon/eta", mu.Eta(), totalWeight, 48, -2.4, 2.4);
                FillHist(basePrefix + "/muon/phi", mu.Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/muon/ptcorr", ptcorr, totalWeight, ptcorr_bins);
                FillHist(basePrefix + "/muon/abseta", abseta, totalWeight, abseta_bins);
                FillHist(basePrefix + "/MT", mT, totalWeight, 600, 0., 300.);
                FillHist(basePrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                FillHist(basePrefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);
                FillHist(basePrefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
                FillHist(basePrefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);

                // Fill flavour-binned histograms in Inclusive/ (MC only)
                if (!IsDATA) {
                    int jetFlavour = (ID == "loose") ? recoObjects.looseMuonJetFlavours[0] : recoObjects.tightMuonJetFlavours[0];
                    TString flavourSubdir = getFlavourSubdir(jetFlavour);
                    TString flavourPrefix = basePrefix + "/" + flavourSubdir;
                    FillHist(flavourPrefix + "/muon/pt", mu.Pt(), totalWeight, 300, 0., 300.);
                    FillHist(flavourPrefix + "/muon/eta", mu.Eta(), totalWeight, 48, -2.4, 2.4);
                    FillHist(flavourPrefix + "/muon/phi", mu.Phi(), totalWeight, 64, -3.2, 3.2);
                    FillHist(flavourPrefix + "/muon/ptcorr", ptcorr, totalWeight, ptcorr_bins);
                    FillHist(flavourPrefix + "/muon/abseta", abseta, totalWeight, abseta_bins);
                    FillHist(flavourPrefix + "/MT", mT, totalWeight, 600, 0., 300.);
                    FillHist(flavourPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                    FillHist(flavourPrefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);
                    FillHist(flavourPrefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
                    FillHist(flavourPrefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);
                }

                // Fill binned histograms: binName/Channel/TrigPrefix/ID/Syst
                TString binnedPrefix = binName + "/" + channelToString(channel) + "/" + trig.prefix + "/" + ID + "/" + syst;
                FillHist(binnedPrefix + "/muon/ptcorr", ptcorr, totalWeight, 200, 0., 200.);
                FillHist(binnedPrefix + "/muon/abseta", abseta, totalWeight, 24, 0., 2.4);
                FillHist(binnedPrefix + "/MT", mT, totalWeight, 500, 0., 500.);
                FillHist(binnedPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                FillHist(binnedPrefix + "/MET", METv.Pt(), totalWeight, 600, 0., 300.);

                // Fill subchannel - determine QCDEnriched or WEnriched
                TString subchannel = "";
                if (mT < 25. && METv.Pt() < 25.) {
                    subchannel = "QCDEnriched";
                } else if (mT > 60.) {
                    subchannel = "WEnriched";
                }

                if (subchannel != "") {
                    TString subchannelPrefix = binName + "/" + subchannel + "/" + trig.prefix + "/" + ID + "/" + syst;
                    FillHist(subchannelPrefix + "/muon/pt", mu.Pt(), totalWeight, 200, 0., 200.);
                    FillHist(subchannelPrefix + "/muon/eta", mu.Eta(), totalWeight, 48, -2.4, 2.4);
                    FillHist(subchannelPrefix + "/muon/ptcorr", ptcorr, totalWeight, 200, 0., 200.);
                    FillHist(subchannelPrefix + "/muon/abseta", abseta, totalWeight, 24, 0., 2.4);
                    FillHist(subchannelPrefix + "/MT", mT, totalWeight, 300, 0., 300.);
                    FillHist(subchannelPrefix + "/MET", METv.Pt(), totalWeight, 300, 0., 300.);

                    // Fill flavour-binned histograms (MC only)
                    if (!IsDATA) {
                        int jetFlavour = (ID == "loose") ? recoObjects.looseMuonJetFlavours[0] : recoObjects.tightMuonJetFlavours[0];
                        TString flavourSubdir = getFlavourSubdir(jetFlavour);

                        // Flavour-binned histograms under binnedPrefix
                        TString flavourBinnedPrefix = binnedPrefix + "/" + flavourSubdir;
                        FillHist(flavourBinnedPrefix + "/muon/ptcorr", ptcorr, totalWeight, 200, 0., 200.);
                        FillHist(flavourBinnedPrefix + "/muon/abseta", abseta, totalWeight, 24, 0., 2.4);
                        FillHist(flavourBinnedPrefix + "/MT", mT, totalWeight, 500, 0., 500.);
                        FillHist(flavourBinnedPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                        FillHist(flavourBinnedPrefix + "/MET", METv.Pt(), totalWeight, 600, 0., 300.);

                        // Flavour-binned histograms under subchannelPrefix
                        TString flavourSubchPrefix = subchannelPrefix + "/" + flavourSubdir;
                        FillHist(flavourSubchPrefix + "/muon/pt", mu.Pt(), totalWeight, 200, 0., 200.);
                        FillHist(flavourSubchPrefix + "/muon/eta", mu.Eta(), totalWeight, 48, -2.4, 2.4);
                        FillHist(flavourSubchPrefix + "/muon/ptcorr", ptcorr, totalWeight, 200, 0., 200.);
                        FillHist(flavourSubchPrefix + "/muon/abseta", abseta, totalWeight, 24, 0., 2.4);
                        FillHist(flavourSubchPrefix + "/MT", mT, totalWeight, 300, 0., 300.);
                        FillHist(flavourSubchPrefix + "/MET", METv.Pt(), totalWeight, 300, 0., 300.);
                    }
                }
            } else { // ELECTRON
                const Electron& el = (ID=="loose") ? recoObjects.looseElectrons[0]: recoObjects.tightElectrons[0];
                float ptcorr = el.Pt()*(1.+max(0., el.MiniPFRelIso()-0.1));

                // Check if lepton passes this trigger's pT threshold
                if (el.Pt() < trig.trigSafePtCut) continue;

                float abseta = fabs(el.scEta());
                float mT = TMath::Sqrt(2.*el.Pt()*METv.Pt()*(1.-TMath::Cos(el.DeltaPhi(METv))));
                float mTfix = TMath::Sqrt(2.*35.*METv.Pt()*(1.-TMath::Cos(el.DeltaPhi(METv))));
                TString binName = getBinPrefix(ptcorr, abseta);

                // Fill inclusive channel histograms
                FillHist(basePrefix + "/electron/pt", el.Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/electron/scEta", el.scEta(), totalWeight, 50, -2.5, 2.5);
                FillHist(basePrefix + "/electron/phi", el.Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/electron/ptcorr", ptcorr, totalWeight, ptcorr_bins);
                FillHist(basePrefix + "/electron/abseta", abseta, totalWeight, abseta_bins);
                FillHist(basePrefix + "/MT", mT, totalWeight, 600, 0., 300.);
                FillHist(basePrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                FillHist(basePrefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);
                FillHist(basePrefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
                FillHist(basePrefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);

                // Fill flavour-binned histograms in Inclusive/ (MC only)
                if (!IsDATA) {
                    int jetFlavour = (ID == "loose") ? recoObjects.looseElectronJetFlavours[0] : recoObjects.tightElectronJetFlavours[0];
                    TString flavourSubdir = getFlavourSubdir(jetFlavour);
                    TString flavourPrefix = basePrefix + "/" + flavourSubdir;
                    FillHist(flavourPrefix + "/electron/pt", el.Pt(), totalWeight, 300, 0., 300.);
                    FillHist(flavourPrefix + "/electron/scEta", el.scEta(), totalWeight, 50, -2.5, 2.5);
                    FillHist(flavourPrefix + "/electron/phi", el.Phi(), totalWeight, 64, -3.2, 3.2);
                    FillHist(flavourPrefix + "/electron/ptcorr", ptcorr, totalWeight, ptcorr_bins);
                    FillHist(flavourPrefix + "/electron/abseta", abseta, totalWeight, abseta_bins);
                    FillHist(flavourPrefix + "/MT", mT, totalWeight, 600, 0., 300.);
                    FillHist(flavourPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                    FillHist(flavourPrefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);
                    FillHist(flavourPrefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
                    FillHist(flavourPrefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);
                }

                // Fill binned histograms: binName/Channel/TrigPrefix/ID/Syst
                TString binnedPrefix = binName + "/" + channelToString(channel) + "/" + trig.prefix + "/" + ID + "/" + syst;
                FillHist(binnedPrefix + "/electron/ptcorr", ptcorr, totalWeight, ptcorr_bins);
                FillHist(binnedPrefix + "/electron/abseta", abseta, totalWeight, abseta_bins);
                FillHist(binnedPrefix + "/MT", mT, totalWeight, 600, 0., 300.);
                FillHist(binnedPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                FillHist(binnedPrefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);

                // Fill subchannel - determine QCDEnriched or WEnriched
                TString subchannel = "";
                if (mT < 25. && METv.Pt() < 25.) {
                    subchannel = "QCDEnriched";
                } else if (mT > 60.) {
                    subchannel = "WEnriched";
                }

                if (subchannel != "") {
                    TString subchannelPrefix = binName + "/" + subchannel + "/" + trig.prefix + "/" + ID + "/" + syst;
                    FillHist(subchannelPrefix + "/electron/pt", el.Pt(), totalWeight, 300, 0., 300.);
                    FillHist(subchannelPrefix + "/electron/scEta", el.scEta(), totalWeight, 50, -2.5, 2.5);
                    FillHist(subchannelPrefix + "/electron/ptcorr", ptcorr, totalWeight, 300, 0., 300.);
                    FillHist(subchannelPrefix + "/electron/abseta", abseta, totalWeight, 25, 0., 2.5);
                    FillHist(subchannelPrefix + "/MT", mT, totalWeight, 300, 0., 300.);
                    FillHist(subchannelPrefix + "/MET", METv.Pt(), totalWeight, 300, 0., 300.);

                    // Fill flavour-binned histograms (MC only)
                    if (!IsDATA) {
                        int jetFlavour = (ID == "loose") ? recoObjects.looseElectronJetFlavours[0] : recoObjects.tightElectronJetFlavours[0];
                        TString flavourSubdir = getFlavourSubdir(jetFlavour);

                        // Flavour-binned histograms under binnedPrefix
                        TString flavourBinnedPrefix = binnedPrefix + "/" + flavourSubdir;
                        FillHist(flavourBinnedPrefix + "/electron/ptcorr", ptcorr, totalWeight, ptcorr_bins);
                        FillHist(flavourBinnedPrefix + "/electron/abseta", abseta, totalWeight, abseta_bins);
                        FillHist(flavourBinnedPrefix + "/MT", mT, totalWeight, 600, 0., 300.);
                        FillHist(flavourBinnedPrefix + "/MTfix", mTfix, totalWeight, 600, 0., 300.);
                        FillHist(flavourBinnedPrefix + "/MET", METv.Pt(), totalWeight, 500, 0., 500.);

                        // Flavour-binned histograms under subchannelPrefix
                        TString flavourSubchPrefix = subchannelPrefix + "/" + flavourSubdir;
                        FillHist(flavourSubchPrefix + "/electron/pt", el.Pt(), totalWeight, 300, 0., 300.);
                        FillHist(flavourSubchPrefix + "/electron/scEta", el.scEta(), totalWeight, 50, -2.5, 2.5);
                        FillHist(flavourSubchPrefix + "/electron/ptcorr", ptcorr, totalWeight, 300, 0., 300.);
                        FillHist(flavourSubchPrefix + "/electron/abseta", abseta, totalWeight, 25, 0., 2.5);
                        FillHist(flavourSubchPrefix + "/MT", mT, totalWeight, 300, 0., 300.);
                        FillHist(flavourSubchPrefix + "/MET", METv.Pt(), totalWeight, 300, 0., 300.);
                    }
                }
            }
        } else if (channel == Channel::ZENRICHED) {
            if (leptonType == LeptonType::MUON) {
                const RVec<Muon>& muons = (ID == "loose") ? recoObjects.looseMuons : recoObjects.tightMuons;

                // For ZENRICHED, check if leading muon passes trigger pT threshold
                float ptcorr_lead = muons[0].Pt()*(1.+max(0., muons[0].MiniPFRelIso()-0.1));
                if (muons[0].Pt() < trig.trigSafePtCut) continue;

                const Particle ZCand = muons[0] + muons[1];
                FillHist(basePrefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
                FillHist(basePrefix + "/ZCand/pt", ZCand.Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/ZCand/eta", ZCand.Eta(), totalWeight, 100, -5., 5.);
                FillHist(basePrefix + "/ZCand/phi", ZCand.Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/muons/1/pt", muons[0].Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/muons/1/eta", muons[0].Eta(), totalWeight, 48, -2.4, 2.4);
                FillHist(basePrefix + "/muons/1/phi", muons[0].Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/muons/2/pt", muons[1].Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/muons/2/eta", muons[1].Eta(), totalWeight, 48, -2.4, 2.4);
                FillHist(basePrefix + "/muons/2/phi", muons[1].Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
                FillHist(basePrefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);

                // Fill flavour-binned histograms for ZENRICHED (MC only)
                if (!IsDATA) {
                    const RVec<int>& jetFlavours = (ID == "loose") ? recoObjects.looseMuonJetFlavours : recoObjects.tightMuonJetFlavours;
                    for (size_t i = 0; i < muons.size() && i < jetFlavours.size(); i++) {
                        TString flavourSubdir = getFlavourSubdir(jetFlavours[i]);
                        TString flavourPrefix = basePrefix + "/" + flavourSubdir;
                        FillHist(flavourPrefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
                        FillHist(flavourPrefix + Form("/muons/%zu/pt", i+1), muons[i].Pt(), totalWeight, 300, 0., 300.);
                        FillHist(flavourPrefix + Form("/muons/%zu/eta", i+1), muons[i].Eta(), totalWeight, 48, -2.4, 2.4);
                    }
                }
            } else { // ELECTRON
                const RVec<Electron>& electrons = (ID == "loose") ? recoObjects.looseElectrons : recoObjects.tightElectrons;

                // For ZENRICHED, check if leading electron passes trigger pT threshold
                float ptcorr_lead = electrons[0].Pt()*(1.+max(0., electrons[0].MiniPFRelIso()-0.1));
                if (electrons[0].Pt() < trig.trigSafePtCut) continue;

                const Particle ZCand = electrons[0] + electrons[1];
                FillHist(basePrefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
                FillHist(basePrefix + "/ZCand/pt", ZCand.Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/ZCand/eta", ZCand.Eta(), totalWeight, 100, -5., 5.);
                FillHist(basePrefix + "/ZCand/phi", ZCand.Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/electrons/1/pt", electrons[0].Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/electrons/1/scEta", electrons[0].scEta(), totalWeight, 50, -2.5, 2.5);
                FillHist(basePrefix + "/electrons/1/phi", electrons[0].Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/electrons/2/pt", electrons[1].Pt(), totalWeight, 300, 0., 300.);
                FillHist(basePrefix + "/electrons/2/scEta", electrons[1].scEta(), totalWeight, 50, -2.5, 2.5);
                FillHist(basePrefix + "/electrons/2/phi", electrons[1].Phi(), totalWeight, 64, -3.2, 3.2);
                FillHist(basePrefix + "/nJets", jets.size(), totalWeight, 10, 0., 10.);
                FillHist(basePrefix + "/nBJets", bjets.size(), totalWeight, 5, 0., 5.);

                // Fill flavour-binned histograms for ZENRICHED (MC only)
                if (!IsDATA) {
                    const RVec<int>& jetFlavours = (ID == "loose") ? recoObjects.looseElectronJetFlavours : recoObjects.tightElectronJetFlavours;
                    for (size_t i = 0; i < electrons.size() && i < jetFlavours.size(); i++) {
                        TString flavourSubdir = getFlavourSubdir(jetFlavours[i]);
                        TString flavourPrefix = basePrefix + "/" + flavourSubdir;
                        FillHist(flavourPrefix + "/ZCand/mass", ZCand.M(), totalWeight, 40, 75., 115.);
                        FillHist(flavourPrefix + Form("/electrons/%zu/pt", i+1), electrons[i].Pt(), totalWeight, 300, 0., 300.);
                        FillHist(flavourPrefix + Form("/electrons/%zu/scEta", i+1), electrons[i].scEta(), totalWeight, 50, -2.5, 2.5);
                    }
                }
            }
        }
    }  // end trigger loop
}

TString fakeV4::getBinPrefix(const double ptcorr, const double abseta) {
    int ptcorr_idx = -1;
    int abseta_idx = -1;
    for (int i = 0; i < ptcorr_bins.size()-1; i++) {
        if (ptcorr_bins[i] <= ptcorr && ptcorr < ptcorr_bins[i+1]) {
            ptcorr_idx = i;
            break;
        }
    }
    if (ptcorr_idx == -1) ptcorr_idx = ptcorr_bins.size()-2;

    for (int i = 0; i < abseta_bins.size()-1; i++) {
        if (abseta_bins[i] <= abseta && abseta < abseta_bins[i+1]) {
            abseta_idx = i;
            break;
        }
    }

    TString etaBin;
    if (abseta_idx == 0) etaBin = "EB1";
    else if (abseta_idx == 1) etaBin = "EB2";
    else if (abseta_idx == 2) etaBin = "EE";
    else etaBin = "EE";

    return TString::Format("ptcorr_%dto%d_%s",
                          static_cast<int>(ptcorr_bins[ptcorr_idx]),
                          static_cast<int>(ptcorr_bins[ptcorr_idx+1]),
                          etaBin.Data());
}

float fakeV4::getJetPtCut(const TString& selection) {
    if (selection.Contains("MotherJetPt_Up"))
        return 60.0;
    else if (selection.Contains("MotherJetPt_Down"))
        return (leptonType == LeptonType::MUON) ? 20.0 : 30.0;
    else
        return 40.0;
}

template<typename T>
int fakeV4::getMotherJetFlavour(const T& lep, const RVec<Jet>& allJets) {
    if (IsDATA) return -999;  // Not applicable for data

    short jetIdx = lep.JetIdx();
    if (jetIdx < 0) return -1;  // No associated jet

    for (const auto& jet : allJets) {
        if (jet.OriginalIndex() == jetIdx) {
            return (jet.genJetIdx() < 0) ? -1 : jet.hadronFlavour();
        }
    }
    return -1;  // Jet not found
}

// Explicit template instantiations
template int fakeV4::getMotherJetFlavour<Muon>(const Muon&, const RVec<Jet>&);
template int fakeV4::getMotherJetFlavour<Electron>(const Electron&, const RVec<Jet>&);

TString fakeV4::getFlavourSubdir(int flavour) {
    if (flavour == 5) return "bjet";
    if (flavour == 4) return "cjet";
    if (flavour == 0) return "ljet";
    return "pujet";  // 0, -1, or any other value
}

void fakeV4::fillCutflow(CutStage stage, const Channel& channel, const TString& ID, float weight, const TString& syst) {
    if (syst != "Central") return;
    TString channelStr = channelToString(channel);
    if (channelStr == "NONE") channelStr = "PreSel";

    int cutIndex = static_cast<int>(stage);
    FillHist(Form("%s/%s/%s/cutflow", channelStr.Data(), ID.Data(), syst.Data()), cutIndex, weight, 9, 0., 9.);
}


//================= WP helpers (official values, verified) =================
bool fakeV4::PassLooseMuon(const Muon &mu) const {
    // Muon::Pass_HcToWALooseRun2 / ...Run3 in the official DataFormats
    if (!mu.isPOGMediumId()) return false;
    if (!(fabs(mu.dZ()) < 0.1)) return false;
    if (!(mu.SIP3D() < (Run == 2 ? 5. : 8.))) return false;
    if (!(mu.TkRelIso() < 0.4)) return false;
    if (!(mu.MiniPFRelIso() < (Run == 2 ? 0.6 : 0.4))) return false;
    return true;
}

bool fakeV4::PassLooseElectron(const Electron &el) const {
    // Electron::Pass_HcToWALooseRun2 / ...Run3 in the official DataFormats
    if (!el.Pass_HcToWABaseline()) return false;
    if (!(el.SIP3D() < (Run == 2 ? 8. : 6.))) return false;
    if (!(el.MiniPFRelIso() < 0.4)) return false;
    const float cutIB = (Run == 2 ? 0.985 : 0.8);
    const float cutOB = (Run == 2 ? 0.96  : 0.5);
    const float cutEC = (Run == 2 ? 0.85  : -0.8);
    float mvaCut = cutEC;
    if (el.etaRegion() == Electron::ETAREGION::IB)      mvaCut = cutIB;
    else if (el.etaRegion() == Electron::ETAREGION::OB) mvaCut = cutOB;
    if (!(el.MvaNoIso() > mvaCut)) return false;
    return true;
}

RVec<Muon> fakeV4::PickMuons(const RVec<Muon> &muons, const TString &wp,
                             float ptmin, float fetamax) const {
    RVec<Muon> out;
    for (const auto &mu : muons) {
        if (!(mu.Pt() > ptmin)) continue;
        if (!(fabs(mu.Eta()) < fetamax)) continue;
        if (wp == "tight") { if (!mu.PassID("HcToWATight")) continue; }
        else               { if (!PassLooseMuon(mu)) continue; }
        out.push_back(mu);
    }
    return out;
}

RVec<Electron> fakeV4::PickElectrons(const RVec<Electron> &els, const TString &wp,
                                     float ptmin, float fetamax, bool vetoHEM) const {
    RVec<Electron> out;
    for (const auto &el : els) {
        if (!(el.Pt() > ptmin)) continue;
        if (!(fabs(el.Eta()) < fetamax)) continue;
        if (wp == "tight") { if (!el.PassID("HcToWATight")) continue; }
        else               { if (!PassLooseElectron(el)) continue; }
        if (vetoHEM && IsHEMElectron(el)) continue;
        out.push_back(el);
    }
    return out;
}

//================= ported from fake.cc =================
TString fakeV4::LeptonTypeToString(int leptonType) const {
    return (leptonType > 0) ? "prompt" : "fake";
}

//==============================================================
// Source jet parton flavor 분류
//   hadronFlavour (ghost B/C hadron matching) 로 b/c 를 먼저 잡고,
//   light (hadronFlavour==0) 는 |partonFlavour| 로 g/s/u/d 세분화한다.
//   partonFlavour: 21=gluon, 1=d, 2=u, 3=s (부호 = quark/antiquark)
//   매칭 실패(partonFlavour==0 등)는 unmatched
//   isPileup (대응하는 gen jet 이 없는 reco jet) 은 flavour 보다 먼저 본다:
//   공식 코드 MeasFakeRateV4::getMotherJetFlavour 가 genJetIdx<0 이면 hadronFlavour
//   를 보지 않고 pujet 을 돌려주므로 같은 우선순위를 쓴다.
//==============================================================
TString fakeV4::FlavorTag(int partonFlavour, int hadronFlavour, bool isPileup) const {
    if (isPileup) return "pileup";
    if (hadronFlavour == 5) return "b";
    if (hadronFlavour == 4) return "c";
    const int ap = abs(partonFlavour);
    if (ap == 21) return "g";
    if (ap == 3)  return "s";
    if (ap == 2)  return "u";
    if (ap == 1)  return "d";
    return "unmatched";
}

//==== muon 이 clustering 된 PF jet 을 NanoAOD 의 jetIdx 포인터로 직접 찾는다.
//==== dR<0.4 최근접 매칭은 muon 이 실제로 들어있지 않은 이웃 jet 을 집을 수 있어서
//==== 공식 코드(MeasFakeRateV4::getMotherJetFlavour)와 같이 jetIdx 를 쓴다.
//==== rawJets 의 OriginalIndex 가 곧 NanoAOD Jet 인덱스라 그대로 대조하면 된다.
TString fakeV4::RecoJetFlavor(const Muon &mu) const {
    const short jetIdx = mu.JetIdx();
    if (jetIdx < 0) return "unmatched";        // 연결된 jet 자체가 없는 muon
    for (const auto &jet : rawJets) {
        if (jet.OriginalIndex() != jetIdx) continue;
        return FlavorTag(jet.partonFlavour(), jet.hadronFlavour(),
                         jet.genJetIdx() < 0);
    }
    return "unmatched";                        // jetIdx 는 있으나 컬렉션에서 못 찾음
}

//==== muon 을 가장 가까운 gen jet 에 dR<0.4 매칭 (source jet 의 진짜 flavor)
TString fakeV4::GenJetFlavor(const Muon &mu) const {
    float bestDR = 0.4f;
    int bestParton = 0, bestHadron = 0;
    bool matched = false;
    for (const auto &gjet : genJets) {
        const float dr = gjet.DeltaR(mu);
        if (dr < bestDR) {
            bestDR = dr;
            bestParton = gjet.partonFlavour();
            bestHadron = gjet.hadronFlavour();
            matched = true;
        }
    }
    return matched ? FlavorTag(bestParton, bestHadron, false) : "unmatched";
}

//==============================================================
// FlavorTag 문자열 -> 5-class 정수 코드 (SetBranch 는 문자열 branch 를 못 만든다)
//   0=b 1=c 2=uds(u+d+s) 3=g 4=pileup, -1=unmatched
//==============================================================
int fakeV4::FlavorCode(const TString &flav) const {
    if (flav == "b") return 0;
    if (flav == "c") return 1;
    if (flav == "s" || flav == "u" || flav == "d") return 2;
    if (flav == "g") return 3;
    if (flav == "pileup") return 4;
    return -1;
}

//==============================================================
// flavor BDT 학습용 flat ntuple (userflag MakeTree, tree 이름 "mu")
//
// loose muon 이 정확히 1개인 측정 영역 이벤트마다 한 줄. FR 을 재는 모집단과
// 같은 population 이어야 하므로 selection 은 measureFakeRate 와 동일하되,
// away jet 컷만 가장 느슨한 30 GeV 로 열어 두고 hasAwayJet30/40/60 을
// 컬럼으로 남긴다 (offline 에서 Central=40 / syst=30,60 를 그대로 복원 가능).
//
// MET / MT 컷도 같은 이유로 **offline** 으로 옮겼다 (2026-08-11). 공식
// (MeasFakeRateV4) MC FR 은 Inclusive 영역 — away jet 만, MET/MT 컷 없음 — 에서
// 재므로, 컷을 트리에 박아 두면 그 맵을 재현할 수 없다. ev_met / ev_mt 컬럼이
// 그대로 있으므로 MeasReg (MET<25 && MT<25) 는 offline 에서 복원된다.
//
// 주의 1: 이 함수는 trigger path loop **밖**에서 불러야 한다. loop 안에서 쓰면
//         Mu8/Mu17 을 모두 통과한 muon 이 두 줄로 중복되어 train/test 폴드에
//         나뉘어 들어가고 검증 성능이 부풀려진다. 경로 정보는 pass_mu8 /
//         pass_mu17 / w_mu8 / w_mu17 컬럼으로 남긴다.
// 주의 2: AnalyzerCore::SetBranch 는 deque 에 값을 push 하고 그 원소의 주소를
//         TBranch 에 넘기는데 FillTrees 가 매번 clear()+shrink_to_fit() 한다.
//         따라서 **모든 branch 를 매 줄마다 빠짐없이 SetBranch 해야 한다.**
//         조건부로 건너뛰면 그 branch 는 해제된 메모리를 가리키게 된다.
//         Run2/Run3 전용 변수와 matched jet 이 없는 경우는 sentinel 로 채운다.

void fakeV4::fillTreeRow(Event &ev, float evtSF) {

    if (!MakeTree) return;

    //==== 측정 영역과 동일: loose muon 1개 + electron 0개
    if (looseMuons.size() != 1) return;
    if (looseElectrons.size() != 0) return;

    const Muon &mu = looseMuons.at(0);

    //==== cone-corrected pT / MT (measureFakeRate 와 같은 정의)
    const float ptCorr = mu.Pt() * (1. + max(0.f, mu.MiniPFRelIso() - treecuts.tight_miniiso_max));
    const float MT = sqrt(2. * mu.Pt() * METv.Pt() * (1. - cos(mu.DeltaPhi(METv))));

    //==== 트리거 경로별 통과 여부 (paths 는 userflag 로 잘릴 수 있어 직접 정의)
    const bool passMu8  = ev.PassTrigger("HLT_Mu8_TrkIsoVVL")
                          && mu.Pt() >= 10. && ptCorr >= 10.;
    const bool passMu17 = ev.PassTrigger("HLT_Mu17_TrkIsoVVL")
                          && mu.Pt() >= 20. && ptCorr >= 30.;
    if (!passMu8 && !passMu17) return;

    //==== away jet: 가장 느슨한 30 GeV 로 열어 두고 40/60 은 플래그로 남긴다
    int nJet25 = 0, nJet30 = 0, nJet40 = 0, nJet60 = 0;
    bool hasAway30 = false, hasAway40 = false, hasAway60 = false;
    float awayJetPt = -999., awayJetDR = -999., awayJetDPhi = -999., HT = 0.;
    for (const auto &jet : jets) {
        HT += jet.Pt();
        nJet25++;
        const bool away = (jet.DeltaR(mu) > treecuts.awayjet_dr);
        if (jet.Pt() > 30.) { nJet30++; if (away) hasAway30 = true; }
        if (jet.Pt() > 40.) {
            nJet40++;
            if (away) {
                hasAway40 = true;
                if (awayJetPt < 0.) {   // jets 는 pT 내림차순 → 첫 번째가 leading
                    awayJetPt   = jet.Pt();
                    awayJetDR   = jet.DeltaR(mu);
                    awayJetDPhi = fabs(jet.DeltaPhi(mu));
                }
            }
        }
        if (jet.Pt() > 60.) { nJet60++; if (away) hasAway60 = true; }
    }
    if (!hasAway30) return;

    //==== 측정 영역의 MET / MT 컷 (W/Z prompt 오염 억제) 은 여기서 걸지 않는다.
    //==== ev_met / ev_mt 를 그대로 남겨 offline 에서 자르는 쪽이 넓다:
    //====   Inclusive (공식 MC FR)      = 컷 없음
    //====   MeasReg   (data FR / bkg)   = ev_met < 25 && ev_mt < 25

    //==== muon 의 source jet: lepton cleaning 이전 컬렉션에서 NanoAOD jetIdx 로 찾는다
    const short jetIdx = mu.JetIdx();
    const Jet *mj = nullptr;
    for (const auto &jet : rawJets) {
        if (jet.OriginalIndex() != jetIdx) continue;
        mj = &jet;
        break;
    }

    //==== matched jet 변수 (없으면 sentinel)
    const float FSENT = -999.;
    const int   ISENT = -99;
    float jetPt = FSENT, jetRawPt = FSENT, jetEta = FSENT, jetMass = FSENT;
    float jetArea = FSENT, jetDR = FSENT, jetPtRatio = FSENT, jetPtRel = FSENT;
    float jetChHEF = FSENT, jetNeHEF = FSENT, jetChEmEF = FSENT, jetNeEmEF = FSENT;
    float jetMuEF = FSENT, jetMuonSubtr = FSENT;
    float jetDeepB = FSENT, jetDeepCvB = FSENT, jetDeepCvL = FSENT, jetDeepQG = FSENT;
    float jetQGL = FSENT, jetPuIdDisc = FSENT;
    int   jetNConst = ISENT, jetNMuons = ISENT, jetNElectrons = ISENT;
    int   jetChMult = ISENT, jetNeMult = ISENT;   // Run3 전용 (아래 참조)
    bool  jetTightId = false, jetLoosePuId = false;
    if (mj) {
        jetPt        = mj->Pt();
        jetRawPt     = mj->GetRawPt();
        jetEta       = mj->Eta();
        jetMass      = mj->M();
        jetDR        = mj->DeltaR(mu);
        jetPtRatio   = (mj->Pt() > 0.) ? mu.Pt() / mj->Pt() : FSENT;
        //==== pTrel: (jet - muon) 축에 대한 muon 운동량의 수직 성분
        const TVector3 axis = mj->Vect() - mu.Vect();
        if (axis.Mag() > 1e-6) jetPtRel = mu.Vect().Perp(axis);
        jetChHEF     = mj->chHEF();
        jetNeHEF     = mj->neHEF();
        jetChEmEF    = mj->chEmEF();
        jetNeEmEF    = mj->neEmEF();
        jetMuEF      = mj->muEF();
        jetNConst    = mj->nConstituents();
        jetNMuons    = mj->nMuons();
        jetNElectrons= mj->nElectrons();
        jetDeepB     = mj->GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        const pair<float,float> cvals = mj->GetCTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        jetDeepCvB   = cvals.first;
        jetDeepCvL   = cvals.second;
        jetDeepQG    = mj->GetQvGTaggerResult(JetTagging::JetFlavTagger::DeepJet);
        jetTightId   = mj->Pass_tightJetID();
        jetLoosePuId = mj->PassID("loosePuId");   // Run3 에서는 항상 참
        //==== charged / neutral multiplicity 는 Run3 (NanoAODv12+) 전용이다.
        //     Run2 v9 에는 Jet_chMultiplicity / Jet_neMultiplicity branch 자체가
        //     없고, AnalyzerCore 의 Run==3 가지에서만 SetHadronMultiplicities 가
        //     불린다. 게다가 Jet 의 생성자는 nConstituents 등과 달리 이 두 멤버를
        //     -999 로 초기화하지 않으므로 Run2 에서 읽으면 쓰레기 값이 나온다.
        //     따라서 반드시 Run==3 로 막고 나머지는 sentinel 로 남긴다
        //     (jet_qgl / jet_puIdDisc 의 Run2 전용 처리와 정확히 대칭).
        if (Run == 3) {
            jetChMult = int(mj->chMultiplicity());
            jetNeMult = int(mj->neMultiplicity());
        }
        if (jetIdx >= 0 && jetIdx < nJet) {
            jetArea      = Jet_area[jetIdx];
            jetMuonSubtr = Jet_muonSubtrFactor[jetIdx];
            if (Run == 2) {
                jetQGL      = mj->qgl();            // Run2 전용
                jetPuIdDisc = Jet_puIdDisc[jetIdx]; // Run2 전용
            }
        }
    }

    //==== source jet 안의 secondary vertex (IVF)
    //   b/c 분리의 고전적인 변수다: b hadron 은 M ~ 5 GeV, c hadron 은 ~ 2 GeV 이고
    //   b 쪽이 track 수도 많고 더 멀리 날아간다. DeepJet 이 내부적으로 SV 를 쓰지만
    //   출력은 3개 숫자로 압축돼 있어 raw SV 변수는 그 위에 정보를 더한다.
    //   jet 축 기준 dR < 0.4 로 매칭하고, 대표 SV 는 dlenSig 가 가장 큰 것
    //   (= 가장 유의미하게 변위된 vertex, HF 판별의 표준 선택) 을 쓴다.
    int   svN = 0, svNTracks = ISENT;
    float svMass = FSENT, svPt = FSENT, svPtRatio = FSENT, svDlenSig = FSENT;
    float svDxySig = FSENT, svPAngle = FSENT, svChi2 = FSENT;
    float svDRJet = FSENT, svDRMu = FSENT, svMassSum = FSENT;
    if (mj) {
        const float jetPhi = mj->Phi();
        float bestDlenSig = -1.;
        int   bestIdx = -1;
        float massSum = 0.;
        for (int i = 0; i < nSV; i++) {
            float dphi = SV_phi[i] - jetPhi;
            while (dphi >  M_PI) dphi -= 2. * M_PI;
            while (dphi < -M_PI) dphi += 2. * M_PI;
            const float deta = SV_eta[i] - jetEta;
            if (sqrt(deta * deta + dphi * dphi) > 0.4) continue;
            svN++;
            massSum += SV_mass[i];
            if (SV_dlenSig[i] > bestDlenSig) { bestDlenSig = SV_dlenSig[i]; bestIdx = i; }
        }
        svMassSum = massSum;
        if (bestIdx >= 0) {
            svMass     = SV_mass[bestIdx];
            svPt       = SV_pt[bestIdx];
            svPtRatio  = (jetPt > 0.) ? SV_pt[bestIdx] / jetPt : FSENT;
            svDlenSig  = SV_dlenSig[bestIdx];
            svDxySig   = SV_dxySig[bestIdx];
            svPAngle   = SV_pAngle[bestIdx];
            svChi2     = SV_chi2[bestIdx];
            svNTracks  = (int)SV_ntracks[bestIdx];
            float dphiJ = SV_phi[bestIdx] - jetPhi;
            while (dphiJ >  M_PI) dphiJ -= 2. * M_PI;
            while (dphiJ < -M_PI) dphiJ += 2. * M_PI;
            const float detaJ = SV_eta[bestIdx] - jetEta;
            svDRJet = sqrt(detaJ * detaJ + dphiJ * dphiJ);
            //==== muon 은 이 SV 에서 나온 것이므로 dR(SV, mu) 도 직접적인 handle
            float dphiM = SV_phi[bestIdx] - mu.Phi();
            while (dphiM >  M_PI) dphiM -= 2. * M_PI;
            while (dphiM < -M_PI) dphiM += 2. * M_PI;
            const float detaM = SV_eta[bestIdx] - mu.Eta();
            svDRMu = sqrt(detaM * detaM + dphiM * dphiM);
        }
    }

    //==== label (MC only). DATA 는 jet flavour 자체가 없으므로 -2
    int flavor = -2, flavorFine = -2, genFlavor = -2, leptonType = ISENT;
    if (!IsDATA) {
        const TString rf = RecoJetFlavor(mu);
        const TString gf = GenJetFlavor(mu);
        flavor    = FlavorCode(rf);
        genFlavor = FlavorCode(gf);
        //==== u/d/s 를 나중에 다시 쪼갤 수 있도록 세분 코드도 남긴다
        flavorFine = (rf == "b") ? 0 : (rf == "c") ? 1 : (rf == "s") ? 2
                   : (rf == "u") ? 3 : (rf == "d") ? 4 : (rf == "g") ? 5
                   : (rf == "pileup") ? 6 : -1;
        leptonType = GetLeptonType(mu, gens);
    }
    const bool isPrompt = (!IsDATA && leptonType > 0);

    //==== weight: 경로마다 trigger lumi 가 달라 하나로 합칠 수 없으므로 따로 남긴다
    float commonW = 1.;
    if (!IsDATA) {
        commonW = MCweight() * GetL1PrefireWeight()
                * myCorr->GetPUWeight(ev.nTrueInt()) * evtSF
                * myCorr->GetMuonRECOSF(mu);
    }
    const float wMu8  = passMu8  ? commonW * (IsDATA ? 1.f : ev.GetTriggerLumi("HLT_Mu8_TrkIsoVVL"))  : 0.f;
    const float wMu17 = passMu17 ? commonW * (IsDATA ? 1.f : ev.GetTriggerLumi("HLT_Mu17_TrkIsoVVL")) : 0.f;

    //==== 아래 SetBranch 블록은 반드시 무조건부 직선 코드여야 한다 (주의 2 참고)
    const TString t = "mu";

    //---- label / bookkeeping
    SetBranch(t, "flavor",      flavor);
    SetBranch(t, "flavor_fine", flavorFine);
    SetBranch(t, "genflavor",   genFlavor);
    SetBranch(t, "isTight",     mu.PassID("HcToWATight"));
    SetBranch(t, "isPrompt",    isPrompt);
    SetBranch(t, "leptonType",  leptonType);
    SetBranch(t, "isData",      (bool)IsDATA);
    SetBranch(t, "pass_mu8",    passMu8);
    SetBranch(t, "pass_mu17",   passMu17);
    SetBranch(t, "w_mu8",       wMu8);
    SetBranch(t, "w_mu17",      wMu17);
    SetBranch(t, "mcweight",    IsDATA ? 1.f : (float)MCweight());
    SetBranch(t, "run",         ev.run());
    SetBranch(t, "lumi",        ev.lumi());
    SetBranch(t, "evt",         ev.event());

    //---- muon kinematics / IP / isolation (전부 data 에서도 얻을 수 있는 양)
    SetBranch(t, "mu_pt",          (float)mu.Pt());
    SetBranch(t, "mu_conept",      ptCorr);
    SetBranch(t, "mu_eta",         (float)mu.Eta());
    SetBranch(t, "mu_abseta",      (float)fabs(mu.Eta()));
    SetBranch(t, "mu_phi",         (float)mu.Phi());
    SetBranch(t, "mu_charge",      (int)mu.Charge());
    SetBranch(t, "mu_dxy",         mu.dXY());
    SetBranch(t, "mu_dxyerr",      mu.dXYerr());
    SetBranch(t, "mu_dxysig",      mu.dXY() / max(mu.dXYerr(), 1e-6f));
    SetBranch(t, "mu_dz",          mu.dZ());
    SetBranch(t, "mu_dzerr",       mu.dZerr());
    SetBranch(t, "mu_dzsig",       mu.dZ() / max(mu.dZerr(), 1e-6f));
    SetBranch(t, "mu_ip3d",        mu.IP3D());
    SetBranch(t, "mu_sip3d",       mu.SIP3D());
    SetBranch(t, "mu_miniiso",     mu.MiniPFRelIso());
    SetBranch(t, "mu_miniiso_chg", mu.MiniPFRelIsoChg());
    SetBranch(t, "mu_miniiso_neu", mu.MiniPFRelIso() - mu.MiniPFRelIsoChg());
    SetBranch(t, "mu_pfiso03",     mu.PfRelIso03());
    SetBranch(t, "mu_pfiso03_chg", mu.PfRelIso03Chg());
    SetBranch(t, "mu_pfiso04",     mu.PfRelIso04());
    SetBranch(t, "mu_tkreliso",    mu.TkRelIso());
    SetBranch(t, "mu_ptErr",       mu.PtErr());
    SetBranch(t, "mu_ptErrRel",    mu.PtErr() / max((float)mu.Pt(), 1e-6f));
    SetBranch(t, "mu_tunepRelPt",  mu.TunepRelPt());
    SetBranch(t, "mu_segmentComp", mu.SegmentComp());
    SetBranch(t, "mu_nStations",   mu.NStations());
    SetBranch(t, "mu_nTrkLayers",  mu.nTrackerLayers());
    SetBranch(t, "mu_mvaTTH",      mu.MvaTTH());
    SetBranch(t, "mu_mvaLowPt",    mu.MvaLowPt());
    SetBranch(t, "mu_softMva",     mu.SoftMva());
    SetBranch(t, "mu_tightCharge", mu.TightCharge());
    SetBranch(t, "mu_isGlobal",    mu.isGlobal());
    SetBranch(t, "mu_isTracker",   mu.isTracker());
    SetBranch(t, "mu_isStandalone",mu.isStandalone());
    SetBranch(t, "mu_isPFcand",    mu.isPFcand());
    SetBranch(t, "mu_highPurity",  mu.highPurity());
    SetBranch(t, "mu_pogMedium",   mu.isPOGMediumId());
    SetBranch(t, "mu_pogTight",    mu.isPOGTightId());

    //---- muon <-> jet 상관 변수 (jetIdx 가 없으면 NanoAOD 가 이미 fallback 값을 준다)
    SetBranch(t, "mu_jetPtRelv2",     mu.JetPtRelv2());
    SetBranch(t, "mu_jetRelIso",      mu.JetRelIso());
    SetBranch(t, "mu_jetNDauCharged", mu.JetNDauCharged());
    SetBranch(t, "mu_hasJet",         (bool)(mj != nullptr));

    //---- source jet 변수 (없으면 sentinel)
    SetBranch(t, "jet_pt",          jetPt);
    SetBranch(t, "jet_rawpt",       jetRawPt);
    SetBranch(t, "jet_eta",         jetEta);
    SetBranch(t, "jet_mass",        jetMass);
    SetBranch(t, "jet_area",        jetArea);
    SetBranch(t, "jet_dr_mu",       jetDR);
    SetBranch(t, "jet_ptratio",     jetPtRatio);
    SetBranch(t, "jet_ptrel",       jetPtRel);
    SetBranch(t, "jet_chHEF",       jetChHEF);
    SetBranch(t, "jet_neHEF",       jetNeHEF);
    SetBranch(t, "jet_chEmEF",      jetChEmEF);
    SetBranch(t, "jet_neEmEF",      jetNeEmEF);
    SetBranch(t, "jet_muEF",        jetMuEF);
    SetBranch(t, "jet_muonSubtr",   jetMuonSubtr);
    SetBranch(t, "jet_nconst",      jetNConst);
    SetBranch(t, "jet_chmult",      jetChMult);     // Run3 전용, Run2 는 sentinel
    SetBranch(t, "jet_nemult",      jetNeMult);     // Run3 전용, Run2 는 sentinel
    SetBranch(t, "jet_nmuons",      jetNMuons);
    SetBranch(t, "jet_nelectrons",  jetNElectrons);
    SetBranch(t, "jet_deepjet_b",   jetDeepB);
    SetBranch(t, "jet_deepjet_cvb", jetDeepCvB);
    SetBranch(t, "jet_deepjet_cvl", jetDeepCvL);
    SetBranch(t, "jet_deepjet_qg",  jetDeepQG);
    SetBranch(t, "jet_qgl",         jetQGL);        // Run2 전용, Run3 는 sentinel
    SetBranch(t, "jet_puIdDisc",    jetPuIdDisc);   // Run2 전용, Run3 는 sentinel
    SetBranch(t, "jet_tightId",     jetTightId);
    SetBranch(t, "jet_loosePuId",   jetLoosePuId);

    //---- source jet 안의 secondary vertex (매칭 SV 가 없으면 sentinel)
    SetBranch(t, "sv_n",         svN);
    SetBranch(t, "sv_mass",      svMass);
    SetBranch(t, "sv_masssum",   svMassSum);
    SetBranch(t, "sv_pt",        svPt);
    SetBranch(t, "sv_ptratio",   svPtRatio);
    SetBranch(t, "sv_dlensig",   svDlenSig);
    SetBranch(t, "sv_dxysig",    svDxySig);
    SetBranch(t, "sv_pangle",    svPAngle);
    SetBranch(t, "sv_chi2",      svChi2);
    SetBranch(t, "sv_ntracks",   svNTracks);
    SetBranch(t, "sv_dr_jet",    svDRJet);
    SetBranch(t, "sv_dr_mu",     svDRMu);

    //---- event / away jet context
    SetBranch(t, "ev_met",          (float)METv.Pt());
    SetBranch(t, "ev_mt",           MT);
    SetBranch(t, "ev_dphi_mu_met",  (float)fabs(mu.DeltaPhi(METv)));
    SetBranch(t, "ev_ht",           HT);
    SetBranch(t, "ev_njet25",       nJet25);
    SetBranch(t, "ev_njet30",       nJet30);
    SetBranch(t, "ev_njet40",       nJet40);
    SetBranch(t, "ev_njet60",       nJet60);
    SetBranch(t, "ev_hasAway30",    hasAway30);
    SetBranch(t, "ev_hasAway40",    hasAway40);
    SetBranch(t, "ev_hasAway60",    hasAway60);
    SetBranch(t, "ev_awayjet_pt",   awayJetPt);
    SetBranch(t, "ev_awayjet_dr",   awayJetDR);
    SetBranch(t, "ev_awayjet_dphi", awayJetDPhi);
    SetBranch(t, "ev_nPV",          ev.nPV());
    SetBranch(t, "ev_nPVsGood",     ev.nPVsGood());
    SetBranch(t, "ev_rho",          fixedGridRhoFastjetAll);
    SetBranch(t, "ev_nTrueInt",     IsDATA ? -1.f : ev.nTrueInt());

    FillTrees(t);
}
