#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2016preVFP -n 40 --tag QCD_1
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2016postVFP -n 40 --tag QCD_1

#SKNano.py -a fake -i DoubleMuon -e 2016preVFP -n 40 --tag data_syst
#SKNano.py -a fake -i DoubleMuon -e 2016postVFP -n 40 --tag data_syst

#SKNano.py -a fake -i DYJets -e 2016preVFP -n 40 --tag data_syst
#SKNano.py -a fake -i DYJets10to50 -e 2016preVFP -n 40 --tag data_syst
#SKNano.py -a fake -i WJets -e 2016preVFP -n 40 --tag data_syst

#SKNano.py -a fake -i DYJets -e 2016postVFP -n 40 --tag data_syst
#SKNano.py -a fake -i DYJets10to50 -e 2016postVFP -n 40 --tag data_syst
#SKNano.py -a fake -i WJets -e 2016postVFP -n 40 --tag data_syst

#SKNano.py -a elecfake -i DoubleEG -e 2016preVFP -n 40 --tag 1
#SKNano.py -a elecfake -i DoubleEG -e 2016postVFP -n 40 --tag 1

#SKNano.py -a elecfake -i DYJets -e 2016preVFP -n 40 --tag 1
#SKNano.py -a elecfake -i DYJets10to50 -e 2016preVFP -n 40 --tag 1
#SKNano.py -a elecfake -i WJets -e 2016preVFP -n 40 --tag 1

#SKNano.py -a elecfake -i DYJets -e 2016postVFP -n 40 --tag 1
#SKNano.py -a elecfake -i DYJets10to50 -e 2016postVFP -n 40 --tag 1
#SKNano.py -a elecfake -i WJets -e 2016postVFP -n 40 --tag 1

#SKNano.py -a elecfake -i 'QCD_Pt-*EMEnriched*' -e 2016preVFP -n 40 --tag QCE_1
#SKNano.py -a elecfake -i 'QCD_Pt-*EMEnriched*' -e 2016postVFP -n 40 --tag QCE_1

# flavor BDT ntuple (userflag MakeTree -> tree "mu")
# 출력: /gv0/Users/snuintern2/SKNanoOutput/fake/MakeTree/{era}_bdtflav/
# MC 먼저, DATA 는 MC 가 빠진 뒤에 (아래 주석 해제)
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2016preVFP  -n 40 --userflags MakeTree --tag bdtflav
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2016postVFP -n 40 --userflags MakeTree --tag bdtflav
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2017        -n 40 --userflags MakeTree --tag bdtflav
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2018        -n 40 --userflags MakeTree --tag bdtflav

#SKNano.py -a fake -i TTLJ_powheg -e 2016preVFP  -n 40 --userflags MakeTree --tag bdtflav
#SKNano.py -a fake -i TTLJ_powheg -e 2016postVFP -n 40 --userflags MakeTree --tag bdtflav
#SKNano.py -a fake -i TTLJ_powheg -e 2017        -n 40 --userflags MakeTree --tag bdtflav
#SKNano.py -a fake -i TTLJ_powheg -e 2018        -n 40 --userflags MakeTree --tag bdtflav

#SKNano.py -a fake -i DoubleMuon -e 2016preVFP  -n 40 --userflags MakeTree --tag bdtflav
#SKNano.py -a fake -i DoubleMuon -e 2016postVFP -n 40 --userflags MakeTree --tag bdtflav
#SKNano.py -a fake -i DoubleMuon -e 2017        -n 40 --userflags MakeTree --tag bdtflav
#SKNano.py -a fake -i DoubleMuon -e 2018        -n 40 --userflags MakeTree --tag bdtflav

# --- bdtflav2: 위와 같은 ntuple 에 sv_* (source jet 안의 IVF secondary vertex)
#     12개 브랜치를 추가한 재생산. b/c 분리용 정보를 DeepJet 출력 너머로 더 준다.
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2016preVFP  -n 40 --userflags MakeTree --tag bdtflav2
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2016postVFP -n 40 --userflags MakeTree --tag bdtflav2
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2017        -n 40 --userflags MakeTree --tag bdtflav2
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2018        -n 40 --userflags MakeTree --tag bdtflav2

#SKNano.py -a fake -i TTLJ_powheg -e 2016preVFP  -n 40 --userflags MakeTree --tag bdtflav2
#SKNano.py -a fake -i TTLJ_powheg -e 2016postVFP -n 40 --userflags MakeTree --tag bdtflav2
#SKNano.py -a fake -i TTLJ_powheg -e 2017        -n 40 --userflags MakeTree --tag bdtflav2
#SKNano.py -a fake -i TTLJ_powheg -e 2018        -n 40 --userflags MakeTree --tag bdtflav2

# --- nff (2026-08-11): neural fake factor (arXiv:2511.06972) 학습용.
#     bdtflav2 와 브랜치는 같고, fillTreeRow 에서 MET<25 / MT<25 컷을 뺐다.
#     공식 MC FR 맵(fr/fakerate_muon_*.json)은 Inclusive 영역(away jet 만)에서
#     재므로 컷이 박힌 bdtflav2 로는 그 맵을 재현할 수 없다. MeasReg 는
#     ev_met / ev_mt 로 offline 복원. MC only — data 는 이 스터디에 안 쓴다.
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2016preVFP  -n 40 --userflags MakeTree --tag nff
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2016postVFP -n 40 --userflags MakeTree --tag nff
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2017        -n 40 --userflags MakeTree --tag nff
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e 2018        -n 40 --userflags MakeTree --tag nff

#SKNano.py -a fake -i TTLJ_powheg -e 2016preVFP  -n 40 --userflags MakeTree --tag nff
#SKNano.py -a fake -i TTLJ_powheg -e 2016postVFP -n 40 --userflags MakeTree --tag nff
#SKNano.py -a fake -i TTLJ_powheg -e 2017        -n 40 --userflags MakeTree --tag nff
#SKNano.py -a fake -i TTLJ_powheg -e 2018        -n 40 --userflags MakeTree --tag nff

# --- elecfake nff (2026-08-18): 위 muon ntuple 의 electron 판 (tree "el").
#     AN Fig. 50 의 eµµ row 를 재현하려면 electron neural FR 이 필요하고,
#     그러려면 elecfake 가 tree 를 뱉어야 한다 (지금까지는 히스토그램만 냈다).
#     fillTreeRow 는 fake.cc 와 같은 규칙: trigger path loop 밖에서 호출,
#     MET/MT 컷 없음 (공식 QCD MC FR 이 Inclusive 영역), away jet 은 30 까지
#     열고 40/60 은 플래그. 축은 (el_conept, el_abssceta).
#     비교 대상 맵: fr/fakerate_elec_*.json 의 fakerate_electron_QCD_EMEnriched.
#     출력: /gv0/Users/snuintern2/SKNanoOutput/elecfake/MakeTree/{era}_nff/
#SKNano.py -a elecfake -i 'QCD_Pt-*EMEnriched*' -e 2016preVFP  -n 40 --userflags MakeTree --tag nff
#SKNano.py -a elecfake -i 'QCD_Pt-*EMEnriched*' -e 2016postVFP -n 40 --userflags MakeTree --tag nff
#SKNano.py -a elecfake -i 'QCD_Pt-*EMEnriched*' -e 2017        -n 40 --userflags MakeTree --tag nff
#SKNano.py -a elecfake -i 'QCD_Pt-*EMEnriched*' -e 2018        -n 40 --userflags MakeTree --tag nff

#SKNano.py -a elecfake -i TTLJ_powheg -e 2016preVFP  -n 40 --userflags MakeTree --tag nff
#SKNano.py -a elecfake -i TTLJ_powheg -e 2016postVFP -n 40 --userflags MakeTree --tag nff
#SKNano.py -a elecfake -i TTLJ_powheg -e 2017        -n 40 --userflags MakeTree --tag nff
#SKNano.py -a elecfake -i TTLJ_powheg -e 2018        -n 40 --userflags MakeTree --tag nff

# --- closnff (2026-08-18 로 갱신): 이제 3Mu + 1E2Mu 두 채널을 한 번에 돈다.
#     lepton multiplicity 가 배타적이라 한 이벤트는 최대 한 채널에만 들어가므로
#     job 을 나눌 이유가 없다 — 아래 명령은 그대로 두고 tree 의 channel 브랜치로
#     (0 = 3Mu, 1 = 1E2Mu) 구분한다. 1E2Mu 는 EMu trigger 를 쓰고 fake 후보가
#     electron 이라 el1_* slot 이 채워진다 (elecfake nff 모델의 입력과 같은 이름).
#     기존 3Mu 결과는 바뀌면 안 된다 — Observed(TTLL) = 389.23 이 regression gate.
#
# --- closnff (2026-08-15): AN-25-154 Fig. 50 의 closure 선택.
#     aa/ClosFakeRate.cc 와 같은 이벤트 선택이지만 fake rate 를 analyzer 안에서
#     계산하지 않고 fail-side muon 의 feature 를 tree "clos" 로 뱉는다
#     (neural FR 은 2D 히스토그램으로 표현할 수 없어 per-muon feature 가 필요).
#     tight/loose WP 은 fake.cc::PassMuonWP 과 동일 — 측정과 적용이 어긋나면 안 된다.
#     TTLL 이 주된 기여다 (dileptonic: prompt muon 2개 + fake 1개).
#     출력: /gv0/Users/snuintern2/SKNanoOutput/closnff/{era}_clos/
#SKNano.py -a closnff -i TTLL_powheg -e 2016preVFP  -n 40 --tag clos
#SKNano.py -a closnff -i TTLL_powheg -e 2016postVFP -n 40 --tag clos
#SKNano.py -a closnff -i TTLL_powheg -e 2017        -n 40 --tag clos
#SKNano.py -a closnff -i TTLL_powheg -e 2018        -n 40 --tag clos

#SKNano.py -a closnff -i TTLJ_powheg -e 2016preVFP  -n 40 --tag clos
#SKNano.py -a closnff -i TTLJ_powheg -e 2016postVFP -n 40 --tag clos
#SKNano.py -a closnff -i TTLJ_powheg -e 2017        -n 40 --tag clos
#SKNano.py -a closnff -i TTLJ_powheg -e 2018        -n 40 --tag clos

# DATA 는 MC 가 빠진 뒤에 (hadd DoubleMuon_*.root -> DoubleMuon.root 는 수동)
#SKNano.py -a fake -i DoubleMuon -e 2016preVFP  -n 40 --userflags MakeTree --tag bdtflav2
#SKNano.py -a fake -i DoubleMuon -e 2016postVFP -n 40 --userflags MakeTree --tag bdtflav2
#SKNano.py -a fake -i DoubleMuon -e 2017        -n 40 --userflags MakeTree --tag bdtflav2
#SKNano.py -a fake -i DoubleMuon -e 2018        -n 40 --userflags MakeTree --tag bdtflav2

# --- Run 3 (2026-08-18): µµµ closure + neural FR 를 2022/2023 으로 확장.
#     목적: sv3 (conePt, |eta|, sv_masssum, jet_nconst) 가 Run 2 우연이 아닌지
#     확인 — 다른 sqrt(s), 다른 pileup, 다른 NanoAOD version 에서 재현되어야 한다.
#     fake.cc 는 이번에 jet_chmult / jet_nemult (Run3 전용 branch) 를 추가했고
#     electron veto 의 loose MVA 컷을 Run3 값으로 고쳤다 -> rebuild 완료.
#     2022 는 ForSNU json 의 sumW/sumsign 이 -1 이라 MCweight 부호가 뒤집혔었다.
#     CommonSampleInfo.json 의 값으로 채워 넣어 해결 (GetEffLumi 재실행 불필요).
#     출력: fake  -> /gv0/.../SKNanoOutput/fake/MakeTree/{era}_nff/
#           closnff-> /gv0/.../SKNanoOutput/closnff/{era}_clos/
#     [2026-08-18 재제출] Gate 4 실패: 우리 Run3 FR 이 공식 map 대비
#     ptcorr 10-12 에서 1.30, 50-100 에서 0.81 로 어긋났다 (4 era, 3 |eta| 모두
#     동일한 monotonic trend). 원인은 muon **loose** WP 이 Run 별로 다른 것:
#       Run2 HcToWALooseRun2 : SIP3D<5, miniiso<0.6
#       Run3 HcToWALooseRun3 : SIP3D<8, miniiso<0.4   (tight 은 동일)
#     공식 aa/Analyzers/src/MeasFakeRateV4.cc:38 참조. fake.cc / closnff.cc 둘 다
#     고쳐서 rebuild 했으므로 아래 16개를 그대로 다시 돌려야 한다.
#     (이전 Run3 출력과 parquet 은 전부 폐기했다. Run2 는 영향 없음.)
#
#     [2026-08-18 2차 재제출] electron loose WP 도 Run 별로 달랐다:
#       Electron Run2 SIP3D<8  ->  Run3 SIP3D<6   (miniiso 는 둘 다 0.4)
#     공식에서 Run 별로 갈리는 WP 은 muon loose / electron loose 딱 둘뿐이고
#     (tight 은 둘 다 Run 무관) 이제 셋 다 맞췄다:
#       fake.cc      muon loose + electron veto
#       closnff.cc   muon loose + electron loose  (측정과 적용이 같아야 한다)
#       elecfake.cc  electron loose + muon veto
#     ** 이전 배치는 condor_rm 으로 죽이고 아래를 다시 돌려야 한다 — 중간에
#        rebuild 가 들어갔으므로 그대로 두면 한 hadd 안에 두 WP 이 섞인다. **
for ERA in 2022 2022EE 2023 2023BPix; do
SKNano.py -a fake    -i 'QCD_Pt-*_MuEnriched' -e $ERA -n 40 --userflags MakeTree --tag nff
SKNano.py -a fake    -i TTLJ_powheg           -e $ERA -n 40 --userflags MakeTree --tag nff
SKNano.py -a closnff -i TTLL_powheg           -e $ERA -n 40 --tag clos
SKNano.py -a closnff -i TTLJ_powheg           -e $ERA -n 40 --tag clos
done
