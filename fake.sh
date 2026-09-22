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
#     (끝남 — RESULTS.md 4h. 다시 돌리지 말 것.)
#for ERA in 2022 2022EE 2023 2023BPix; do
#SKNano.py -a fake    -i 'QCD_Pt-*_MuEnriched' -e $ERA -n 40 --userflags MakeTree --tag nff
#SKNano.py -a fake    -i TTLJ_powheg           -e $ERA -n 40 --userflags MakeTree --tag nff
#SKNano.py -a closnff -i TTLL_powheg           -e $ERA -n 40 --tag clos
#SKNano.py -a closnff -i TTLJ_powheg           -e $ERA -n 40 --tag clos
#done

# --- nffdata (2026-08-27): neural fake factor 를 DATA 에서 재기 위한 ntuple.
#     지금까지의 nff 스터디는 전부 MC (QCD -> TTLJ) 라서, sv3 가 데이터에서도
#     되는지는 확인된 적이 없다. RESULTS.md 4c/4d/4f 참조.
#
#     tag 만 다르고 analyzer 는 nff 와 완전히 같다 (fake.cc 수정 없음):
#     fillTreeRow 는 이미 DATA 에서 돌고 (w=1, isPrompt=false, flavor=-2),
#     fillZEnriched 도 DATA/MC 양쪽에서 돌기 때문에, 이 job 한 세트로
#       - tree "mu"                  -> neural FR 학습용
#       - {path}/MeasReg/...         -> binned data FR (Gate 1 의 비교 대상)
#       - {path}/ZEnriched/mll       -> prescale normalization scale
#     이 한꺼번에 나온다.
#
#     prompt 제거는 Tools/fakerate_data_pre.py 와 같은 방식:
#       scale = N_data(Z) / N_MC(Z)  (trigger path 별),
#       FR = [T_data - scale*T_prompt] / [L_data - scale*L_prompt]
#     DYJets / DYJets10to50 / WJets 가 그 prompt MC 다. MC 의 fake 성분은
#     빼면 안 된다 (데이터의 fake 에 이미 들어있다).
#
#     2016 두 era 만 먼저: 그 두 개만 histogram level 데이터 측정이 이미
#     돌아가 있어서 fr/fakerate_muon_{2016pre,2016post}.json 의
#     `fakerate_muon_Central` 로 정확한 cross-check 가 된다.
#
#     ** DoubleMuon 은 period 별로 나오므로 hadd 로 DoubleMuon.root 하나로
#        합쳐야 한다. hadd -j 8 은 8 GB cgroup 에서 죽는다 (elecfake 때 그랬다,
#        RESULTS.md 4e) — 직렬로, condor 밖에서. 그리고 skim 의 glob 은 반드시
#        "DoubleMuon.root" 여야 한다; "DoubleMuon*.root" 로 하면 합친 파일과
#        period 파일이 둘 다 잡혀서 모든 row 가 두 번 세어진다. **
#     prompt MC 목록은 공식과 맞췄다 (aa/scripts/submission/MeasFakeRate.sh:27-29:
#     DYJets, WJets, TTLJ_powheg, TTLL_powheg + SampleLists/Run2NanoV9/DiLepton.txt).
#     Tools/fakerate_data_pre.py 는 DY+WJets 만 뺐고, 그것이 우리 data FR 이
#     공식 맵보다 높은 이유다: Mu17 영역에서 prompt 를 ~3% 덜 빼면 prompt 가
#     tight data 의 73-82% 인 마지막 conePt bin 에서 25% 차이로 증폭된다.
#     RESULTS.md 4j 참조.
#     2016preVFP / 2016postVFP 는 이미 돌았다 (RESULTS.md 4j). 남은 두 era 만.
#     주의 1: 2018 의 DoubleMuon_E / _F 는 ForSNU json 의 path 가 비어 있다
#             (2018 실제 period 는 A-D). 빈 출력이 나오는 건 정상이고 hadd 에서
#             무시하면 된다.
#     주의 2: 2017/2018 DoubleMuon 은 2016 보다 훨씬 크다 (2018 period D 하나가
#             파일 2059개). hadd 는 반드시 직렬로, condor 밖에서.
#for ERA in 2017 2018; do
#SKNano.py -a fake -i DoubleMuon   -e $ERA -n 40 --userflags MakeTree --tag nffdata
# prompt MC: 공식과 동일한 13개
#SKNano.py -a fake -i DYJets       -e $ERA -n 40 --userflags MakeTree --tag nffdata
#SKNano.py -a fake -i DYJets10to50 -e $ERA -n 40 --userflags MakeTree --tag nffdata
#SKNano.py -a fake -i WJets        -e $ERA -n 40 --userflags MakeTree --tag nffdata
#SKNano.py -a fake -i TTLJ_powheg  -e $ERA -n 40 --userflags MakeTree --tag nffdata
#SKNano.py -a fake -i TTLL_powheg  -e $ERA -n 40 --userflags MakeTree --tag nffdata
#SKNano.py -a fake -i WW_pythia    -e $ERA -n 40 --userflags MakeTree --tag nffdata
#SKNano.py -a fake -i WZ_pythia    -e $ERA -n 40 --userflags MakeTree --tag nffdata
#SKNano.py -a fake -i ZZ_pythia    -e $ERA -n 40 --userflags MakeTree --tag nffdata
#SKNano.py -a fake -i SingleTop_sch_Lep            -e $ERA -n 40 --userflags MakeTree --tag nffdata
#SKNano.py -a fake -i SingleTop_tch_top_Incl       -e $ERA -n 40 --userflags MakeTree --tag nffdata
#SKNano.py -a fake -i SingleTop_tch_antitop_Incl   -e $ERA -n 40 --userflags MakeTree --tag nffdata
#SKNano.py -a fake -i SingleTop_tW_top_NoFullyHad     -e $ERA -n 40 --userflags MakeTree --tag nffdata
#SKNano.py -a fake -i SingleTop_tW_antitop_NoFullyHad -e $ERA -n 40 --userflags MakeTree --tag nffdata
#done

# --- nffb (2026-09-02): DY / W / DY+b 로 closure test 를 넓히기 위한 ntuple.
#     지금까지 closure 는 QCD -> TTLJ (per-muon) 와 TTLL (trilepton) 뿐이라
#     binned map 의 실패와 sv3/sv3m 의 우세가 ttbar 고유의 사실인지 알 수 없다.
#     DY = 낮은 HT / light-gluon jet 환경, W = 진짜 neutrino 하나,
#     DY+b = top decay kinematics 없는 b flavour fake.
#
#     tag 를 nff 와 나눈 이유: 이 배치부터 fake.cc 가 away jet 관련 브랜치
#     두 개를 더 뱉는다 (Central 정의 pT>40 && dR(mu,jet)>0.7 기준):
#       ev_nawayjet40   away jet 개수
#       ev_nawayjetb40  그 중 DeepJet Medium b-tag 개수  (closnff.cc 와 같은 WP)
#     DY+b 영역은 ev_nawayjetb40 == 1 로 offline 에서 자른다.
#     기존 nff (QCD/TTLJ) 트리에는 이 브랜치가 없으므로 한 디렉토리에 섞지 않는다.
#     bdtflav -> bdtflav2 -> nff 와 같은 이유.
#
#     Run2 DY/WJets 는 {era}_nffdata 에 이미 있지만 그건 prompt subtraction 용
#     (skim 이 isPrompt 만 남기고 weight 부호를 뒤집는다) + away b 브랜치가 없다.
#     8 era x 3 sample = 24 개, /gv0 에 약 10 GB.
#     출력: /gv0/Users/snuintern2/SKNanoOutput/fake/MakeTree/{era}_nffb/
#for ERA in 2016preVFP 2016postVFP 2017 2018 2022 2022EE 2023 2023BPix; do
#SKNano.py -a fake -i DYJets       -e $ERA -n 40 --userflags MakeTree --tag nffb
#SKNano.py -a fake -i DYJets10to50 -e $ERA -n 40 --userflags MakeTree --tag nffb
#SKNano.py -a fake -i WJets        -e $ERA -n 40 --userflags MakeTree --tag nffb
#done

# --- nffs (2026-09-02): ttbar 의 muR/muF, ISR/FSR 를 흔들었을 때에도 ML fake
#     rate 가 닫히는지 보기 위한 배치. 지금까지의 closure (RESULTS.md 4b/4c/4f)
#     는 전부 TTLJ_powheg 의 nominal scale 하나에서만 잰 것이라, sv3 의 우세가
#     방법의 성질인지 ttbar tune 하나의 성질인지 구분이 안 된다.
#
#     tag 를 nff/nffb 와 나눈 이유: fillTreeRow 가 branch 13 개를 더 뱉는다.
#       w_scale0..8   LHEScaleWeight raw 9-entry (nominal 은 [4], 값 1)
#       w_isr_up/dn, w_fsr_up/dn   PSWeight
#     전부 w_var / w_nominal 이라 offline 에서 weight 에 그냥 곱하면 된다.
#     AnalyzerCore::GetScaleVariation() 은 nominal 이 빠진 8-entry 규약을
#     가정해서 여기 NanoAOD (8 era 전부 nLHEScaleWeight == 9) 에서는 "up" 이
#     한 칸씩 밀린다. 그래서 fake.cc 가 raw 벡터를 직접 남긴다.
#
#     QCD 는 넣지 않는다: pythia8 이라 LHEScaleWeight 브랜치 자체가 없고
#     (nPSWeight == 1), fake rate 를 재는 쪽은 어차피 흔들 수 없다. 그게 오히려
#     data 에서 잰 fake rate 와 같은 상황이라 이 테스트에 맞다.
#
#     8 era x 1 sample. 출력: /gv0/Users/snuintern2/SKNanoOutput/fake/MakeTree/{era}_nffs/
#for ERA in 2016preVFP 2016postVFP 2017 2018 2022 2022EE 2023 2023BPix; do
#SKNano.py -a fake -i TTLJ_powheg -e $ERA -n 40 --userflags MakeTree --tag nffs
#done

# --- closdy (2026-09-03): DY 의 closure 를 FR 측정 영역이 아니라 **trilepton
#     영역** 에서 한다. per-muon closure (tag nffb) 는 fake.cc 의 측정 영역을
#     그대로 쓰는데, 그 영역은 loose muon 이 정확히 1개라 평범한 Z->mumu 를
#     통째로 veto 하고, 게다가 fake muon 본인이 prescale 된 Mu8/Mu17 을 직접
#     쏴야 한다. 그래서 DY 의 유효 통계가 ttbar 의 1/540 (N_eff 1036 vs 559k)
#     까지 떨어진다.
#
#     closnff 의 3Mu 영역에서는 그 두 제약이 동시에 사라진다:
#       - Z(mumu) + fake mu 가 veto 대상이 아니라 바로 그 신호 배치다
#       - DblMu trigger 는 prompt Z muon 두 개가 쏘므로 fake muon 에는
#         trigger 요구가 아예 없다
#       - away jet 요구도 없다
#     dilepton 관측량(ZCand_mass/pt, pair_lowM/highM_mass, m3l)은 이미 clos
#     tree 의 브랜치다.
#
#     --userflags OpenJet 이 jet>=2 && b-jet>=1 을 끈다 (closnff.cc:66-72).
#     그 컷은 ttbar 에 맞춘 것이라 DY 에 그대로 쓰면 b 가 없어 통계가 안 남는다.
#     ev_njet / ev_nbjet 는 브랜치로 나가므로 offline --select 로 재현한다:
#       DY      = 컷 없음
#       DY+b    = ev_nbjet >= 1
#       AN 영역 = ev_njet >= 2, ev_nbjet >= 1   (ttbar 와 같은 조건에서 비교)
#
#     tag 를 clos 와 나눈 이유: 이 배치는 selection 이 열려 있어 같은 디렉토리에
#     섞으면 --closnff-glob 이 어느 쪽인지 구분할 수 없다. 기존 {era}_clos/
#     의 TTLL/TTLJ 는 건드리지 않는다 (Observed(TTLL) = 389.23 regression gate).
#
#     16 job. 출력: /gv0/Users/snuintern2/SKNanoOutput/closnff/{era}_closdy/
#for ERA in 2016preVFP 2016postVFP 2017 2018 2022 2022EE 2023 2023BPix; do
#SKNano.py -a closnff -i DYJets       -e $ERA -n 40 --userflags OpenJet --tag closdy
#SKNano.py -a closnff -i DYJets10to50 -e $ERA -n 40 --userflags OpenJet --tag closdy
#done
#
#     그리고 Stage 1 편집이 무해함을 증명하는 job 하나 (userflag 없이).
#     row 수가 2018_clos/TTLL_powheg.root 과 같아야 한다.
#SKNano.py -a closnff -i TTLL_powheg -e 2018 -n 40 --tag closchk

# --- closdimu (2026-09-03): DY 의 closure 를 **OS dimuon** 영역에서. bkg.cc 가
#     matrix method 를 실제로 배포하는 영역이고 (fake/Tools 의 fake background
#     숫자가 거기서 나온다), README.md:155 가 FR universality — 이 방법의 주된
#     계통오차 — 를 재는 방법으로 바로 그 MC closure 를 지목한다.
#     지금 그 closure 는 binned map 에만 있다: bkg.cc 는 히스토그램만 쓰고
#     per-lepton feature tree 가 없어 neural FR 을 아예 평가할 수 없다.
#
#     그래서 closnff 에 channel = 2 (2Mu) 를 추가했다:
#       loose muon 2개 + veto electron 0개
#       trigger 는 **FR 을 잰 경로** Mu8/Mu17_TrkIsoVVL (leading pT > 10 / 20)
#       charge / mll / Z veto / MET / jet 컷 **없음** — 전부 offline
#     bkg.cc 영역(OS + mll>20 + |mll-91.2|>15 + njet>=2 + MET>40)은 그 중
#     한 slice 이고, DY 용 Z window 도 같은 트리에서 나온다.
#
#     주의 1: Mu8/Mu17 은 prescale 되어 있는데 closnff 의 weight 는
#             GetTriggerLumi("Full") (상수) 이다. 따라서 2Mu 의 절대 수율은
#             물리적 이벤트 수가 아니다. Exp/Obs 에서는 상쇄되지만 3Mu 채널의
#             수율과 나란히 놓으면 안 된다.
#     주의 2: 2Mu 에서는 SR(전부 tight)이 진짜 prompt pair 로 뒤덮이므로
#             offline 에서 --truth-fake 가 **필수**다 (bkg.cc:231 과 같은 분할).
#
#     16 job. 출력: /gv0/Users/snuintern2/SKNanoOutput/closnff/OpenJet/{era}_closdimu/
#             (--userflags 가 디렉토리 한 겹을 끼워 넣는다)
#for ERA in 2016preVFP 2016postVFP 2017 2018 2022 2022EE 2023 2023BPix; do
#SKNano.py -a closnff -i DYJets       -e $ERA -n 40 --userflags OpenJet --tag closdimu
#SKNano.py -a closnff -i DYJets10to50 -e $ERA -n 40 --userflags OpenJet --tag closdimu
#done
#
#     그리고 2Mu 추가가 기존 채널을 안 건드렸다는 증명 job 하나 (userflag 없이).
#     2018_clos/TTLL_powheg.root 과 같아야 한다: 77430 rows,
#     3Mu 34718 / SR 4667 / sum w 167.29244, 1E2Mu 42712 / SR 6787 / 243.65506
#SKNano.py -a closnff -i TTLL_powheg -e 2018 -n 40 --tag closchk2

# --- closos (2026-09-03): 위 closdimu 를 **대체**한다. closdimu 는 prescale 된
#     Mu8/Mu17 (FR 을 잰 경로) 로 trigger 했는데, 경로만 같고 정작 중요한 bias
#     는 같지 않았다 — 측정에서는 fake muon 본인이 trigger object 인데 dimuon
#     에서는 prompt muon 이 쏜다. 게다가 weight 가 GetTriggerLumi("Full")
#     (상수) 라 prescale 된 경로에서는 절대 수율이 물리적이지 않았다.
#
#     이제 unprescaled single muon trigger 를 쓴다 (bkg.cc:14-22 와 동일):
#       2016 IsoMu24 + IsoTkMu24 / 2017 IsoMu27 / 그 외 IsoMu24
#       leading pT > 26 (2017 은 29), subleading > 10
#     하나만 쏘면 되고 HLT leg 이 tight isolation 이라 사실상 prompt muon 이
#     쏜다 -> **fake muon 에 trigger 요구가 없다** (3Mu 영역이 깨끗한 이유와 동일).
#     prescale 이 없으니 수율도 물리적이다.
#
#     analyzer 에서 자르는 것: OS (q1+q2==0) + M(mumu) > 12 (AN 의 12, bkg.cc 의
#     mll>20 이 아니다). Z veto / MET / jet 은 여전히 offline:
#       bkg.cc 영역 = mll>20 && |mll-91.2|>15 && ev_njet>=2 && ev_met>40
#       DY+b        = ev_nbjet>=1
#
#     주의: offline 에서 --truth-fake 가 **필수**다. 2Mu 의 SR(전부 tight)은
#           진짜 prompt pair 로 뒤덮이는데 matrix method 는 그걸 예측하겠다고
#           한 적이 없다 (bkg.cc:231 과 같은 분할).
#
#     16 job. 출력: /gv0/Users/snuintern2/SKNanoOutput/closnff/OpenJet/{era}_closos/
#for ERA in 2016preVFP 2016postVFP 2017 2018 2022 2022EE 2023 2023BPix; do
#SKNano.py -a closnff -i DYJets       -e $ERA -n 40 --userflags OpenJet --tag closos
#SKNano.py -a closnff -i DYJets10to50 -e $ERA -n 40 --userflags OpenJet --tag closos
#done
#
#     3Mu / 1E2Mu 가 안 변했다는 증명 job (userflag 없이). **총 row 수는 늘어난다**
#     (2Mu row 가 추가되므로) — 채널별로 비교해야 한다:
#       3Mu 34718 / SR 4667 / sum w 167.29244
#       1E2Mu 42712 / SR 6787 / sum w 243.65506
#SKNano.py -a closnff -i TTLL_powheg -e 2018 -n 40 --tag closchk3

# --- clostune (2026-09-11): ttbar modelling 변형이 trilepton closure 를 흔드는가.
#     RESULTS 4u. 4l 은 muR/muF/ISR/FSR 를 **reweighting** 으로 했다
#     (LHEScaleWeight / PSWeight = w_var / w_nominal, 같은 row, 한 번의 dump).
#     여기서는 재생성된 **별도 샘플** 이다. 따라서
#       - nominal 과 event 를 공유하지 않는다 -> shift 오차가 4l 처럼 상쇄되지
#         않고 통째로 남는다. Run2 3Mu nominal 이 Observed 389.23 +- 3.82
#         (0.98%) 이고 변형 샘플은 그 36-43% 뿐이라, 변형 하나당 closure ratio
#         오차 ~1.85%, nominal 대비 shift 오차 ~2.2%. **~2% 수준의 null 이
#         예상 결과다** — 실패가 아니라 이 방법의 한계다. 먼저 말하고 숫자를 읽자.
#       - TTLJ 에는 이 변형들이 없다 (TTToSemiLeptonic_TuneCP5up 는 존재하지
#         않는다). 그래서 4l(TTLJ per-muon) 과 직접 비교할 수 없고, baseline 은
#         TTLL nominal 이어야 한다.
#
#     nominal 은 **다시 돌릴 필요가 없다**: closnff/OpenJet_Skip2Mu/{era}_closttop/
#     TTLL_powheg.root (2026-09-08) 가 지금 HEAD 의 closnff 로, 같은 userflag 로
#     이미 만들어져 있다. 확인했다 — offline 로 njet>=2 && nbjet>=1 을 걸면
#     기존 {era}_clos (analyzer 안에서 자른 것) 와 정확히 같다:
#       Observed 389.23 +- 3.82 / flat 1.1168 / sv3 1.0127
#     closnff.cc 의 uncommitted diff 는 branch 추가뿐이라 rebuild 해도 무효화되지
#     않는다 (selection/weight 불변).
#
#     --userflags OpenJet,Skip2Mu 인 이유:
#       OpenJet  jet>=2 && b-jet>=1 을 offline 으로 넘긴다. AN 영역은 --select 가
#                되고, 통계가 모자라면 컷 없는 open 영역(약 3배)을 같은 트리에서
#                쓸 수 있다. 샘플이 nominal 의 40% 뿐이라 이 탈출구가 중요하다.
#       Skip2Mu  OS dimuon 벌크(트리의 99.9%, 예전에 hadd 를 8GB 에서 죽인 것)만
#                버린다. 3Mu / 1E2Mu / 2Mu-SS 는 남는다.
#
#     샘플 8종은 2026-09-11 에 4 era 모두 sample DB 에 등록했다 (CommonSampleInfo
#     + ForSNU). sumW/nmc 는 Runs tree 에서, sumsign 은 sumW/|genWeight| 로
#     구했다 — powheg 라 |genWeight| 가 상수다. nominal 로 gate 했고 nmc 는 정확히,
#     sumW/sumsign 은 4e-8 로 일치했다. GetEffLumi 제출은 필요 없었다.
#
#     32 job. 출력: /gv0/Users/snuintern2/SKNanoOutput/closnff/OpenJet_Skip2Mu/{era}_clostune/
#     (디렉토리는 --tag 가 아니라 **userflags** 를 인코딩한다. 두 세션이 여기 걸렸다.)
#     node 0 개짜리 DAG 가 나오면 alias 등록이 안 된 것이다.
#for ERA in 2016preVFP 2016postVFP 2017 2018; do
#for S in TTLL_powheg_CP5up TTLL_powheg_CP5down TTLL_powheg_CR1 TTLL_powheg_CR2 \
#         TTLL_powheg_hdampUP TTLL_powheg_hdampDOWN \
#         TTLL_powheg_mtop171p5 TTLL_powheg_mtop173p5; do
#SKNano.py -a closnff -i $S -e $ERA -n 40 --userflags OpenJet,Skip2Mu --tag clostune
#done; done
#
#     2026-09-12: 두 개가 계속 죽었다. 원인은 framework 이 아니라 **중앙
#     NanoAOD 두 파일이 깨진 것**이다 (uproot 로 branch 11 개가 안 읽히고,
#     ROOT 는 "TBasket::Streamer: fObjlen is incorrect (-270467072)" 를 내고
#     zombie buffer 를 돌려주는데 SKNanoLoader 가 그걸 dereference 한다):
#       2016postVFP TTLL_powheg_mtop173p5  NANOAOD_259.root
#       2018        TTLL_powheg_hdampUP    NANOAOD_1052.root
#     결정론적이다 — 재제출해도 같은 entry 에서 죽고, -n 을 바꾸면 그 파일을
#     맡는 job 번호만 바뀐다. 주의: crash 는 **앞 파일의 끝**에서 보고된다
#     (chain 이 거기서 다음 파일을 연다). 나머지 파일은 전부 정상임을 확인했다.
#
#     조치: 두 파일을 ForSNU path 목록에서 빼고 nmc/sumW/sumsign 을 남은
#     파일로 **다시 계산**했다 (각각 0.246% / 0.088% 손실). 생성 총량을 그대로
#     두면 그 비율만큼 수율이 낮게 나온다. 그 다음 이 둘만 재제출:
#SKNano.py -a closnff -i TTLL_powheg_mtop173p5 -e 2016postVFP -n 40 --userflags OpenJet,Skip2Mu --tag clostune
#SKNano.py -a closnff -i TTLL_powheg_hdampUP   -e 2018        -n 40 --userflags OpenJet,Skip2Mu --tag clostune

# =========================================================================
# --- closnff, MATCHED LOOSE WP (2026-09-13):  --userflags MatchWP
#
#     왜: muon loose WP 이 Run 별로 다르다 (fake.cc:47 / closnff.cc:71).
#       Run2  |SIP3D| < 5 && miniiso < 0.6
#       Run3  |SIP3D| < 8 && miniiso < 0.4      tight 은 둘 다 |SIP3D| < 3
#     Run3 loose 안에만 never-tight (SIP3D 5-8) 띠가 들어가고, ttbar fake 는
#     QCD 보다 displaced 해서 그 띠가 ttbar 를 더 희석한다. 그 결과 binned map
#     의 over-prediction 이 Run2 +10.7% / Run3 +24.7% 로 갈라진다.
#     per-muon 쪽에서는 parquet 을 offline 으로 잘라 이미 확인했다: 두 run 의
#     **교집합** |SIP3D| < 5 && miniiso < 0.4 로 맞추면 8 era 전부 1.138-1.164
#     한 밴드로 모이고 run 간 격차가 +12.6% -> +0.4% 가 된다
#     (neuralff/RESULTS.md 4ab, fake/Tools/JETLESS.md).
#
#     왜 job 이 필요한가: trilepton 쪽은 offline 으로 못 한다.
#     closnff.cc:230 의 is3Mu 는 loose muon 이 **정확히 3개**일 것을 요구하므로,
#     예전 WP 에서 4개였다가 새 WP 에서 3개가 되는 event 는 애초에 tree 에
#     없다. 그래서 analyzer 단에서 WP 을 바꿔 다시 돌려야 한다.
#
#     --userflags MatchWP 는 fake.cc / closnff.cc 둘 다에 있고 값이 같다.
#     flag 를 안 주면 공식 WP 그대로다 (기존 tag 들은 영향 없음).
#     주의: MatchWP 는 muon 만 건드린다. electron veto WP 은 그대로다.
#
#     새 tag 를 쓴다 — 기존 clos / closdy 와 절대 섞지 말 것. 같은 디렉토리에
#     두 WP 이 섞이면 --closnff-glob 이 구분할 수 없다.
#
#     A) ttbar trilepton (AN-25-154 Fig. 50).  16 job.
#        출력: /gv0/.../closnff/MatchWP/{era}_closwp/   <- 한 단계 깊다
#        (userflag 는 SKNano.py:426 에서 디렉토리 level 을 만든다. 기존
#         {era}_clos/ 와 달리 --closnff 를 한 단계 더 내려 줘야 한다.)
#for ERA in 2016preVFP 2016postVFP 2017 2018 2022 2022EE 2023 2023BPix; do
#SKNano.py -a closnff -i TTLL_powheg -e $ERA -n 40 --userflags MatchWP --tag closwp
#SKNano.py -a closnff -i TTLJ_powheg -e $ERA -n 40 --userflags MatchWP --tag closwp
#done
#
#     B) DY trilepton.  OpenJet 이 jet>=2 && b-jet>=1 을 끈다 (DY 는 그 컷을
#        견딜 통계가 없다). userflag 가 디렉토리 level 을 하나 더 만든다:
#        출력: /gv0/.../closnff/OpenJet_MatchWP/{era}_closdywp/   <- 한 단계 깊다
#        16 job.
#for ERA in 2016preVFP 2016postVFP 2017 2018 2022 2022EE 2023 2023BPix; do
#SKNano.py -a closnff -i DYJets       -e $ERA -n 40 --userflags OpenJet,MatchWP --tag closdywp
#SKNano.py -a closnff -i DYJets10to50 -e $ERA -n 40 --userflags OpenJet,MatchWP --tag closdywp
#done
#
#     제출 전 확인: rebuild 는 이미 되어 있다 (2026-09-13, BUILD_EXIT=0).
#     제출 후 gate: Run2 4 era 의 Observed 합이 aggregate 와 같아야 한다
#     (closure_root.py --era-table 가 assert 한다). 공식 WP 의 389.23 과는
#     다른 값이 나오는 것이 정상이다 — loose 가 좁아졌으니 fail 쪽이 준다.

# =========================================================================
# --- Run 3 를 Run 2 의 loose WP 으로 (2026-09-14):  --userflags Run2WP
#
#     MatchWP 와 무엇이 다른가:
#       MatchWP : 두 WP 의 **교집합** |SIP3D|<5 && miniiso<0.4 를 양쪽에.
#                 Run 2 의 miniiso 까지 0.6 -> 0.4 로 조이므로 **어느 run 도
#                 자기 analysis WP 에 있지 않다**. Run 2 는 1.1074 -> 1.1484
#                 로 움직였고 (+4.1 pt), 그 숫자는 배포된 값이 아니다.
#       Run2WP  : Run 3 만 Run 2 것 그대로 |SIP3D|<5 && miniiso<0.6 으로.
#                 Run 2 는 손도 대지 않는다 (이미 그 값이라 no-op).
#                 기준계가 Run 2 의 실제 배포 WP 가 된다.
#     기대값 (§4ab 항 분해: SIP3D 8->5 = -9.4 pt, miniiso 0.6->0.4 = +4.1 pt):
#       Run 3 flat ~ 1.11  vs  Run 2 의 1.1074 — 격차가 다시 닫혀야 한다.
#
#     왜 job 이 필요한가 (MatchWP 때와 이유가 다르다): parquet 에서 loose 는
#     **좁히는 것만** 된다. 여기는 Run 3 miniiso 를 0.4 -> 0.6 으로 **넓혀야**
#     하는데 0.4-0.6 띠는 애초에 tree 에 쓰이지 않았다. per-muon 쪽까지
#     offline 으로 못 하므로 fake MakeTree 자체를 다시 돌린다.
#
#     Run 3 4 era 만. Run 2 는 재제출도 재학습도 없다.
#     Run2WP 와 MatchWP 를 **같이 주지 말 것** (뒤의 것이 이긴다).
#
#     A) per-muon arm — fake MakeTree.  8 submission.
#        출력: /gv0/.../fake/MakeTree_Run2WP/{era}_nffr2/   <- 한 단계 깊다
#        (userflag 가 SKNano.py:426 에서 디렉토리 level 을 만든다.)
#        QCD 는 map 과 model 의 소스, TTLJ 는 closure 의 target 이다.
#for ERA in 2022 2022EE 2023 2023BPix; do
#SKNano.py -a fake -i 'QCD_Pt-*MuEnriched*' -e $ERA -n 40 --userflags MakeTree,Run2WP --tag nffr2
#SKNano.py -a fake -i TTLJ_powheg           -e $ERA -n 40 --userflags MakeTree,Run2WP --tag nffr2
#done
#
#     B) trilepton arm — closnff (AN-25-154 Fig. 50).  4 submission.
#        출력: /gv0/.../closnff/Run2WP/{era}_closr2wp/      <- 한 단계 깊다
#        Run 2 쪽은 기존 {era}_clos/ 를 그대로 쓴다 (WP 이 안 바뀌었으므로).
#for ERA in 2022 2022EE 2023 2023BPix; do
#SKNano.py -a closnff -i TTLL_powheg -e $ERA -n 40 --userflags Run2WP --tag closr2wp
#done
#
#     rebuild: 2026-09-14, BUILD_EXIT=0 (Run2WP 가 libAnalyzers.so 에 들어있음).
#     디스크: /gv0 에 ~30 GB (Run 3 nff twin 이 26.5 GB 였다; miniiso 를 넓혔으니
#     조금 더 크다). 85 T 남아 있다.
#
#     제출 후 gate (가장 강한 것부터):
#       1. tight 은 불변이어야 한다. tight (SIP3D<3 && miniiso<0.1) 은 두 loose
#          안에 완전히 들어가므로, 새 parquet 의 tight row 수와 per-era
#          Observed 가 기존 ntuples_nff3 / plots/permuon 과 **정확히** 같아야
#          한다 (Run 3 합 11225.0). 다르면 loose WP 말고 다른 것이 바뀐 것이다.
#       2. tree 에 WP 이 실제로 들어갔는지: max mu_sip3d < 5.0, max mu_miniiso
#          < 0.6, 그리고 miniiso 0.4-0.6 띠에 row 가 **있어야** 한다.
#          (flag 없이 만든 tree 는 1,2 의 상한은 통과하고 이것만 실패한다.)
