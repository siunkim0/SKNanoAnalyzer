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
for ERA in 2017 2018; do
SKNano.py -a fake -i DoubleMuon   -e $ERA -n 40 --userflags MakeTree --tag nffdata
# prompt MC: 공식과 동일한 13개
SKNano.py -a fake -i DYJets       -e $ERA -n 40 --userflags MakeTree --tag nffdata
SKNano.py -a fake -i DYJets10to50 -e $ERA -n 40 --userflags MakeTree --tag nffdata
SKNano.py -a fake -i WJets        -e $ERA -n 40 --userflags MakeTree --tag nffdata
SKNano.py -a fake -i TTLJ_powheg  -e $ERA -n 40 --userflags MakeTree --tag nffdata
SKNano.py -a fake -i TTLL_powheg  -e $ERA -n 40 --userflags MakeTree --tag nffdata
SKNano.py -a fake -i WW_pythia    -e $ERA -n 40 --userflags MakeTree --tag nffdata
SKNano.py -a fake -i WZ_pythia    -e $ERA -n 40 --userflags MakeTree --tag nffdata
SKNano.py -a fake -i ZZ_pythia    -e $ERA -n 40 --userflags MakeTree --tag nffdata
SKNano.py -a fake -i SingleTop_sch_Lep            -e $ERA -n 40 --userflags MakeTree --tag nffdata
SKNano.py -a fake -i SingleTop_tch_top_Incl       -e $ERA -n 40 --userflags MakeTree --tag nffdata
SKNano.py -a fake -i SingleTop_tch_antitop_Incl   -e $ERA -n 40 --userflags MakeTree --tag nffdata
SKNano.py -a fake -i SingleTop_tW_top_NoFullyHad     -e $ERA -n 40 --userflags MakeTree --tag nffdata
SKNano.py -a fake -i SingleTop_tW_antitop_NoFullyHad -e $ERA -n 40 --userflags MakeTree --tag nffdata
done

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
