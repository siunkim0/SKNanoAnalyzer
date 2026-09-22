# muid: muon loose WP 이 fake muon 을 얼마나 통과시키는지, 다섯 컷을 하나씩 분해한다.
#
# 왜: trilepton 의 `== 3 loose mu` 단계가 전체의 0.00107% 만 남긴다 (DY Run2,
#     3.385e9 -> 36,285). 그게 (a) 세 번째 muon **후보**가 애초에 없어서인지
#     (b) loose ID 가 fake 에게 가혹해서인지를 디스크의 트리로는 알 수 없다.
#     fake/MakeTree 도 closnff 도 **이미 loose 를 통과한** muon 만 저장한다.
#     그래서 GetAllMuons() 에서 시작하는 새 pass 가 필요하다.
#
# muid1: 첫 제출. Run2 4 era x (QCD MuEnriched 12 pT bin + TTLJ + DYJets +
#        DYJets10to50) = 60 job. tree 없음, 히스토그램만 -> 파일 ~100 kB,
#        hadd 메모리 문제 없음. 비용은 muon 마다 도는 GetLeptonType 인데
#        track.cc 가 이미 같은 일을 -n 40 으로 하고 있으므로 들어간다.
#
#        출력 /gv0/Users/snuintern2/SKNanoOutput/muid/{era}_muid1/
#        디렉토리 All/ (모든 muon) + Fake/ (GetLeptonType<=0). prompt 는 All-Fake.
#
#        대조 게이트 (plan Verification): CutFlow bin 0-2 가 closnff 의
#        CutFlow_3Mu bin 0-2 와, nLooseMu 가 CutFlow_nLooseMu 와 일치해야 한다.
#        preamble 과 weight 을 closnff.cc:210-225 와 같게 맞춘 이유가 그것이다.
#
# Run3 는 같은 레시피로 60 job 더 (WP 분기는 코드에 있으나 아직 제출 안 함).
# Quote wildcards!
for ERA in 2016preVFP 2016postVFP 2017 2018; do
  SKNano.py -a muid -i 'QCD_Pt-*MuEnriched*' -e $ERA -n 40 --tag muid1
  SKNano.py -a muid -i TTLJ_powheg           -e $ERA -n 40 --tag muid1
  SKNano.py -a muid -i DYJets                -e $ERA -n 40 --tag muid1
  SKNano.py -a muid -i DYJets10to50          -e $ERA -n 40 --tag muid1
done
