# Wtag v2b: adds the Branch/{W,QCD}_* raw-FatJet-branch histograms for the
#   merged-W vs QCD shape comparison (PLAN_branch_compare.md). The `jets` tree
#   and every v2 histogram are unchanged, so v2 stays reproducible.
#   Outputs → /gv0/Users/snuintern2/SKNanoOutput/Wtag/<era>_v2b/
# (v2 [tag v2]: per-AK8-jet ntuple production for the ML W tagger — PLAN_v2.md.
#  v1 [tag v1]: TTLJ/WW/WZ only, histogram study — results in PLAN.md)
#   signal W jets: TTLJ/TTJJ (powheg) + WW/WZ (pythia); background: inclusive
#   QCD_PT bins (registered 2026-07-17 with exact Runs-tree sums).
for ERA in 2022 2022EE; do
  SKNano.py -a Wtag -i "TTLJ_powheg" -e $ERA -n 40 --tag v2b
  SKNano.py -a Wtag -i "TTJJ_powheg" -e $ERA -n 40 --tag v2b
  SKNano.py -a Wtag -i "WW_pythia"   -e $ERA -n 40 --tag v2b
  SKNano.py -a Wtag -i "WZ_pythia"   -e $ERA -n 40 --tag v2b
  SKNano.py -a Wtag -i 'QCD_Pt-170to300'   -e $ERA -n 40 --tag v2b
  SKNano.py -a Wtag -i 'QCD_Pt-300to470'   -e $ERA -n 40 --tag v2b
  SKNano.py -a Wtag -i 'QCD_Pt-470to600'   -e $ERA -n 40 --tag v2b
  SKNano.py -a Wtag -i 'QCD_Pt-600to800'   -e $ERA -n 40 --tag v2b
  SKNano.py -a Wtag -i 'QCD_Pt-800to1000'  -e $ERA -n 40 --tag v2b
  SKNano.py -a Wtag -i 'QCD_Pt-1000to1400' -e $ERA -n 40 --tag v2b
  SKNano.py -a Wtag -i 'QCD_Pt-1400to1800' -e $ERA -n 40 --tag v2b
  SKNano.py -a Wtag -i 'QCD_Pt-1800to2400' -e $ERA -n 40 --tag v2b
  SKNano.py -a Wtag -i 'QCD_Pt-2400to3200' -e $ERA -n 40 --tag v2b
done
