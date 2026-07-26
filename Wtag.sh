# Wtag v2d: v2c plus the branches needed to benchmark the *mass-decorrelated*
#   ParticleNet fairly (PLAN_pnet_benchmark.md). Everything else is unchanged,
#   so v2d reproduces v2c's histograms exactly.
#   - `jets` tree gains pnet_qcd, pnet_xccvsqcd, pnet_xbbvsqcd. Benchmarking our
#     MD tagger against XqqVsQCD alone is unfair to ParticleNet: Xqq is the
#     light-quark head and ~half of hadronic W decays contain charm (W->cs), and
#     measured, xcc out-rejects xqq on these W jets (6.5 vs 5.4 at eps_S=0.5).
#     The fair discriminant is (Xqq+Xcc+Xbb)/(...+QCD), which needs the raw QCD
#     probability to invert the stored ratios: X = QCD*r/(1-r).
#   - Branch/ + BranchNoMSD/ also gain tau31/tau41/tau42; v2c got those three
#     offline from the tree instead (agreed to 0.0000 in S), so this only makes
#     a future rerun self-contained.
#   Outputs -> /gv0/Users/snuintern2/SKNanoOutput/Wtag/<era>_v2d/
# (v2c [tag v2c]: n3b1 axis + padded score axis + BranchNoMSD/ — the numbers in
#  PLAN_branch_compare.md.
#  v2b [tag v2b]: first Branch/ pass — its n3b1/PNet_x* and
#  all-conditional-on-m_SD numbers are superseded by v2c.
#  v2 [tag v2]: per-AK8-jet ntuple production for the ML W tagger — PLAN_v2.md.
#  v1 [tag v1]: TTLJ/WW/WZ only, histogram study — results in PLAN.md)
#   signal W jets: TTLJ/TTJJ (powheg) + WW/WZ (pythia); background: inclusive
#   QCD_PT bins (registered 2026-07-17 with exact Runs-tree sums).
for ERA in 2022 2022EE; do
  SKNano.py -a Wtag -i "TTLJ_powheg" -e $ERA -n 40 --tag v2d
  SKNano.py -a Wtag -i "TTJJ_powheg" -e $ERA -n 40 --tag v2d
  SKNano.py -a Wtag -i "WW_pythia"   -e $ERA -n 40 --tag v2d
  SKNano.py -a Wtag -i "WZ_pythia"   -e $ERA -n 40 --tag v2d
  SKNano.py -a Wtag -i 'QCD_Pt-170to300'   -e $ERA -n 40 --tag v2d
  SKNano.py -a Wtag -i 'QCD_Pt-300to470'   -e $ERA -n 40 --tag v2d
  SKNano.py -a Wtag -i 'QCD_Pt-470to600'   -e $ERA -n 40 --tag v2d
  SKNano.py -a Wtag -i 'QCD_Pt-600to800'   -e $ERA -n 40 --tag v2d
  SKNano.py -a Wtag -i 'QCD_Pt-800to1000'  -e $ERA -n 40 --tag v2d
  SKNano.py -a Wtag -i 'QCD_Pt-1000to1400' -e $ERA -n 40 --tag v2d
  SKNano.py -a Wtag -i 'QCD_Pt-1400to1800' -e $ERA -n 40 --tag v2d
  SKNano.py -a Wtag -i 'QCD_Pt-1800to2400' -e $ERA -n 40 --tag v2d
  SKNano.py -a Wtag -i 'QCD_Pt-2400to3200' -e $ERA -n 40 --tag v2d
done
