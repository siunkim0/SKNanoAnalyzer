#!/bin/bash
# Wtag v2 post-processing: hadd condor outputs (2022 + 2022EE) into
# /data6/.../wtag/SKNanoAnalyzer/wtag/2022_v2/. Trees hadd cleanly; `weight` is already
# xsec/sumsign-normalized per sample so QCD bins stitch by weight.
# Run: rot && bash merge_v2.sh
set -e
SRC=/gv0/Users/snuintern2/SKNanoOutput/Wtag
DST=/data6/Users/snuintern2/wtag/SKNanoAnalyzer/wtag/2022_v2
ERAS="2022_v2 2022EE_v2"
mkdir -p $DST

paths() { for e in $ERAS; do for p in "$@"; do echo $SRC/$e/$p; done; done; }

hadd -f $DST/TTLJ.root $(paths 'TTLJ_powheg.root')
hadd -f $DST/TTJJ.root $(paths 'TTJJ_powheg.root')
hadd -f $DST/WW.root   $(paths 'WW_pythia.root')
hadd -f $DST/WZ.root   $(paths 'WZ_pythia.root')
for BIN in 170to300 300to470 470to600 600to800 800to1000 \
           1000to1400 1400to1800 1800to2400 2400to3200; do
  hadd -f $DST/QCD_Pt-$BIN.root $(paths "QCD_Pt-$BIN.root")
done

echo "=== merged into $DST ==="
ls -lh $DST
