#!/bin/bash
# Wtag v1 post-processing: hadd condor outputs (2022 + 2022EE) into
# /data6/.../wtag/SKNanoAnalyzer/wtag/2022_v1/{TTLJ,WW,WZ}.root. Run: rot && bash merge_v1.sh
set -e
SRC=/gv0/Users/snuintern2/SKNanoOutput/Wtag
DST=/data6/Users/snuintern2/wtag/SKNanoAnalyzer/wtag/2022_v1
ERAS="2022_v1 2022EE_v1"
mkdir -p $DST

paths() { for e in $ERAS; do for p in "$@"; do echo $SRC/$e/$p; done; done; }

hadd -f $DST/TTLJ.root $(paths 'TTLJ_powheg.root')
hadd -f $DST/WW.root   $(paths 'WW_pythia.root')
hadd -f $DST/WZ.root   $(paths 'WZ_pythia.root')

echo "=== merged into $DST ==="
ls -lh $DST
