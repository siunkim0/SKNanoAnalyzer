#!/bin/bash
# Submit the HiggsBDT analyzer once per BDT model.
#
# The model is chosen at RUN TIME via the `model=<stem>` userflag, so you do
# NOT recompile to switch models — build once, then just rerun this script.
# Each run is tagged (--tag) so its output lands in a distinct
#   $SKNANO_OUTPUT/HiggsBDT/<era>_<tag>/
# folder and never overwrites another model's result.
#
# Add/remove a model by editing the MODELS list below. A stem here must have a
# matching <stem>.onnx + <stem>_features.txt under
#   data/Run3_v13_Run2_v9/<era>/BDT/HZZ4mu/

ERA=2018
NJOBS=40

MODELS=(
    bdt_v1_new
    bdt_v2_new
    bdt_v3_new
    bdt_v5_new
    bdt_v6_new
)

for MODEL in "${MODELS[@]}"; do
    echo "==== submitting HiggsBDT with model=${MODEL} (tag=${MODEL}) ===="
    SKNano.py -a HiggsBDT -i SingleMuon    -n $NJOBS -e $ERA --userflags model=${MODEL} --tag ${MODEL}
    SKNano.py -a HiggsBDT -i HiggsList.txt -n $NJOBS -e $ERA --userflags model=${MODEL} --tag ${MODEL}
done
