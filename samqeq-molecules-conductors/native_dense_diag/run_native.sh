#!/bin/bash
# Native water (Section sec:water): run 0 on the NPT frame at three settings, then reconstruct each.
#   a6p  : pppm 1e-6, g_ewald pinned at the production value 0.33957131   (the operator production applies)
#   a8p  : pppm 1e-8, same pin                                            (particle-mesh limit of the validation)
#   a6u  : pppm 1e-6, g_ewald retuned by PPPM on this frame               (what a fresh run 0 would apply)
# Usage: bash run_native.sh <frame.data>
set -u
cd "$(dirname "$0")"
FRAME=${1:?frame data file}
LMP=/home/tpascal/codes/lammps/lammps-dev2/build/lmp_parallel
export OMP_NUM_THREADS=1
run1() {  # tag acc gew mesh
  local tag=$1 acc=$2 gew=$3 mesh=$4
  $LMP -in in.nat_run0 -var start file -var data "$FRAME" -var acc $acc -var gew $gew -var mesh $mesh -var tag $tag -log log.$tag > out.$tag 2>&1
  echo "rc=$? $tag"; grep -h "G vector\|grid = \|SOLVER-DIAG\|RUN0_DONE\|ERROR\|WARNING: samqeq\|did not converge" out.$tag | cut -c1-160
}
run1 nat_a6p 1.0e-6 0.33957131 36
run1 nat_a8p 1.0e-8 0.33957131 0
run1 nat_a6u 1.0e-6 0 0
export OMP_NUM_THREADS=4
for tag in nat_a6p nat_a8p nat_a6u; do
  G=$(grep -m1 "G vector" out.$tag | awk '{print $NF}')
  echo "=== reconstruct $tag  g_ewald=$G  rs=0.37933"
  python3 hd_reconstruct.py --dump frame_$tag.dump --param spcfq_expt.param --kernel cbrt --gewald $G --rs 0.37933 \
      --swb 10.0 --solve-types 1 2 --out rec_$tag.json --label "native NPT frame $tag" 2>&1 | grep -v "^closest"
done
