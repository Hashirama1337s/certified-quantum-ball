#!/bin/bash
cd <LAB>/quantum/seats/build-4
export QB2B_NOPATCH=1
k=$(py -3.11 <MEMCAP> 8 py -3.11 -c "import json;print(json.load(open('pass_ckpt.json'))['k'])" 2>/dev/null | head -1)
[ -z "$k" ] && k=0
while true; do
  k1=$((k+240))
  py -3.11 <MEMCAP> 8 py -3.11 pass_b.py $k $k1 >> pass.log 2>&1
  rc=$?
  echo "chunk $k-$k1 rc=$rc $(date +%T)" >> pass.log
  if [ $rc -eq 0 ]; then echo PASS_DONE >> pass.log; break; fi
  if [ $rc -ne 3 ]; then echo PASS_FAILED >> pass.log; break; fi
  k=$k1
done
