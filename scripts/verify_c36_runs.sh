#!/usr/bin/env bash
# Verify force field, temperature, length and replica independence from actual run files.
set -u
GMX=${GMX:-gmx}
B=/home/raza/p53_thermal_xai/data/structures
OUT=/home/raza/p53_thermal_xai/analysis/verify
mkdir -p "$OUT"; TMP=$(mktemp -d)
TSV="$OUT/c36_run_verification.tsv"

lasttime() { awk '{for(i=1;i<=NF;i++) if($i=="time") v=$(i+1)} END{print v}'; }
ps2ns()    { awk -v t="$1" 'BEGIN{ if(t=="") print "NA"; else printf "%.2f", t/1000 }'; }

ff_class() {
  local d=$1 types cmap
  types=$(grep -o 'name="[^"]*",nameB' "$d" | sed 's/name="//;s/",nameB//' | sort -u | tr '\n' ' ')
  cmap=$(grep -ci 'cmap' "$d")
  if   grep -q 'opls_' <<<"$types"; then echo "OPLS-AA"
  elif grep -qw 'CT1' <<<"$types" && [ "$cmap" -gt 0 ]; then echo "CHARMM36-family(CMAP)"
  elif grep -qw 'CT1' <<<"$types"; then echo "CHARMM-types-NO-CMAP"
  else echo "OTHER:$(echo $types | cut -c1-30 | tr ' ' ',')"; fi
}
water() {
  local d=$1
  if   grep -q 'name="OT",nameB' "$d"; then echo "TIP3P(CHARMM-OT/HT)"
  elif grep -q 'name="OW' "$d"; then echo "OW-type"
  elif grep -q 'opls_11[12]' "$d"; then echo "TIP3P(OPLS)"
  else echo "unknown"; fi
}
tprinfo() {  # prints: FF Water Tcoupl RefT PlannedNs Tinit LdSeed VelMd5
  local d=$1 nsteps dt tinit reft tc ld vel
  nsteps=$(awk '$1=="nsteps"{print $3; exit}' "$d")
  dt=$(awk '($1=="delta-t"||$1=="delta_t"){print $3; exit}' "$d")
  tinit=$(awk '$1=="tinit"{print $3; exit}' "$d")
  reft=$(awk '$1=="ref-t:"||$1=="ref_t:"{$1=""; print; exit}' "$d" | xargs | tr ' ' '/')
  tc=$(awk '$1=="tcoupl"{print $3; exit}' "$d")
  ld=$(awk '($1=="ld-seed"||$1=="ld_seed"){print $3; exit}' "$d")
  vel=$(awk '/^ +v \(/{f=1;next} f&&/^ +v\[/{print;next} f{exit}' "$d")
  if [ -z "$vel" ]; then vel="NO_VELOCITIES"; else vel=$(md5sum <<<"$vel" | cut -c1-10); fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" "$(ff_class "$d")" "$(water "$d")" "${tc:-NA}" "${reft:-NA}" \
    "$(awk -v n="$nsteps" -v t="$dt" 'BEGIN{printf "%.2f", n*t/1000}')" "${tinit:-NA}" "${ld:-NA}" "$vel"
}

echo "################ PART A: 21-system CHARMM36m production runs ################"
printf "System\tRep\tFF_types\tWater\tTcoupl\tRef_T_K\tPlanned_ns\ttinit_ps\tld_seed\tVel_md5\tEDR_end_ns\tXTC_pbcfix_end_ns\tSource\n" > "$TSV"
systems=$(ls "$B"/c36_prod/*_r[12].tpr "$B"/c36_build/*_r[12].tpr 2>/dev/null | xargs -n1 basename | sed -E 's/_r[12]\.tpr$//' | sort -u)
for s in $systems; do
  for r in r1 r2; do
    tpr=""; src=""
    for dd in c36_prod c36_build; do [ -f "$B/$dd/${s}_${r}.tpr" ] && { tpr="$B/$dd/${s}_${r}.tpr"; src=$dd; break; }; done
    if [ -z "$tpr" ]; then printf "%s\t%s\tMISSING_TPR\n" "$s" "$r" >> "$TSV"; continue; fi
    dump="$TMP/${s}_${r}.txt"; $GMX dump -s "$tpr" > "$dump" 2>/dev/null
    edr="${tpr%.tpr}.edr"; xtc="$B/c36_pbcfix/${s}_${r}.xtc"
    e=""; [ -f "$edr" ] && e=$($GMX check -e "$edr" 2>&1 | grep -i "last energy frame" | lasttime)
    x=""; [ -f "$xtc" ] && x=$($GMX check -f "$xtc" 2>&1 | grep -i "^last frame" | lasttime)
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$s" "$r" "$(tprinfo "$dump")" "$(ps2ns "$e")" "$(ps2ns "$x")" "$src" >> "$TSV"
  done
done
column -t -s $'\t' "$TSV"

echo; echo "---- Replica independence (starting-velocity hash and ld-seed, r1 vs r2) ----"
awk -F'\t' 'NR>1{v[$1,$2]=$10; l[$1,$2]=$9; s[$1]=1}
  END{for(k in s){st=(v[k,"r1"]==v[k,"r2"])?"IDENTICAL_VELOCITIES(!)":"different";
      ls=(l[k,"r1"]==l[k,"r2"])?"same_ld_seed":"diff_ld_seed"; print k, st, ls}}' "$TSV" | sort | column -t

echo; echo "---- Consistency: distinct values across all 42 runs ----"
for c in 3:FF_types 4:Water 5:Tcoupl 6:Ref_T_K 7:Planned_ns 11:EDR_end_ns 12:XTC_pbcfix_end_ns; do
  i=${c%%:*}; n=${c#*:}; echo "$n: $(awk -F'\t' -v i=$i 'NR>1{print $i}' "$TSV" | sort | uniq -c | xargs)"; done

echo; echo "---- Force-field version (topology includes) ----"
find "$B" -maxdepth 2 -name "*.top" -print0 | xargs -0 grep -H -m1 'forcefield.itp' 2>/dev/null | sed "s|$B/||" | sort
echo; echo "---- Local .ff directories and their descriptions ----"
find /home/raza -maxdepth 6 -type d -name "*.ff" 2>/dev/null | grep -i -E 'charmm|opls|amber' | while read f; do
  echo "$f :: $(head -3 "$f/forcefield.doc" 2>/dev/null | tr '\n' ' ')"; done

echo; echo "---- Production .mdp files found (gen_vel / gen_seed / continuation / ref_t / nsteps) ----"
find "$B/c36_prod" "$B/c36_build" "$B/c36_temps" -maxdepth 1 -name "*.mdp" 2>/dev/null | sort | while read m; do
  echo "$(basename "$m"): $(grep -E '^\s*(gen[-_]vel|gen[-_]seed|continuation|ref[-_]t|nsteps|tcoupl)\s*=' "$m" | sed 's/;.*//' | xargs)"; done

echo; echo "################ PART B: old six-temperature long_simulations ################"
printf "Run\tFF_types\tWater\tTcoupl\tRef_T_K\tPlanned_ns\ttinit_ps\tld_seed\tVel_md5\tEDR_end_ns\n" > "$OUT/long_sim_verification.tsv"
for T in 248 273 310 350 400 450; do
  tpr="$B/long_simulations/md_100ns_${T}K.tpr"; [ -f "$tpr" ] || { echo "missing $tpr"; continue; }
  dump="$TMP/long_${T}.txt"; $GMX dump -s "$tpr" > "$dump" 2>/dev/null
  e=$($GMX check -e "${tpr%.tpr}.edr" 2>&1 | grep -i "last energy frame" | lasttime)
  printf "md_100ns_%sK\t%s\t%s\n" "$T" "$(tprinfo "$dump")" "$(ps2ns "$e")" >> "$OUT/long_sim_verification.tsv"
done
column -t -s $'\t' "$OUT/long_sim_verification.tsv"
echo "Top files near long_simulations:"
find "$B" "$B/long_simulations" -maxdepth 1 -name "*.top" -print0 | xargs -0 grep -H -m1 'forcefield.itp' 2>/dev/null
rm -rf "$TMP"
echo; echo "Saved: $TSV and $OUT/long_sim_verification.tsv"
