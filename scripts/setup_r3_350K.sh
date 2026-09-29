#!/usr/bin/env bash
# Build r3 production inputs for 5 systems: identical to r2 except a new gen_seed; verify via tpr diff.
set -u
GMX=${GMX:-gmx}
B=/home/raza/p53_thermal_xai/data/structures
LOG=/home/raza/p53_thermal_xai/analysis/r3_setup; mkdir -p "$LOG"
SEEDLOG=$LOG/r3_seeds.tsv; printf "System\tDir\tr1_seed\tr2_seed\tr3_seed\n" > "$SEEDLOG"
used=$(grep -h -E '^\s*gen[-_]seed\s*=' "$B"/c36_build/*.mdp "$B"/c36_prod/*.mdp "$B"/c36_temps/*.mdp 2>/dev/null | sed 's/;.*//' | awk -F= '{gsub(/ /,"",$2); print $2}' | sort -u)
for s in P151S M237I G245S R249S R248W; do
  echo "================ $s"
  dir=""; for d in c36_prod c36_build; do [ -f "$B/$d/${s}_r2.tpr" ] && { dir="$B/$d"; break; }; done
  mdp2=$(ls "$B"/c36_prod/md_${s}_r2.mdp "$B"/c36_build/md_${s}_r2.mdp 2>/dev/null | head -1)
  top="$B/c36_build/${s}.top"
  [ -z "$dir" ] && { echo "ERROR: ${s}_r2.tpr not found"; continue; }
  [ -z "$mdp2" ] && { echo "ERROR: md_${s}_r2.mdp not found"; continue; }
  [ -f "$top" ] || { echo "ERROR: $top not found"; continue; }
  [ -f "$dir/${s}_r3.tpr" ] && { echo "SKIP: $dir/${s}_r3.tpr already exists"; continue; }
  while :; do seed=$(shuf -i 100000-999999 -n 1); grep -qx "$seed" <<<"$used" || break; done
  used=$(printf "%s\n%s" "$used" "$seed")
  mdp3="$(dirname "$mdp2")/md_${s}_r3.mdp"
  sed -E "s/^([[:space:]]*gen[-_]seed[[:space:]]*=[[:space:]]*)[-0-9]+/\1${seed}/" "$mdp2" > "$mdp3"
  echo "-- mdp diff r2 -> r3 (must be gen_seed only):"; diff "$mdp2" "$mdp3"
  grep -qE '^\s*gen[-_]vel\s*=\s*yes' "$mdp3" && grep -qE '^\s*continuation\s*=\s*no' "$mdp3" || { echo "ERROR: gen_vel/continuation not as expected in $mdp3"; continue; }
  s1=$(grep -E '^\s*gen[-_]seed' "$(dirname "$mdp2")/md_${s}_r1.mdp" 2>/dev/null | sed 's/;.*//' | awk -F= '{gsub(/ /,"",$2);print $2}')
  s2=$(grep -E '^\s*gen[-_]seed' "$mdp2" | sed 's/;.*//' | awk -F= '{gsub(/ /,"",$2);print $2}')
  printf "%s\t%s\t%s\t%s\t%s\n" "$s" "$(basename "$dir")" "${s1:-NA}" "${s2:-NA}" "$seed" >> "$SEEDLOG"
  ndx=""; for n in "$B/c36_build/${s}.ndx" "$B/c36_build/index_${s}.ndx" "$dir/${s}.ndx"; do [ -f "$n" ] && { ndx="-n $n"; break; }; done
  ( cd "$B/c36_build" && $GMX grompp -f "$mdp3" -c "$dir/${s}_r2.tpr" -p "$top" $ndx \
      -o "$dir/${s}_r3.tpr" -po "$LOG/${s}_r3_mdout.mdp" > "$LOG/${s}_r3_grompp.log" 2>&1 )
  if [ ! -f "$dir/${s}_r3.tpr" ]; then echo "ERROR: grompp failed for $s:"; grep -E -A4 "ERROR|Fatal" "$LOG/${s}_r3_grompp.log" | head -30; continue; fi
  grep -E "^WARNING" -A3 "$LOG/${s}_r3_grompp.log" | head -20
  $GMX check -s1 "$dir/${s}_r2.tpr" -s2 "$dir/${s}_r3.tpr" > "$LOG/${s}_r2_vs_r3_tprdiff.txt" 2>&1
  nv=$(grep -cE '^\s*v\[' "$LOG/${s}_r2_vs_r3_tprdiff.txt")
  echo "-- tpr r2 vs r3: velocity lines differing = $nv (expected >0); non-velocity differences:"
  grep -vE '^\s*v\[' "$LOG/${s}_r2_vs_r3_tprdiff.txt" | grep -E ' - |differ|!=|<|>' | grep -viE 'GROMACS|reading|note|command line|:-\)' | head -20
  echo "-- OK: $dir/${s}_r3.tpr"
done
echo; column -t "$SEEDLOG"
echo; ls -la "$B"/c36_prod/*_r3.tpr "$B"/c36_build/*_r3.tpr 2>/dev/null
