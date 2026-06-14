#!/usr/bin/env bash
in="$1"
out="$2"

awk -v OFS="\t" '{
  id = sprintf("peak_%07d", NR);
  strand = $5;
  if (strand=="+" || strand=="0") strand_out="0";
  else if (strand=="-" || strand=="1") strand_out="1";
  else strand_out="0";
  print $1, $2, $3, id, 0, strand_out, $6, $4
}' "$in" > "$out"

