# Extract SV breakpoints from the merged SV VCF (CN1 coordinates).
# Output BED6: chrom, start(0-based), end(0-based), id, svtype, side
#   DEL/DUP/INV -> 2 breakpoints (POS=L, END=R)
#   INS         -> 1 breakpoint (POS=S)
#   TRA         -> 2 breakpoints (POS=L on CHROM, END=R on CHR2)
BEGIN { FS = OFS = "\t"; idx = 0 }
/^#/ { next }
{
  svtype = ""; end = 0; chr2 = ""
  n = split($8, f, ";")
  for (i = 1; i <= n; i++) {
    if (f[i] ~ /^SVTYPE=/) { svtype = substr(f[i], 8) }
    else if (f[i] ~ /^END=/) { end = substr(f[i], 5) + 0 }
    else if (f[i] ~ /^CHR2=/) { chr2 = substr(f[i], 6) }
  }
  pos = $2 + 0
  id = "SV" idx
  if (svtype == "INS") {
    print $1, pos-1, pos, id ":S", svtype, "S"
  } else if (svtype == "TRA") {
    print $1, pos-1, pos, id ":L", svtype, "L"
    print chr2, end-1, end, id ":R", svtype, "R"
  } else {
    print $1, pos-1, pos, id ":L", svtype, "L"
    if (end > pos) print $1, end-1, end, id ":R", svtype, "R"
  }
  idx++
}
