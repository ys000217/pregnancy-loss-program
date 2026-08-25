/^>/ {
  acc = $1; sub(/^>/, "", acc)
  name = ""
  if (acc ~ /^NC_0000/) {
    num = substr(acc, 4, 6) + 0
    if (num >= 1 && num <= 22) name = "chr" num
    else if (num == 23) name = "chrX"
    else if (num == 24) name = "chrY"
  } else if (acc ~ /^NC_012920/) {
    name = "chrM"
  }
  keep = (name != "")
  if (keep) print ">" name
  next
}
{ if (keep) print }
