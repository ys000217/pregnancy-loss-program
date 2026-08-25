/^>/ {
  name = substr($1, 2)
  for (i = 1; i <= NF; i++) {
    if ($i ~ /^OriSeqID=/) { ori = $i; sub(/^OriSeqID=/, "", ori); name = ori }
  }
  print ">" name
  next
}
{ print }
