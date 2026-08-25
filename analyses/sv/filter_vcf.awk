BEGIN { FS = OFS = "\t" }
FNR == NR { excl[$1] = 1; next }
/^##/ { print; next }
/^#CHROM/ {
    for (i = 1; i <= 9; i++) keep[i] = 1
    for (i = 10; i <= NF; i++) keep[i] = (($i) in excl) ? 0 : 1
    first = 1; line = ""
    for (i = 1; i <= NF; i++) {
        if (keep[i]) { if (!first) line = line OFS; line = line $i; first = 0 }
    }
    print line
    next
}
{
    first = 1; line = ""
    for (i = 1; i <= NF; i++) {
        if (keep[i]) { if (!first) line = line OFS; line = line $i; first = 0 }
    }
    print line
}
