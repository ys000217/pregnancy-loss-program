# 移除 INFO 中按样本数计算、筛选后已失真的字段: SUPP, SUPP_VEC, AC, AN
BEGIN { FS = OFS = "\t" }
/^##INFO=<ID=(SUPP|SUPP_VEC|AC|AN),/ { next }
/^#/ { print; next }
{
    n = split($8, f, ";")
    out = ""; first = 1
    for (i = 1; i <= n; i++) {
        key = f[i]; sub(/=.*/, "", key)
        if (key == "SUPP" || key == "SUPP_VEC" || key == "AC" || key == "AN") continue
        if (!first) out = out ";"
        out = out f[i]; first = 0
    }
    $8 = out
    print
}
