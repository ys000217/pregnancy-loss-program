import sys
out = []
def check(name):
    try:
        m = __import__(name)
        out.append("%s: %s" % (name, getattr(m, '__version__', 'ok')))
    except Exception as e:
        out.append("%s: MISSING (%s)" % (name, e))
for n in ["numpy", "pandas", "scipy", "statsmodels", "polars", "pyarrow", "patsy"]:
    check(n)
out.append("python: " + sys.version.split()[0])
open(r"D:\ONT\_env_check.txt", "w", encoding="utf-8").write("\n".join(out))
print("\n".join(out))
