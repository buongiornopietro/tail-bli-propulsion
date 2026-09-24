import re, sys, json, numpy as np
def leggi(fn, vec=False):
    s = open(fn).read(); i = s.index("internalField"); seg = s[i:s.index("boundaryField")]
    m = re.search(r"(\d+)\s*\(", seg); n = int(m.group(1)); body = seg[m.end():]
    if vec:
        return np.array([l.split() for l in re.findall(r"\(([^()]*)\)", body)[:n]], float)
    return np.array(body.split(")")[0].split()[:n], float)
t = sys.argv[1]
C = leggi(f"{t}/C", True); nut = leggi(f"{t}/nut")
x = C[:, 0]; r = np.hypot(C[:, 1], C[:, 2])
for a in np.arange(0.05, 1.0, 0.1):
    m = (np.abs(x - a) < 0.02) & (r < 0.12)
    print(f"x={a:.2f}  nut_max={nut[m].max():.2e}")
