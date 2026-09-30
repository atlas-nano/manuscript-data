import numpy as np, re, glob, os
KB=8.617333262e-5; CAP=0.5*384*KB
num=re.compile(r'^-?\d+(\.\d+)?([eE][-+]?\d+)?$')
def load(fn):
    d=[]
    for l in open(fn):
        s=l.split()
        if len(s)==5 and all(num.match(x) for x in s): d.append([float(x) for x in s])
    return np.array(d)
res={}
for rt,lab in (("50","pair 5 A (shipped)"),("80","pair 8 A")):
    ks=[]
    for s in ("731","4127","9013"):
        f=f"out.A.s{s}_r{rt}"
        if not os.path.exists(f): print(f"  MISSING {f}"); continue
        a=load(f)
        if len(a)<100: print(f"  SHORT {f}: {len(a)} rows"); continue
        step,T,pe,et,es=a.T
        t=step*1e-6; m=t>=0.02
        k=np.polyfit(t[m],et[m],1)[0]/CAP
        h=m.sum()//2
        k1=np.polyfit(t[m][:h],et[m][:h],1)[0]/CAP
        k2=np.polyfit(t[m][h:],et[m][h:],1)[0]/CAP
        ks.append(k)
        print(f"  {lab:20s} seed {s:>5s}: drift {k:8.1f} K/ns  (halves {k1:7.1f}/{k2:7.1f})  "
              f"sd(E) {et[m].std():.4f} eV  <T> {T[m].mean():6.1f} K")
    res[lab]=np.array(ks)
print()
for lab,v in res.items():
    if len(v): print(f"{lab:20s} drift = {v.mean():8.1f} +/- {v.std(ddof=1) if len(v)>1 else float('nan'):.1f} K/ns   (n={len(v)} seeds)")
