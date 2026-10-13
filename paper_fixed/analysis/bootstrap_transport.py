import numpy as np, sys
sys.path.insert(0,'/home/IITB/multiscale-mechanics/ashwani12/polymer/paper_fixed/analysis')
import traj_io as T
from scipy.stats import linregress

d=T.load(); pos,mol,typ=d['pos'],d['mol'],d['typ']
sel_li=(typ==T.T_LI); sel_tf=np.isin(typ,T.T_TFSI)
sel=sel_li|sel_tf|np.isin(typ,(T.T_PEO_C,T.T_PEO_H,T.T_PEO_O))
pu=T.unwrap_xy(pos)
w=T.masses_of(typ).copy(); w[~sel]=0.0
pu=T.remove_com_drift(pu,w)
com_tf,mids=T.molecule_com(pu,mol,typ,sel_tf)
nf=pos.shape[0]; maxlag=nf//2; lag_ns=np.arange(maxlag)*T.DT_FRAME_NS
m=(lag_ns>=3.0)&(lag_ns<=10.0); x=lag_ns[m]*1e-9

def perion(a):
    n=a.shape[1]; per=np.zeros((maxlag,n))
    for lg in range(1,maxlag):
        dd=a[lg:]-a[:-lg]; per[lg]=(dd**2).sum(-1).mean(0)
    return per
pl=perion(pu[:,sel_li,:2]); pt=perion(com_tf[:,:,:2])
print(f"Li ions={pl.shape[1]}  TFSI molecules={pt.shape[1]}")
D=lambda p: linregress(x,p[m]*1e-16).slope/4.0

Dl,Dt=D(pl.mean(1)),D(pt.mean(1))
rng=np.random.default_rng(1); N=2000
bl=np.array([D(pl[:,rng.integers(0,pl.shape[1],pl.shape[1])].mean(1)) for _ in range(N)])
bt=np.array([D(pt[:,rng.integers(0,pt.shape[1],pt.shape[1])].mean(1)) for _ in range(N)])
rat=bt/bl; tp=bl/(bl+bt)
# sigma scales linearly with (D_Li*z^2*n_Li + D_TFSI*z^2*n_TFSI); reuse paper's sigma/D scaling
sig0=6.029e-03
sig=sig0*( (bl*62+bt*pt.shape[1]) / (Dl*62+Dt*pt.shape[1]) )

def rep(lab,v,b,paper,fmt="{:.3e}"):
    print(f"{lab:22s} {fmt.format(v)}  bootstrap +/- {fmt.format(b.std())}"
          f"   95% [{fmt.format(np.percentile(b,2.5))}, {fmt.format(np.percentile(b,97.5))}]"
          f"   | paper: {paper}")
print()
rep("D_Li  (in-plane)",Dl,bl,"4.401e-07 +/- 1.2e-09")
rep("D_TFSI(in-plane)",Dt,bt,"1.167e-06 +/- 2.6e-09")
rep("D_TFSI/D_Li",Dt/Dl,rat,"2.65 (no error quoted)","{:.3f}")
rep("t+ (ideal NE)",Dl/(Dl+Dt),tp,"0.274 +/- 0.001","{:.4f}")
rep("sigma NE (S/cm)",sig0,sig,"6.029e-03")
print(f"\nP(D_TFSI > D_Li) over bootstrap = {(bt>bl).mean()*100:.1f}%   <- ordering robustness")
