#!/usr/bin/env python3
"""Turn <prefix>_history.csv into a stress-strain curve, quality checks and a summary.

usage: python analyse_run.py T0p01_C2_history.csv --name T0p01_C2 --gaps 0.12 \
           --H 60.5 --W 66.706 --D 10 --mass-g 8.61154
Outputs: <name>_curve.csv, <name>_summary.txt, <name>_plots.png
"""
import argparse
import numpy as np
import pandas as pd
_trapz = getattr(np, 'trapezoid', None) or np.trapz
np.seterr(divide='ignore', invalid='ignore')
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("csv")
p.add_argument("--name", default="run")
p.add_argument("--gaps", type=float, default=0.12, help="top+bottom gap [mm], subtracted from |U2|")
p.add_argument("--H", type=float, default=60.5, help="patch height along loading [mm]")
p.add_argument("--W", type=float, default=66.706, help="width perpendicular to loading, in plane [mm]")
p.add_argument("--D", type=float, default=10.0, help="depth [mm]")
p.add_argument("--mass-g", type=float, default=8.61154, help="patch mass [g] for specific energy")
p.add_argument("--slope-window", type=float, nargs=2, default=[0.002, 0.004], help="strain window for the initial slope")
p.add_argument("--sea-strain", type=float, default=0.03, help="strain for energy absorbed")
p.add_argument("--axis", type=int, default=2, choices=[1, 2, 3], help="loaded axis: 1, 2 or 3 (U2/RF2 is axis 2)")
a = p.parse_args()

df = pd.read_csv(a.csv)

# platen reference point: among regions with both U and RF on the loaded axis,
# take the one that moves most (the driven top platen; the bottom one is fixed)
UV, RV = "U%d" % a.axis, "RF%d" % a.axis
cand = [r for r, g in df.groupby("region") if {UV, RV} <= set(g["variable"])]
if not cand:
    raise SystemExit("no region with both %s and %s found" % (UV, RV))
def _travel(r):
    g = df[(df.region == r) & (df.variable == UV)]
    return float(np.abs(g["value"]).max()) if len(g) else 0.0
reg = max(cand, key=_travel)
print("platen region used: %s  (max |%s| = %.3f mm; %d candidate region(s))" % (reg, UV, _travel(reg), len(cand)))
u = df[(df.region == reg) & (df.variable == UV)].sort_values("time")
f = df[(df.region == reg) & (df.variable == RV)].sort_values("time")
t = f["time"].to_numpy()
U = np.abs(np.interp(t, u["time"], u["value"]))
F = np.abs(f["value"].to_numpy())

def energy(name):
    g = df[df.variable == name].sort_values("time")
    return np.interp(t, g["time"], g["value"]) if len(g) else np.full_like(t, np.nan)

KE, IE, AE, ET = (energy(n) for n in ("ALLKE", "ALLIE", "ALLAE", "ETOTAL"))

eps = np.clip(U - a.gaps, 0, None) / a.H
sig = F / (a.W * a.D)
out = pd.DataFrame(dict(time=t, U2_abs=U, strain=eps, stress_MPa=sig, force_N=F,
                        ALLKE=KE, ALLIE=IE, ALLAE=AE, ETOTAL=ET))
out.to_csv(f"{a.name}_curve.csv", index=False)

# --- quality gates
mask = IE > 0.05 * np.nanmax(IE)               # ignore the start, where IE is tiny
ke_ie = np.where(mask, KE / IE, np.nan)
ae_ie = np.where(mask, AE / IE, np.nan)
etot_ie = np.where(mask, np.abs(ET) / IE, np.nan)

# --- metrics
ip = int(np.nanargmax(sig))
lo, hi = a.slope_window
w = (eps >= lo) & (eps <= hi)
slope = np.polyfit(eps[w], sig[w], 1)[0] if w.sum() >= 3 else np.nan
sea_mask = eps <= a.sea_strain
work = _trapz(F[sea_mask], U[sea_mask]) if sea_mask.sum() > 2 else np.nan   # N*mm = mJ
lines = [
    f"run                         : {a.name}",
    f"max ALLKE/ALLIE (IE>5% max) : {np.nanmax(ke_ie):.3f}   (limit 0.05)",
    f"fraction of points <= 0.05  : {np.mean((KE / IE)[mask] <= 0.05):.2f}",
    f"max ALLAE/ALLIE             : {np.nanmax(ae_ie):.3f}   (limit 0.05)",
    f"max |ETOTAL|/ALLIE          : {np.nanmax(etot_ie):.4f}",
    f"peak stress [MPa]           : {sig[ip]:.5f}  at strain {eps[ip]:.4f}",
    f"initial slope [MPa] ({lo}-{hi}) : {slope:.3f}",
    f"energy to strain {a.sea_strain} [mJ] : {work:.3f}",
    f"specific energy [J/kg]      : {work / a.mass_g:.3f}",
    f"final strain                : {eps[-1]:.4f}",
]
open(f"{a.name}_summary.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(eps, sig); ax[0].set_xlabel("strain"); ax[0].set_ylabel("nominal stress [MPa]")
ax[0].set_title(a.name)
ax[1].plot(eps, ke_ie, label="ALLKE/ALLIE"); ax[1].plot(eps, ae_ie, label="ALLAE/ALLIE")
ax[1].axhline(0.05, color="k", ls="--", lw=0.8); ax[1].set_xlabel("strain"); ax[1].set_ylim(0, 0.5)
ax[1].legend()
fig.tight_layout(); fig.savefig(f"{a.name}_plots.png", dpi=150)
