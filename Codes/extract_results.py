# -*- coding: utf-8 -*-
"""
extract_results.py - pull the force-displacement, energy and lateral-strain data out of an
Abaqus/Standard .odb written by build_planestrain_run.py. Works with Abaqus Python 2.7.

Run (from the folder holding the .odb):
    abaqus python extract_results.py PS_v5_C6_dqs_m025.odb
or  abaqus python extract_results.py PS_v5_C6_dqs_m025.odb my_prefix
In Abaqus/CAE (File > Run Script) set ODB_PATH below instead.

Writes <prefix>_history.csv and <prefix>_lateral.csv, and prints a summary.
Plot them with plot_results.py in an ordinary Python with matplotlib.
"""
from __future__ import print_function, division
import sys
import os
import math
from odbAccess import openOdb

# ----------------------------------------------------------------- settings (match the run)
ODB_PATH = 'PS_v5_C6_dqs_m025.odb'   # used only when no command-line argument is given
STEP = 'Step-1'
HC = 'HC-1'               # honeycomb instance name
AXIS = 2                  # 2: load along Y (U2, RF2); 1: load along X
H = 57.1577               # mm, patch height along the load axis (RP top - gap, from the .dat)
W = 56.0                  # mm, patch width (5 * 11 + 1)
DEPTH = 1.0               # mm, plane strain depth
GAP_TOP = 0.05            # mm, initial gap below the top platen
MASS_T = 6.9541739e-07    # tonne, honeycomb mass from the .dat (for DEPTH = 1 mm)
STRAIN_LEVELS = (0.05, 0.10, 0.20, 0.30)   # energy absorption is reported up to these strains
BAD = 1.0e20              # history values above this are treated as garbage


def interp(xs, ys, x):
    """Linear interpolation of ys(xs) at x, clamped at the ends."""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    lo, hi = 0, len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    d = xs[hi] - xs[lo]
    f = (x - xs[lo]) / d if d != 0.0 else 0.0
    return ys[lo] + f * (ys[hi] - ys[lo])


def clean(v):
    return v if abs(v) < BAD else float('nan')


def ratio(a, b):
    if b is None or b != b or a != a or b <= 0.0:
        return float('nan')
    return a / b


def get_series(reg, name):
    data = reg.historyOutputs[name].data
    return [float(p[0]) for p in data], [float(p[1]) for p in data]


def main():
    odb_path = sys.argv[1] if len(sys.argv) > 1 else ODB_PATH
    prefix = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(odb_path)[0]
    odb = openOdb(path=odb_path, readOnly=True)
    step = odb.steps[STEP]
    ax = str(AXIS)

    # --- the two platen reference points: the driven one has the larger displacement
    cands = []
    for key, reg in step.historyRegions.items():
        ho = reg.historyOutputs
        if ('U' + ax) in ho and ('RF' + ax) in ho:
            tt, uu = get_series(reg, 'U' + ax)
            cands.append((abs(uu[-1]), key, reg))
    if len(cands) < 2:
        raise Exception('expected two reference points with U%s and RF%s, found %d' % (ax, ax, len(cands)))
    cands.sort(key=lambda c: -c[0])
    top_key, top = cands[0][1], cands[0][2]
    bot_key, bot = cands[1][1], cands[1][2]
    print('top platen history region   : %s' % top_key)
    print('bottom platen history region: %s' % bot_key)

    t, u_top = get_series(top, 'U' + ax)
    tf, f_top = get_series(top, 'RF' + ax)
    tb, f_bot = get_series(bot, 'RF' + ax)

    # --- energies (whole model)
    en = {}
    akeys = [k for k in step.historyRegions.keys() if k.startswith('Assembly')]
    if akeys:
        reg = step.historyRegions[akeys[0]]
        for name in ('ALLIE', 'ALLKE', 'ALLAE', 'ALLSD', 'ALLWK', 'ETOTAL'):
            if name in reg.historyOutputs:
                en[name] = get_series(reg, name)
    n_bad = 0
    for name in en:
        n_bad += len([v for v in en[name][1] if abs(v) >= BAD])
    if n_bad:
        print('NOTE: %d energy values above %g were ignored (garbage in the history output)' % (n_bad, BAD))

    def e_at(name, time):
        if name not in en:
            return float('nan')
        xs, ys = en[name]
        good = [(x, y) for x, y in zip(xs, ys) if abs(y) < BAD]
        if not good:
            return float('nan')
        return interp([g[0] for g in good], [g[1] for g in good], time)

    # --- history table
    rows = []
    cum = 0.0
    prev_d = 0.0
    prev_f = 0.0
    e_at_strain = {}
    for i in range(len(t)):
        tm = t[i]
        comp = max(0.0, abs(u_top[i]) - GAP_TOP)        # mm of patch compression
        eps = comp / H
        ft = f_top[i]
        fb = interp(tb, f_bot, tm)
        stress = abs(ft) / (W * DEPTH)                  # MPa (force for DEPTH is per mm of depth)
        cum += 0.5 * (abs(ft) + prev_f) * (comp - prev_d) if i > 0 else 0.0
        prev_d, prev_f = comp, abs(ft)
        ie, ke, ae = e_at('ALLIE', tm), e_at('ALLKE', tm), e_at('ALLAE', tm)
        sd, wk, et = e_at('ALLSD', tm), e_at('ALLWK', tm), e_at('ETOTAL', tm)
        for lv in STRAIN_LEVELS:
            if lv not in e_at_strain and eps >= lv:
                e_at_strain[lv] = cum
        rows.append((tm, u_top[i], eps, ft, fb, ft + fb, stress, ie, ke, ae, sd, wk, et,
                     ratio(ke, ie), ratio(ae, ie), ratio(sd, ie),
                     cum, cum / (W * H * DEPTH), cum / (MASS_T * 1.0e6) / 1000.0))
    head = ('time_s,U_top_mm,strain,RF_top_N,RF_bot_N,RF_sum_N,stress_MPa,ALLIE_mJ,ALLKE_mJ,ALLAE_mJ,'
            'ALLSD_mJ,ALLWK_mJ,ETOTAL_mJ,KE_over_IE,AE_over_IE,SD_over_IE,'
            'energy_mJ,energy_per_volume_MJ_m3,energy_per_mass_J_g')
    path_h = prefix + '_history.csv'
    f = open(path_h, 'w')
    f.write(head + '\n')
    for r in rows:
        f.write(','.join('%.8g' % v for v in r) + '\n')
    f.close()

    # --- lateral strain from the honeycomb edge nodes, frame by frame
    inst = odb.rootAssembly.instances[HC]
    lat = 1 - (AXIS - 1)            # transverse axis index: 1 (X) when loading along Y, 0 (Y) otherwise
    coords = [(n.label, n.coordinates[lat]) for n in inst.nodes]
    cmin = min(c for l, c in coords)
    cmax = max(c for l, c in coords)
    tol = 1.0e-6
    left = set(l for l, c in coords if abs(c - cmin) < tol)
    right = set(l for l, c in coords if abs(c - cmax) < tol)
    f = open(prefix + '_lateral.csv', 'w')
    f.write('time_s,strain_load,strain_transverse,poisson_ratio_eff\n')
    nu_last = float('nan')
    for fr in step.frames:
        tm = fr.frameValue
        if 'U' not in fr.fieldOutputs:
            continue
        uf = fr.fieldOutputs['U'].getSubset(region=inst)
        sl = sr = 0.0
        nl = nr = 0
        for v in uf.values:
            if v.nodeLabel in left:
                sl += v.data[lat]
                nl += 1
            elif v.nodeLabel in right:
                sr += v.data[lat]
                nr += 1
        if nl == 0 or nr == 0:
            continue
        width = (cmax + sr / nr) - (cmin + sl / nl)
        ex = width / (cmax - cmin) - 1.0
        ey = max(0.0, abs(interp(t, u_top, tm)) - GAP_TOP) / H
        nu = ex / ey if ey > 1.0e-4 else float('nan')    # ey is the compressive strain magnitude
        nu_last = nu
        f.write('%.8g,%.8g,%.8g,%.8g\n' % (tm, ey, ex, nu))
    f.close()

    # --- summary
    stress = [r[6] for r in rows]
    smax = max(stress)
    imax = stress.index(smax)
    print('')
    print('==================== SUMMARY ====================')
    print('history points        : %d, step time to %.2f s' % (len(rows), t[-1]))
    print('end strain            : %.2f %%' % (100.0 * rows[-1][2]))
    print('peak stress           : %.4f MPa at strain %.2f %%' % (smax, 100.0 * rows[imax][2]))
    fmax = max(abs(x) for x in f_top)
    bal = max(abs(r[5]) for r in rows)
    print('platen force balance  : max |RF_top + RF_bot| = %.3g N (peak force %.3g N)' % (bal, fmax))
    for name, col in (('KE/IE', 13), ('AE/IE', 14), ('SD/IE', 15)):
        vals = [r[col] for r in rows if r[col] == r[col] and r[2] > 0.005]
        if vals:
            print('max %-6s after 0.5%%  : %.4f' % (name, max(vals)))
    for lv in STRAIN_LEVELS:
        if lv in e_at_strain:
            print('energy to %4.0f %% strain : %.3f mJ  (%.4f MJ/m3, %.4f J/g)'
                  % (100 * lv, e_at_strain[lv], e_at_strain[lv] / (W * H * DEPTH),
                     e_at_strain[lv] / (MASS_T * 1.0e6) / 1000.0))
        else:
            print('energy to %4.0f %% strain : not reached' % (100 * lv))
    print('effective Poisson ratio at the last frame: %.3f' % nu_last)
    print('wrote %s and %s' % (path_h, prefix + '_lateral.csv'))
    odb.close()


main()
