# -*- coding: utf-8 -*-
"""
Regular hexagonal honeycomb (RHC) patch, 2D planar, for Abaqus/CAE 2020 (Python 2.7).

Run:  File > Run Script...   or   abaqus cae script=build_rhc_patch.py

Geometry (pointy-top hexagons, load along y, origin at bottom-left):
  - D  : clear (hole) distance across flats, mm
  - T  : wall thickness between neighbouring holes, mm
  - pitch = D + T
  - NX : holes in the even rows (odd rows have NX-1), NY : number of rows (odd)
  - Even rows start at x = pitch/2; odd rows start at x = pitch.
  - Left/right edges follow the wall centrelines (boundary walls are T/2 thick).
  - Top and bottom are cut flat through the hole centres of the first and last row,
    leaving half-hole notches open to the platens.
  - Hole corners are filleted with radius R_FILLET. The two cut corners of each
    notch are left sharp.
"""
from __future__ import division, print_function
import math

# ----------------------------------------------------------------- parameters
D = 10.0
T = 1.0
R_FILLET = 0.5
NX = 5
NY = 7
MODEL_NAME = 'RHC_Patch'
PART_NAME = 'Core'
SHEET_SIZE = 200.0


# ------------------------------------------------------------ pure geometry
def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _norm(a):
    return math.hypot(a[0], a[1])


def _unit(a):
    n = _norm(a)
    return (a[0] / n, a[1] / n)


def _fillet_corner(prev_pt, v, next_pt, r):
    """Return (tangent point on prev edge, tangent point on next edge,
    arc centre, 'CCW' or 'CW')."""
    u1 = _unit(_sub(prev_pt, v))
    u2 = _unit(_sub(next_pt, v))
    cos_phi = max(-1.0, min(1.0, u1[0] * u2[0] + u1[1] * u2[1]))
    half = math.acos(cos_phi) / 2.0
    tan_len = r / math.tan(half)
    t1 = (v[0] + u1[0] * tan_len, v[1] + u1[1] * tan_len)
    t2 = (v[0] + u2[0] * tan_len, v[1] + u2[1] * tan_len)
    b = _unit((u1[0] + u2[0], u1[1] + u2[1]))
    dist = r / math.sin(half)
    c = (v[0] + b[0] * dist, v[1] + b[1] * dist)
    d_in = _sub(v, prev_pt)
    d_out = _sub(next_pt, v)
    cross = d_in[0] * d_out[1] - d_in[1] * d_out[0]
    return t1, t2, c, ('CCW' if cross > 0 else 'CW')


def _chain(points, flags, r):
    """Closed chain of (x, y) points. flags[i] True -> fillet vertex i.
    Returns list of segments:
      ('L', p1, p2)  or  ('A', centre, p1, p2, 'CCW'|'CW')"""
    n = len(points)
    entry = []   # point where the chain arrives at vertex i
    exit_ = []   # point where the chain leaves vertex i
    arcs = []
    for i in range(n):
        v = points[i]
        if flags[i]:
            t1, t2, c, d = _fillet_corner(points[i - 1], v, points[(i + 1) % n], r)
            entry.append(t1)
            exit_.append(t2)
            arcs.append((c, t1, t2, d))
        else:
            entry.append(v)
            exit_.append(v)
            arcs.append(None)
    segs = []
    for i in range(n):
        if arcs[i] is not None:
            c, t1, t2, d = arcs[i]
            segs.append(('A', c, t1, t2, d))
        j = (i + 1) % n
        if _norm(_sub(exit_[i], entry[j])) > 1e-9:
            segs.append(('L', exit_[i], entry[j]))
    return segs


def build_loops(d=D, t=T, r=R_FILLET, nx=NX, ny=NY):
    """Return (loops, width, height). loops[0] is the outer loop,
    the others are interior holes."""
    if ny % 2 != 1:
        raise ValueError('NY must be odd so that the first and last rows are both even-type rows')
    pitch = d + t
    big_r = pitch / math.sqrt(3.0)      # wall-centreline hexagon circumradius
    rh = d / math.sqrt(3.0)             # hole circumradius
    hx = d / 2.0                        # hole half-width across flats
    row_pitch = 1.5 * big_r
    width = nx * pitch
    height = (ny - 1) * row_pitch
    even_x = [pitch / 2.0 + i * pitch for i in range(nx)]
    odd_x = [pitch + i * pitch for i in range(nx - 1)]

    # left edge, going up
    left = [(0.0, 0.0)]
    for j in range(0, ny - 1, 2):
        cy = j * row_pitch
        left += [(0.0, cy + big_r / 2.0),
                 (pitch / 2.0, cy + big_r),
                 (pitch / 2.0, cy + 2.0 * big_r),
                 (0.0, cy + 2.5 * big_r)]
    left.append((0.0, height))
    right = [(width - x, y) for (x, y) in left]

    pts = []
    flg = []

    def add(p, f):
        pts.append(p)
        flg.append(f)

    add((0.0, 0.0), False)
    for cx in even_x:                                   # bottom, left to right
        add((cx - hx, 0.0), False)
        add((cx - hx, rh / 2.0), True)
        add((cx, rh), True)
        add((cx + hx, rh / 2.0), True)
        add((cx + hx, 0.0), False)
    add((width, 0.0), False)
    for p in right[1:]:                                 # right side, up
        add(p, False)
    for cx in reversed(even_x):                         # top, right to left
        add((cx + hx, height), False)
        add((cx + hx, height - rh / 2.0), True)
        add((cx, height - rh), True)
        add((cx - hx, height - rh / 2.0), True)
        add((cx - hx, height), False)
    for p in list(reversed(left))[:-1]:                 # left side, down (start point excluded)
        add(p, False)

    out_p, out_f = _dedupe(pts, flg)
    loops = [_chain(out_p, out_f, r)]

    # interior holes, rows 1 .. ny-2
    for j in range(1, ny - 1):
        cy = j * row_pitch
        xs = even_x if j % 2 == 0 else odd_x
        for cx in xs:
            vs = [(cx + rh * math.cos(math.radians(a)), cy + rh * math.sin(math.radians(a)))
                  for a in (90, 150, 210, 270, 330, 30)]
            loops.append(_chain(vs, [True] * 6, r))
    return loops, width, height


def _dedupe(pts, flg):
    out_p, out_f = [], []
    for p, f in zip(pts, flg):
        if out_p and _norm(_sub(p, out_p[-1])) < 1e-9:
            out_f[-1] = out_f[-1] or f
            continue
        out_p.append(p)
        out_f.append(f)
    if _norm(_sub(out_p[0], out_p[-1])) < 1e-9:
        out_p.pop()
        out_f.pop()
    return out_p, out_f


def loop_area(segs):
    """Signed area by Green's theorem (CCW positive)."""
    a = 0.0
    for s in segs:
        if s[0] == 'L':
            (x1, y1), (x2, y2) = s[1], s[2]
            a += 0.5 * (x1 * y2 - x2 * y1)
        else:
            c, p1, p2, d = s[1], s[2], s[3], s[4]
            rho = _norm(_sub(p1, c))
            a1 = math.atan2(p1[1] - c[1], p1[0] - c[0])
            a2 = math.atan2(p2[1] - c[1], p2[0] - c[0])
            if d == 'CCW':
                sweep = (a2 - a1) % (2.0 * math.pi)
            else:
                sweep = -((a1 - a2) % (2.0 * math.pi))
            e = a1 + sweep
            integral = rho * (c[0] * (math.sin(e) - math.sin(a1))
                              - c[1] * (math.cos(e) - math.cos(a1))) + rho * rho * sweep
            a += 0.5 * integral
    return a


def expected_area(loops):
    return abs(loop_area(loops[0])) - sum(abs(loop_area(l)) for l in loops[1:])


# --------------------------------------------------------------- Abaqus part
def build_in_abaqus():
    from abaqus import mdb
    from abaqusConstants import (TWO_D_PLANAR, DEFORMABLE_BODY,
                                 CLOCKWISE, COUNTERCLOCKWISE)

    if MODEL_NAME in mdb.models.keys():
        raise RuntimeError("Model '%s' already exists. Rename it or change MODEL_NAME." % MODEL_NAME)

    loops, width, height = build_loops()
    model = mdb.Model(name=MODEL_NAME)
    sk = model.ConstrainedSketch(name='coreProfile', sheetSize=SHEET_SIZE)

    n_line = n_arc = 0
    for loop in loops:
        for s in loop:
            if s[0] == 'L':
                sk.Line(point1=s[1], point2=s[2])
                n_line += 1
            else:
                sk.ArcByCenterEnds(center=s[1], point1=s[2], point2=s[3],
                                   direction=(COUNTERCLOCKWISE if s[4] == 'CCW' else CLOCKWISE))
                n_arc += 1

    part = model.Part(name=PART_NAME, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    part.BaseShell(sketch=sk)

    xs = [v.pointOn[0][0] for v in part.vertices]
    ys = [v.pointOn[0][1] for v in part.vertices]
    area = part.getArea(part.faces)
    print('--- RHC patch check ---')
    print('holes (interior loops)    : %d' % (len(loops) - 1))
    print('sketch lines / arcs       : %d / %d' % (n_line, n_arc))
    print('faces in part             : %d' % len(part.faces))
    print('bounding box x            : %.4f to %.4f  (target 0 to %.4f)' % (min(xs), max(xs), width))
    print('bounding box y            : %.4f to %.4f  (target 0 to %.4f)' % (min(ys), max(ys), height))
    print('Abaqus face area  [mm^2]  : %.4f' % area)
    print('script area       [mm^2]  : %.4f' % expected_area(loops))
    print('-----------------------')


if __name__ == '__main__':
    build_in_abaqus()
