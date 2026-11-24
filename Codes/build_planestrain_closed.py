# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# build_planestrain_closed.py - plane-strain quasi-static compression of the honeycomb patch CLOSED by two
# 1 mm plates, between two rigid platens. Abaqus/CAE 2020 (Python 2.7).  Run:  File > Run Script.
#
# Self-contained: no other model, part or script is needed. The geometry is built from scratch (same
# geometry rules as Build_2D_Geometry.py) in a NEW model (MODEL_NAME, then _2, _3, ...).
#
# Assembly:
#   - ONE deformable part "ClosedCore": the honeycomb with a solid plate of PLATE_T on the top and on
#     the bottom, one face (the plates and the ribs share nodes, the half-hole notches become closed).
#   - TWO rigid platens (discrete rigid wires, each tied to a reference point). The bottom one is fixed,
#     the top one moves at SPEED_MM_MIN. Lateral motion of both platens is fixed.
#   - Contact: as in build_planestrain_run.py: honeycomb self-contact and platen contact, penalty friction.
# Run it in an empty database (File > New Model Database) if you want a .cae with nothing else in it.
# NOT tested in Abaqus by the author: read the printout, then read the data check.
# ---------------------------------------------------------------------------
from __future__ import division
from abaqus import *
from abaqusConstants import *
import regionToolset
import math, os, traceback
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *

# ======================= SETTINGS (edit here only) =========================
MODEL_NAME = 'PlaneStrain_CLOSED'   # base name of the NEW model. An existing model is never replaced:
                                    # if the name is taken, _2, _3, ... is used
PART_NAME  = 'ClosedCore'

# lattice geometry (mm), pointy-top hexagons, load along Y, origin at the bottom-left of the LATTICE
D        = 9.393                 # hole distance across flats
T        = 1.608                  # wall thickness
R_FILLET = 0.5                  # hole corner fillet; the two corners where a notch meets a plate stay sharp
NX       = 5                    # holes in the even rows (odd rows have NX-1)
NY       = 7                    # rows, must be odd
PLATE_T  = 1.0                  # mm, thickness of the closing plate on the top and on the bottom
FILLET_PLATE_CORNERS = False    # True: also round the two corners where each end notch meets the plate
DEPTH    = 1.0                  # mm, plane strain thickness (forces are per this depth)
MAKE_3D_COPY = True             # also create part 'ClosedCore_3D' (the same shape extruded by DEPTH) for the STEP
                                # export. It is not used by the analysis.

# material (placeholder TPU, neo-Hookean, rate independent)
E_MPA = 67.0
NU    = 0.49
RHO   = 1.22e-9                 # tonne/mm^3

# contact: the same setup as build_planestrain_run.py (one property for all pairs: penalty friction, hard normal)
FRICTION = 0.3                  # honeycomb self-contact and platen contact (0.3 as in the completed 0.5 and 0.25 mm runs)
PLATEN_LINK = 'TIE'         # how the platens act on the plates:
                                # 'FINGERS': friction contact, but ONLY over the 6 finger segments of each plate face
                                #            (the skin over the holes carries no force and is left out of the contact)
                                #            GAP_TOP / GAP_BOT may be > 0 (0.05 as in the open-core runs).
                                # 'TIE'    : the plate outer face is tied to the rigid platen (bonded, no slip, NO contact
                                #            search). Needs GAP_TOP = GAP_BOT = 0. This is what the professor describes:
                                #            a prescribed displacement on the plate.
                                # 'CONTACT': surface-to-surface contact with friction FRICTION. With the flat plates this
                                #            chatters in increment 1 (see the .msg of the first closed-core runs).

# loading (quasi-static)
SPEED_MM_MIN = 5.0              # mm/min platen speed  (5 mm/min = 0.08333 mm/s)
COMPRESS = 30.0                 # mm of compression after contact
GAP_TOP  = 0.0                # mm, initial clearance between the top plate and the top platen (as in the open-core
                                # runs). PLATEN_LINK = 'TIE' needs 0.
GAP_BOT  = 0.0                  # mm, 0 = bottom platen starts touching the bottom plate

# platens: discrete rigid wires
PLATEN_MARGIN = 5.0             # mm, platen overhang on each side of the specimen
PLATEN_MESH   = 2.0             # mm
FLIP_PLATEN   = False           # True: reverse both wires (the .dat scan at the end reports orientation warnings)

# mesh and element: 'CPE4' (plane strain, full integration), 'CPE4H', 'CPE4I' or 'CPE4R'
INPLANE_SIZE = 0.25              # mm. 0.5 for a quick test; 0.25 for the final run
ELEMENT      = 'CPE4H'          # hybrid: recommended for hyperelastic materials (Guide 22.5.1)

# step. 'DYNAMIC_QS': implicit dynamic, application = quasi-static. 'STATIC': static (with stabilization).
STEP_TYPE   = 'DYNAMIC_QS'
STABILIZE   = True              # STATIC only: automatic stabilization (damping)
STAB_FACTOR = 2.0e-4            # STATIC only: damping factor magnitude (Abaqus default)
INIT_FRAC   = 1.0e-3            # initial increment as a fraction of the step time
MAX_FRAC    = 1.0e-3            # maximum increment as a fraction of the step time
MIN_FRAC    = 1.0e-9
MAX_NUM_INC = 20000

# output
FIELD_INTERVALS = 100
HIST_INTERVALS  = 400

# job
TAG            = 'tie_t161_ctrl_12'
CPUS           = 4
DATACHECK_ONLY = True           # True: data check only. False: full analysis.
SUBMIT         = False           # False: only write the .inp
SAVE_CAE       = True           # save the database as MODEL_NAME in the work directory
# ===========================================================================


def run(label, fn):
    try:
        out = fn()
        print('OK    : ' + label)
        return out
    except Exception as e:
        print('FAILED: %s  ->  %s' % (label, str(e)))
        traceback.print_exc()
        raise Exception('STOPPED at: ' + label + '. Send the lines above.')


ELEMS = {'CPE4':  (CPE4,  CPE3),
         'CPE4R': (CPE4R, CPE3),
         'CPE4I': (CPE4I, CPE3),
         'CPE4H': (CPE4H, CPE3H)}
lax = 1                                          # the load axis is Y
SPEED = SPEED_MM_MIN / 60.0                      # mm/s
TRAVEL = COMPRESS + GAP_TOP + GAP_BOT            # mm of top platen travel
T_STEP = TRAVEL / SPEED                          # s


# ---- 0. pure geometry (no Abaqus calls) --------------------------------------------------------
def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _norm(a):
    return math.hypot(a[0], a[1])


def _unit(a):
    n = _norm(a)
    return (a[0] / n, a[1] / n)


def _fillet_corner(prev_pt, v, next_pt, r):
    """(tangent point on the previous edge, tangent point on the next edge, arc centre, 'CCW' or 'CW')."""
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
    """Closed chain of points; flags[i] True = fillet vertex i.
    Returns segments ('L', p1, p2) or ('A', centre, p1, p2, 'CCW'|'CW')."""
    n = len(points)
    entry, exit_, arcs = [], [], []
    for i in range(n):
        v = points[i]
        if flags[i] and r > 0.0:
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


def build_loops(d, t, r, nx, ny, plate_t, fillet_plate):
    """loops[0] = outer loop (plates included), the others = holes (end notches closed by the plates).
    The lattice occupies 0 <= y <= height; the plates lie in y < 0 and y > height."""
    if ny % 2 != 1:
        raise ValueError('NY must be odd so that the first and last rows are both even-type rows')
    pitch = d + t
    big_r = pitch / math.sqrt(3.0)          # wall-centreline hexagon circumradius
    rh = d / math.sqrt(3.0)                 # hole circumradius
    hx = d / 2.0                            # hole half-width across flats
    row_pitch = 1.5 * big_r
    o = t / 2.0
    delta = o / math.sqrt(3.0)
    width = nx * pitch + t
    height = (ny - 1) * row_pitch
    even_x = [pitch / 2.0 + o + i * pitch for i in range(nx)]
    odd_x = [pitch + o + i * pitch for i in range(nx - 1)]

    # left edge, going up (zigzag between the plates), x = 0 at the even rows
    left = [(0.0, -plate_t)]
    for j in range(0, ny - 1, 2):
        cy = j * row_pitch
        left += [(0.0, cy + big_r / 2.0 + delta),
                 (pitch / 2.0, cy + big_r + delta),
                 (pitch / 2.0, cy + 2.0 * big_r - delta),
                 (0.0, cy + 2.5 * big_r - delta)]
    left.append((0.0, height + plate_t))
    right = [(width - x, y) for (x, y) in left]

    pts, flg = [], []
    pts.append((0.0, -plate_t)); flg.append(False)
    pts.append((width, -plate_t)); flg.append(False)
    for p in right[1:]:                                  # right side, up to the top plate
        pts.append(p); flg.append(False)
    for p in list(reversed(left))[:-1]:                  # top edge corner, then the left side down
        pts.append(p); flg.append(False)
    out_p, out_f = _dedupe(pts, flg)
    loops = [_chain(out_p, out_f, r)]

    # bottom row: closed half holes above y = 0 (the plate closes them)
    for cx in even_x:
        vs = [(cx - hx, 0.0), (cx - hx, rh / 2.0), (cx, rh), (cx + hx, rh / 2.0), (cx + hx, 0.0)]
        fl = [fillet_plate, True, True, True, fillet_plate]
        loops.append(_chain(vs, fl, r))
    # interior rows: full hexagons
    for j in range(1, ny - 1):
        cy = j * row_pitch
        xs = even_x if j % 2 == 0 else odd_x
        for cx in xs:
            vs = [(cx + rh * math.cos(math.radians(a)), cy + rh * math.sin(math.radians(a)))
                  for a in (90, 150, 210, 270, 330, 30)]
            loops.append(_chain(vs, [True] * 6, r))
    # top row: closed half holes below y = height
    for cx in even_x:
        vs = [(cx + hx, height), (cx + hx, height - rh / 2.0), (cx, height - rh),
              (cx - hx, height - rh / 2.0), (cx - hx, height)]
        fl = [fillet_plate, True, True, True, fillet_plate]
        loops.append(_chain(vs, fl, r))
    return loops, width, height


def loop_area(segs):
    """Signed area by Green's theorem (CCW positive)."""
    a = 0.0
    for s in segs:
        if s[0] == 'L':
            (x1, y1), (x2, y2) = s[1], s[2]
            a += 0.5 * (x1 * y2 - x2 * y1)
        else:
            c, p1, p2, dr = s[1], s[2], s[3], s[4]
            rho = _norm(_sub(p1, c))
            a1 = math.atan2(p1[1] - c[1], p1[0] - c[0])
            a2 = math.atan2(p2[1] - c[1], p2[0] - c[0])
            sweep = ((a2 - a1) % (2.0 * math.pi)) if dr == 'CCW' else -((a1 - a2) % (2.0 * math.pi))
            e = a1 + sweep
            integral = rho * (c[0] * (math.sin(e) - math.sin(a1))
                              - c[1] * (math.cos(e) - math.cos(a1))) + rho * rho * sweep
            a += 0.5 * integral
    return a


def expected_area(loops):
    holes = 0.0
    for l in loops[1:]:
        holes += abs(loop_area(l))
    return abs(loop_area(loops[0])) - holes


if PLATEN_LINK not in ('FINGERS', 'TIE', 'CONTACT'):
    raise Exception("PLATEN_LINK must be 'FINGERS', 'TIE' or 'CONTACT'")
if PLATEN_LINK == 'TIE' and (GAP_TOP != 0.0 or GAP_BOT != 0.0):
    raise Exception("PLATEN_LINK = 'TIE' needs GAP_TOP = 0 and GAP_BOT = 0 (the platens must touch the plates)")
if STEP_TYPE not in ('DYNAMIC_QS', 'STATIC'):
    raise Exception("STEP_TYPE must be 'DYNAMIC_QS' or 'STATIC'")
if ELEMENT not in ELEMS:
    raise Exception('ELEMENT must be one of %s' % ', '.join(sorted(ELEMS.keys())))
if PLATE_T <= 0.0:
    raise Exception('PLATE_T must be positive')
LOOPS, W_TRANS, H_LAT = build_loops(D, T, R_FILLET, NX, NY, PLATE_T, FILLET_PLATE_CORNERS)
Y_LO = -PLATE_T                                  # bottom face of the bottom plate
Y_HI = H_LAT + PLATE_T                           # top face of the top plate
cen_t = 0.5 * W_TRANS
# x ranges of the six finger segments on the plate faces (between the end holes of the first and last rows)
_pitch = D + T
_hole_spans = [(_pitch / 2.0 + T / 2.0 + i * _pitch - D / 2.0, _pitch / 2.0 + T / 2.0 + i * _pitch + D / 2.0)
               for i in range(NX)]
FINGERS = [(0.0, _hole_spans[0][0])]
for _i in range(NX - 1):
    FINGERS.append((_hole_spans[_i][1], _hole_spans[_i + 1][0]))
FINGERS.append((_hole_spans[-1][1], W_TRANS))
SPLIT_X = []
for _a, _b in _hole_spans:
    SPLIT_X += [_a, _b]


# ---- 1. new empty model ---------------------------------------------------------------------------
def b_model():
    global MODEL_NAME
    base = MODEL_NAME
    k = 1
    while MODEL_NAME in mdb.models.keys() or os.path.exists(os.path.join(os.getcwd(), MODEL_NAME + '.cae')):
        k += 1
        MODEL_NAME = '%s_%d' % (base, k)
    if MODEL_NAME != base:
        print("      '%s' already exists; the new model is called '%s'" % (base, MODEL_NAME))
    mdb.Model(name=MODEL_NAME)
    others = [n for n in mdb.models.keys() if n != MODEL_NAME]
    if others:
        print('      NOTE: this database also holds %s. The new model uses none of them.' % ', '.join(others))
run('new empty model %s' % MODEL_NAME, b_model)
m = mdb.models[MODEL_NAME]


# ---- 2. geometry: the closed core, one part ----------------------------------------------------------
def b_geometry():
    global p
    sk = m.ConstrainedSketch(name='coreProfile', sheetSize=200.0)
    n_line = n_arc = 0
    for loop in LOOPS:
        for s in loop:
            if s[0] == 'L':
                sk.Line(point1=s[1], point2=s[2])
                n_line += 1
            else:
                sk.ArcByCenterEnds(center=s[1], point1=s[2], point2=s[3],
                                   direction=(COUNTERCLOCKWISE if s[4] == 'CCW' else CLOCKWISE))
                n_arc += 1
    p = m.Part(name=PART_NAME, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    p.BaseShell(sketch=sk)
    xs = [v.pointOn[0][0] for v in p.vertices]
    ys = [v.pointOn[0][1] for v in p.vertices]
    area = p.getArea(p.faces)
    print('      holes (closed loops)  : %d' % (len(LOOPS) - 1))
    print('      sketch lines / arcs   : %d / %d' % (n_line, n_arc))
    print('      faces in part         : %d' % len(p.faces))
    print('      bounding box x        : %.4f to %.4f  (target 0 to %.4f)' % (min(xs), max(xs), W_TRANS))
    print('      bounding box y        : %.4f to %.4f  (target %.4f to %.4f)' % (min(ys), max(ys), Y_LO, Y_HI))
    print('      Abaqus face area      : %.4f mm^2' % area)
    print('      script area           : %.4f mm^2' % expected_area(LOOPS))
    if len(p.faces) != 1:
        raise Exception('the part has %d faces instead of 1: the sketch loops are not closed or intersect' % len(p.faces))
    if PLATEN_LINK == 'FINGERS':
        # split the two plate outer edges at the finger boundaries so the finger segments can be selected
        for plane in (Y_HI, Y_LO):
            for xs_ in SPLIT_X:
                e = p.edges.findAt(((xs_, plane, 0.0),))
                dp = p.DatumPointByCoordinate(coords=(xs_, plane, 0.0))
                p.PartitionEdgeByPoint(edge=e[0], point=p.datums[dp.id])
        print('      plate outer edges split into %d segments each (%d finger segments)'
              % (2 * NX + 1, NX + 1))
    if abs(area - expected_area(LOOPS)) > 1.0e-3 * expected_area(LOOPS):
        raise Exception('Abaqus area differs from the script area by more than 0.1 %')
    print('      lattice height H = %.4f mm (0 to H); plates %g mm on each side; total height %.4f mm'
          % (H_LAT, PLATE_T, Y_HI - Y_LO))
run('closed core geometry (lattice + %g mm plates, one face)' % PLATE_T, b_geometry)


def b_geometry_3d():
    if not MAKE_3D_COPY:
        return
    try:
        p3 = m.Part(name=PART_NAME + '_3D', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p3.BaseSolidExtrude(sketch=m.sketches['coreProfile'], depth=DEPTH)
        print('      %s_3D: solid extruded %g mm, %d cell(s). Export it with File > Export > Part.'
              % (PART_NAME, DEPTH, len(p3.cells)))
    except Exception as e:
        print('      (3D copy for the STEP export was not created: %s; the analysis is not affected)' % str(e))
b_geometry_3d()


# ---- 3. material and section ---------------------------------------------------------------------------
def b_material():
    C10 = E_MPA / (4.0 * (1.0 + NU))
    D1 = 6.0 * (1.0 - 2.0 * NU) / E_MPA
    print('      neo-Hooke: C10 = %.4f MPa, D1 = %.6f 1/MPa (E = %g, nu = %g)' % (C10, D1, E_MPA, NU))
    mat = m.Material(name='TPU')
    mat.Density(table=((RHO,),))
    mat.Hyperelastic(materialType=ISOTROPIC, table=((C10, D1),), testData=OFF,
                     type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)
    m.HomogeneousSolidSection(name='Section-TPU', material='TPU', thickness=DEPTH)
    p.Set(faces=p.faces, name='HC_FACES')
    p.SectionAssignment(region=p.sets['HC_FACES'], sectionName='Section-TPU', offset=0.0,
                        offsetType=MIDDLE_SURFACE, offsetField='',
                        thicknessAssignment=FROM_SECTION)
run('TPU material, plane strain section (depth %g mm), lattice and plates' % DEPTH, b_material)

if NU > 0.475 and ELEMENT in ('CPE4', 'CPE4R', 'CPE4I'):
    print('NOTE: nu = %g is nearly incompressible. Use ELEMENT = CPE4H for hyperelastic material.' % NU)


# ---- 4. mesh (part level) ----------------------------------------------------------------------------------
def b_mesh():
    quad, tri = ELEMS[ELEMENT]
    p.setElementType(regions=(p.faces,), elemTypes=(
        ElemType(elemCode=quad, elemLibrary=STANDARD),
        ElemType(elemCode=tri, elemLibrary=STANDARD)))
    p.seedPart(size=INPLANE_SIZE, deviationFactor=0.1, minSizeFactor=0.1)
    p.setMeshControls(regions=p.faces, elemShape=QUAD, technique=FREE, algorithm=MEDIAL_AXIS)
    p.generateMesh()
    unmeshed = None
    try:
        unmeshed = p.getUnmeshedRegions()
    except AttributeError:
        pass
    if unmeshed:
        print('      quad-only medial-axis mesh failed, retrying quad-dominated')
        p.setMeshControls(regions=p.faces, elemShape=QUAD_DOMINATED, technique=FREE)
        p.generateMesh()
        unmeshed = p.getUnmeshedRegions()
        if unmeshed:
            raise Exception('part has unmeshed regions: check the geometry or INPLANE_SIZE')
    n_tri = len([e for e in p.elements if len(e.connectivity) == 3])
    print('      closed core: %d elements (%d triangles), %d nodes, size %.3f mm, element %s'
          % (len(p.elements), n_tri, len(p.nodes), INPLANE_SIZE, ELEMENT))
    if n_tri:
        print('      WARNING: %d triangular elements: CPE3 is stiff; consider a finer seed near the fillets' % n_tri)
    try:
        q = p.verifyMeshQuality(criterion=ANALYSIS_CHECKS, regions=(p.faces,))
        print('      analysis checks: %s' % str(q))
    except Exception as e:
        print('      (mesh verify skipped: %s)' % str(e))
run('mesh, %s plane strain (Abaqus/Standard)' % ELEMENT, b_mesh)


# ---- 5. platen parts (discrete rigid wires in absolute coordinates) ---------------------------------
PLANE_TOP = Y_HI + GAP_TOP
PLANE_BOT = Y_LO - GAP_BOT
T0 = -PLATEN_MARGIN
T1 = W_TRANS + PLATEN_MARGIN


def wire_ends(plane, normal_sign):
    """End points of a platen wire at y = plane. The positive normal (SPOS) of an R2D2 element is the
    direction node 1 -> node 2 rotated 90 degrees COUNTERCLOCKWISE (Analysis User's Guide, Rigid elements).
    normal_sign = +1: normal along +Y (bottom platen), -1: along -Y (top platen)."""
    sgn = -normal_sign if FLIP_PLATEN else normal_sign
    n = (0.0, float(sgn))
    tang = (n[1], -n[0])
    a_pt, b_pt = (T0, plane), (T1, plane)
    d = (b_pt[0] - a_pt[0], b_pt[1] - a_pt[1])
    if d[0] * tang[0] + d[1] * tang[1] < 0.0:
        a_pt, b_pt = b_pt, a_pt
    return a_pt, b_pt


def make_platen(name, plane, normal_sign):
    p1, p2 = wire_ends(plane, normal_sign)
    s = m.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    s.Line(point1=p1, point2=p2)
    pp = m.Part(name=name, dimensionality=TWO_D_PLANAR, type=DISCRETE_RIGID_SURFACE)
    pp.BaseWire(sketch=s)
    del m.sketches['__profile__']
    pp.seedPart(size=PLATEN_MESH, deviationFactor=0.1, minSizeFactor=0.1)
    pp.generateMesh()
    print('      %s: wire from (%.3f, %.3f) to (%.3f, %.3f), %d elements'
          % (name, p1[0], p1[1], p2[0], p2[1], len(pp.elements)))


def b_platens():
    make_platen('Platen-Top', PLANE_TOP, -1)
    make_platen('Platen-Bot', PLANE_BOT, +1)
run('platen parts (discrete rigid wires, %g mm overhang)' % PLATEN_MARGIN, b_platens)


# ---- 6. assembly, reference points, rigid bodies, contact surfaces -----------------------------------
a = m.rootAssembly


def b_assembly():
    a.DatumCsysByDefault(CARTESIAN)
    a.Instance(name='HC-1', part=p, dependent=ON)
    for name, plane, tag in (('Platen-Top', PLANE_TOP, 'TOP'), ('Platen-Bot', PLANE_BOT, 'BOT')):
        a.Instance(name=name, part=m.parts[name], dependent=ON)
        rpf = a.ReferencePoint(point=(cen_t, plane, 0.0))
        rp = a.referencePoints[rpf.id]
        a.Set(referencePoints=(rp,), name='RP_' + tag)
        m.RigidBody(name='RB-' + tag, refPointRegion=regionToolset.Region(referencePoints=(rp,)),
                    bodyRegion=regionToolset.Region(edges=a.instances[name].edges))
    a.Surface(side1Edges=a.instances['HC-1'].edges, name='S_HC')
    a.Surface(side1Edges=a.instances['Platen-Top'].edges, name='S_PTOP')
    a.Surface(side1Edges=a.instances['Platen-Bot'].edges, name='S_PBOT')
    if PLATEN_LINK == 'FINGERS':
        for sname, plane in (('S_FING_TOP', Y_HI), ('S_FING_BOT', Y_LO)):
            mids = [(0.5 * (x0 + x1), plane, 0.0) for (x0, x1) in FINGERS]
            ed = a.instances['HC-1'].edges.findAt(*[(mp,) for mp in mids])
            if len(ed) != len(FINGERS):
                raise Exception('found %d finger edges at y = %.4f instead of %d' % (len(ed), plane, len(FINGERS)))
            a.Surface(side1Edges=ed, name=sname)
        print('      contact surfaces: %d finger segments on each plate face (skin over the holes excluded)' % len(FINGERS))
    if PLATEN_LINK == 'TIE':
        for sname, plane in (('S_PLATE_TOP', Y_HI), ('S_PLATE_BOT', Y_LO)):
            ed = a.instances['HC-1'].edges.findAt(((cen_t, plane, 0.0),))
            if len(ed) != 1:
                raise Exception('could not find the plate outer edge at y = %.4f (found %d)' % (plane, len(ed)))
            a.Surface(side1Edges=ed, name=sname)
    print('      CHECK in CAE: Assembly > Surfaces > S_PTOP / S_PBOT must show the highlighted side')
    print('      facing the specimen. If not, set FLIP_PLATEN = True and rerun.')
    print('      top platen contact face at y = %.4f (top plate face %.4f, gap %.3f)' % (PLANE_TOP, Y_HI, GAP_TOP))
    print('      bottom platen contact face at y = %.4f (bottom plate face %.4f, gap %.3f)' % (PLANE_BOT, Y_LO, GAP_BOT))
run('assembly, platens as rigid bodies, contact surfaces', b_assembly)


# ---- 7. step -----------------------------------------------------------------------------------------------
def b_step():
    inc = dict(maxNumInc=MAX_NUM_INC, initialInc=T_STEP * INIT_FRAC,
               minInc=T_STEP * MIN_FRAC, maxInc=T_STEP * MAX_FRAC)
    if STEP_TYPE == 'DYNAMIC_QS':
        # amplitude=RAMP: the platen displacement is ramped over the step (the dynamic default is STEP)
        m.ImplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=T_STEP,
                               application=QUASI_STATIC, nlgeom=ON, amplitude=RAMP, **inc)
    else:
        kw = {}
        if STABILIZE:
            kw = dict(stabilizationMethod=DAMPING_FACTOR, stabilizationMagnitude=STAB_FACTOR)
        m.StaticStep(name='Step-1', previous='Initial', timePeriod=T_STEP, nlgeom=ON, **dict(inc, **kw))
run('%s step %.1f s, nlgeom on' % ('implicit dynamic quasi-static' if STEP_TYPE == 'DYNAMIC_QS' else 'static',
                                   T_STEP), b_step)


# ---- 8. boundary conditions ------------------------------------------------------------------------------
def b_bc():
    # a 2D reference point has U1, U2, UR3: all three fixed in the Initial step (no lateral motion)
    m.DisplacementBC(name='BC-Bot', createStepName='Initial', region=a.sets['RP_BOT'],
                     u1=SET, u2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM,
                     fieldName='', localCsys=None)
    m.DisplacementBC(name='BC-Top', createStepName='Initial', region=a.sets['RP_TOP'],
                     u1=SET, u2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM,
                     fieldName='', localCsys=None)
    m.boundaryConditions['BC-Top'].setValuesInStep(stepName='Step-1', u2=-TRAVEL)
    print('      top platen: U2 = %g mm over %.2f s (constant %.5f mm/s = %g mm/min), U1 and UR3 fixed'
          % (-TRAVEL, T_STEP, SPEED, SPEED_MM_MIN))
run('boundary conditions (bottom platen fixed, top platen driven along Y)', b_bc)


# ---- 9. contact ---------------------------------------------------------------------------------------------
def b_contact():
    prop = m.ContactProperty('IntProp-1')
    prop.TangentialBehavior(formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF,
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0,
        table=((FRICTION,),), shearStressLimit=None, maximumElasticSlip=FRACTION,
        fraction=0.005, elasticSlipStiffness=None)
    prop.NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON,
                        constraintEnforcementMethod=DEFAULT)
    m.SelfContactStd(name='Int-Self', createStepName='Initial', surface=a.surfaces['S_HC'],
                     interactionProperty='IntProp-1', thickness=ON)
    if PLATEN_LINK in ('CONTACT', 'FINGERS'):
        for tag, sname in (('TOP', 'S_PTOP'), ('BOT', 'S_PBOT')):
            slave_name = 'S_HC' if PLATEN_LINK == 'CONTACT' else 'S_FING_' + tag
            m.SurfaceToSurfaceContactStd(name='Int-' + tag, createStepName='Initial',
                master=a.surfaces[sname], slave=a.surfaces[slave_name], sliding=FINITE, thickness=ON,
                interactionProperty='IntProp-1', adjustMethod=NONE, initialClearance=OMIT,
                datumAxis=None, clearanceRegion=None)
    else:
        for tag, sname, pname in (('TOP', 'S_PTOP', 'S_PLATE_TOP'), ('BOT', 'S_PBOT', 'S_PLATE_BOT')):
            m.Tie(name='Tie-' + tag, master=a.surfaces[sname], slave=a.surfaces[pname],
                  adjust=OFF, positionToleranceMethod=COMPUTED, tieRotations=OFF, thickness=ON)
run('contact: honeycomb self-contact friction %g; platens %s' % (
    FRICTION, 'tied to the plates (no contact search)' if PLATEN_LINK == 'TIE'
    else ('friction contact over the finger segments only' if PLATEN_LINK == 'FINGERS'
          else 'surface-to-surface contact, friction %g' % FRICTION)), b_contact)


# ---- 10. output requests ----------------------------------------------------------------------------------
def b_output():
    # S and LE do not exist for the R2D2 platen elements, so they are requested on the core only
    m.fieldOutputRequests['F-Output-1'].setValues(variables=('U', 'RF'),
                                                  numIntervals=FIELD_INTERVALS, timeMarks=OFF)
    m.FieldOutputRequest(name='F-Output-HC', createStepName='Step-1',
                         region=a.instances['HC-1'].sets['HC_FACES'], variables=('S', 'LE'),
                         numIntervals=FIELD_INTERVALS, timeMarks=OFF,
                         sectionPoints=DEFAULT, rebar=EXCLUDE)
    m.historyOutputRequests['H-Output-1'].setValues(
        variables=('ALLAE', 'ALLIE', 'ALLKE', 'ALLSD', 'ALLWK', 'ALLFD', 'ETOTAL'),
        numIntervals=HIST_INTERVALS, timeMarks=OFF)
    for tag in ('TOP', 'BOT'):
        m.HistoryOutputRequest(name='H-Output-' + tag, createStepName='Step-1',
            region=a.sets['RP_' + tag], variables=('U1', 'U2', 'RF1', 'RF2'),
            numIntervals=HIST_INTERVALS, timeMarks=OFF, sectionPoints=DEFAULT, rebar=EXCLUDE)
run('output requests (energies incl. ALLSD and ALLFD, platen U and RF, field output)', b_output)


# ---- 11. job and data check ------------------------------------------------------------------------------
mass_t = None
try:
    mp = a.getMassProperties(regions=(a.instances['HC-1'],))
    mass_t = mp['mass']
except Exception as e:
    print('      (mass properties not available: %s)' % str(e))

jname = 'PS_v%s_C%s_%s' % (('%g' % SPEED_MM_MIN).replace('.', 'p'), ('%g' % COMPRESS).replace('.', 'p'), TAG)
_k = 1
_base = jname
while jname in mdb.jobs.keys() or os.path.exists(os.path.join(os.getcwd(), jname + '.odb')):
    _k += 1
    jname = '%s_%d' % (_base, _k)


def scan_dat(name):
    """Print the contact orientation warnings found in the .dat file, if any."""
    path = os.path.join(os.getcwd(), name + '.dat')
    if not os.path.exists(path):
        print('      (no %s.dat found to scan)' % name)
        return
    f = open(path)
    try:
        lines = f.read().splitlines()
    finally:
        f.close()
    hits = [ln.strip() for ln in lines if 'wrong direction' in ln.lower()]
    if hits:
        print('      ORIENTATION WARNINGS STILL PRESENT in %s.dat:' % name)
        for ln in hits[:5]:
            print('        ' + ln)
        print('      -> set FLIP_PLATEN = %s and rerun.' % (not FLIP_PLATEN))
    else:
        print('      no "wrong direction" facet warnings in %s.dat' % name)


def b_job():
    j = mdb.Job(name=jname, model=MODEL_NAME, description='', type=ANALYSIS,
                memory=90, memoryUnits=PERCENTAGE, numCpus=CPUS, numDomains=CPUS,
                multiprocessingMode=DEFAULT, resultsFormat=ODB,
                echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF)
    j.writeInput(consistencyChecking=OFF)
    print('      wrote %s.inp in %s' % (jname, os.getcwd()))
    if SAVE_CAE:
        mdb.saveAs(pathName=os.path.join(os.getcwd(), MODEL_NAME))
        print('      saved %s.cae' % MODEL_NAME)
    if SUBMIT:
        j.submit(consistencyChecking=OFF, datacheckJob=DATACHECK_ONLY)
        j.waitForCompletion()
        scan_dat(jname)
run('job %s (%s)' % (jname, 'DATA CHECK ONLY' if DATACHECK_ONLY else 'FULL RUN'), b_job)

print('')
print('==================== SUMMARY ====================')
print('job                : %s' % jname)
print('model              : %s' % MODEL_NAME)
if STEP_TYPE == 'DYNAMIC_QS':
    print('procedure          : Abaqus/Standard implicit dynamic, application quasi-static, nlgeom ON, ramp')
    print('                     check ALLKE/ALLIE after the run; keep the ratio small (Guide 6.3.3: typically < 10%)')
else:
    print('procedure          : Abaqus/Standard static, nlgeom ON, stabilization %s (factor %g)'
          % ('ON' if STABILIZE else 'OFF', STAB_FACTOR))
print('element            : %s, in-plane size %.3f mm' % (ELEMENT, INPLANE_SIZE))
print('lattice height H   : %.4f mm;  plates %g mm each;  total height %.4f mm;  width %.4f mm;  depth %g mm'
      % (H_LAT, PLATE_T, Y_HI - Y_LO, W_TRANS, DEPTH))
if mass_t is not None:
    print('mass (core+plates) : %.5e tonne  (%.4f g)' % (mass_t, mass_t * 1.0e6))
print('platen travel      : %.3f mm (compression %.3f + gaps %.3f)' % (TRAVEL, COMPRESS, GAP_TOP + GAP_BOT))
print('platen speed       : %g mm/min = %.5f mm/s;  step time %.2f s' % (SPEED_MM_MIN, SPEED, T_STEP))
print('nominal strain     : %.1f %% of the lattice height H  (%.1f %% of the total height with plates)'
      % (100.0 * COMPRESS / H_LAT, 100.0 * COMPRESS / (Y_HI - Y_LO)))
print('platen link        : %s;  self-contact friction %g (penalty)' % (PLATEN_LINK, FRICTION))
if PLATEN_LINK == 'FINGERS':
    print('CHECK AFTER RUN    : the skin over the holes is outside the contact. In the odb, confirm that no plate node')
    print('                     passes the platen plane (U2 of the top plate face nodes must not exceed the platen U2).')
print('EXTRACTION         : in extract_results.py keep H = %.4f, GAP_TOP = %g, and set MASS_T from the .dat.' % (H_LAT, GAP_TOP))
print('NEXT: open %s.dat and %s.msg; check element counts, warnings, contact status and any' % (jname, jname))
print('      "zero pivot" or "numerical singularity" lines. Then set DATACHECK_ONLY = False.')