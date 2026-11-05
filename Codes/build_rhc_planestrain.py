# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# build_rhc_planestrain.py - ONE script: honeycomb geometry + plane strain quasi-static setup.
# Abaqus/CAE 2020 (Python 2.7).  Run:  File > Run Script.
#
# Creates a brand new model (MODEL_NAME) and builds everything in it, nothing is read from or
# copied from any other model:
#   2D planar regular-hexagon honeycomb patch (pointy-top, load along y, every wall T thick,
#   half-cell notches cut through the first and last row, filleted hole corners) ->
#   TPU neo-Hooke -> plane strain section -> CPE4 mesh (Abaqus/Standard) ->
#   two discrete rigid platens (R2D2 wires tied to reference points) ->
#   contact pairs (honeycomb self-contact + platen contact, friction FRICTION) ->
#   static step, nlgeom, automatic stabilization -> top platen driven at SPEED_MM_MIN ->
#   output requests -> job (.inp written, data check).
# NOT tested in Abaqus by the author: read the printout, then read the data check.
# ---------------------------------------------------------------------------
from __future__ import division, print_function
import math
from abaqus import *
from abaqusConstants import *
import regionToolset
import os, traceback
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
MODEL_NAME = 'RHC_PlaneStrain'  # new model; an existing model with this name is replaced
PART_NAME  = 'Core'

# geometry (pointy-top hexagons, load along y, origin at bottom-left)
D = 10.0                        # mm, clear (hole) distance across flats
T = 1.0                         # mm, wall thickness, every wall including the boundary walls
R_FILLET = 0.5                  # mm, hole corner fillet radius
NX = 5                          # holes in the even rows (odd rows have NX-1)
NY = 7                          # number of rows, must be odd
DEPTH = 1.0                     # mm, plane strain thickness (forces are per this depth)

# material (placeholder TPU, neo-Hookean, rate independent)
E_MPA = 67.0
NU    = 0.49
RHO   = 1.22e-9                 # tonne/mm^3
FRICTION = 0.2

# loading (quasi-static)
SPEED_MM_MIN = 5.0              # mm/min platen speed  (5 mm/min = 0.08333 mm/s)
COMPRESS = 20.0                 # mm of patch compression after contact
GAP_TOP  = 0.05                 # mm, initial clearance below the top platen
GAP_BOT  = 0.0                  # mm, 0 = bottom platen starts touching (keeps the patch supported)

# platens: discrete rigid wires. The wire direction fixes the contact normal (see wire_ends).
PLATEN_MARGIN = 5.0             # mm, platen overhang on each side of the patch
PLATEN_MESH   = 2.0             # mm
FLIP_PLATEN   = True            # reverses both wires. True because the data check reported the platen
                                # normals facing away from the honeycomb with the other setting.
                                # If a facet-orientation warning appears again, try False.

# mesh and element: 'CPE4' (plane strain, full integration), 'CPE4H', 'CPE4I' or 'CPE4R'
INPLANE_SIZE = 0.25             # mm
ELEMENT      = 'CPE4'

# static step
STABILIZE   = True              # automatic stabilization (damping) for contact and collapse
STAB_FACTOR = 2.0e-4            # damping factor magnitude (Abaqus default)
INIT_FRAC   = 1.0e-3            # initial increment as a fraction of the step time
MAX_FRAC    = 5.0e-3            # maximum increment as a fraction of the step time
MIN_FRAC    = 1.0e-9
MAX_NUM_INC = 20000

# output
FIELD_INTERVALS = 100
HIST_INTERVALS  = 400

# job
TAG            = 'qs_m025'
CPUS           = 4
DATACHECK_ONLY = True           # True: data check only. False: full analysis.
SUBMIT         = True           # False: only write the .inp
SAVE_CAE       = True           # save the model database as MODEL_NAME in the work directory
SHEET_SIZE     = 200.0
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
SPEED = SPEED_MM_MIN / 60.0                      # mm/s
TRAVEL = COMPRESS + GAP_TOP + GAP_BOT            # mm of platen travel
T_STEP = TRAVEL / SPEED                          # s


# ------------------------------------------------------------ pure geometry
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
    o = t / 2.0                         # outward offset of the boundary from the wall centreline
    delta = o / math.sqrt(3.0)          # mitre shift along y of the offset zigzag (o * tan 30)
    width = nx * pitch + t
    height = (ny - 1) * row_pitch
    even_x = [pitch / 2.0 + o + i * pitch for i in range(nx)]
    odd_x = [pitch + o + i * pitch for i in range(nx - 1)]

    # left edge, going up (wall centreline zigzag offset outward by o, then shifted by +o in x)
    left = [(0.0, 0.0)]
    for j in range(0, ny - 1, 2):
        cy = j * row_pitch
        left += [(0.0, cy + big_r / 2.0 + delta),
                 (pitch / 2.0, cy + big_r + delta),
                 (pitch / 2.0, cy + 2.0 * big_r - delta),
                 (0.0, cy + 2.5 * big_r - delta)]
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




# ======================= ABAQUS MODEL BUILD =======================================

# ---- 1. new model and honeycomb geometry --------------------------------------------
LOOPS, W_TRANS, H_LOAD = build_loops()
LO = (0.0, 0.0)
HI = (W_TRANS, H_LOAD)


def b_model():
    if ELEMENT not in ELEMS:
        raise Exception('ELEMENT must be one of %s' % ', '.join(sorted(ELEMS.keys())))
    if MODEL_NAME in mdb.models.keys():
        del mdb.models[MODEL_NAME]
    mdb.Model(name=MODEL_NAME)
    mm = mdb.models[MODEL_NAME]
    sk = mm.ConstrainedSketch(name='coreProfile', sheetSize=SHEET_SIZE)
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
    pt = mm.Part(name=PART_NAME, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    pt.BaseShell(sketch=sk)
    print('      sketch: %d interior holes, %d lines, %d arcs' % (len(LOOPS) - 1, n_line, n_arc))
run('new model %s and honeycomb geometry' % MODEL_NAME, b_model)

m = mdb.models[MODEL_NAME]
p = m.parts[PART_NAME]


def b_geometry_check():
    xs = [v.pointOn[0][0] for v in p.vertices]
    ys = [v.pointOn[0][1] for v in p.vertices]
    area = p.getArea(p.faces)
    print('      faces in part      : %d' % len(p.faces))
    print('      bounding box x     : %.4f to %.4f  (target 0 to %.4f)' % (min(xs), max(xs), W_TRANS))
    print('      bounding box y     : %.4f to %.4f  (target 0 to %.4f)' % (min(ys), max(ys), H_LOAD))
    print('      Abaqus face area   : %.4f mm^2' % area)
    print('      script area        : %.4f mm^2' % expected_area(LOOPS))
    if len(p.faces) != 1:
        raise Exception('the part has %d faces, expected one connected face' % len(p.faces))
    if abs(area - expected_area(LOOPS)) > 1.0e-3 * expected_area(LOOPS):
        raise Exception('Abaqus face area differs from the script area by more than 0.1 percent')
run('geometry check (faces, bounding box, area)', b_geometry_check)


# ---- 2. material and section --------------------------------------------------------------
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
run('TPU material, plane strain section (depth %g mm)' % DEPTH, b_material)

if NU > 0.475 and ELEMENT in ('CPE4', 'CPE4R', 'CPE4I'):
    print('NOTE: nu = %g is nearly incompressible (K/mu about %.0f). Full-integration CPE4 can lock'
          % (NU, (E_MPA / (3.0 * (1.0 - 2.0 * NU))) / (E_MPA / (2.0 * (1.0 + NU)))))
    print('      volumetrically and Abaqus will warn about it. Use ELEMENT = CPE4H for production runs.')


# ---- 3. honeycomb mesh (part level) -----------------------------------------------------
def b_mesh_hc():
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
            raise Exception('part has unmeshed regions: check INPLANE_SIZE')
    n_tri = len([e for e in p.elements if len(e.connectivity) == 3])
    print('      honeycomb: %d elements (%d triangles), %d nodes, size %.3f mm, element %s'
          % (len(p.elements), n_tri, len(p.nodes), INPLANE_SIZE, ELEMENT))
    if n_tri:
        print('      WARNING: %d triangular elements: CPE3 is stiff; consider a finer seed near the fillets' % n_tri)
    try:
        q = p.verifyMeshQuality(criterion=ANALYSIS_CHECKS, regions=(p.faces,))
        print('      analysis checks: %s' % str(q))
    except Exception as e:
        print('      (mesh verify skipped: %s)' % str(e))
run('honeycomb mesh, %s plane strain (Abaqus/Standard)' % ELEMENT, b_mesh_hc)


# ---- 4. platen parts (discrete rigid wires in absolute coordinates, no instance transform) ----
PLANE_TOP = HI[1] + GAP_TOP
PLANE_BOT = LO[1] - GAP_BOT
X0 = LO[0] - PLATEN_MARGIN
X1 = HI[0] + PLATEN_MARGIN


def wire_ends(plane, normal_sign):
    """End points (start, end) of a horizontal platen wire at y = plane.
    Convention used here: the contact normal is the wire tangent rotated 90 degrees clockwise
    (to the right of the drawing direction), as for a CCW deformable element face.
    normal_sign = +1: normal along +y (bottom platen), -1: along -y (top platen).
    FLIP_PLATEN reverses the result."""
    sgn = -normal_sign if FLIP_PLATEN else normal_sign
    tx = -float(sgn)                              # tangent (tx, 0) has clockwise normal (0, -tx) = (0, sgn)
    if tx > 0.0:
        return (X0, plane), (X1, plane)
    return (X1, plane), (X0, plane)


def make_platen(name, plane, normal_sign):
    p1, p2 = wire_ends(plane, normal_sign)
    s = m.ConstrainedSketch(name='__profile__', sheetSize=SHEET_SIZE)
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


# ---- 5. assembly, reference points, rigid bodies, contact surfaces -------------------------
a = m.rootAssembly


def b_assembly():
    a.DatumCsysByDefault(CARTESIAN)
    a.Instance(name='HC-1', part=p, dependent=ON)
    for name, plane, tag in (('Platen-Top', PLANE_TOP, 'TOP'), ('Platen-Bot', PLANE_BOT, 'BOT')):
        a.Instance(name=name, part=m.parts[name], dependent=ON)
        rpf = a.ReferencePoint(point=(0.5 * (LO[0] + HI[0]), plane, 0.0))
        rp = a.referencePoints[rpf.id]
        a.Set(referencePoints=(rp,), name='RP_' + tag)
        m.RigidBody(name='RB-' + tag, refPointRegion=regionToolset.Region(referencePoints=(rp,)),
                    bodyRegion=regionToolset.Region(edges=a.instances[name].edges))
    a.Surface(side1Edges=a.instances['HC-1'].edges, name='S_HC')
    a.Surface(side1Edges=a.instances['Platen-Top'].edges, name='S_PTOP')
    a.Surface(side1Edges=a.instances['Platen-Bot'].edges, name='S_PBOT')
    a.Set(faces=a.instances['HC-1'].faces, name='HC_ALL')
    print('      CHECK in CAE: Assembly > Surfaces > S_PTOP / S_PBOT must show the highlighted side')
    print('      facing the honeycomb. If not, change FLIP_PLATEN and rerun.')
    print('      top platen contact line at y = %.4f (patch top %.4f, gap %.3f)' % (PLANE_TOP, HI[1], GAP_TOP))
    print('      bottom platen contact line at y = %.4f (patch bottom %.4f, gap %.3f)' % (PLANE_BOT, LO[1], GAP_BOT))
run('assembly, platens as rigid bodies, contact surfaces', b_assembly)


# ---- 6. static step ---------------------------------------------------------------------------
def b_step():
    kw = {}
    if STABILIZE:
        kw = dict(stabilizationMethod=DAMPING_FACTOR, stabilizationMagnitude=STAB_FACTOR)
    m.StaticStep(name='Step-1', previous='Initial', timePeriod=T_STEP, nlgeom=ON,
                 maxNumInc=MAX_NUM_INC, initialInc=T_STEP * INIT_FRAC,
                 minInc=T_STEP * MIN_FRAC, maxInc=T_STEP * MAX_FRAC, **kw)
run('static step %.1f s, nlgeom on, stabilization %s' % (T_STEP, 'on' if STABILIZE else 'off'), b_step)


# ---- 7. boundary conditions --------------------------------------------------------------------
def b_bc():
    # a 2D reference point has U1, U2, UR3: all three fixed in the Initial step
    m.DisplacementBC(name='BC-Bot', createStepName='Initial', region=a.sets['RP_BOT'],
                     u1=SET, u2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM,
                     fieldName='', localCsys=None)
    m.DisplacementBC(name='BC-Top', createStepName='Initial', region=a.sets['RP_TOP'],
                     u1=SET, u2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM,
                     fieldName='', localCsys=None)
    m.boundaryConditions['BC-Top'].setValuesInStep(stepName='Step-1', u2=-TRAVEL)
    print('      top platen: U2 = %g mm over %.2f s (constant %.5f mm/s = %g mm/min)'
          % (-TRAVEL, T_STEP, SPEED, SPEED_MM_MIN))
run('boundary conditions (bottom fixed, top driven along y)', b_bc)


# ---- 8. contact ----------------------------------------------------------------------------------
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
    for tag, sname in (('TOP', 'S_PTOP'), ('BOT', 'S_PBOT')):
        m.SurfaceToSurfaceContactStd(name='Int-' + tag, createStepName='Initial',
            master=a.surfaces[sname], slave=a.surfaces['S_HC'], sliding=FINITE, thickness=ON,
            interactionProperty='IntProp-1', adjustMethod=NONE, initialClearance=OMIT,
            datumAxis=None, clearanceRegion=None)
run('contact: honeycomb self-contact and two platen pairs, friction %g, hard normal' % FRICTION, b_contact)


# ---- 9. output requests ------------------------------------------------------------------------
def b_output():
    # U and RF on the whole model (platens move in the animation); S and LE only on the
    # honeycomb, because R2D2 rigid elements have no stress or strain output.
    # timeMarks=OFF: do not force increments to land on exact output times.
    m.fieldOutputRequests['F-Output-1'].setValues(variables=('U', 'RF'),
        numIntervals=FIELD_INTERVALS, timeMarks=OFF)
    m.FieldOutputRequest(name='F-Output-HC', createStepName='Step-1', region=a.sets['HC_ALL'],
        variables=('S', 'LE'), numIntervals=FIELD_INTERVALS, timeMarks=OFF)
    m.historyOutputRequests['H-Output-1'].setValues(
        variables=('ALLAE', 'ALLIE', 'ALLKE', 'ALLSD', 'ALLWK', 'ETOTAL'),
        numIntervals=HIST_INTERVALS, timeMarks=OFF)
    for tag in ('TOP', 'BOT'):
        m.HistoryOutputRequest(name='H-Output-' + tag, createStepName='Step-1',
            region=a.sets['RP_' + tag], variables=('U1', 'U2', 'RF1', 'RF2'),
            numIntervals=HIST_INTERVALS, timeMarks=OFF, sectionPoints=DEFAULT, rebar=EXCLUDE)
run('output requests (energies incl. ALLSD, platen U and RF, honeycomb S and LE)', b_output)


# ---- 10. job and data check ----------------------------------------------------------------------
mass_t = None
try:
    mp = a.getMassProperties(regions=(a.instances['HC-1'],))
    mass_t = mp['mass']
except Exception as e:
    print('      (mass properties not available: %s)' % str(e))

jname = 'PS_v%s_C%s_%s' % (('%g' % SPEED_MM_MIN).replace('.', 'p'), ('%g' % COMPRESS).replace('.', 'p'), TAG)


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
run('job %s (%s)' % (jname, 'DATA CHECK ONLY' if DATACHECK_ONLY else 'FULL RUN'), b_job)

print('')
print('==================== SUMMARY ====================')
print('job                : %s' % jname)
print('geometry           : %g x %g mm patch, %d holes, cell %g mm, wall %g mm, fillet %g mm'
      % (W_TRANS, H_LOAD, len(LOOPS) - 1 + 2 * NX, D, T, R_FILLET))
print('procedure          : Abaqus/Standard static, nlgeom ON, stabilization %s (factor %g)'
      % ('ON' if STABILIZE else 'OFF', STAB_FACTOR))
print('element            : %s, in-plane size %.3f mm, depth %g mm' % (ELEMENT, INPLANE_SIZE, DEPTH))
if mass_t is not None:
    print('honeycomb mass     : %.5e tonne  (%.4f g)' % (mass_t, mass_t * 1.0e6))
print('platen travel      : %.3f mm (compression %.3f + gaps %.3f)' % (TRAVEL, COMPRESS, GAP_TOP + GAP_BOT))
print('platen speed       : %g mm/min = %.5f mm/s;  step time %.2f s' % (SPEED_MM_MIN, SPEED, T_STEP))
print('end strain         : %.1f %%  (after subtracting the gaps)' % (100.0 * COMPRESS / H_LOAD))
print('friction           : %g' % FRICTION)
print('NEXT: open %s.dat and %s.msg; check element counts, warnings, contact status and any' % (jname, jname))
print('      "zero pivot" or "numerical singularity" lines. Then set DATACHECK_ONLY = False.')