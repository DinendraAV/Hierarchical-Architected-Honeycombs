# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# build_planestrain_run.py - general plane-strain, quasi-static compression setup.
# Abaqus/CAE 2020 (Python 2.7).  Run:  File > Run Script.
#
# Starts from a model that already holds the 2D planar honeycomb PART (for example
# RHC_Patch / Core from build_rhc_patch.py). Builds in a NEW model (MODEL_NAME), a copy
# of the source model, so the source model is not touched:
#   TPU neo-Hooke material -> plane strain section -> CPE4 mesh (Abaqus/Standard) ->
#   two discrete rigid platens (R2D2 wires tied to reference points) ->

#   contact pairs (honeycomb self-contact + platen contact, friction FRICTION) ->
#   static step, nlgeom, optional automatic stabilization -> displacement-driven top
#   platen at SPEED_MM_MIN -> output requests -> job (.inp written, data check).
# Works for any single-face planar honeycomb part: the platens are sized from the part's
# bounding box. Compression axis X or Y (LOAD_AXIS).
# NOT tested in Abaqus by the author: read the printout, then read the data check.
# ---------------------------------------------------------------------------
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
# source model and part (must already exist in the open model database, or in CAE_PATH)
CAE_PATH   = ''                 # '' = use the database that is open in CAE; else full path of the .cae
SRC_MODEL  = 'RHC_Patch'
PART_NAME  = 'Core'
MODEL_NAME = 'PlaneStrain_QS'   # new model created by this script (replaced if it already exists)

LOAD_AXIS  = 'Y'                # 'X' or 'Y': axis along which the platens press
DEPTH      = 1.0                # mm, plane strain thickness (forces are per this depth)

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

# platens: discrete rigid wires. The wire direction fixes the contact normal (see make_platen).
PLATEN_MARGIN = 5.0             # mm, platen overhang on each side of the patch
PLATEN_MESH   = 2.0             # mm
FLIP_PLATEN   = False           # True: reverse both wires if the platens face away from the patch

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
SAVE_CAE       = True           # save the new model database as MODEL_NAME in the work directory
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
AX = 'XY'
SPEED = SPEED_MM_MIN / 60.0                      # mm/s
TRAVEL = COMPRESS + GAP_TOP + GAP_BOT            # mm of platen travel
T_STEP = TRAVEL / SPEED                          # s


# ---- 1. new model from the source model --------------------------------------------
def b_model():
    global mdb
    if SRC_MODEL not in mdb.models.keys():
        if CAE_PATH and os.path.exists(CAE_PATH):
            mdb = openMdb(pathName=CAE_PATH)
        else:
            raise Exception("model '%s' is not in the open database; open the .cae or set CAE_PATH" % SRC_MODEL)
    if SRC_MODEL not in mdb.models.keys():
        raise Exception("model '%s' not found in %s" % (SRC_MODEL, CAE_PATH))
    if MODEL_NAME == SRC_MODEL:
        raise Exception('MODEL_NAME must differ from SRC_MODEL: the source model is never overwritten')
    if PART_NAME not in mdb.models[SRC_MODEL].parts.keys():
        raise Exception("part '%s' not found in model '%s'" % (PART_NAME, SRC_MODEL))
    if LOAD_AXIS not in ('X', 'Y'):
        raise Exception("LOAD_AXIS must be 'X' or 'Y'")
    if ELEMENT not in ELEMS:
        raise Exception('ELEMENT must be one of %s' % ', '.join(sorted(ELEMS.keys())))
    if MODEL_NAME in mdb.models.keys():
        del mdb.models[MODEL_NAME]
    mdb.Model(name=MODEL_NAME, objectToCopy=mdb.models[SRC_MODEL])
run('new model %s copied from %s' % (MODEL_NAME, SRC_MODEL), b_model)

m = mdb.models[MODEL_NAME]
p = m.parts[PART_NAME]

pts = [v.pointOn[0] for v in p.vertices]
lo = [min(q[i] for q in pts) for i in range(2)]
hi = [max(q[i] for q in pts) for i in range(2)]
lax = AX.index(LOAD_AXIS)
tax = 1 - lax
H_LOAD = hi[lax] - lo[lax]
W_TRANS = hi[tax] - lo[tax]
cen_t = 0.5 * (lo[tax] + hi[tax])
print('      part bounding box (from vertices): X %.4f to %.4f, Y %.4f to %.4f' % (lo[0], hi[0], lo[1], hi[1]))
print('      loaded axis %s, H = %.4f mm, transverse width = %.4f mm, faces in part: %d'
      % (LOAD_AXIS, H_LOAD, W_TRANS, len(p.faces)))
if len(p.faces) != 1:
    print('WARNING: the part has %d faces; the script expects one connected face' % len(p.faces))


# ---- 2. materials and sections --------------------------------------------------------
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
    print('NOTE: nu = %g is nearly incompressible. Full-integration CPE4 can lock volumetrically;'
          % NU)
    print('      if stiffness looks too high, repeat with ELEMENT = CPE4H or CPE4I.')


# ---- 3. honeycomb mesh (part level) --------------------------------------------------
def b_mesh_hc():
    quad, tri = ELEMS[ELEMENT]
    try:
        p.deleteMesh()
    except Exception:
        pass
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
def xy(t, l):
    """(transverse, load-axis) coordinates -> (x, y)."""
    return (t, l) if lax == 1 else (l, t)


PLANE_TOP = hi[lax] + GAP_TOP
PLANE_BOT = lo[lax] - GAP_BOT
T0 = lo[tax] - PLATEN_MARGIN
T1 = hi[tax] + PLATEN_MARGIN


def wire_ends(plane, normal_sign):
    """End points (start, end) of a platen wire at load-axis coordinate `plane`.
    Convention used here: the contact normal is the wire tangent rotated 90 degrees clockwise
    (to the right of the drawing direction), as for a CCW deformable element face.
    normal_sign = +1: normal along +load axis (bottom platen), -1: along -load axis (top platen)."""
    sgn = -normal_sign if FLIP_PLATEN else normal_sign
    n = [0.0, 0.0]
    n[lax] = float(sgn)
    tang = (-n[1], n[0])                          # tangent whose clockwise normal is n
    a_pt, b_pt = xy(T0, plane), xy(T1, plane)
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
    print('      %s: wire from (%.3f, %.3f) to (%.3f, %.3f), %d elements, normal %s%s'
          % (name, p1[0], p1[1], p2[0], p2[1], len(pp.elements),
             '+' if (normal_sign > 0) != FLIP_PLATEN else '-', LOAD_AXIS))


def b_platens():
    make_platen('Platen-Top', PLANE_TOP, -1)
    make_platen('Platen-Bot', PLANE_BOT, +1)
run('platen parts (discrete rigid wires, %g mm overhang)' % PLATEN_MARGIN, b_platens)


# ---- 5. assembly, reference points, rigid bodies, contact surfaces ----------------------
a = m.rootAssembly


def b_assembly():
    a.DatumCsysByDefault(CARTESIAN)
    a.Instance(name='HC-1', part=p, dependent=ON)
    for name, plane, tag in (('Platen-Top', PLANE_TOP, 'TOP'), ('Platen-Bot', PLANE_BOT, 'BOT')):
        a.Instance(name=name, part=m.parts[name], dependent=ON)
        c = xy(cen_t, plane)
        rpf = a.ReferencePoint(point=(c[0], c[1], 0.0))
        rp = a.referencePoints[rpf.id]
        a.Set(referencePoints=(rp,), name='RP_' + tag)
        m.RigidBody(name='RB-' + tag, refPointRegion=regionToolset.Region(referencePoints=(rp,)),
                    bodyRegion=regionToolset.Region(edges=a.instances[name].edges))
    a.Surface(side1Edges=a.instances['HC-1'].edges, name='S_HC')
    a.Surface(side1Edges=a.instances['Platen-Top'].edges, name='S_PTOP')
    a.Surface(side1Edges=a.instances['Platen-Bot'].edges, name='S_PBOT')
    print('      CHECK in CAE: Assembly > Surfaces > S_PTOP / S_PBOT must show the highlighted side')
    print('      facing the honeycomb. If not, set FLIP_PLATEN = True and rerun.')
    print('      top platen contact face at %s = %.4f (patch top %.4f, gap %.3f)'
          % (LOAD_AXIS, PLANE_TOP, hi[lax], GAP_TOP))
    print('      bottom platen contact face at %s = %.4f (patch bottom %.4f, gap %.3f)'
          % (LOAD_AXIS, PLANE_BOT, lo[lax], GAP_BOT))
run('assembly, platens as rigid bodies, contact surfaces', b_assembly)


# ---- 6. static step ----------------------------------------------------------------------
def b_step():
    kw = {}
    if STABILIZE:
        kw = dict(stabilizationMethod=DAMPING_FACTOR, stabilizationMagnitude=STAB_FACTOR)
    m.StaticStep(name='Step-1', previous='Initial', timePeriod=T_STEP, nlgeom=ON,
                 maxNumInc=MAX_NUM_INC, initialInc=T_STEP * INIT_FRAC,
                 minInc=T_STEP * MIN_FRAC, maxInc=T_STEP * MAX_FRAC, **kw)
run('static step %.1f s, nlgeom on, stabilization %s' % (T_STEP, 'on' if STABILIZE else 'off'), b_step)


# ---- 7. boundary conditions ---------------------------------------------------------------
def b_bc():
    # a 2D reference point has U1, U2, UR3: all three fixed in the Initial step
    m.DisplacementBC(name='BC-Bot', createStepName='Initial', region=a.sets['RP_BOT'],
                     u1=SET, u2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM,
                     fieldName='', localCsys=None)
    m.DisplacementBC(name='BC-Top', createStepName='Initial', region=a.sets['RP_TOP'],
                     u1=SET, u2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM,
                     fieldName='', localCsys=None)
    kw = {'u%d' % (lax + 1): -TRAVEL}
    m.boundaryConditions['BC-Top'].setValuesInStep(stepName='Step-1', **kw)
    print('      top platen: U%d = %g mm over %.2f s (constant %.5f mm/s = %g mm/min)'
          % (lax + 1, -TRAVEL, T_STEP, SPEED, SPEED_MM_MIN))
run('boundary conditions (bottom fixed, top driven along %s)' % LOAD_AXIS, b_bc)


# ---- 8. contact ----------------------------------------------------------------------------
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
    m.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'LE', 'U', 'RF'),
                                                  numIntervals=FIELD_INTERVALS)
    m.historyOutputRequests['H-Output-1'].setValues(
        variables=('ALLAE', 'ALLIE', 'ALLKE', 'ALLSD', 'ALLWK', 'ETOTAL'),
        numIntervals=HIST_INTERVALS)
    for tag in ('TOP', 'BOT'):
        m.HistoryOutputRequest(name='H-Output-' + tag, createStepName='Step-1',
            region=a.sets['RP_' + tag], variables=('U1', 'U2', 'RF1', 'RF2'),
            numIntervals=HIST_INTERVALS, sectionPoints=DEFAULT, rebar=EXCLUDE)
run('output requests (energies incl. ALLSD, platen U and RF, field output)', b_output)


# ---- 10. job and data check -----------------------------------------------------------------------
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
print('procedure          : Abaqus/Standard static, nlgeom ON, stabilization %s (factor %g)'
      % ('ON' if STABILIZE else 'OFF', STAB_FACTOR))
print('element            : %s, in-plane size %.3f mm' % (ELEMENT, INPLANE_SIZE))
print('loaded axis        : %s,  H = %.4f mm;  width %.4f mm;  depth %g mm' % (LOAD_AXIS, H_LOAD, W_TRANS, DEPTH))
if mass_t is not None:
    print('honeycomb mass     : %.5e tonne  (%.4f g)' % (mass_t, mass_t * 1.0e6))
print('platen travel      : %.3f mm (compression %.3f + gaps %.3f)' % (TRAVEL, COMPRESS, GAP_TOP + GAP_BOT))
print('platen speed       : %g mm/min = %.5f mm/s;  step time %.2f s' % (SPEED_MM_MIN, SPEED, T_STEP))
print('end strain         : %.1f %%  (after subtracting the gaps)' % (100.0 * COMPRESS / H_LOAD))
print('friction           : %g' % FRICTION)
print('NEXT: open %s.dat and %s.msg; check element counts, warnings, contact status and any' % (jname, jname))
print('      "zero pivot" or "numerical singularity" lines. Then set DATACHECK_ONLY = False.')
