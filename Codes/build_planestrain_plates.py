# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# build_planestrain_plates.py - plane-strain quasi-static compression of the honeycomb patch with
# two END PLATES merged into the lattice (Abaqus/Standard, implicit dynamic quasi-static).
# Abaqus/CAE 2020 (Python 2.7).  Run:  File > Run Script.
#
# Difference from build_planestrain_run.py: there are NO rigid platens and NO platen contact.
#   - A solid plate (PLATE_T mm thick, same TPU material) is added on top and on the bottom of the
#     core and MERGED with it (one part, shared nodes), so the ends of the ribs are bonded to the
#     plates and the half-hole notches at the two ends become closed cells.
#   - The outer face of each plate is coupled (kinematic, translation only) to a reference point.
#     The bottom reference point is fixed. The top one is given the prescribed displacement, with
#     the lateral motion fixed. The plates therefore cannot slide sideways: this is a bonded,
#     no-slip end condition, equivalent to infinite friction.
#   - Only the wall-to-wall self-contact of the lattice remains.
#
# A NEW model is created every time (MODEL_NAME, then _2, _3, ...). The source model is not touched.
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
CAE_PATH   = '"D:\Temporary\PlaneStrain_QS.cae"'                 # '' = use the database that is open in CAE; else full path of the .cae
SRC_MODEL  = 'PlaneStrain_QS'  # source model (must exist in the open database or in CAE_PATH)
PART_NAME  = 'Core'
MODEL_NAME = 'PlaneStrain_PLT'  # base name of the NEW model. An existing model is never replaced:
                                # if the name is taken, _2, _3, ... is used

LOAD_AXIS  = 'Y'                # 'X' or 'Y': axis along which the plates are pressed
DEPTH      = 1.0                # mm, plane strain thickness (forces are per this depth)

# end plates, merged with the core
PLATE_T    = 2.0                # mm, thickness of each plate. PLACEHOLDER: set the real value of the specimen.
                                # The plates have the same width as the core and the same material.

# material (placeholder TPU, neo-Hookean, rate independent)
E_MPA = 67.0
NU    = 0.49
RHO   = 1.22e-9                 # tonne/mm^3
FRICTION = 0.3                  # wall-to-wall self-contact only (there is no platen contact any more)

# loading (quasi-static)
SPEED_MM_MIN = 5.0              # mm/min speed of the top plate (5 mm/min = 0.08333 mm/s)
COMPRESS = 20.0                 # mm of plate travel (no gap: the plates are bonded to the core)

# mesh and element: 'CPE4' (plane strain, full integration), 'CPE4H', 'CPE4I' or 'CPE4R'
INPLANE_SIZE = 0.5              # mm. 0.5 for a quick test; use 0.25 for the final run (see the mesh check)
ELEMENT      = 'CPE4H'          # hybrid: recommended for hyperelastic materials (Guide 22.5.1)

# step. 'DYNAMIC_QS': implicit dynamic, application = quasi-static. 'STATIC': static (with stabilization).
STEP_TYPE   = 'DYNAMIC_QS'
STABILIZE   = True              # STATIC only: automatic stabilization (damping)
STAB_FACTOR = 2.0e-4            # STATIC only: damping factor magnitude (Abaqus default)
INIT_FRAC   = 1.0e-3            # initial increment as a fraction of the step time
MAX_FRAC    = 5.0e-3            # maximum increment as a fraction of the step time
MIN_FRAC    = 1.0e-9
MAX_NUM_INC = 20000

# output
FIELD_INTERVALS = 100
HIST_INTERVALS  = 400

# job
TAG            = 'plt_m05'
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
MERGED = 'CorePlated'
SPEED = SPEED_MM_MIN / 60.0                      # mm/s
TRAVEL = COMPRESS                                # mm of top plate travel (no gap)
T_STEP = TRAVEL / SPEED                          # s


# ---- 1. new model from the source model --------------------------------------------
def b_model():
    global mdb, MODEL_NAME
    if STEP_TYPE not in ('DYNAMIC_QS', 'STATIC'):
        raise Exception("STEP_TYPE must be 'DYNAMIC_QS' or 'STATIC'")
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
    if PLATE_T <= 0.0:
        raise Exception('PLATE_T must be positive')
    base = MODEL_NAME
    k = 1
    while MODEL_NAME in mdb.models.keys() or os.path.exists(os.path.join(os.getcwd(), MODEL_NAME + '.cae')):
        k += 1
        MODEL_NAME = '%s_%d' % (base, k)
    if MODEL_NAME != base:
        print("      '%s' already exists; the new model is called '%s'" % (base, MODEL_NAME))
    mdb.Model(name=MODEL_NAME, objectToCopy=mdb.models[SRC_MODEL])
run('new model %s copied from %s' % (MODEL_NAME, SRC_MODEL), b_model)

m = mdb.models[MODEL_NAME]
p0 = m.parts[PART_NAME]

pts = [v.pointOn[0] for v in p0.vertices]
lo = [min(q[i] for q in pts) for i in range(2)]
hi = [max(q[i] for q in pts) for i in range(2)]
lax = AX.index(LOAD_AXIS)
tax = 1 - lax
H_LOAD = hi[lax] - lo[lax]
W_TRANS = hi[tax] - lo[tax]
cen_t = 0.5 * (lo[tax] + hi[tax])
print('      core bounding box (from vertices): X %.4f to %.4f, Y %.4f to %.4f' % (lo[0], hi[0], lo[1], hi[1]))
print('      loaded axis %s, H = %.4f mm, transverse width = %.4f mm, faces in part: %d'
      % (LOAD_AXIS, H_LOAD, W_TRANS, len(p0.faces)))
if len(p0.faces) != 1:
    print('WARNING: the part has %d faces; the script expects one connected face' % len(p0.faces))


def xy(t, l):
    """(transverse, load-axis) coordinates -> (x, y)."""
    return (t, l) if lax == 1 else (l, t)


PLANE_TOP = hi[lax] + PLATE_T                    # outer face of the top plate (load-axis coordinate)
PLANE_BOT = lo[lax] - PLATE_T                    # outer face of the bottom plate


# ---- 2. end plates merged into the core (one part, shared nodes) ---------------------------
def make_plate(name, l0, l1):
    s = m.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    s.rectangle(point1=xy(lo[tax], l0), point2=xy(hi[tax], l1))
    pp = m.Part(name=name, dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
    pp.BaseShell(sketch=s)
    del m.sketches['__profile__']
    return pp


def purge():
    """Strip the copied model down to the bare Core part (the source model is a finished model)."""
    a0 = m.rootAssembly
    gone = []

    def drop(label, repo, keys=None):
        if repo is None:
            return
        for k in (keys if keys is not None else list(repo.keys())):
            try:
                del repo[k]
                gone.append('%s %s' % (label, k))
            except Exception:
                pass
    drop('interaction', m.interactions)
    drop('interaction property', getattr(m, 'interactionProperties', None))
    drop('constraint', m.constraints)
    drop('BC', m.boundaryConditions)
    drop('predefined field', getattr(m, 'predefinedFields', None))
    drop('load', m.loads)
    # steps: last first, never 'Initial' (its default output requests go with the steps)
    drop('step', m.steps, [k for k in reversed(m.steps.keys()) if k != 'Initial'])
    drop('assembly set', a0.sets)
    drop('assembly surface', a0.surfaces)
    for k in list(a0.instances.keys()):
        try:
            a0.deleteFeatures((k,))
            gone.append('instance ' + k)
        except Exception:
            pass
    for k in list(a0.features.keys()):
        try:
            a0.deleteFeatures((k,))
            gone.append('feature ' + k)
        except Exception:
            pass
    drop('part', m.parts, [k for k in m.parts.keys() if k != PART_NAME])
    drop('section', m.sections)
    drop('material', m.materials)
    drop('amplitude', getattr(m, 'amplitudes', None))
    print('      purged from the copy: %d items (%s)' % (len(gone), ', '.join(gone[:12]) + (' ...' if len(gone) > 12 else '')))
    left = len(a0.instances.keys())
    if left:
        raise Exception('%d instance(s) are still in the assembly: delete them by hand in CAE and rerun' % left)


def b_merge():
    global p
    purge()
    try:
        p0.deleteMesh()
    except Exception:
        pass
    make_plate('PlateTop', hi[lax], PLANE_TOP)
    make_plate('PlateBot', PLANE_BOT, lo[lax])
    a0 = m.rootAssembly
    a0.DatumCsysByDefault(CARTESIAN)
    i1 = a0.Instance(name='MRG-Core', part=p0, dependent=OFF)
    i2 = a0.Instance(name='MRG-Top', part=m.parts['PlateTop'], dependent=OFF)
    i3 = a0.Instance(name='MRG-Bot', part=m.parts['PlateBot'], dependent=OFF)
    a0.InstanceFromBooleanMerge(name=MERGED, instances=(i1, i2, i3), keepIntersections=False,
                                originalInstances=DELETE, domain=GEOMETRY)
    # the temporary merged instance is replaced later by a dependent instance of the merged part
    tmp = [k for k in a0.instances.keys() if a0.instances[k].partName == MERGED]
    if tmp:
        a0.deleteFeatures(tuple(tmp))
    p = m.parts[MERGED]
    nf = len(p.faces)
    print('      merged part %s: %d face(s)' % (MERGED, nf))
    if nf != 1:
        raise Exception('the plates did not merge into one face (%d faces): the ribs would not be bonded. '
                        'Check the plate and core edges in CAE (Part > %s)' % (nf, MERGED))
    print('      plates: %.2f mm thick, %.2f mm wide; plate outer faces at %s = %.4f (bottom) and %.4f (top)'
          % (PLATE_T, W_TRANS, LOAD_AXIS, PLANE_BOT, PLANE_TOP))
run('end plates (%g mm) merged into the core' % PLATE_T, b_merge)


# ---- 3. materials and sections (on the merged part) -----------------------------------------
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
run('TPU material, plane strain section (depth %g mm), core and plates' % DEPTH, b_material)

if NU > 0.475 and ELEMENT in ('CPE4', 'CPE4R', 'CPE4I'):
    print('NOTE: nu = %g is nearly incompressible. The Guide (22.5.1) does not recommend pure' % NU)
    print('      displacement elements for hyperelastic materials; use ELEMENT = CPE4H.')


# ---- 4. mesh (part level) ------------------------------------------------------------------------
def b_mesh():
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
    print('      core + plates: %d elements (%d triangles), %d nodes, size %.3f mm, element %s'
          % (len(p.elements), n_tri, len(p.nodes), INPLANE_SIZE, ELEMENT))
    if n_tri:
        print('      WARNING: %d triangular elements: CPE3 is stiff; consider a finer seed near the fillets' % n_tri)
    try:
        q = p.verifyMeshQuality(criterion=ANALYSIS_CHECKS, regions=(p.faces,))
        print('      analysis checks: %s' % str(q))
    except Exception as e:
        print('      (mesh verify skipped: %s)' % str(e))
run('mesh, %s plane strain (Abaqus/Standard)' % ELEMENT, b_mesh)


# ---- 5. assembly, reference points, couplings, contact surface ----------------------------------
a = m.rootAssembly


def outer_edge(plane):
    """The straight outer edge of a plate, found by a point on it."""
    pt = xy(cen_t, plane)
    edges = a.instances['HC-1'].edges.findAt(((pt[0], pt[1], 0.0),))
    if len(edges) != 1:
        raise Exception('could not find the plate outer edge at %s = %.4f (found %d)'
                        % (LOAD_AXIS, plane, len(edges)))
    return edges


def b_assembly():
    a.Instance(name='HC-1', part=p, dependent=ON)
    for plane, tag in ((PLANE_TOP, 'TOP'), (PLANE_BOT, 'BOT')):
        c = xy(cen_t, plane)
        rpf = a.ReferencePoint(point=(c[0], c[1], 0.0))
        rp = a.referencePoints[rpf.id]
        a.Set(referencePoints=(rp,), name='RP_' + tag)
        edges = outer_edge(plane)
        # translation-only kinematic coupling: the whole outer face moves with the reference point
        m.Coupling(name='C-' + tag, controlPoint=regionToolset.Region(referencePoints=(rp,)),
                   surface=regionToolset.Region(edges=edges), influenceRadius=WHOLE_SURFACE,
                   couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, ur3=OFF)
    a.Surface(side1Edges=a.instances['HC-1'].edges, name='S_HC')
    print('      reference points: RP_TOP at %s = %.4f, RP_BOT at %s = %.4f'
          % (LOAD_AXIS, PLANE_TOP, LOAD_AXIS, PLANE_BOT))
run('assembly, reference points, coupling of the plate outer faces', b_assembly)


# ---- 6. step -----------------------------------------------------------------------------
def b_step():
    inc = dict(maxNumInc=MAX_NUM_INC, initialInc=T_STEP * INIT_FRAC,
               minInc=T_STEP * MIN_FRAC, maxInc=T_STEP * MAX_FRAC)
    if STEP_TYPE == 'DYNAMIC_QS':
        # amplitude=RAMP: the displacement is ramped over the step (the dynamic default is STEP).
        m.ImplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=T_STEP,
                               application=QUASI_STATIC, nlgeom=ON, amplitude=RAMP, **inc)
    else:
        kw = {}
        if STABILIZE:
            kw = dict(stabilizationMethod=DAMPING_FACTOR, stabilizationMagnitude=STAB_FACTOR)
        m.StaticStep(name='Step-1', previous='Initial', timePeriod=T_STEP, nlgeom=ON, **dict(inc, **kw))
run('%s step %.1f s, nlgeom on' % ('implicit dynamic quasi-static' if STEP_TYPE == 'DYNAMIC_QS' else 'static',
                                   T_STEP), b_step)


# ---- 7. boundary conditions ---------------------------------------------------------------
def b_bc():
    # a 2D reference point has U1, U2, UR3: all three fixed in the Initial step (no lateral motion)
    m.DisplacementBC(name='BC-Bot', createStepName='Initial', region=a.sets['RP_BOT'],
                     u1=SET, u2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM,
                     fieldName='', localCsys=None)
    m.DisplacementBC(name='BC-Top', createStepName='Initial', region=a.sets['RP_TOP'],
                     u1=SET, u2=SET, ur3=SET, amplitude=UNSET, distributionType=UNIFORM,
                     fieldName='', localCsys=None)
    kw = {'u%d' % (lax + 1): -TRAVEL}
    m.boundaryConditions['BC-Top'].setValuesInStep(stepName='Step-1', **kw)
    print('      top plate: U%d = %g mm over %.2f s (constant %.5f mm/s = %g mm/min), lateral motion fixed'
          % (lax + 1, -TRAVEL, T_STEP, SPEED, SPEED_MM_MIN))
run('boundary conditions (bottom plate fixed, top plate driven along %s)' % LOAD_AXIS, b_bc)


# ---- 8. contact: wall-to-wall self-contact only ------------------------------------------------
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
run('contact: lattice self-contact only, friction %g, hard normal' % FRICTION, b_contact)


# ---- 9. output requests ------------------------------------------------------------------------
def b_output():
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
run('output requests (energies incl. ALLSD and ALLFD, plate reference point U and RF, field output)', b_output)


# ---- 10. job and data check -----------------------------------------------------------------------
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


def scan_dat(name):
    """Print any 'unconnected regions' lines found in the .dat file."""
    path = os.path.join(os.getcwd(), name + '.dat')
    if not os.path.exists(path):
        print('      (no %s.dat found to scan)' % name)
        return
    f = open(path)
    try:
        lines = f.read().splitlines()
    finally:
        f.close()
    hits = [ln.strip() for ln in lines if 'unconnected region' in ln.lower()]
    if hits:
        for ln in hits[:3]:
            print('      ' + ln)
        print('      -> the plate faces must be coupled to the reference points; an unconnected region is a problem here')
    else:
        print('      no "unconnected regions" message in %s.dat' % name)
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
print('loaded axis        : %s,  core H = %.4f mm;  width %.4f mm;  depth %g mm' % (LOAD_AXIS, H_LOAD, W_TRANS, DEPTH))
print('end plates         : %g mm thick each, merged (bonded) to the core, lateral motion fixed at both' % PLATE_T)
if mass_t is not None:
    print('mass (core+plates) : %.5e tonne  (%.4f g)' % (mass_t, mass_t * 1.0e6))
print('plate travel       : %.3f mm (no gap)' % TRAVEL)
print('plate speed        : %g mm/min = %.5f mm/s;  step time %.2f s' % (SPEED_MM_MIN, SPEED, T_STEP))
print('nominal strain     : %.1f %% at the end (travel / core height, the plates are not counted in H)'
      % (100.0 * COMPRESS / H_LOAD))
print('friction           : %g (lattice self-contact only)' % FRICTION)
print('EXTRACTION         : in extract_results.py set GAP_TOP = 0.0 and MASS_T from the .dat; H and W stay as they are.')
print('NEXT: open %s.dat and %s.msg; check element counts, warnings and any' % (jname, jname))
print('      "zero pivot" or "numerical singularity" lines. Then set DATACHECK_ONLY = False.')
