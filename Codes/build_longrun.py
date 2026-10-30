# -*- coding: mbcs -*-
# ---------------------------------------------------------------------------
# build_longrun.py  -  complete in-plane compression model from a STEP file.
# Abaqus/CAE 2020 (Python 2.7).  Run:  File > Run Script.
#
# Builds in a NEW model (MODEL_NAME); your old 'Model-1' is not touched:
#   STEP import -> TPU (neo-Hookean) -> C3D8R explicit, enhanced hourglass ->
#   two discrete-rigid platens with reference points and rigid-body constraints ->
#   general contact (friction 0.3, hard) -> explicit step, smooth-step amplitude,
#   semi-automatic mass scaling -> output requests -> job (.inp written).
# API calls follow your own recorded journal (Oct 21-22) wherever possible.
# NOT tested in Abaqus by the author: read the printout, then run a data check first.
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
STEP_PATH   = r'C:\path\to\your_patch.STEP'   # <-- EDIT: full path of the STEP file of the current patch
MODEL_NAME  = 'LongRun'
PART_NAME   = 'HC'

# Current patch: 57.18 (loaded) x 55.0 (transverse) x 50.0 (depth) mm.
# LOAD_AXIS : axis of the STEP file along which the platens press: 'X', 'Y', 'Z' or 'AUTO'.
#             AUTO picks the in-plane axis whose extent is closest to EXPECT_LOAD_EXTENT.
# DEPTH_AXIS: extrusion axis: 'X', 'Y', 'Z' or 'AUTO' (= the axis with the smallest extent).
#             With depth 50 and width 55 these are close, so set it explicitly if AUTO looks wrong.
LOAD_AXIS          = 'AUTO'
DEPTH_AXIS         = 'AUTO'
EXPECT_LOAD_EXTENT = 57.18     # mm, loaded dimension
EXPECT_DEPTH       = 50.0      # mm

# material (placeholder TPU, neo-Hookean, rate independent)
E_MPA = 67.0
NU    = 0.49
RHO   = 1.22e-9                # tonne/mm^3
FRICTION = 0.3

# loading
# Diagnostic run (hourglass tests): COMPRESS = 5, T = 0.1  (about 51 mm/s average).
# Long run to see the collapse:     COMPRESS = 20, T = 0.3 (about 67 mm/s average).
COMPRESS = 5.0                 # mm of patch compression after contact
GAP_TOP  = 0.05                # mm
GAP_BOT  = 0.05                # mm
T        = 0.1                 # s, step time (average speed = travel / T)

# platens
PLATEN_SIZE = 75.0             # mm, square
PLATEN_MESH = 7.5              # mm

# mesh
INPLANE_SIZE = 0.5             # mm
DEPTH_ELEMENT = 2.5            # mm, target element size through the depth (50 mm -> 20 layers, aspect ratio 5)
DEPTH_LAYERS  = None           # None = depth / DEPTH_ELEMENT; or give an integer
HOURGLASS     = 'ENHANCED'     # ENHANCED, RELAX_STIFFNESS, COMBINED, STIFFNESS or VISCOUS
DT_TARGET    = 2.0e-7          # s, mass-scaling target
CPUS         = 4

# output
FIELD_INTERVALS = 100
HIST_INTERVALS  = 400
WRITE_ELASE     = True         # element artificial strain energy (locates hourglass energy)
PLATENS_IN_ODB  = True         # True: field output on the whole model so the platens move in the animation
                               # (Abaqus then prints harmless 'LE/S not available for R3D4' warnings)

# job
TAG            = 'd50_m050'    # goes into the job name
DATACHECK_ONLY = True          # True: data check only. False: full analysis.
SUBMIT         = True          # False: only write the .inp
OVERWRITE      = False         # False: refuse to run if the ODB already exists
SAVE_CAE       = True          # save the .cae as MODEL_NAME in the work directory
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


def axis_info(points, load_axis, expect, depth_axis='AUTO'):
    """points: list of (x, y, z). Returns lo, hi, ext, depth axis, load axis, transverse axis."""
    lo = [min(q[i] for q in points) for i in range(3)]
    hi = [max(q[i] for q in points) for i in range(3)]
    ext = [hi[i] - lo[i] for i in range(3)]
    dax = ext.index(min(ext)) if depth_axis == 'AUTO' else 'XYZ'.index(depth_axis)
    inpl = [i for i in range(3) if i != dax]
    if load_axis == 'AUTO':
        if abs(ext[inpl[0]] - ext[inpl[1]]) < 0.5:
            raise Exception('in-plane extents %.3f and %.3f are too close for AUTO: set LOAD_AXIS to X, Y or Z'
                            % (ext[inpl[0]], ext[inpl[1]]))
        lax = min(inpl, key=lambda i: abs(ext[i] - expect))
    else:
        lax = 'XYZ'.index(load_axis)
        if lax == dax:
            raise Exception('LOAD_AXIS %s is the depth (extrusion) axis' % load_axis)
    tax = [i for i in inpl if i != lax][0]
    return lo, hi, ext, dax, lax, tax


# ---- 1. model and geometry ------------------------------------------------
def b_model():
    if not os.path.exists(STEP_PATH):
        raise Exception('STEP file not found: ' + STEP_PATH)
    if MODEL_NAME in mdb.models.keys():
        del mdb.models[MODEL_NAME]
    mdb.Model(name=MODEL_NAME)
    mdb.openStep(STEP_PATH, scaleFromFile=OFF)
    mdb.models[MODEL_NAME].PartFromGeometryFile(name=PART_NAME, geometryFile=mdb.acis,
        combine=False, dimensionality=THREE_D, type=DEFORMABLE_BODY)
run('model %s and STEP import' % MODEL_NAME, b_model)

m = mdb.models[MODEL_NAME]
p = m.parts[PART_NAME]

pts = [v.pointOn[0] for v in p.vertices]
lo, hi, ext, dax, lax, tax = axis_info(pts, LOAD_AXIS, EXPECT_LOAD_EXTENT, DEPTH_AXIS)
AX = 'XYZ'
H_LOAD, W_TRANS, D_DEPTH = ext[lax], ext[tax], ext[dax]
cen = [0.5 * (lo[i] + hi[i]) for i in range(3)]
print('      extents  X %.3f   Y %.3f   Z %.3f' % tuple(ext))
print('      depth axis = %s (%.3f)   LOADED axis = %s (%.3f)   transverse = %s (%.3f)'
      % (AX[dax], D_DEPTH, AX[lax], H_LOAD, AX[tax], W_TRANS))
if abs(H_LOAD - EXPECT_LOAD_EXTENT) > 0.5:
    print('WARNING: loaded extent %.3f differs from EXPECT_LOAD_EXTENT %.3f: check the geometry' % (H_LOAD, EXPECT_LOAD_EXTENT))
if abs(D_DEPTH - EXPECT_DEPTH) > 0.5:
    print('WARNING: depth is %.3f, expected %.3f: check DEPTH_AXIS and the geometry' % (D_DEPTH, EXPECT_DEPTH))
N_LAYERS = DEPTH_LAYERS if DEPTH_LAYERS else max(1, int(round(D_DEPTH / DEPTH_ELEMENT)))
print('      depth layers: %d  (%.2f mm each, aspect ratio %.1f against %.2f mm in-plane)'
      % (N_LAYERS, D_DEPTH / N_LAYERS, (D_DEPTH / N_LAYERS) / INPLANE_SIZE, INPLANE_SIZE))
if (D_DEPTH / N_LAYERS) / INPLANE_SIZE > 6.0:
    print('WARNING: element aspect ratio above 6: expect high hourglass energy. Raise DEPTH_LAYERS.')
if max(W_TRANS, D_DEPTH) > PLATEN_SIZE - 5.0:
    print('WARNING: PLATEN_SIZE %.1f is too small for the %.1f x %.1f footprint' % (PLATEN_SIZE, W_TRANS, D_DEPTH))


# ---- 2. material and section ----------------------------------------------
def b_material():
    C10 = E_MPA / (4.0 * (1.0 + NU))
    D1 = 6.0 * (1.0 - 2.0 * NU) / E_MPA
    print('      neo-Hooke: C10 = %.4f MPa, D1 = %.6f 1/MPa (E = %g, nu = %g)' % (C10, D1, E_MPA, NU))
    mat = m.Material(name='TPU')
    mat.Density(table=((RHO,),))
    mat.Hyperelastic(materialType=ISOTROPIC, table=((C10, D1),), testData=OFF,
                     type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)
    m.HomogeneousSolidSection(name='Section-1', material='TPU', thickness=None)
    p.Set(cells=p.cells, name='HC_CELLS')
    p.SectionAssignment(region=p.sets['HC_CELLS'], sectionName='Section-1', offset=0.0,
                        offsetType=MIDDLE_SURFACE, offsetField='',
                        thicknessAssignment=FROM_SECTION)
run('TPU material and solid section', b_material)


# ---- 3. honeycomb mesh (part level) -----------------------------------------
HG = {'ENHANCED': ENHANCED, 'RELAX_STIFFNESS': RELAX_STIFFNESS, 'COMBINED': COMBINED,
      'STIFFNESS': STIFFNESS, 'VISCOUS': VISCOUS}


def b_mesh_hc():
    p.setElementType(regions=(p.cells,), elemTypes=(
        ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, hourglassControl=HG[HOURGLASS],
                 distortionControl=DEFAULT),
        ElemType(elemCode=C3D6, elemLibrary=EXPLICIT),
        ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)))
    vp = dict((v.index, v.pointOn[0]) for v in p.vertices)
    sel = []
    for e in p.edges:
        ids = e.getVertices()
        if len(ids) != 2 or ids[0] == ids[1]:
            continue
        a_, b_ = vp[ids[0]], vp[ids[1]]
        d = [abs(a_[i] - b_[i]) for i in range(3)]
        straight = True
        for i in range(3):
            if i != dax and d[i] > 1e-3:
                straight = False
        if straight and abs(d[dax] - D_DEPTH) < 1e-3:
            sel.append(e)
    print('      straight depth edges found: %d' % len(sel))
    if len(sel) == 0:
        raise Exception('no depth edges found')
    p.seedPart(size=INPLANE_SIZE, deviationFactor=0.1, minSizeFactor=0.1)
    p.seedEdgeByNumber(edges=tuple(sel), number=N_LAYERS, constraint=FIXED)
    p.setMeshControls(regions=p.cells, elemShape=HEX, technique=SWEEP)
    p.generateMesh()
    print('      honeycomb elements: %d   nodes: %d   (in-plane %.2f mm, %d depth layers, hourglass %s)'
          % (len(p.elements), len(p.nodes), INPLANE_SIZE, N_LAYERS, HOURGLASS))
    try:
        um = p.getUnmeshedRegions()
        if um:
            raise Exception('part has unmeshed regions')
    except AttributeError:
        pass
    try:
        q = p.verifyMeshQuality(criterion=ANALYSIS_CHECKS, regions=(p.cells,))
        print('      analysis checks: %s' % str(q))
    except Exception as e:
        print('      (mesh verify skipped: %s)' % str(e))
run('honeycomb mesh: C3D8R explicit, %s hourglass, swept hex' % HOURGLASS, b_mesh_hc)


# ---- 4. platen part ---------------------------------------------------------
def b_platen():
    h = 0.5 * PLATEN_SIZE
    s = m.ConstrainedSketch(name='__profile__', sheetSize=200.0)
    s.rectangle(point1=(-h, -h), point2=(h, h))
    pp = m.Part(name='Platen', dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
    pp.BaseShell(sketch=s)
    del m.sketches['__profile__']
    pp.seedPart(size=PLATEN_MESH, deviationFactor=0.1, minSizeFactor=0.1)
    pp.generateMesh()
    print('      platen elements: %d' % len(pp.elements))
run('platen part (discrete rigid, %g mm square)' % PLATEN_SIZE, b_platen)


# ---- 5. assembly, platen placement, reference points, rigid bodies -------------
a = m.rootAssembly


def b_assembly():
    a.DatumCsysByDefault(CARTESIAN)
    a.Instance(name='HC-1', part=p, dependent=ON)
    plane_top = hi[lax] + GAP_TOP
    plane_bot = lo[lax] - GAP_BOT
    for name, plane in (('Platen-Top', plane_top), ('Platen-Bot', plane_bot)):
        a.Instance(name=name, part=m.parts['Platen'], dependent=ON)
        # platen is built in its XY plane (normal Z): rotate the normal onto the load axis
        if lax == 0:
            a.rotate(instanceList=(name,), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 10.0, 0.0), angle=90.0)
        elif lax == 1:
            a.rotate(instanceList=(name,), axisPoint=(0.0, 0.0, 0.0), axisDirection=(10.0, 0.0, 0.0), angle=-90.0)
        v = [cen[0], cen[1], cen[2]]
        v[lax] = plane
        a.translate(instanceList=(name,), vector=tuple(v))
        rpf = a.ReferencePoint(point=tuple(v))
        rp = a.referencePoints[rpf.id]
        tag = 'TOP' if name == 'Platen-Top' else 'BOT'
        a.Set(referencePoints=(rp,), name='RP_' + tag)
        m.RigidBody(name='RB-' + tag, refPointRegion=regionToolset.Region(referencePoints=(rp,)),
                    bodyRegion=regionToolset.Region(faces=a.instances[name].faces))
    a.Set(cells=a.instances['HC-1'].cells, name='HC_ALL')
    print('      %s platen plane at %.4f (patch top %.4f, gap %.3f)' % (AX[lax], plane_top, hi[lax], GAP_TOP))
    print('      %s platen plane at %.4f (patch bottom %.4f, gap %.3f)' % (AX[lax], plane_bot, lo[lax], GAP_BOT))
run('assembly, platens, reference points, rigid bodies', b_assembly)


# ---- 6. step, amplitude, mass scaling --------------------------------------------
TRAVEL = COMPRESS + GAP_TOP + GAP_BOT


def b_step():
    m.ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=T, improvedDtMethod=ON)
    m.steps['Step-1'].setValues(massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, 0.0,
        DT_TARGET, BELOW_MIN, 0, 0, 0.0, 0.0, 0, None),))
    m.SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((0.0, 0.0), (T, 1.0)))
run('explicit step %g s, smooth-step amplitude, mass scaling dt = %g' % (T, DT_TARGET), b_step)


# ---- 7. boundary conditions ------------------------------------------------------------
def b_bc():
    m.EncastreBC(name='BC-Bot', createStepName='Initial', localCsys=None,
                 region=a.sets['RP_BOT'])
    # all six DOF fixed in Initial (a rigid body needs this), then the load-axis DOF is driven in Step-1
    m.DisplacementBC(name='BC-Top', createStepName='Initial',
                     region=a.sets['RP_TOP'],
                     u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET,
                     distributionType=UNIFORM, fieldName='', localCsys=None)
    kw = {'u%d' % (lax + 1): -TRAVEL}
    m.boundaryConditions['BC-Top'].setValuesInStep(stepName='Step-1', amplitude='Amp-1', **kw)
    print('      top platen: U%d = %g mm over %g s  (average %.1f mm/s, peak of smooth step %.1f mm/s)'
          % (lax + 1, -TRAVEL, T, TRAVEL / T, 1.875 * TRAVEL / T))
run('boundary conditions (bottom fixed, top driven along %s)' % AX[lax], b_bc)


# ---- 8. contact ----------------------------------------------------------------------------
def b_contact():
    prop = m.ContactProperty('IntProp-1')
    prop.TangentialBehavior(formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF,
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0,
        table=((FRICTION,),), shearStressLimit=None, maximumElasticSlip=FRACTION,
        fraction=0.005, elasticSlipStiffness=None)
    prop.NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON,
                        constraintEnforcementMethod=DEFAULT)
    m.ContactExp(name='Int-1', createStepName='Step-1')
    m.interactions['Int-1'].includedPairs.setValuesInStep(stepName='Step-1', useAllstar=ON)
    m.interactions['Int-1'].contactPropertyAssignments.appendInStep(stepName='Step-1',
        assignments=((GLOBAL, SELF, 'IntProp-1'),))
run('general contact, all exterior, friction %g, hard normal' % FRICTION, b_contact)


# ---- 9. output requests ----------------------------------------------------------------------
def b_output():
    fo = m.fieldOutputRequests['F-Output-1']
    fo.setValues(variables=('S', 'LE', 'U', 'RF'), numIntervals=FIELD_INTERVALS)
    if not PLATENS_IN_ODB:
        fo.setValues(region=a.sets['HC_ALL'])      # honeycomb only: platens will NOT move in the animation
    if WRITE_ELASE:
        try:
            fo.setValues(variables=('S', 'LE', 'U', 'RF', 'ELASE'))
            print('      ELASE added to field output')
        except Exception as e:
            print('      ELASE not accepted (%s): continuing without it' % str(e))
    m.historyOutputRequests['H-Output-1'].setValues(
        variables=('ALLAE', 'ALLFD', 'ALLIE', 'ALLKE', 'ALLWK', 'ETOTAL'),
        numIntervals=HIST_INTERVALS)
    for tag in ('TOP', 'BOT'):
        m.HistoryOutputRequest(name='H-Output-' + tag, createStepName='Step-1',
            region=a.sets['RP_' + tag], variables=('U1', 'U2', 'U3', 'RF1', 'RF2', 'RF3'),
            numIntervals=HIST_INTERVALS, sectionPoints=DEFAULT, rebar=EXCLUDE)
run('output requests (energies, platen U and RF, honeycomb field output)', b_output)


# ---- 10. report, job ----------------------------------------------------------------------------
mass_t = None
try:
    mp = a.getMassProperties(regions=(a.instances['HC-1'],))
    mass_t = mp['mass']
except Exception as e:
    print('      (mass properties not available: %s)' % str(e))

jname = 'LR_T%s_C%s_%s' % (('%g' % T).replace('.', 'p'), ('%g' % COMPRESS).replace('.', 'p'), TAG)


def b_job():
    if (not DATACHECK_ONLY) and SUBMIT and (not OVERWRITE) and os.path.exists(os.path.join(os.getcwd(), jname + '.odb')):
        raise Exception('%s.odb exists in %s: change TAG or set OVERWRITE = True' % (jname, os.getcwd()))
    j = mdb.Job(name=jname, model=MODEL_NAME, description='', type=ANALYSIS,
                explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, numCpus=CPUS,
                numDomains=CPUS, multiprocessingMode=DEFAULT, parallelizationMethodExplicit=DOMAIN,
                memory=90, memoryUnits=PERCENTAGE, resultsFormat=ODB, activateLoadBalancing=False,
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
print('loaded axis        : %s,  H = %.3f mm;  transverse %.3f mm;  depth %.3f mm' % (AX[lax], H_LOAD, W_TRANS, D_DEPTH))
print('area (W x D)       : %.2f mm2' % (W_TRANS * D_DEPTH))
if mass_t is not None:
    print('honeycomb mass     : %.5e tonne  (%.3f g)' % (mass_t, mass_t * 1.0e6))
print('platen travel      : %.3f mm (compression %.3f + gaps %.3f)' % (TRAVEL, COMPRESS, GAP_TOP + GAP_BOT))
print('end strain         : %.1f %%' % (100.0 * COMPRESS / H_LOAD))
print('increments (approx): %.2e  at dt = %g s' % (T / DT_TARGET, DT_TARGET))
print('analyse_run.py args: --H %.3f --W %.3f --D %.3f --gaps %.3f%s'
      % (H_LOAD, W_TRANS, D_DEPTH, GAP_TOP + GAP_BOT,
         (' --mass-g %.3f' % (mass_t * 1.0e6)) if mass_t is not None else ''))
print('NEXT: open %s.dat, check NUMBER OF ELEMENTS, warnings, and "surface intersection" lines.' % jname)
