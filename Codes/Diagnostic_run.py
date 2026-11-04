# -*- coding: mbcs -*-
# ---------------------------------------------------------------------------
# build_strain35.py - plane-strain slice of the honeycomb patch, set up to reproduce the
# 35 % run: 20 mm of compression in 0.4 s (platen travel 20.1 mm with the gaps), 0.25 mm mesh,
# 1 layer, plane strain, mass-scaling target 1e-7 s. Same script as build_strain50.py, with
# END_STRAIN, AVG_SPEED, SAVE_NAME, the model name and the batch file name changed.
# Abaqus/CAE 2020 (Python 2.7).  Run:  File > Run Script.
#
# One model, one job. The loading is defined by END_STRAIN and AVG_SPEED:
#   nominal strain = (platen displacement - gaps) / loaded height
#   platen travel  = END_STRAIN * loaded height + gaps
#   step time      = platen travel / AVG_SPEED
# The loaded height is read from the STEP geometry, so no number is typed twice.
#
# The script writes the .inp, sets the restart request, saves the .cae and writes
# run_strain35.bat. The analysis is NOT run from CAE: start the batch file from an
# Abaqus Command window (see the printout).
# Built from build_slice.py, which has been run in Abaqus. The loading and job-name
# logic are new and not yet tested in Abaqus: read the printout and run the data
# check first.
# ---------------------------------------------------------------------------
from abaqus import *
from abaqusConstants import *
import regionToolset
import os, re, traceback
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
STEP_PATH  = r'C:\Users\Dino\Documents\GitHub\Hierarchical-Architected-Honeycombs\geometry\rhc\test1.STEP'
PART_NAME  = 'HC'
SAVE_NAME  = 'Strain35'              # the .cae written to the work directory

# Axes of the STEP file.
LOAD_AXIS          = 'AUTO'          # 'X', 'Y', 'Z' or 'AUTO'
DEPTH_AXIS         = 'AUTO'          # 'X', 'Y', 'Z' or 'AUTO' (= smallest extent)
EXPECT_LOAD_EXTENT = 57.18           # mm, loaded dimension
EXPECT_DEPTH       = 50.0            # mm

# loading
END_STRAIN = 20.013 / 57.18            # nominal strain at the end of the step (20.013 mm of 57.18 mm = 35.0 %)
AVG_SPEED  = 50.25                   # mm/s, average platen speed: 20.1 mm in 0.4 s (the smooth step peaks at 1.875 times this)
GAP_TOP    = 0.05                    # mm
GAP_BOT    = 0.05                    # mm

# mesh
INPLANE_SIZE = 0.25                  # mm
DEPTH_LAYERS = 1
PLANE_STRAIN = True                  # fixes the depth displacement of every honeycomb node
HOURGLASS    = 'ENHANCED'            # ENHANCED, RELAX_STIFFNESS, COMBINED, STIFFNESS or VISCOUS

# material: neo-Hookean, rate independent. Enter your own coefficients.
# C10 and D1 below are placeholders (the values that E = 67 MPa, nu = 0.49 give): replace them.
C10_MPA    = 11.2416                 # MPa, C10 = initial shear modulus / 2
D1_PER_MPA = 0.001791                # 1/MPa, D1 = 2 / initial bulk modulus (must be above 0 in Abaqus/Explicit)
RHO        = 1.22e-9                 # tonne/mm^3
FRICTION   = 0.3

# platens
PLATEN_SIZE = 75.0                   # mm, square
PLATEN_MESH = 7.5                    # mm

# solver
DT_TARGET         = 5.0e-8           # s, mass-scaling target. None = no mass scaling
CPUS              = 4
RESTART_INTERVALS = 40               # restart states written per step. 0 = leave the restart request unchanged

# output
FIELD_INTERVALS = 100
HIST_INTERVALS  = 400
WRITE_ELASE     = True               # element artificial strain energy (locates hourglass energy)
PLATENS_IN_ODB  = True

# job
SUBMIT_DATACHECK = True              # True: run a data check from CAE (takes seconds)
# ===========================================================================

AX = 'XYZ'
HG = {'ENHANCED': ENHANCED, 'RELAX_STIFFNESS': RELAX_STIFFNESS, 'COMBINED': COMBINED,
      'STIFFNESS': STIFFNESS, 'VISCOUS': VISCOUS}
RATE_PER_10K = 10150.0               # increments per minute for 10,000 elements (measured: gate run, 4 CPUs, double precision)


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


def patch_restart(jname):
    """Set the number of restart intervals in the written .inp file (Abaqus/Explicit)."""
    if not RESTART_INTERVALS:
        return
    path = os.path.join(os.getcwd(), jname + '.inp')
    if not os.path.exists(path):
        print('      (restart request not changed: %s not found)' % path)
        return
    f = open(path, 'r')
    txt = f.read()
    f.close()
    pat = r'(?im)^\*Restart, write, number interval=\d+'
    if re.search(pat, txt):
        new = re.sub(pat, '*Restart, write, number interval=%d' % RESTART_INTERVALS, txt)
        print('      restart request set to %d intervals' % RESTART_INTERVALS)
    else:
        new = re.sub(r'(?im)^\*End Step',
                     '*Restart, write, number interval=%d, time marks=NO\n*End Step' % RESTART_INTERVALS,
                     txt, count=1)
        if new == txt:
            print('WARNING: no restart request written to %s.inp: no *Restart and no *End Step line found' % jname)
            return
        print('      restart request added: %d intervals' % RESTART_INTERVALS)
    f = open(path, 'w')
    f.write(new)
    f.close()


# ---- 1. model and geometry -------------------------------------------------
MODEL_NAME = 'S35_m%s_L%d' % (('%g' % INPLANE_SIZE).replace('.', 'p'), DEPTH_LAYERS)


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
cen = [0.5 * (lo[i] + hi[i]) for i in range(3)]
H_LOAD, W_TRANS, D_DEPTH = ext[lax], ext[tax], ext[dax]
print('      extents  X %.3f   Y %.3f   Z %.3f' % tuple(ext))
print('      depth axis = %s (%.3f)   LOADED axis = %s (%.3f)   transverse = %s (%.3f)'
      % (AX[dax], D_DEPTH, AX[lax], H_LOAD, AX[tax], W_TRANS))
if abs(H_LOAD - EXPECT_LOAD_EXTENT) > 0.5:
    print('WARNING: loaded extent %.3f differs from EXPECT_LOAD_EXTENT %.3f: check the geometry' % (H_LOAD, EXPECT_LOAD_EXTENT))
if abs(D_DEPTH - EXPECT_DEPTH) > 0.5:
    print('WARNING: depth is %.3f, expected %.3f: check DEPTH_AXIS and the geometry' % (D_DEPTH, EXPECT_DEPTH))
if max(W_TRANS, D_DEPTH) > PLATEN_SIZE - 5.0:
    print('WARNING: PLATEN_SIZE %.1f is too small for the %.1f x %.1f footprint' % (PLATEN_SIZE, W_TRANS, D_DEPTH))

# loading derived from the geometry
COMPRESS = END_STRAIN * H_LOAD
TRAVEL = COMPRESS + GAP_TOP + GAP_BOT
T = TRAVEL / AVG_SPEED
print('      loading: %.1f %% strain = %.3f mm of %.3f mm, platen travel %.3f mm incl. gaps'
      % (100.0 * END_STRAIN, COMPRESS, H_LOAD, TRAVEL))
print('      step time %.4f s, average %.1f mm/s, peak of smooth step %.1f mm/s'
      % (T, AVG_SPEED, 1.875 * AVG_SPEED))


# ---- 2. material and section -----------------------------------------------
def b_material():
    print('      neo-Hooke: C10 = %.4f MPa, D1 = %.6f 1/MPa' % (C10_MPA, D1_PER_MPA))
    mat = m.Material(name='TPU')
    mat.Density(table=((RHO,),))
    mat.Hyperelastic(materialType=ISOTROPIC, table=((C10_MPA, D1_PER_MPA),), testData=OFF,
                     type=NEO_HOOKE, volumetricResponse=VOLUMETRIC_DATA)
    m.HomogeneousSolidSection(name='Section-1', material='TPU', thickness=None)
    p.Set(cells=p.cells, name='HC_CELLS')
    p.SectionAssignment(region=p.sets['HC_CELLS'], sectionName='Section-1', offset=0.0,
                        offsetType=MIDDLE_SURFACE, offsetField='',
                        thicknessAssignment=FROM_SECTION)
run('TPU material and solid section', b_material)

n_elements = [0]


# ---- 3. honeycomb mesh -------------------------------------------------------
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
    p.seedEdgeByNumber(edges=tuple(sel), number=DEPTH_LAYERS, constraint=FIXED)
    p.setMeshControls(regions=p.cells, elemShape=HEX, technique=SWEEP)
    p.generateMesh()
    n_elements[0] = len(p.elements)
    print('      honeycomb elements: %d   nodes: %d   (in-plane %.2f mm, %d depth layer(s), hourglass %s)'
          % (len(p.elements), len(p.nodes), INPLANE_SIZE, DEPTH_LAYERS, HOURGLASS))
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


# ---- 4. platen part ----------------------------------------------------------
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


# ---- 5. assembly, platens, reference points, rigid bodies ---------------------
a = m.rootAssembly


def b_assembly():
    a.DatumCsysByDefault(CARTESIAN)
    a.Instance(name='HC-1', part=p, dependent=ON)
    plane_top = hi[lax] + GAP_TOP
    plane_bot = lo[lax] - GAP_BOT
    for name, plane in (('Platen-Top', plane_top), ('Platen-Bot', plane_bot)):
        a.Instance(name=name, part=m.parts['Platen'], dependent=ON)
        # the platen is built in its XY plane (normal Z): rotate the normal onto the load axis
        if lax == 0:
            a.rotate(instanceList=(name,), axisPoint=(0.0, 0.0, 0.0), axisDirection=(0.0, 10.0, 0.0), angle=90.0)
        elif lax == 1:
            a.rotate(instanceList=(name,), axisPoint=(0.0, 0.0, 0.0), axisDirection=(10.0, 0.0, 0.0), angle=-90.0)
        v = [cen[0], cen[1], cen[2]]
        v[lax] = plane
        a.translate(instanceList=(name,), vector=tuple(v))
        rpf = a.ReferencePoint(point=tuple(v))
        rp = a.referencePoints[rpf.id]
        rtag = 'TOP' if name == 'Platen-Top' else 'BOT'
        a.Set(referencePoints=(rp,), name='RP_' + rtag)
        m.RigidBody(name='RB-' + rtag, refPointRegion=regionToolset.Region(referencePoints=(rp,)),
                    bodyRegion=regionToolset.Region(faces=a.instances[name].faces))
    a.Set(cells=a.instances['HC-1'].cells, name='HC_ALL')
    print('      %s platen plane at %.4f (patch top %.4f, gap %.3f)' % (AX[lax], plane_top, hi[lax], GAP_TOP))
    print('      %s platen plane at %.4f (patch bottom %.4f, gap %.3f)' % (AX[lax], plane_bot, lo[lax], GAP_BOT))
run('assembly, platens, reference points, rigid bodies', b_assembly)


# ---- 6. step, amplitude, mass scaling -------------------------------------------
def b_step():
    m.ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=T, improvedDtMethod=ON)
    if DT_TARGET:
        m.steps['Step-1'].setValues(massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, 0.0,
            DT_TARGET, BELOW_MIN, 0, 0, 0.0, 0.0, 0, None),))
    else:
        print('      no mass scaling')
    m.SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((0.0, 0.0), (T, 1.0)))
run('explicit step %.4f s, smooth-step amplitude, mass scaling dt = %s' % (T, str(DT_TARGET)), b_step)


# ---- 7. boundary conditions -------------------------------------------------------
def b_bc():
    m.EncastreBC(name='BC-Bot', createStepName='Initial', localCsys=None,
                 region=a.sets['RP_BOT'])
    m.DisplacementBC(name='BC-Top', createStepName='Initial',
                     region=a.sets['RP_TOP'],
                     u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET,
                     distributionType=UNIFORM, fieldName='', localCsys=None)
    kw = {'u%d' % (lax + 1): -TRAVEL}
    m.boundaryConditions['BC-Top'].setValuesInStep(stepName='Step-1', amplitude='Amp-1', **kw)
    print('      top platen: U%d = %g mm over %.4f s' % (lax + 1, -TRAVEL, T))
run('boundary conditions (bottom fixed, top driven along %s)' % AX[lax], b_bc)


def b_planestrain():
    a.Set(nodes=a.instances['HC-1'].nodes, name='HC_NODES')
    kw = {'u%d' % (dax + 1): SET}
    m.DisplacementBC(name='BC-PlaneStrain', createStepName='Initial',
                     region=a.sets['HC_NODES'], amplitude=UNSET,
                     distributionType=UNIFORM, fieldName='', localCsys=None, **kw)
if PLANE_STRAIN:
    run('plane strain: U%d = 0 on all honeycomb nodes' % (dax + 1), b_planestrain)


# ---- 8. contact ----------------------------------------------------------------------
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


# ---- 9. output requests ------------------------------------------------------------------
def b_output():
    fo = m.fieldOutputRequests['F-Output-1']
    fo.setValues(variables=('S', 'LE', 'U', 'RF'), numIntervals=FIELD_INTERVALS)
    if not PLATENS_IN_ODB:
        fo.setValues(region=a.sets['HC_ALL'])
    if WRITE_ELASE:
        try:
            fo.setValues(variables=('S', 'LE', 'U', 'RF', 'ELASE'))
            print('      ELASE added to field output')
        except Exception as e:
            print('      ELASE not accepted (%s): continuing without it' % str(e))
    m.historyOutputRequests['H-Output-1'].setValues(
        variables=('ALLAE', 'ALLFD', 'ALLIE', 'ALLKE', 'ALLWK', 'ETOTAL'),
        numIntervals=HIST_INTERVALS)
    for rtag in ('TOP', 'BOT'):
        m.HistoryOutputRequest(name='H-Output-' + rtag, createStepName='Step-1',
            region=a.sets['RP_' + rtag], variables=('U1', 'U2', 'U3', 'RF1', 'RF2', 'RF3'),
            numIntervals=HIST_INTERVALS, sectionPoints=DEFAULT, rebar=EXCLUDE)
run('output requests (energies, platen U and RF, honeycomb field output)', b_output)


# ---- 10. job ---------------------------------------------------------------------------------
mass_t = None
try:
    mp = a.getMassProperties(regions=(a.instances['HC-1'],))
    mass_t = mp['mass']
except Exception as e:
    print('      (mass properties not available: %s)' % str(e))

if DT_TARGET:
    dtstr = 'dt' + ('%g' % DT_TARGET).replace('.', 'p').replace('-', 'm')
else:
    dtstr = 'dtnone'
jname = 'PS_m%s_L%d_E%d_T%s_%s' % (('%g' % INPLANE_SIZE).replace('.', 'p'), DEPTH_LAYERS,
                                   int(round(100.0 * END_STRAIN)), ('%.2f' % T).replace('.', 'p'), dtstr)


def b_job():
    j = mdb.Job(name=jname, model=MODEL_NAME, description='', type=ANALYSIS,
                explicitPrecision=DOUBLE_PLUS_PACK, nodalOutputPrecision=FULL, numCpus=CPUS,
                numDomains=CPUS, multiprocessingMode=DEFAULT, parallelizationMethodExplicit=DOMAIN,
                memory=90, memoryUnits=PERCENTAGE, resultsFormat=ODB, activateLoadBalancing=False,
                echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF)
    j.writeInput(consistencyChecking=OFF)
    print('      wrote %s.inp in %s' % (jname, os.getcwd()))
    try:
        mdb.saveAs(pathName=os.path.join(os.getcwd(), SAVE_NAME))
        print('      saved %s.cae' % SAVE_NAME)
    except Exception as e:
        print('      (could not save the .cae: %s)' % str(e))
    if SUBMIT_DATACHECK:
        j.submit(consistencyChecking=OFF, datacheckJob=True)
        j.waitForCompletion()
    # CAE rewrites the .inp on every submit, so the restart request is set last
    patch_restart(jname)
run('job %s (input file%s)' % (jname, ' and data check' if SUBMIT_DATACHECK else ''), b_job)

bat = os.path.join(os.getcwd(), 'run_strain35.bat')
f = open(bat, 'w')
f.write('@echo off\n')
f.write('cd /d "%s"\n' % os.getcwd())
f.write('call abaqus job=%s input=%s.inp cpus=%d double=both output_precision=full ask_delete=OFF interactive\n'
        % (jname, jname, CPUS))
f.write('echo job finished\n')
f.write('pause\n')
f.close()

print('')
print('==================== SUMMARY ====================')
print('job                : %s' % jname)
print('loaded axis        : %s (U%d),  H = %.3f mm;  transverse %.3f mm;  depth %.3f mm'
      % (AX[lax], lax + 1, H_LOAD, W_TRANS, D_DEPTH))
print('end strain         : %.1f %%  (%.3f mm of platen travel after the gaps)' % (100.0 * END_STRAIN, COMPRESS))
print('step time          : %.4f s at %.1f mm/s average' % (T, AVG_SPEED))
if mass_t is not None:
    print('honeycomb mass     : %.5e tonne  (%.3f g)' % (mass_t, mass_t * 1.0e6))
if DT_TARGET:
    n_inc = T / DT_TARGET
    rate = RATE_PER_10K * 10000.0 / max(1, n_elements[0])
    print('increments (approx): %.2e at dt = %g s' % (n_inc, DT_TARGET))
    print('run time           : at least %.1f h (the increment falls below the target as the walls deform and touch)'
          % (n_inc / rate / 60.0))
print('analyse_run.py args: --H %.3f --W %.3f --D %.3f --gaps %.3f --axis %d%s --sea-strain %.2f'
      % (H_LOAD, W_TRANS, D_DEPTH, GAP_TOP + GAP_BOT, lax + 1,
         (' --mass-g %.3f' % (mass_t * 1.0e6)) if mass_t is not None else '', 0.9 * END_STRAIN))
print('batch file         : %s' % bat)
print('NEXT: check the data-check .sta and .dat (PERCENT CHNG MASS, warnings), then in an Abaqus Command window:')
print('      start "abaqus strain35" /min cmd /c run_strain35.bat')
