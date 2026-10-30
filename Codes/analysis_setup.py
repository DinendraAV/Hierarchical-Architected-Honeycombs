# -*- coding: mbcs -*-
# ---------------------------------------------------------------------------
# continue_from_assembly.py
# Continues an EXISTING model whose Part, Property and Assembly work is finished:
# honeycomb instance + two platen instances in place, TPU material and section assigned.
# Abaqus/CAE 2020 (Python 2.7).  Run:  File > Run Script.   SAVE A BACKUP .cae FIRST.
#
# What it does (nothing is assumed about axes):
#   - finds the loading axis from the platens (the axis in which both are flat)
#   - sets equal gaps, centres the platens on the patch
#   - meshes the honeycomb (C3D8R explicit, enhanced hourglass, swept hex) and the platens
#   - replaces reference points: one at each platen centre, with rigid-body constraints
#   - explicit step, smooth-step amplitude, mass scaling, BCs, general contact, outputs
#   - writes the .inp and runs a data check (or the full run)
# It does NOT touch the geometry, the material or the section assignment.
# It DELETES existing steps, BCs, interactions, constraints, amplitudes and reference
# points in this model, then recreates them (so it is safe to run again).
# NOT tested in Abaqus by the author: it stops at the first failing block and says which.
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
MODEL_NAME   = 'Model-1'
HC_INSTANCE  = 'test1-1'
PLATEN_INSTANCES = None        # None = every other instance (must be exactly two)

COMPRESS = 20.0                # mm of patch compression after contact
GAP_TOP  = 0.05                # mm
GAP_BOT  = 0.05                # mm
T        = 0.3                 # s; average speed = (COMPRESS + gaps) / T
FRICTION = 0.3
CENTRE_PLATENS = True          # align platen centres with the patch centre

INPLANE_SIZE = 0.5             # mm
DEPTH_LAYERS = 6
PLATEN_MESH  = 7.5             # mm
DT_TARGET    = 2.0e-7          # s
CPUS         = 4

FIELD_INTERVALS = 100
HIST_INTERVALS  = 400
WRITE_ELASE     = True

TAG            = 'm050'        # goes into the job name
DATACHECK_ONLY = True          # True: data check only. False: full analysis.
SUBMIT         = True          # False: only write the .inp
OVERWRITE      = False         # False: refuse to run if the ODB already exists
SAVE_AS_NAME   = 'HC_longrun_setup'   # saves the .cae under this name in the work directory ('' = do not save)
# ===========================================================================

AX = 'XYZ'


def run(label, fn):
    try:
        out = fn()
        print('OK    : ' + label)
        return out
    except Exception as e:
        print('FAILED: %s  ->  %s' % (label, str(e)))
        traceback.print_exc()
        raise Exception('STOPPED at: ' + label + '. Send the lines above.')


def mean_of(vals):
    # Abaqus replaces the built-in sum() with its own function, so do not use sum() in this script
    t = 0.0
    n = 0
    for v in vals:
        t += v
        n += 1
    return t / n


def coords(inst):
    return [v.pointOn[0] for v in inst.vertices]


def extents(points):
    lo = [min(q[i] for q in points) for i in range(3)]
    hi = [max(q[i] for q in points) for i in range(3)]
    return lo, hi, [hi[i] - lo[i] for i in range(3)]


def flat_axis(points, tol=1e-6):
    """Index of the axis in which all points share one coordinate (a planar platen)."""
    lo, hi, ext = extents(points)
    flat = [i for i in range(3) if ext[i] <= tol]
    if len(flat) != 1:
        raise Exception('platen is not a plane normal to one axis (extents %s)' % str(ext))
    return flat[0]


def plan_platen_moves(patch_lo, patch_hi, lax, plane_top, plane_bot, gap_top, gap_bot,
                      pl_centre_top, pl_centre_bot, patch_centre, centre):
    """Translation vectors that give the requested gaps (and optionally centre the platens)."""
    moves = []
    for plane, target, c in ((plane_top, patch_hi[lax] + gap_top, pl_centre_top),
                             (plane_bot, patch_lo[lax] - gap_bot, pl_centre_bot)):
        v = [0.0, 0.0, 0.0]
        v[lax] = target - plane
        if centre:
            for i in range(3):
                if i != lax:
                    v[i] = patch_centre[i] - c[i]
        moves.append(tuple(v))
    return moves


# ---- 0. look at what is there ------------------------------------------------------------------
if MODEL_NAME not in mdb.models.keys():
    raise Exception('model %s not found. Models in this database: %s' % (MODEL_NAME, str(mdb.models.keys())))
m = mdb.models[MODEL_NAME]
a = m.rootAssembly
print('      instances: %s' % str(a.instances.keys()))
if HC_INSTANCE not in a.instances.keys():
    raise Exception('instance %s not found: set HC_INSTANCE' % HC_INSTANCE)
if PLATEN_INSTANCES is None:
    pnames = [k for k in a.instances.keys() if k != HC_INSTANCE]
else:
    pnames = list(PLATEN_INSTANCES)
if len(pnames) != 2:
    raise Exception('expected exactly two platen instances, found %s: set PLATEN_INSTANCES' % str(pnames))
hc = a.instances[HC_INSTANCE]
print('      materials: %s' % str(m.materials.keys()))
for mk in m.materials.keys():
    try:
        print('      %s hyperelastic table: %s' % (mk, str(m.materials[mk].hyperelastic.table)))
        print('      %s density table    : %s' % (mk, str(m.materials[mk].density.table)))
    except Exception:
        pass
print('      sections: %s' % str(m.sections.keys()))

pl_pts = dict((n, coords(a.instances[n])) for n in pnames)
laxes = [flat_axis(pl_pts[n]) for n in pnames]
if laxes[0] != laxes[1]:
    raise Exception('platens are flat in different axes %s: check the assembly' % str(laxes))
lax = laxes[0]
hc_lo, hc_hi, hc_ext = extents(coords(hc))
dax = [i for i in range(3) if i != lax and hc_ext[i] == min(hc_ext[j] for j in range(3) if j != lax)][0]
tax = [i for i in range(3) if i not in (lax, dax)][0]
H_LOAD, W_TRANS, D_DEPTH = hc_ext[lax], hc_ext[tax], hc_ext[dax]
patch_centre = [0.5 * (hc_lo[i] + hc_hi[i]) for i in range(3)]
print('      honeycomb extents  X %.3f  Y %.3f  Z %.3f' % tuple(hc_ext))
print('      LOADED axis = %s (H = %.3f mm), transverse = %s (%.3f mm), depth = %s (%.3f mm)'
      % (AX[lax], H_LOAD, AX[tax], W_TRANS, AX[dax], D_DEPTH))
print('      (loaded axis taken from the platens: both are flat in %s)' % AX[lax])
if abs(D_DEPTH - 10.0) > 0.1:
    print('WARNING: depth axis extent is %.3f, not the expected 10 mm' % D_DEPTH)

# platen order: top = larger coordinate along the loaded axis
planes = dict((n, mean_of([q[lax] for q in pl_pts[n]])) for n in pnames)
top_name = max(pnames, key=lambda n: planes[n])
bot_name = min(pnames, key=lambda n: planes[n])
if top_name == bot_name or planes[top_name] == planes[bot_name]:
    raise Exception('both platens are in the same plane')
if not (planes[bot_name] < patch_centre[lax] < planes[top_name]):
    raise Exception('the patch is not between the platens along %s: check the assembly' % AX[lax])
print('      top platen = %s, bottom platen = %s' % (top_name, bot_name))
for n in (top_name, bot_name):
    lo_p, hi_p, ext_p = extents(pl_pts[n])
    other = [i for i in range(3) if i != lax]
    print('      %s size %.2f x %.2f mm' % (n, ext_p[other[0]], ext_p[other[1]]))


# ---- 1. clear old definitions (so the script can be run again) -----------------------------------
def b_wipe():
    for repo in (m.boundaryConditions, m.interactions, m.constraints):
        for k in list(repo.keys()):
            del repo[k]
    for k in list(m.interactionProperties.keys()):
        del m.interactionProperties[k]
    for k in list(m.steps.keys()):
        if k != 'Initial':
            del m.steps[k]
    for k in list(m.amplitudes.keys()):
        del m.amplitudes[k]
    rps = [k for k in a.features.keys() if k.startswith('RP')]
    if rps:
        a.deleteFeatures(tuple(rps))
    for nm in ('HC_ALL', 'RP_TOP', 'RP_BOT'):
        if nm in a.sets.keys():
            del a.sets[nm]
    print('      cleared steps, BCs, interactions, constraints, amplitudes and %d reference point(s)' % len(rps))
run('clear old step/BC/contact/constraint/RP definitions', b_wipe)


# ---- 2. gaps and platen centring ------------------------------------------------------------------
def b_gaps():
    ctop = [mean_of([q[i] for q in pl_pts[top_name]]) for i in range(3)]
    cbot = [mean_of([q[i] for q in pl_pts[bot_name]]) for i in range(3)]
    moves = plan_platen_moves(hc_lo, hc_hi, lax, planes[top_name], planes[bot_name],
                              GAP_TOP, GAP_BOT, ctop, cbot, patch_centre, CENTRE_PLATENS)
    for name, v in ((top_name, moves[0]), (bot_name, moves[1])):
        if max(abs(x) for x in v) > 1e-9:
            a.translate(instanceList=(name,), vector=v)
        print('      %s moved by (%.4f, %.4f, %.4f)' % (name, v[0], v[1], v[2]))
    a.regenerate()
    t = [q[lax] for q in coords(a.instances[top_name])]
    b = [q[lax] for q in coords(a.instances[bot_name])]
    print('      gap top    = %.4f mm (target %.3f)' % (min(t) - hc_hi[lax], GAP_TOP))
    print('      gap bottom = %.4f mm (target %.3f)' % (hc_lo[lax] - max(b), GAP_BOT))
run('gaps %.3f / %.3f mm along %s' % (GAP_TOP, GAP_BOT, AX[lax]), b_gaps)


# ---- 3. mesh ---------------------------------------------------------------------------------------------
def mesh_honeycomb():
    inst = a.instances[HC_INSTANCE]
    dep = inst.dependent
    if dep:
        owner = m.parts[inst.partName]
        obj = owner
        cells = owner.cells
    else:
        owner = a
        obj = inst
        cells = inst.cells
    vp = dict((v.index, v.pointOn[0]) for v in obj.vertices)
    olo, ohi, oext = extents(list(vp.values()))
    odax = oext.index(min(oext))            # extrusion axis in the object's own frame
    depth = oext[odax]
    sel = []
    for e in obj.edges:
        ids = e.getVertices()
        if len(ids) != 2 or ids[0] == ids[1]:
            continue
        p1, p2 = vp[ids[0]], vp[ids[1]]
        d = [abs(p1[i] - p2[i]) for i in range(3)]
        straight = True
        for i in range(3):
            if i != odax and d[i] > 1e-3:
                straight = False
        if straight and abs(d[odax] - depth) < 1e-3:
            sel.append(e)
    print('      %s instance; depth %.3f mm; straight depth edges found: %d'
          % ('dependent' if dep else 'independent', depth, len(sel)))
    if len(sel) == 0:
        raise Exception('no depth edges found')
    elemtypes = (ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, hourglassControl=ENHANCED,
                          distortionControl=DEFAULT),
                 ElemType(elemCode=C3D6, elemLibrary=EXPLICIT),
                 ElemType(elemCode=C3D4, elemLibrary=EXPLICIT))
    owner.setElementType(regions=(cells,), elemTypes=elemtypes)
    if dep:
        owner.seedPart(size=INPLANE_SIZE, deviationFactor=0.1, minSizeFactor=0.1)
    else:
        a.seedPartInstance(regions=(inst,), size=INPLANE_SIZE, deviationFactor=0.1, minSizeFactor=0.1)
    owner.seedEdgeByNumber(edges=tuple(sel), number=DEPTH_LAYERS, constraint=FIXED)
    owner.setMeshControls(regions=cells, elemShape=HEX, technique=SWEEP)
    if dep:
        owner.generateMesh()
    else:
        a.generateMesh(regions=(inst,))
    n_el = len(obj.elements) if dep else len(inst.elements)
    print('      honeycomb elements: %d  (in-plane %.2f mm, %d depth layers)' % (n_el, INPLANE_SIZE, DEPTH_LAYERS))
    try:
        q = obj.verifyMeshQuality(criterion=ANALYSIS_CHECKS, regions=(cells,))
        print('      analysis checks: %s' % str(q))
    except Exception as e:
        print('      (mesh verify skipped: %s)' % str(e))
run('honeycomb mesh: C3D8R explicit, enhanced hourglass, swept hex', mesh_honeycomb)


def mesh_platens():
    done_parts = []
    for n in (top_name, bot_name):
        inst = a.instances[n]
        if inst.dependent:
            pt = m.parts[inst.partName]
            if inst.partName not in done_parts:
                pt.seedPart(size=PLATEN_MESH, deviationFactor=0.1, minSizeFactor=0.1)
                pt.generateMesh()
                done_parts.append(inst.partName)
            print('      %s: %d elements (part mesh)' % (n, len(pt.elements)))
        else:
            a.seedPartInstance(regions=(inst,), size=PLATEN_MESH, deviationFactor=0.1, minSizeFactor=0.1)
            a.generateMesh(regions=(inst,))
            print('      %s: %d elements' % (n, len(inst.elements)))
run('platen meshes (R3D4)', mesh_platens)


# ---- 4. reference points, rigid bodies -----------------------------------------------------------------------
def b_rigid():
    for tag, name in (('TOP', top_name), ('BOT', bot_name)):
        pts = coords(a.instances[name])
        c = [mean_of([q[i] for q in pts]) for i in range(3)]
        rpf = a.ReferencePoint(point=tuple(c))
        rp = a.referencePoints[rpf.id]
        a.Set(referencePoints=(rp,), name='RP_' + tag)
        m.RigidBody(name='RB-' + tag, refPointRegion=regionToolset.Region(referencePoints=(rp,)),
                    bodyRegion=regionToolset.Region(faces=a.instances[name].faces))
    a.Set(cells=a.instances[HC_INSTANCE].cells, name='HC_ALL')
run('reference points at the platen centres and rigid-body constraints', b_rigid)


# ---- 5. step, amplitude, mass scaling ---------------------------------------------------------------------------
TRAVEL = COMPRESS + GAP_TOP + GAP_BOT


def b_step():
    m.ExplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=T, improvedDtMethod=ON)
    m.steps['Step-1'].setValues(massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, 0.0,
        DT_TARGET, BELOW_MIN, 0, 0, 0.0, 0.0, 0, None),))
    m.SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((0.0, 0.0), (T, 1.0)))
run('explicit step %g s, smooth-step amplitude, mass scaling dt = %g' % (T, DT_TARGET), b_step)


# ---- 6. boundary conditions -----------------------------------------------------------------------------------------
def b_bc():
    m.EncastreBC(name='BC-Bot', createStepName='Initial', localCsys=None, region=a.sets['RP_BOT'])
    m.DisplacementBC(name='BC-Top', createStepName='Initial', region=a.sets['RP_TOP'],
                     u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET,
                     distributionType=UNIFORM, fieldName='', localCsys=None)
    kw = {'u%d' % (lax + 1): -TRAVEL}
    m.boundaryConditions['BC-Top'].setValuesInStep(stepName='Step-1', amplitude='Amp-1', **kw)
    print('      top platen: U%d = %g mm over %g s (average %.1f mm/s, smooth-step peak %.1f mm/s)'
          % (lax + 1, -TRAVEL, T, TRAVEL / T, 1.875 * TRAVEL / T))
run('boundary conditions (bottom fixed, top driven along %s)' % AX[lax], b_bc)


# ---- 7. contact ------------------------------------------------------------------------------------------------------------
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


# ---- 8. output requests ---------------------------------------------------------------------------------------------------------
def b_output():
    fo = m.fieldOutputRequests['F-Output-1']
    fo.setValues(variables=('S', 'LE', 'U', 'RF'), numIntervals=FIELD_INTERVALS, region=a.sets['HC_ALL'])
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


# ---- 9. report and job ---------------------------------------------------------------------------------------------------------------
mass_t = None
try:
    mass_t = a.getMassProperties(regions=(a.instances[HC_INSTANCE],))['mass']
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
    if SAVE_AS_NAME:
        mdb.saveAs(pathName=os.path.join(os.getcwd(), SAVE_AS_NAME))
        print('      saved %s.cae' % SAVE_AS_NAME)
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
print('NEXT: open %s.dat: check NUMBER OF ELEMENTS, warnings and any "surface intersection" lines.' % jname)
