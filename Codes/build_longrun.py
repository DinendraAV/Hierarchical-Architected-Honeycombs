# -*- coding: mbcs -*-
# ---------------------------------------------------------------------------
# build_slice.py - plane-strain slice study for the hourglass energy gate.
# Abaqus/CAE 2020 (Python 2.7).  Run:  File > Run Script.
#
# Based on build_longrun.py.  Differences:
#   1. A list of cases (in-plane element size, number of depth layers) is built,
#      each in its own model, so the whole mesh study needs one script run.
#   2. PLANE_STRAIN = True fixes the depth displacement of every honeycomb node.
#      One element layer then represents the full depth of the specimen, and
#      forces, masses and energies stay comparable with the full 3D model.
#   3. The jobs are written for double precision.
#   4. The analyses are NOT run from CAE. The script writes the .inp files and
#      run_study.bat. Run the batch file from an Abaqus Command window, so the
#      solver does not depend on the CAE session.
# NOT tested in Abaqus by the author: read the printout and run the data check first.
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
STEP_PATH   = r'C:\Users\Dino\Documents\GitHub\Hierarchical-Architected-Honeycombs\geometry\rhc\test1.STEP'
PART_NAME   = 'HC'
SAVE_NAME   = 'LongRun_m025'         # the .cae holding all cases

# Axes of the STEP file, as in build_longrun.py.
LOAD_AXIS          = 'AUTO'          # 'X', 'Y', 'Z' or 'AUTO'
DEPTH_AXIS         = 'AUTO'          # 'X', 'Y', 'Z' or 'AUTO' (= smallest extent)
EXPECT_LOAD_EXTENT = 57.18           # mm, loaded dimension
EXPECT_DEPTH       = 50.0            # mm

# Cases: (in-plane element size in mm, number of element layers through the depth).
# With PLANE_STRAIN = True one layer is enough. The last case repeats one mesh
# with two layers: if its energies and forces match the one-layer case, the
# depth modelling is consistent.
CASES        = [(0.25, 1)]
PLANE_STRAIN = True

# Run type. GATE: short check of the energy ratios (about 52 mm/s average).
# LONG: compression to about 35 % strain (about 67 mm/s average).
RUN = 'LONG'
PRESETS = {'GATE': (2.0, 0.04), 'LONG': (20.0, 0.3)}     # (COMPRESS in mm, T in s)
COMPRESS, T = PRESETS[RUN]
GAP_TOP = 0.05                       # mm
GAP_BOT = 0.05                       # mm

# material (placeholder TPU, neo-Hookean, rate independent)
E_MPA    = 67.0
NU       = 0.49
RHO      = 1.22e-9                   # tonne/mm^3
FRICTION = 0.3

# platens
PLATEN_SIZE = 75.0                   # mm, square
PLATEN_MESH = 7.5                    # mm

# element and solver
HOURGLASS = 'ENHANCED'               # ENHANCED, RELAX_STIFFNESS, COMBINED, STIFFNESS or VISCOUS
DT_TARGET = 1.0e-7                   # s, mass-scaling target. None = no mass scaling
CPUS      = 4
RESTART_INTERVALS = 40           # restart states written per step. 0 = leave the restart request unchanged

# output
FIELD_INTERVALS = 100
HIST_INTERVALS  = 400
WRITE_ELASE     = True               # element artificial strain energy (locates hourglass energy)
PLATENS_IN_ODB  = True

# job
SUBMIT_DATACHECK = True              # True: data check of every case from CAE (takes seconds)
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


HG = {'ENHANCED': ENHANCED, 'RELAX_STIFFNESS': RELAX_STIFFNESS, 'COMBINED': COMBINED,
      'STIFFNESS': STIFFNESS, 'VISCOUS': VISCOUS}
AX = 'XYZ'
TRAVEL = COMPRESS + GAP_TOP + GAP_BOT


def build_case(size, layers):
    ssz = ('%g' % size).replace('.', 'p')
    tag = 'm%s_L%d' % (ssz, layers)
    model_name = ('PS_' if PLANE_STRAIN else 'S3D_') + tag
    res = {}
    print('')
    print('=== case %s : in-plane %.3f mm, %d depth layer(s), plane strain %s ==='
          % (model_name, size, layers, str(PLANE_STRAIN)))

    # ---- 1. model and geometry ------------------------------------------------
    def b_model():
        if not os.path.exists(STEP_PATH):
            raise Exception('STEP file not found: ' + STEP_PATH)
        if model_name in mdb.models.keys():
            del mdb.models[model_name]
        mdb.Model(name=model_name)
        mdb.openStep(STEP_PATH, scaleFromFile=OFF)
        mdb.models[model_name].PartFromGeometryFile(name=PART_NAME, geometryFile=mdb.acis,
            combine=False, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    run('model %s and STEP import' % model_name, b_model)

    m = mdb.models[model_name]
    p = m.parts[PART_NAME]

    pts = [v.pointOn[0] for v in p.vertices]
    lo, hi, ext, dax, lax, tax = axis_info(pts, LOAD_AXIS, EXPECT_LOAD_EXTENT, DEPTH_AXIS)
    cen = [0.5 * (lo[i] + hi[i]) for i in range(3)]
    H_LOAD, W_TRANS, D_DEPTH = ext[lax], ext[tax], ext[dax]
    print('      extents  X %.3f   Y %.3f   Z %.3f' % tuple(ext))
    print('      depth axis = %s (%.3f)   LOADED axis = %s (%.3f)   transverse = %s (%.3f)'
          % (AX[dax], D_DEPTH, AX[lax], H_LOAD, AX[tax], W_TRANS))
    if abs(H_LOAD - EXPECT_LOAD_EXTENT) > 0.5:
        print('WARNING: loaded extent %.3f differs from EXPECT_LOAD_EXTENT %.3f: check the geometry'
              % (H_LOAD, EXPECT_LOAD_EXTENT))
    if abs(D_DEPTH - EXPECT_DEPTH) > 0.5:
        print('WARNING: depth is %.3f, expected %.3f: check DEPTH_AXIS and the geometry' % (D_DEPTH, EXPECT_DEPTH))
    if max(W_TRANS, D_DEPTH) > PLATEN_SIZE - 5.0:
        print('WARNING: PLATEN_SIZE %.1f is too small for the %.1f x %.1f footprint' % (PLATEN_SIZE, W_TRANS, D_DEPTH))
    print('      depth layers: %d  (%.2f mm each, aspect ratio %.1f against %.2f mm in-plane)'
          % (layers, D_DEPTH / layers, (D_DEPTH / layers) / size, size))
    if PLANE_STRAIN:
        print('      plane strain: depth displacement fixed on all honeycomb nodes. '
              'The aspect ratio should not affect the in-plane solution: compare the 1-layer and 2-layer cases')
    elif (D_DEPTH / layers) / size > 6.0:
        print('WARNING: element aspect ratio above 6 without plane strain: expect high hourglass energy')

    # ---- 2. material and section ------------------------------------------------
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

    # ---- 3. honeycomb mesh --------------------------------------------------------
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
        p.seedPart(size=size, deviationFactor=0.1, minSizeFactor=0.1)
        p.seedEdgeByNumber(edges=tuple(sel), number=layers, constraint=FIXED)
        p.setMeshControls(regions=p.cells, elemShape=HEX, technique=SWEEP)
        p.generateMesh()
        res['elements'] = len(p.elements)
        print('      honeycomb elements: %d   nodes: %d   (in-plane %.2f mm, %d depth layers, hourglass %s)'
              % (len(p.elements), len(p.nodes), size, layers, HOURGLASS))
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

    # ---- 4. platen part -------------------------------------------------------------
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

    # ---- 5. assembly, platens, reference points, rigid bodies -------------------
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
    run('explicit step %g s, smooth-step amplitude, mass scaling dt = %s' % (T, str(DT_TARGET)), b_step)

    # ---- 7. boundary conditions ------------------------------------------------------
    def b_bc():
        m.EncastreBC(name='BC-Bot', createStepName='Initial', localCsys=None,
                     region=a.sets['RP_BOT'])
        m.DisplacementBC(name='BC-Top', createStepName='Initial',
                         region=a.sets['RP_TOP'],
                         u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, amplitude=UNSET,
                         distributionType=UNIFORM, fieldName='', localCsys=None)
        kw = {'u%d' % (lax + 1): -TRAVEL}
        m.boundaryConditions['BC-Top'].setValuesInStep(stepName='Step-1', amplitude='Amp-1', **kw)
        print('      top platen: U%d = %g mm over %g s  (average %.1f mm/s, peak of smooth step %.1f mm/s)'
              % (lax + 1, -TRAVEL, T, TRAVEL / T, 1.875 * TRAVEL / T))
    run('boundary conditions (bottom fixed, top driven along %s)' % AX[lax], b_bc)

    def b_planestrain():
        a.Set(nodes=a.instances['HC-1'].nodes, name='HC_NODES')
        kw = {'u%d' % (dax + 1): SET}
        m.DisplacementBC(name='BC-PlaneStrain', createStepName='Initial',
                         region=a.sets['HC_NODES'], amplitude=UNSET,
                         distributionType=UNIFORM, fieldName='', localCsys=None, **kw)
    if PLANE_STRAIN:
        run('plane strain: U%d = 0 on all honeycomb nodes' % (dax + 1), b_planestrain)

    # ---- 8. contact --------------------------------------------------------------------
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

    # ---- 9. output requests ---------------------------------------------------------------
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

    # ---- 10. job ---------------------------------------------------------------------------
    mass_t = None
    try:
        mp = a.getMassProperties(regions=(a.instances['HC-1'],))
        mass_t = mp['mass']
    except Exception as e:
        print('      (mass properties not available: %s)' % str(e))

    if DT_TARGET:
        dtstr = '_dt' + ('%g' % DT_TARGET).replace('.', 'p').replace('-', 'm')
    else:
        dtstr = '_dtnone'
    jname = '%s_T%s_C%s%s' % (model_name, ('%g' % T).replace('.', 'p'), ('%g' % COMPRESS).replace('.', 'p'), dtstr)
    res['job'] = jname

    def b_job():
        j = mdb.Job(name=jname, model=model_name, description='', type=ANALYSIS,
                    explicitPrecision=DOUBLE_PLUS_PACK, nodalOutputPrecision=FULL, numCpus=CPUS,
                    numDomains=CPUS, multiprocessingMode=DEFAULT, parallelizationMethodExplicit=DOMAIN,
                    memory=90, memoryUnits=PERCENTAGE, resultsFormat=ODB, activateLoadBalancing=False,
                    echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF)
        j.writeInput(consistencyChecking=OFF)
        print('      wrote %s.inp in %s' % (jname, os.getcwd()))
        patch_restart(jname)
        if SUBMIT_DATACHECK:
            j.submit(consistencyChecking=OFF, datacheckJob=True)
            j.waitForCompletion()
    run('job %s (input file%s)' % (jname, ' and data check' if SUBMIT_DATACHECK else ''), b_job)

    res['H'] = H_LOAD
    res['W'] = W_TRANS
    res['D'] = D_DEPTH
    res['mass_g'] = (mass_t * 1.0e6) if mass_t is not None else None
    return res


def write_batch(jobs):
    path = os.path.join(os.getcwd(), 'run_study.bat')
    f = open(path, 'w')
    f.write('@echo off\n')
    f.write('cd /d "%s"\n' % os.getcwd())
    for j in jobs:
        f.write('call abaqus job=%s input=%s.inp cpus=%d double=both output_precision=full ask_delete=OFF interactive\n'
                % (j, j, CPUS))
    f.write('echo all jobs finished\n')
    f.write('pause\n')
    f.close()
    return path


results = []
for (size, layers) in CASES:
    results.append(build_case(size, layers))

try:
    mdb.saveAs(pathName=os.path.join(os.getcwd(), SAVE_NAME))
    print('saved %s.cae' % SAVE_NAME)
except Exception as e:
    print('(could not save the .cae: %s)' % str(e))

bat = write_batch([r['job'] for r in results])

print('')
print('==================== SUMMARY ====================')
print('run type           : %s  (compression %.2f mm + gaps %.2f mm, T = %g s)' % (RUN, COMPRESS, GAP_TOP + GAP_BOT, T))
print('increments (approx): %.2e at dt = %s' % ((T / DT_TARGET) if DT_TARGET else 0.0, str(DT_TARGET)))
for r in results:
    print('job %-34s elements %7d' % (r['job'], r['elements']))
    if DT_TARGET:
        # rate measured on the 0.25 mm slice: about 5700 increments/min per 10,174 elements (double precision, 4 CPUs)
        est_h = (T / DT_TARGET) / (5711.0 * 10174.0 / r['elements']) / 60.0
        print('   estimated run time at least %.1f h (the increment falls below the target as the walls deform)' % est_h)
r0 = results[0]
print('analyse_run.py args: --H %.3f --W %.3f --D %.3f --gaps %.3f%s'
      % (r0['H'], r0['W'], r0['D'], GAP_TOP + GAP_BOT,
         (' --mass-g %.3f' % r0['mass_g']) if r0['mass_g'] is not None else ''))
print('batch file         : %s' % bat)
print('NEXT: check each .dat (NUMBER OF ELEMENTS, warnings), then run the batch file from an Abaqus Command window.')