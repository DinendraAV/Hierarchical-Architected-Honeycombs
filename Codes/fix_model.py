# -*- coding: mbcs -*-
# Patches the open model (Abaqus/CAE 2020, Python 2.7).
# Run: File > Run Script. SAVE THE .cae FIRST (File > Save As a backup).
from abaqus import *
from abaqusConstants import *
from mesh import *

# ---------------- user settings ----------------
MODEL_NAME = 'Model-1'
JOB      = 'Job-1'
T        = 0.002      # step time [s]            (use 0.002 for a quick test)
DELTA    = 2.0      # platen travel [mm]       (use 2.0 with T = 0.002)
DT_TARGET = 2.0e-7   # mass-scaling target dt [s]
DEPTH_LAYERS = 6     # elements through the 10 mm depth
CPUS     = 4
RESEED   = False      # False = skip the mesh change
ONLY_MESH = False     # True = run only the mesh step (the other blocks already ran OK)
# -----------------------------------------------

m = mdb.models[MODEL_NAME]
a = m.rootAssembly


def step(label, fn, mesh=False):
    if ONLY_MESH and not mesh:
        return
    try:
        fn()
        print('OK    : ' + label)
    except Exception as e:
        print('FAILED: ' + label + '  ->  ' + str(e))


# 1. step time, amplitude, platen displacement ------------------------------
def f1():
    m.steps['Step-1'].setValues(timePeriod=T)
step('step time = %g s' % T, f1)


def f2():
    m.amplitudes['Amp-1'].setValues(timeSpan=STEP, data=((0.0, 0.0), (T, 1.0)))
step('Amp-1 now ramps 0 -> 1 over %g s' % T, f2)


def f3():
    m.boundaryConditions['BC-3'].setValuesInStep(stepName='Step-1', u2=-DELTA)
step('BC-3: U2 = %g mm' % -DELTA, f3)


# 2. history output on the reference point (not on the face set) -------------
def f4():
    rp = a.referencePoints[15]          # top platen RP, as in your journal
    a.Set(referencePoints=(rp,), name='RP_TOP')
    m.historyOutputRequests['H-Output-2'].setValues(region=a.sets['RP_TOP'])
step('H-Output-2 (U2, RF2) now on RP_TOP', f4)


# 3. field output only on the honeycomb (removes the R3D4 S/LE warnings) -----
def f5():
    inst = a.instances['HC-1']
    a.Set(cells=inst.cells, name='HC_ALL')
    m.fieldOutputRequests['F-Output-1'].setValues(region=a.sets['HC_ALL'])
step('F-Output-1 restricted to honeycomb set', f5)


# 4. Explicit element library, enhanced hourglass control -------------------
def f6():
    inst = a.instances['HC-1']
    a.setElementType(regions=(inst.cells,), elemTypes=(
        ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT,
                 hourglassControl=ENHANCED, distortionControl=DEFAULT),
        ElemType(elemCode=C3D6, elemLibrary=EXPLICIT),
        ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)))
step('HC-1 elements: C3D8R Explicit, enhanced hourglass', f6)


# 5. semi-automatic mass scaling (tuple copied from the form CAE records) ----
def f7():
    m.steps['Step-1'].setValues(massScaling=((SEMI_AUTOMATIC, MODEL,
        AT_BEGINNING, 0.0, DT_TARGET, BELOW_MIN, 0, 0, 0.0, 0.0, 0, None), ))
step('mass scaling, target dt = %g s' % DT_TARGET, f7)


# 6. parallel run ------------------------------------------------------------
def f8():
    mdb.jobs[JOB].setValues(numCpus=CPUS, numDomains=CPUS)
step('%s: %d CPUs' % (JOB, CPUS), f8)


# 7. fewer layers through the depth -------------------------------------------
def f9():
    inst = a.instances['HC-1']
    vp = dict((v.index, v.pointOn[0]) for v in inst.vertices)
    ext = [max(p[i] for p in vp.values()) - min(p[i] for p in vp.values())
           for i in range(3)]
    ax = ext.index(min(ext))            # the thin direction = extrusion axis
    depth = ext[ax]
    print('      bounding box extents (x, y, z) = %s' % str(ext))
    print('      extrusion axis = %s, depth = %.3f mm' % ('XYZ'[ax], depth))
    sel = []
    for e in inst.edges:
        ids = e.getVertices()
        if len(ids) != 2 or ids[0] == ids[1]:
            continue
        p, q = vp[ids[0]], vp[ids[1]]
        d = [abs(p[i] - q[i]) for i in range(3)]
        others_zero = True
        for i in range(3):
            if i != ax and d[i] > 1e-3:
                others_zero = False
        if others_zero and abs(d[ax] - depth) < 1e-3:
            sel.append(e)
    print('      straight depth edges found: %d' % len(sel))
    if len(sel) == 0:
        raise Exception('no depth edges found - mesh left unchanged')
    a.deleteMesh(regions=(inst,))
    a.seedPartInstance(regions=(inst,), size=0.5, deviationFactor=0.1,
                       minSizeFactor=0.1)
    a.seedEdgeByNumber(edges=tuple(sel), number=DEPTH_LAYERS, constraint=FIXED)
    a.setMeshControls(regions=inst.cells, elemShape=HEX, technique=SWEEP)
    a.generateMesh(regions=(inst,))
    print('      elements in HC-1: %d' % len(inst.elements))
if RESEED:
    step('reseed depth to %d layers and remesh' % DEPTH_LAYERS, f9, mesh=True)

print('Done. Now run Job > Data Check, then read the .dat/.sta file.')