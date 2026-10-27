# -*- coding: mbcs -*-
# One script for a clean case: gaps + step time + amplitude + platen travel + job.
# Abaqus/CAE 2020 (Python 2.7). SAVE A BACKUP .cae FIRST. Safe to run repeatedly.
from abaqus import *
from abaqusConstants import *

# ---------------- user settings ----------------
MODEL_NAME = 'Model-1'
BASE_JOB   = 'Job-1'     # job to copy settings from (CPUs, mass scaling is in the step)
T          = 0.04        # step time [s]
COMPRESS   = 2.0         # compression of the patch after contact [mm]
GAP_TOP    = 0.1         # mm, top platen to patch
GAP_BOT    = 0.02        # mm, bottom platen to patch
DATACHECK_ONLY = True    # True: only data check (seconds). False: run the analysis.
WRITE_PLATENS = False    # True: also write platen motion to the ODB (untested; crashed CAE once)
TAG = 'm033'              # added to the job name, e.g. 'm033' for the 0.33 mm mesh, so old results are not overwritten
OVERWRITE = False        # False: refuse to run if an ODB with this job name already exists
# -----------------------------------------------

m = mdb.models[MODEL_NAME]
a = m.rootAssembly
U2_TOTAL = COMPRESS + GAP_TOP + GAP_BOT     # platen travel so that the patch is compressed by COMPRESS


def say(msg):
    print(msg)


# 1. gaps (idempotent: it measures the current distance and sets it to the gap)
hc = a.instances['HC-1']
ys = [v.pointOn[0][1] for v in hc.vertices]
ymin, ymax = min(ys), max(ys)
ymid = 0.5 * (ymin + ymax)
for name, gap in (('Part-2-1', GAP_TOP), ('Part-2-2', GAP_BOT)):
    p = a.instances[name]
    py = [v.pointOn[0][1] for v in p.vertices]
    if max(py) - min(py) > 1e-6:
        raise Exception(name + ' is not flat in Y - stop and check orientation')
    y = sum(py) / len(py)
    top = y > ymid
    dist = (y - ymax) if top else (ymin - y)      # + gap, - penetration
    dy = (gap - dist) if top else -(gap - dist)
    if abs(dy) > 1e-9:
        a.translate(instanceList=(name,), vector=(0.0, dy, 0.0))
    say('%s: %s platen was %.6f mm from the patch, moved %.6f mm -> gap %.4f mm'
        % (name, 'TOP' if top else 'BOTTOM', dist, dy, gap))
if not (a.instances['Part-2-1'].vertices[0].pointOn[0][1] > ymid):
    say('WARNING: Part-2-1 is not the top platen. BC-3 acts on referencePoints[15]; check which platen it drives.')

# 2. step time, amplitude, platen displacement
m.steps['Step-1'].setValues(timePeriod=T)
m.amplitudes['Amp-1'].setValues(timeSpan=STEP, data=((0.0, 0.0), (T, 1.0)))
m.boundaryConditions['BC-3'].setValuesInStep(stepName='Step-1', u2=-U2_TOTAL)

# 2b. write platen motion to the ODB (field output was restricted to the honeycomb)
if WRITE_PLATENS:
    for tag, inst_name in (('TOP', 'Part-2-1'), ('BOT', 'Part-2-2')):
        a.Set(name='PLATEN_' + tag, nodes=a.instances[inst_name].nodes)
        if 'F-Output-' + tag not in m.fieldOutputRequests.keys():
            m.FieldOutputRequest(name='F-Output-' + tag, createStepName='Step-1',
                                 variables=('U', 'UR'), region=a.sets['PLATEN_' + tag],
                                 numIntervals=40)
    say('field output for both platens: U, UR')
else:
    say('platen field output skipped (WRITE_PLATENS = False)')

# 3. read back
say('step time      : %g s' % m.steps['Step-1'].timePeriod)
say('Amp-1 data     : %s' % str(m.amplitudes['Amp-1'].data))
say('platen travel  : U2 = %g mm (compression after contact = %g mm)' % (-U2_TOTAL, COMPRESS))
say('average speed  : %g mm/s' % (U2_TOTAL / T))

# 4. job
jname = 'T%s_C%s' % (('%g' % T).replace('.', 'p'), ('%g' % COMPRESS).replace('.', 'p'))
if TAG:
    jname = jname + '_' + TAG
import os
say('honeycomb elements in the model: %d' % len(a.instances['HC-1'].elements))
if (not DATACHECK_ONLY) and (not OVERWRITE) and os.path.exists(os.path.join(os.getcwd(), jname + '.odb')):
    raise Exception('%s.odb already exists in %s - set TAG to a new name (or OVERWRITE = True)' % (jname, os.getcwd()))
if jname not in mdb.jobs.keys():
    mdb.Job(name=jname, objectToCopy=mdb.jobs[BASE_JOB])
say('job: ' + jname + ('  (DATA CHECK ONLY)' if DATACHECK_ONLY else '  (FULL RUN)'))
mdb.jobs[jname].submit(consistencyChecking=OFF, datacheckJob=DATACHECK_ONLY)
mdb.jobs[jname].waitForCompletion()
say('finished. Read ' + jname + '.dat / .msg for the surface-intersection warning.')