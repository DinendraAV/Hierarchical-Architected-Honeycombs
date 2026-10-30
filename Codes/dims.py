# -*- coding: mbcs -*-
# Prints the size of the honeycomb patch in the open model and the arguments for analyse_run.py.
# Abaqus/CAE: File > Run Script. Change the three names below if yours differ.
from abaqus import *
from abaqusConstants import *

MODEL_NAME  = 'Model-1'
HC_INSTANCE = 'test1-1'

a = mdb.models[MODEL_NAME].rootAssembly
print('instances: %s' % str(a.instances.keys()))
hc = a.instances[HC_INSTANCE]
pts = [v.pointOn[0] for v in hc.vertices]
lo = [min([q[i] for q in pts]) for i in range(3)]
hi = [max([q[i] for q in pts]) for i in range(3)]
ext = [hi[i] - lo[i] for i in range(3)]
for i in range(3):
    print('%s: %.4f to %.4f   size %.4f mm' % ('XYZ'[i], lo[i], hi[i], ext[i]))

# loaded axis = the axis in which the platens are flat
lax = None
for n in a.instances.keys():
    if n == HC_INSTANCE:
        continue
    pp = [v.pointOn[0] for v in a.instances[n].vertices]
    e = [max([q[i] for q in pp]) - min([q[i] for q in pp]) for i in range(3)]
    flat = [i for i in range(3) if e[i] < 1e-6]
    print('platen %s: extents %.2f %.2f %.2f   flat in %s' % (n, e[0], e[1], e[2], ''.join(['XYZ'[i] for i in flat])))
    if len(flat) == 1:
        lax = flat[0]
dax = ext.index(min(ext))
if lax is not None:
    tax = [i for i in range(3) if i not in (lax, dax)][0]
    mass_g = None
    try:
        mass_g = a.getMassProperties(regions=(hc,))['mass'] * 1.0e6
    except Exception as e:
        print('mass not available: %s' % str(e))
    print('')
    print('loaded axis %s: H = %.3f   transverse W = %.3f   depth D = %.3f' % ('XYZ'[lax], ext[lax], ext[tax], ext[dax]))
    print('analyse_run.py: --axis %d --H %.3f --W %.3f --D %.3f%s' % (lax + 1, ext[lax], ext[tax], ext[dax],
          (' --mass-g %.3f' % mass_g) if mass_g is not None else ''))
