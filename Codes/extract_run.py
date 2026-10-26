# -*- coding: mbcs -*-
# Dump every history output of an Abaqus ODB to one CSV (long format).
# Run from the Abaqus command prompt (Python 2.7, no GUI needed):
#     abaqus python extract_run.py D:\Temporary\T0p01_C2.odb T0p01_C2
# Writes T0p01_C2_history.csv with columns: step, region, variable, time, value
import sys
import csv
from odbAccess import openOdb


def main(odb_path, prefix):
    odb = openOdb(odb_path, readOnly=True)
    out = open(prefix + '_history.csv', 'wb')
    w = csv.writer(out)
    w.writerow(['step', 'region', 'variable', 'time', 'value'])
    n = 0
    for step_name in odb.steps.keys():
        step = odb.steps[step_name]
        for reg_name in step.historyRegions.keys():
            reg = step.historyRegions[reg_name]
            for var_name in reg.historyOutputs.keys():
                for t, v in reg.historyOutputs[var_name].data:
                    w.writerow([step_name, reg_name, var_name, repr(t), repr(v)])
                    n += 1
    out.close()
    odb.close()
    print('wrote %d rows to %s_history.csv' % (n, prefix))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('usage: abaqus python extract_run.py <file.odb> <output_prefix>')
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
