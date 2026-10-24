"""CHEM-N6, part 1: layered-front statistics on the saved CHEM-N1 seed-0 runs
(no re-run). Each run kept its colormode-6 frame at t=6000 and the middle row
over steps 5400-6000. Fields: alive (red channel > 0) and a waste proxy (the
blue channel: environmental W on dead sites, private w on alive ones).
  ring_*  ring_spectrum of the final frame (b=4): one frame, so no moving fraction
  kw_*    kw_spectrum of the space-time plot (moving fronts against stationary)
Joined to the N1 results row (seed 0) for the parameters and the defect excess."""
import csv
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
while not os.path.exists(os.path.join(ROOT, 'C', 'evoca.h')):
    ROOT = os.path.dirname(ROOT)
sys.path.insert(0, os.path.join(ROOT, 'python'))
import evoca_waves as ew                                     # noqa: E402

N1 = os.path.join(HERE, '..', '2026-10-23_N1_twoway_waves')
PAT = re.compile(r'(\w+)_N(\d+)_ws([\d.]+)_wd([\d.]+)_tox([\d.]+)_x(\d)\.npz')


def fields(rgb):
    return (rgb[..., 0] > 0).astype(float), rgb[..., 2].astype(float) / 255


def main():
    rows = {}
    for r in csv.DictReader(open(os.path.join(N1, 'results.csv'))):
        if r['seed'] == '0':
            key = (r['group'], int(r['N']), float(r['param_waste_sink']), float(r['param_wdiff']),
                   float(r['param_waste_tox']), int(float(r['param_waste_exchange'])))
            rows[key] = r
    out = []
    for f in sorted(glob.glob(os.path.join(N1, 'kymo', '*.npz'))):
        g, n, ws, wd, tox, x = PAT.match(os.path.basename(f)).groups()
        key = (g, int(n), float(ws), float(wd), float(tox), int(x))
        d = np.load(f)
        a_last, w_last = fields(d['last'])
        a_ky, w_ky = fields(d['kymo'])
        row = dict(zip(('group', 'N', 'waste_sink', 'wdiff', 'waste_tox', 'waste_exchange'), key))
        for name, fr, ky in (('alive', a_last, a_ky), ('waste', w_last, w_ky)):
            r = ew.ring_spectrum(fr, b=4)
            k = ew.kw_spectrum(ky, b=4)
            row.update({f'ring_{name}_{s}': r[s] for s in ('wavelength', 'prominence', 'contrast', 'layers', 'isotropy')})
            row.update({f'kw_{name}_{s}': k[s] for s in ('moving', 'wavelength', 'speed', 'peak_share')})
            row[f'period_{name}'], row[f'period_{name}_peak'] = ew.site_period(ky)
        n1 = rows.get(key, {})
        row.update({'n1_surv_480_excess': n1.get('surv_480_excess', ''), 'n1_R_mean': n1.get('R_mean', ''),
                    'n1_alive_mean': n1.get('alive_mean', '')})
        out.append(row)
    with open(os.path.join(HERE, 'n1_fronts.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f'wrote {len(out)} rows')


if __name__ == '__main__':
    main()
