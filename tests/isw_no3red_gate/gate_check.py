#!/usr/bin/env python3
"""Decision of the isw_no3red port gate (run by gate.pjsub). jsasaki 2026-11-25: new file.

  gate_check.py make_yaml SRC DST SWITCH   SWITCH = none (copy unchanged) | 0 | 1 (isw_no3red for every ersem/primary_producer)
  gate_check.py check RUNS                 exit 0 = PASS, 1 = FAIL, 2 = could not decide (fail-closed)

Checks (skill bit-identity-gates):
  G0  every run finished: output.nc readable, last time = the stop time of its run.nml
  G1  floor: base_r1 == base_r2, every variable, every value, bit for bit
  G2  switch off: new_off == base_r1 and new_off0 == base_r1, same variable set, bit for bit
  G3  switch on, whole year: P<n>_ledger_no3red_O2o == 2 * P<n>_fN3PIn at every record (exact), not all zero;
      the diagnostic is absent in every run with the switch off
  G4  one Euler step from the same state (step_off, step_on): at the first step every variable except O2_o is identical
      and O2_o(on) - O2_o(off) = dt * sum_n P<n>_ledger_no3red_O2o / 86400 (relative 1e-12): O2 gained = 2 x nitrate uptake
"""
import re
import sys
from pathlib import Path

import numpy as np


def make_yaml(src, dst, switch):
    text = Path(src).read_text()
    if switch == 'none':
        Path(dst).write_text(text)
        return 0
    if switch not in ('0', '1'):
        print(f'bad switch {switch}')
        return 2
    lines = text.splitlines(keepends=True)
    out, in_pp, n = [], False, 0
    for line in lines:
        out.append(line)
        if re.match(r'^  \S', line):           # a new instance
            in_pp = False
        if re.match(r'^    model:\s*ersem/primary_producer\s*$', line):
            in_pp = True
        if in_pp and re.match(r'^    parameters:\s*$', line):
            out.append(f'      isw_no3red: {switch}\n')
            in_pp, n = False, n + 1
    if n != 4:
        print(f'expected 4 producer instances, found {n}')
        return 2
    Path(dst).write_text(''.join(out))
    return 0


def load(path):
    import netCDF4
    data = {}
    with netCDF4.Dataset(path) as nc:
        nc.set_auto_mask(False)
        for name, var in nc.variables.items():
            data[name] = np.array(var[:])
        t = nc.variables['time']
        units = t.units
    return data, units


def stop_seconds(run):
    nml = (run / 'run.nml').read_text()
    import datetime as dt
    get = lambda k: dt.datetime.strptime(re.search(rf'^{k}="([^"]+)"', nml, re.M).group(1), '%Y-%m-%d %H:%M:%S')
    return (get('stop') - get('start')).total_seconds(), float(re.search(r'^dt=([0-9.]+)', nml, re.M).group(1))


def identical(a, b):
    """0 identical, 1 differs; raises on anything that cannot be compared."""
    if set(a) != set(b):
        raise ValueError(f'variable sets differ: {sorted(set(a) ^ set(b))}')
    diff = []
    for k in a:
        x, y = a[k], b[k]
        if x.shape != y.shape or x.dtype != y.dtype:
            raise ValueError(f'{k}: shape/dtype {x.shape} {x.dtype} vs {y.shape} {y.dtype}')
        if x.tobytes() != y.tobytes():
            diff.append(k)
    return diff


def check(runs):
    runs = Path(runs)
    names = ['base_r1', 'base_r2', 'new_off', 'new_off0', 'new_on', 'step_off', 'step_on']
    D, ok = {}, True
    for n in names:
        d, units = load(runs / n / 'output.nc')
        span, dt = stop_seconds(runs / n)
        tsec = d['time'].astype(float) * (1.0 if units.startswith('seconds') else float('nan'))
        good = np.isfinite(tsec[-1]) and abs(tsec[-1] - span) < 1e-6
        print(f'G0 {n}: {len(tsec)} records, last time {tsec[-1]:.0f} s, expected {span:.0f} s ({units}): '
              f'{"PASS" if good else "FAIL"}')
        ok &= bool(good)
        D[n] = (d, tsec, dt)
    for label, a, b in [('G1 floor base_r1 == base_r2', 'base_r1', 'base_r2'),
                        ('G2 new_off == base_r1', 'new_off', 'base_r1'),
                        ('G2 new_off0 == base_r1', 'new_off0', 'base_r1')]:
        diff = identical(D[a][0], D[b][0])
        print(f'{label}: {len(D[a][0])} variables, {len(diff)} differ {diff[:10]}: {"PASS" if not diff else "FAIL"}')
        ok &= not diff
    # G3
    on = D['new_on'][0]
    for n in names:
        if n not in ('new_on', 'step_on'):
            present = [k for k in D[n][0] if 'no3red' in k]
            print(f'G3 {n}: no3red variables {present}: {"PASS" if not present else "FAIL"}')
            ok &= not present
    tot = 0.0
    for p in ('P1', 'P2', 'P3', 'P4'):
        led, up = on[f'{p}_ledger_no3red_O2o'], on[f'{p}_fN3PIn']
        exact = led.tobytes() == (2.0 * up).tobytes()
        tot += float(np.abs(led).sum())
        print(f'G3 {p}: ledger_no3red_O2o == 2 fN3PIn at {led.size} records: {"PASS" if exact else "FAIL"}; '
              f'max ledger {led.max():.6g} mmol O2/m3/d, mean {led.mean():.6g}')
        ok &= exact
    print(f'G3 the term is active (sum |ledger| > 0): {"PASS" if tot > 0 else "FAIL"}')
    ok &= tot > 0
    off_y = D['new_off'][0]
    dO2 = on['O2_o'] - off_y['O2_o']
    print(f'    year: O2_o(on) - O2_o(off): final {dO2.ravel()[-1]:.6g}, max {dO2.max():.6g} mmol O2/m3 (orientation only)')
    # G4
    s0, t0, dt = D['step_off']
    s1, t1, _ = D['step_on']
    if not np.array_equal(t0, t1):
        raise ValueError('step runs have different time axes')
    i = int(np.flatnonzero(np.abs(t0 - dt) < 1e-6)[0])
    common = sorted(set(s0) & set(s1))
    differ = [k for k in common if k != 'time' and s0[k].ndim and s0[k][i].tobytes() != s1[k][i].tobytes()]
    init_same = [k for k in common if s0[k].ndim and s0[k][0].tobytes() != s1[k][0].tobytes() and 'no3red' not in k]
    print(f'G4 record 0 (t = {t0[0]:.0f} s) identical except the new diagnostic: {"PASS" if not init_same else "FAIL " + str(init_same)}')
    ok &= not init_same
    print(f'G4 record {i} (t = {t0[i]:.0f} s): variables that differ: {differ}: {"PASS" if differ == ["O2_o"] else "FAIL"}')
    ok &= differ == ['O2_o']
    gained = float(s1['O2_o'][i].ravel()[0] - s0['O2_o'][i].ravel()[0])
    rows = {}
    for j in (i - 1, i):
        rows[j] = sum(float(s1[f'{p}_ledger_no3red_O2o'][j].ravel()[0]) for p in ('P1', 'P2', 'P3', 'P4')) * dt / 86400.0
    up = sum(float(s1[f'{p}_fN3PIn'][i].ravel()[0]) for p in ('P1', 'P2', 'P3', 'P4')) * dt / 86400.0
    exp = rows[i]
    rel = abs(gained - exp) / abs(exp) if exp != 0 else float('inf')
    print(f'G4 O2 gained in the first step {gained:.17g} mmol O2/m3; dt x sum ledger (record {i}) {exp:.17g} '
          f'(record {i - 1}: {rows[i - 1]:.17g}); 2 x dt x sum fN3PIn {2 * up:.17g}; relative difference {rel:.3g}: '
          f'{"PASS" if rel < 1e-12 else "FAIL"}')
    ok &= rel < 1e-12
    print('GATE', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    try:
        if sys.argv[1] == 'make_yaml':
            sys.exit(make_yaml(*sys.argv[2:5]))
        if sys.argv[1] == 'check':
            sys.exit(check(sys.argv[2]))
        print('usage: see the module docstring')
        sys.exit(2)
    except SystemExit:
        raise
    except Exception as exc:   # anything that prevents the decision is "could not decide", never "differs"
        print(f'COULD NOT DECIDE: {type(exc).__name__}: {exc}')
        sys.exit(2)
