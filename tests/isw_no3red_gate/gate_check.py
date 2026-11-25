#!/usr/bin/env python3
"""Decision of the isw_no3red port gate (run by gate.pjsub). jsasaki 2026-11-25: new file.

  gate_check.py make_yaml SRC DST SWITCH   SWITCH = none (copy unchanged) | 0 | 1 (isw_no3red for every ersem/primary_producer)
  gate_check.py make_env SRC DST TIME     constant environment (or wind) file around TIME (step runs)
  gate_check.py pick_state NC              the step start time (and +30 min, +2 h) from the new_off year run
  gate_check.py set_state YAML NC TIME     every initial value of YAML := the state in NC at TIME
  gate_check.py check RUNS                 exit 0 = PASS, 1 = FAIL, 2 = could not decide (fail-closed)

Checks (skill bit-identity-gates):
  G0  every run finished: output.nc readable, last time = the stop time of its run.nml
  G1  floor: base_r1 == base_r2, every variable, every value, bit for bit
  G2  switch off: new_off == base_r1 and new_off0 == base_r1, same variable set, bit for bit
  G3  switch on, whole year: P<n>_ledger_no3red_O2o == 2 * P<n>_fN3PIn at every record (exact), not all zero;
      the diagnostic is absent in every run with the switch off
  G4  one Euler step (step_off, step_on) from the new_off year-run state at the record where all of P1-P4 take up nitrate
      (pick_state, set_state), forcing held constant (make_env for PAR / temperature / salinity and for wind): after the
      first step every STATE variable except O2_o is identical (diagnostics evaluated from the new state may differ, e.g.
      O2_AOU, O2_fair, pel_nit_fN4N3n); every P<n>_fN3PIn of the step's rates is > 0; and O2_o(on) - O2_o(off) =
      dt * sum_n P<n>_ledger_no3red_O2o / 86400 to the round-off of the O2 states (double-precision output), where the
      rates are read from record 0 of step_ref (same state, started one step later: fabm0d evaluates a step's rates at the
      time of the step's end): O2 gained = 2 x nitrate uptake
  jsasaki 2026-11-25: G4 docstring after review r1 #3; pick_state / set_state after review r1 #1
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


def make_env(src, dst, at):
    """The environment line at time `at`, repeated hourly from 2 h before to 4 h after it: constant forcing over the step
    runs, so that the rates of the first step are those written in record 0 (fabm0d evaluates the rates of a step with the
    environment of the step's end time, job 7685488)."""
    import datetime as dt
    t = dt.datetime.strptime(at, '%Y-%m-%d %H:%M:%S')
    line = next((l for l in Path(src).read_text().splitlines() if l.startswith(at)), None)
    if line is None:
        print(f'no environment line at {at}')
        return 2
    values = line[len(at):]
    Path(dst).write_text(''.join(f'{(t + dt.timedelta(hours=h)).strftime("%Y-%m-%d %H:%M:%S")}{values}\n' for h in range(-2, 5)))
    return 0


def states_of(yaml_path):
    """Output names of the state variables (<instance>_<variable> of every initialization entry)."""
    import yaml
    d = yaml.safe_load(Path(yaml_path).read_text())
    return {f'{i}_{k}' for i, v in d['instances'].items() for k in (v.get('initialization') or {})}


#---> jsasaki 2026-11-25: review r1 #1: a step state in which all four producers take up nitrate
PRODUCERS = ('P1', 'P2', 'P3', 'P4')


def pick_state(nc_path):
    """The record of the year run where the smallest nitrate uptake of P1-P4 is largest (all > 0); prints its time,
    +30 min and +2 h as ISO strings (the step runs' start, reference start and stop)."""
    import datetime as dt
    import netCDF4
    with netCDF4.Dataset(nc_path) as d:
        up = np.array([np.asarray(d[f'{p}_fN3PIn'][:], float).reshape(len(d['time']), -1)[:, 0] for p in PRODUCERS])
        t = netCDF4.num2date(d['time'][:], d['time'].units, only_use_cftime_datetimes=False)
    low = up.min(axis=0); r = int(np.argmax(low))
    if not low[r] > 0:
        print('no record where all producers take up nitrate', file=sys.stderr)
        return 2
    t0 = dt.datetime(t[r].year, t[r].month, t[r].day, t[r].hour, t[r].minute, t[r].second)
    print(*(x.strftime('%Y-%m-%dT%H:%M:%S') for x in (t0, t0 + dt.timedelta(minutes=30), t0 + dt.timedelta(hours=2))))
    return 0


def set_state(yaml_path, nc_path, at):
    """Replace every initialization value in YAML by the state of NC at time `at` (exact repr of the double value)."""
    import datetime as dt
    import netCDF4
    states = states_of(yaml_path)
    with netCDF4.Dataset(nc_path) as d:
        t = netCDF4.num2date(d['time'][:], d['time'].units, only_use_cftime_datetimes=False)
        target = dt.datetime.strptime(at, '%Y-%m-%d %H:%M:%S')
        hit = [i for i, x in enumerate(t) if (x.year, x.month, x.day, x.hour, x.minute, x.second) ==
               (target.year, target.month, target.day, target.hour, target.minute, target.second)]
        if len(hit) != 1:
            print(f'time {at} not in {nc_path}', file=sys.stderr)
            return 2
        if d['O2_o'].dtype != np.float64:
            print('the year run must write double precision', file=sys.stderr)
            return 2
        value = {k: float(np.asarray(d[k][hit[0]]).ravel()[0]) for k in states}
    out, inst, in_init, n = [], None, False, 0
    for line in Path(yaml_path).read_text().splitlines(keepends=True):
        m = re.match(r'^  (\S[^:]*):\s*$', line)
        if m:
            inst, in_init = m.group(1), False
        elif re.match(r'^    \S', line):
            in_init = bool(re.match(r'^    initialization:\s*$', line))
        elif in_init:
            m = re.match(r'^(      )([^:#\s]+):\s*([^#\n]*?)(\s*#.*)?$', line.rstrip('\n'))
            if m and f'{inst}_{m.group(2)}' in value:
                line = f'{m.group(1)}{m.group(2)}: {value[f"{inst}_{m.group(2)}"]!r}{m.group(4) or ""}\n'; n += 1
        out.append(line)
    if n != len(states):
        print(f'replaced {n} of {len(states)} initial values', file=sys.stderr)
        return 2
    Path(yaml_path).write_text(''.join(out))
    return 0
#<--- jsasaki 2026-11-25


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
    names = ['base_r1', 'base_r2', 'new_off', 'new_off0', 'new_on', 'step_off', 'step_on', 'step_ref']
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
        if n not in ('new_on', 'step_on', 'step_ref'):
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
    nt = len(t0)
    common = sorted(k for k in set(s0) & set(s1) if k != 'time' and s0[k].ndim and s0[k].shape[0] == nt)   # time series only
    fixed = sorted(k for k in set(s0) & set(s1) if k not in common and k != 'time')
    if any(s0[k].tobytes() != s1[k].tobytes() for k in fixed):
        raise ValueError('a variable without the time dimension differs between the step runs')
    differ = [k for k in common if s0[k][i].tobytes() != s1[k][i].tobytes()]
    init_same = [k for k in common if s0[k][0].tobytes() != s1[k][0].tobytes()]
    print(f'G4 record 0 (t = {t0[0]:.0f} s) identical except the new diagnostic: {"PASS" if not init_same else "FAIL " + str(init_same)}')
    ok &= not init_same
    # diagnostics written in a record are evaluated from that record's state, so the O2-dependent ones follow O2_o;
    # among the STATE variables only O2_o may differ after the first step
    states = states_of(runs / 'step_on' / 'fabm.yaml')
    if not {'O2_o', 'N3_n', 'P1_c'} <= states:
        raise ValueError('state list from fabm.yaml incomplete')
    dstates = sorted(set(differ) & states)
    good4 = dstates == ['O2_o']
    print(f'G4 record {i} (t = {t0[i]:.0f} s): state variables that differ: {dstates}; diagnostics that differ: '
          f'{sorted(set(differ) - states)}: {"PASS" if good4 else "FAIL"}')
    ok &= good4
    for k in ('O2_o', 'P1_ledger_no3red_O2o', 'P1_fN3PIn'):
        if s1[k].dtype != np.float64:
            raise ValueError(f'{k} is {s1[k].dtype}: the step runs must write double precision (output.yaml xtype: double)')
    # Euler: the step from record 0 to record 1 uses the record-0 state with the time of record 1 (zenith angle); step_ref
    # starts from the same initial state at that time, so its record 0 holds exactly those rates
    sr = D['step_ref'][0]
    if i != 1:
        raise ValueError('the first step must end at record 1')
    for k in states:
        if k in sr and k in s1 and sr[k][0].tobytes() != s1[k][0].tobytes():
            raise ValueError(f'step_ref does not start from the initial state of step_on: {k}')
    # jsasaki 2026-11-25: review r1 #1: every producer must contribute to the step, or its O2 source is not tested
    ups = {p: float(sr[f'{p}_fN3PIn'][0].ravel()[0]) for p in PRODUCERS}
    allup = all(v > 0 for v in ups.values())
    print(f'G4 nitrate uptake of the step (step_ref record 0, mmol N/m3/d): {ups}: {"PASS" if allup else "FAIL"}')
    ok &= allup
    o2_on, o2_off = float(s1['O2_o'][i].ravel()[0]), float(s0['O2_o'][i].ravel()[0])
    gained = o2_on - o2_off
    exp = sum(float(sr[f'{p}_ledger_no3red_O2o'][0].ravel()[0]) for p in ('P1', 'P2', 'P3', 'P4')) * dt / 86400.0
    up = sum(float(sr[f'{p}_fN3PIn'][0].ravel()[0]) for p in ('P1', 'P2', 'P3', 'P4')) * dt / 86400.0
    prev = sum(float(s1[f'{p}_ledger_no3red_O2o'][0].ravel()[0]) for p in ('P1', 'P2', 'P3', 'P4')) * dt / 86400.0
    print(f'    (record 0 of step_on, evaluated at its own start time: dt x sum ledger {prev:.17g}; orientation only)')
    tol = 4 * np.finfo(float).eps * max(abs(o2_on), abs(o2_off))   # round-off of the two O2 states
    good = exp > 0 and abs(gained - exp) <= tol and abs(2 * up - exp) <= 4 * np.finfo(float).eps * exp
    print(f'G4 O2 gained in the step {i - 1} -> {i}: {gained:.17g} mmol O2/m3; dt x sum ledger_no3red_O2o (step_ref record 0) {exp:.17g}; '
          f'2 x dt x sum fN3PIn {2 * up:.17g}; |difference| {abs(gained - exp):.3g} (round-off bound {tol:.3g}): {"PASS" if good else "FAIL"}')
    ok &= good
    print('GATE', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    try:
        if sys.argv[1] == 'make_yaml':
            sys.exit(make_yaml(*sys.argv[2:5]))
        if sys.argv[1] == 'pick_state':
            sys.exit(pick_state(sys.argv[2]))
        if sys.argv[1] == 'set_state':
            sys.exit(set_state(*sys.argv[2:5]))
        if sys.argv[1] == 'make_env':
            sys.exit(make_env(*sys.argv[2:5]))
        if sys.argv[1] == 'check':
            sys.exit(check(sys.argv[2]))
        print('usage: see the module docstring')
        sys.exit(2)
    except SystemExit:
        raise
    except Exception as exc:   # anything that prevents the decision is "could not decide", never "differs"
        print(f'COULD NOT DECIDE: {type(exc).__name__}: {exc}')
        sys.exit(2)
