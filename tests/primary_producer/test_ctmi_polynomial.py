#!/usr/bin/env python3
"""
Regression tests for the polynomial CTMI temperature response
and Michaelis-Menten nutrient uptake in primary_producer.F90.

Tests verify mathematical properties without requiring PyFABM.
Run with: python3 test_ctmi_polynomial.py
    or:   pytest test_ctmi_polynomial.py -v
"""

import sys
import math


# ===================================================================
# Polynomial CTMI implementation (mirrors primary_producer.F90)
# ===================================================================

def ctmi_coefficients(Tmin, Topt, Tmax):
    """Precompute polynomial CTMI coefficients (as in initialize)."""
    a = Topt - Tmin
    b = Topt - Tmax
    ab2 = (a * b) ** 2
    ctmi_a = -(a + b) / ab2
    ctmi_b = (a * b + (a + b) * Topt) / ab2
    return ctmi_a, ctmi_b


def ctmi_eval(T, Tmin, Topt, Tmax):
    """Evaluate polynomial CTMI (as in do subroutine)."""
    if T <= Tmin or T >= Tmax:
        return 0.0
    ctmi_a, ctmi_b = ctmi_coefficients(Tmin, Topt, Tmax)
    et = (T - Tmin) * (T - Tmax) * (ctmi_a * T + ctmi_b)
    return max(0.0, min(1.0, et))


def michaelis_menten(S, Vmax, Ks):
    """Michaelis-Menten saturation: V = Vmax * S / (Ks + S)."""
    return Vmax * S / (Ks + S)


# ===================================================================
# Test parameter sets
# ===================================================================

# Tokyo Bay phytoplankton groups (from doc/phytoplankton_temperature_ctmi.md)
GROUPS_DOC = {
    'P1_diatoms':   (2.0, 15.0, 30.0),
    'P2_nano':      (5.0, 20.0, 33.0),
    'P3_pico':      (8.0, 25.0, 35.0),
    'P4_dino':      (10.0, 25.0, 35.0),
}

# TB-GOTM production parameters
GROUPS_GOTM = {
    'P1_GOTM': (0.0, 12.0, 25.0),
    'P2_GOTM': (10.0, 25.0, 35.0),
    'P3_GOTM': (5.0, 20.0, 33.0),
    'P4_GOTM': (8.0, 25.0, 35.0),
    'P5_GOTM': (10.0, 25.0, 35.0),
}

# Edge cases: Topt at or near midpoint
GROUPS_EDGE = {
    'symmetric':     (0.0, 15.0, 30.0),   # Topt exactly at midpoint
    'near_midpoint': (0.0, 14.9, 30.0),   # Topt just below midpoint
    'wide_range':    (0.0, 10.0, 40.0),   # Very wide range, Topt far from midpoint
    'narrow_range':  (10.0, 15.0, 20.0),  # Narrow range
}

ALL_GROUPS = {**GROUPS_DOC, **GROUPS_GOTM, **GROUPS_EDGE}


# ===================================================================
# Tests
# ===================================================================

def test_f_at_topt_equals_one():
    """f(Topt) must equal 1.0 exactly for all parameter sets."""
    for name, (Tmin, Topt, Tmax) in ALL_GROUPS.items():
        val = ctmi_eval(Topt, Tmin, Topt, Tmax)
        assert abs(val - 1.0) < 1e-12, \
            f"{name}: f(Topt={Topt}) = {val}, expected 1.0"


def test_f_at_tmin_equals_zero():
    """f(Tmin) must equal 0.0."""
    for name, (Tmin, Topt, Tmax) in ALL_GROUPS.items():
        val = ctmi_eval(Tmin, Tmin, Topt, Tmax)
        assert val == 0.0, f"{name}: f(Tmin={Tmin}) = {val}, expected 0.0"


def test_f_at_tmax_equals_zero():
    """f(Tmax) must equal 0.0."""
    for name, (Tmin, Topt, Tmax) in ALL_GROUPS.items():
        val = ctmi_eval(Tmax, Tmin, Topt, Tmax)
        assert val == 0.0, f"{name}: f(Tmax={Tmax}) = {val}, expected 0.0"


def test_f_bounded_zero_one():
    """f(T) must be in [0, 1] for all T in [Tmin, Tmax]."""
    for name, (Tmin, Topt, Tmax) in ALL_GROUPS.items():
        n_points = 1000
        for i in range(n_points + 1):
            T = Tmin + (Tmax - Tmin) * i / n_points
            val = ctmi_eval(T, Tmin, Topt, Tmax)
            assert 0.0 <= val <= 1.0, \
                f"{name}: f({T:.3f}) = {val}, outside [0, 1]"


def test_f_outside_range_is_zero():
    """f(T) must be 0 for T <= Tmin or T >= Tmax."""
    for name, (Tmin, Topt, Tmax) in ALL_GROUPS.items():
        for T in [Tmin - 10.0, Tmin - 1.0, Tmin, Tmax, Tmax + 1.0, Tmax + 10.0]:
            val = ctmi_eval(T, Tmin, Topt, Tmax)
            assert val == 0.0, \
                f"{name}: f({T}) = {val}, expected 0.0 outside range"


def test_f_prime_at_topt_is_zero():
    """f'(Topt) must be approximately 0 (numerical derivative)."""
    h = 1e-6
    for name, (Tmin, Topt, Tmax) in ALL_GROUPS.items():
        fp = ctmi_eval(Topt + h, Tmin, Topt, Tmax)
        fm = ctmi_eval(Topt - h, Tmin, Topt, Tmax)
        deriv = (fp - fm) / (2 * h)
        assert abs(deriv) < 1e-4, \
            f"{name}: f'(Topt={Topt}) = {deriv:.6e}, expected ~0"


def test_no_singularity_below_midpoint():
    """No singularity when Topt < (Tmin+Tmax)/2 (Rosso formula would fail)."""
    # These cases have Topt below the midpoint
    cases = {
        'P1_diatoms': (2.0, 15.0, 30.0),    # midpoint=16, Topt=15
        'P1_GOTM':    (0.0, 12.0, 25.0),    # midpoint=12.5, Topt=12
        'wide_low':   (0.0, 5.0, 40.0),     # midpoint=20, Topt=5
    }
    for name, (Tmin, Topt, Tmax) in cases.items():
        midpoint = (Tmin + Tmax) / 2.0
        assert Topt < midpoint, f"{name}: Topt not below midpoint"
        # Sweep at fine resolution through the entire range
        n_points = 10000
        for i in range(n_points + 1):
            T = Tmin + (Tmax - Tmin) * i / n_points
            val = ctmi_eval(T, Tmin, Topt, Tmax)
            assert math.isfinite(val), \
                f"{name}: f({T:.4f}) is not finite ({val})"
            assert 0.0 <= val <= 1.0, \
                f"{name}: f({T:.4f}) = {val}, outside [0, 1]"


def test_maximum_at_topt():
    """f(Topt) must be the global maximum in [Tmin, Tmax]."""
    for name, (Tmin, Topt, Tmax) in ALL_GROUPS.items():
        n_points = 1000
        max_val = 0.0
        max_T = Tmin
        for i in range(n_points + 1):
            T = Tmin + (Tmax - Tmin) * i / n_points
            val = ctmi_eval(T, Tmin, Topt, Tmax)
            if val > max_val:
                max_val = val
                max_T = T
        assert abs(max_T - Topt) < (Tmax - Tmin) / n_points * 2, \
            f"{name}: maximum at T={max_T:.3f}, expected Topt={Topt}"
        assert abs(max_val - 1.0) < 1e-6, \
            f"{name}: max value = {max_val}, expected 1.0"


def test_unimodal():
    """f(T) must be monotonically increasing below Topt and decreasing above."""
    for name, (Tmin, Topt, Tmax) in ALL_GROUPS.items():
        n_points = 500
        # Below Topt: monotonically increasing
        prev = 0.0
        for i in range(1, n_points):
            T = Tmin + (Topt - Tmin) * i / n_points
            val = ctmi_eval(T, Tmin, Topt, Tmax)
            assert val >= prev - 1e-10, \
                f"{name}: not monotonically increasing at T={T:.3f} (prev={prev:.6f}, val={val:.6f})"
            prev = val
        # Above Topt: monotonically decreasing
        prev = 1.0
        for i in range(1, n_points):
            T = Topt + (Tmax - Topt) * i / n_points
            val = ctmi_eval(T, Tmin, Topt, Tmax)
            assert val <= prev + 1e-10, \
                f"{name}: not monotonically decreasing at T={T:.3f} (prev={prev:.6f}, val={val:.6f})"
            prev = val


# ===================================================================
# Michaelis-Menten tests
# ===================================================================

def test_mm_at_ks_equals_half_vmax():
    """At S = Ks, V must equal Vmax/2."""
    for Ks in [0.1, 0.5, 1.0, 5.0, 10.0]:
        Vmax = 1.0
        V = michaelis_menten(Ks, Vmax, Ks)
        assert abs(V - 0.5) < 1e-12, \
            f"Ks={Ks}: V(Ks) = {V}, expected 0.5"


def test_mm_saturation():
    """V must approach Vmax as S >> Ks."""
    Vmax, Ks = 2.0, 1.0
    V = michaelis_menten(1000.0, Vmax, Ks)
    assert abs(V - Vmax) < 0.01, \
        f"V(S=1000) = {V}, expected ~{Vmax}"


def test_mm_zero_substrate():
    """V(0) must be 0."""
    V = michaelis_menten(0.0, 1.0, 1.0)
    assert V == 0.0, f"V(0) = {V}, expected 0.0"


def test_mm_positive_ks_prevents_divzero():
    """With Ks > 0, V is finite even at S = 0."""
    for Ks in [1e-6, 0.001, 0.5, 1.0]:
        V = michaelis_menten(0.0, 1.0, Ks)
        assert math.isfinite(V), f"Ks={Ks}: V(0) not finite"


# ===================================================================
# Expected values regression (from computed polynomial CTMI table)
# ===================================================================

def test_expected_values_p1():
    """P1 (Tmin=2, Topt=15, Tmax=30) polynomial CTMI expected values."""
    expected = {5: 0.42, 10: 0.86, 15: 1.00, 20: 0.88, 25: 0.53, 28: 0.23}
    Tmin, Topt, Tmax = 2.0, 15.0, 30.0
    for T, exp_val in expected.items():
        val = ctmi_eval(float(T), Tmin, Topt, Tmax)
        assert abs(val - exp_val) < 0.01, \
            f"P1 f({T}) = {val:.4f}, expected ~{exp_val}"


def test_expected_values_p2():
    """P2 (Tmin=5, Topt=20, Tmax=33) polynomial CTMI expected values."""
    expected = {10: 0.53, 15: 0.88, 20: 1.00, 25: 0.86, 30: 0.42}
    Tmin, Topt, Tmax = 5.0, 20.0, 33.0
    for T, exp_val in expected.items():
        val = ctmi_eval(float(T), Tmin, Topt, Tmax)
        assert abs(val - exp_val) < 0.01, \
            f"P2 f({T}) = {val:.4f}, expected ~{exp_val}"


# ===================================================================
# Runner
# ===================================================================

def run_all_tests():
    """Run all tests and report results."""
    test_functions = [v for k, v in globals().items()
                      if k.startswith('test_') and callable(v)]
    passed = 0
    failed = 0
    for test_fn in test_functions:
        try:
            test_fn()
            passed += 1
            print(f"  PASS: {test_fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL: {test_fn.__name__}: {e}")
        except Exception as e:
            failed += 1
            print(f"  ERROR: {test_fn.__name__}: {e}")

    print(f"\n{'=' * 50}")
    print(f"Total: {passed + failed}, Passed: {passed}, Failed: {failed}")
    print(f"{'=' * 50}")
    return failed == 0


if __name__ == '__main__':
    print("Polynomial CTMI & Michaelis-Menten regression tests")
    print("=" * 50)
    success = run_all_tests()
    sys.exit(0 if success else 1)
