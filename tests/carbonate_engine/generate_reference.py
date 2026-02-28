#!/usr/bin/env python3
"""
generate_reference.py

Generate reference test cases for carbonate_engine validation using PyCO2SYS.

This script is for DEVELOPER USE ONLY and is NOT used at runtime by ERSEM.
The Fortran test driver reads the generated CSV file directly.

Requirements:
    pip install PyCO2SYS

Usage:
    python generate_reference.py

Output:
    reference_cases.csv - CSV file with test cases and expected results

PyCO2SYS settings used (to match carbonate_engine with opt_k_carbonic=10):
    - K1K2: 10 (Lueker et al. 2000, total pH scale)
    - KSO4: 1 (Dickson 1990)
    - KF: 2 (Perez & Fraga 1987)
    - Total boron: 2 (Lee et al. 2010)
    - pH scale: 1 (Total)
"""

import sys

try:
    import PyCO2SYS as pyco2
except ImportError:
    print("ERROR: PyCO2SYS is required to generate reference cases.")
    print("Install with: pip install PyCO2SYS")
    sys.exit(1)

# Define test cases: (T, S, DIC_mol/kg, TA_mol/kg)
# Cover a range of typical marine conditions
TEST_CASES = [
    # Standard seawater conditions
    (25.0, 35.0, 0.002100, 0.002300),
    (25.0, 35.0, 0.002000, 0.002200),
    (25.0, 35.0, 0.002200, 0.002400),

    # Cooler water
    (15.0, 35.0, 0.002100, 0.002300),
    (15.0, 35.0, 0.002000, 0.002200),

    # Lower salinity (coastal)
    (10.0, 30.0, 0.002050, 0.002250),
    (20.0, 32.0, 0.002100, 0.002350),

    # Varying conditions
    (25.0, 34.0, 0.001900, 0.002100),
    (18.0, 33.0, 0.002150, 0.002380),
    (22.0, 36.0, 0.002080, 0.002280),

    # Edge cases (still valid marine)
    (5.0, 28.0, 0.002000, 0.002200),
    (30.0, 38.0, 0.002100, 0.002300),

    # Tokyo Bay-like conditions
    (20.0, 25.0, 0.001800, 0.002000),
    (25.0, 28.0, 0.001900, 0.002100),
    (15.0, 30.0, 0.002000, 0.002150),

    # Low salinity cases for estuarine conditions
    (20.0, 5.0, 0.002000, 0.002200),
    (20.0, 10.0, 0.002000, 0.002200),
    (20.0, 15.0, 0.002000, 0.002200),
    (25.0, 10.0, 0.002000, 0.002200),
    (15.0, 10.0, 0.002000, 0.002200),

    # Sub-zero temperature (polar waters)
    (-1.5, 34.0, 0.002200, 0.002350),
    (-1.0, 35.0, 0.002100, 0.002300),
]

# PyCO2SYS settings to match carbonate_engine (opt_k_carbonic=10 in PyCO2SYS)
PYCO2_KWARGS = {
    'opt_k_carbonic': 10,      # Lueker et al. 2000 (Total scale)
    'opt_k_bisulfate': 1,      # Dickson 1990
    'opt_k_fluoride': 2,       # Perez & Fraga 1987
    'opt_total_borate': 2,     # Lee et al. 2010
    'opt_pH_scale': 1,         # Total scale
}


def generate_reference_cases(output_file='reference_cases.csv'):
    """Generate reference cases and write to CSV file."""

    print("Generating reference cases using PyCO2SYS...")
    print(f"Settings: {PYCO2_KWARGS}")
    print()

    results = []

    for T, S, DIC_molkg, TA_molkg in TEST_CASES:
        # Convert to units expected by PyCO2SYS
        # PyCO2SYS expects DIC and TA in umol/kg
        DIC_umolkg = DIC_molkg * 1e6
        TA_umolkg = TA_molkg * 1e6

        # Run PyCO2SYS
        result = pyco2.sys(
            par1=TA_umolkg,
            par2=DIC_umolkg,
            par1_type=1,       # 1 = Total Alkalinity
            par2_type=2,       # 2 = DIC
            temperature=T,
            salinity=S,
            pressure=0,        # Surface pressure
            **PYCO2_KWARGS
        )

        # Extract results
        pH_total = result['pH_total']
        pCO2_uatm = result['pCO2']
        pCO2_atm = pCO2_uatm * 1e-6  # Convert uatm to atm

        results.append({
            'T': T,
            'S': S,
            'DIC_molkg': DIC_molkg,
            'TA_molkg': TA_molkg,
            'pH_total': pH_total,
            'pCO2_atm': pCO2_atm,
        })

        print(f"T={T:5.1f}C, S={S:5.1f}, DIC={DIC_molkg:.6f}, "
              f"TA={TA_molkg:.6f} -> pH={pH_total:.6f}, "
              f"pCO2={pCO2_atm:.6e} atm")

    # Write to CSV
    with open(output_file, 'w', newline='') as f:
        f.write("# Reference test cases for carbonate_engine\n")
        f.write("# Uses: K0 Weiss 1974, K1/K2 Lueker 2000, KB Dickson 1990, "
                "KW Millero 1995,\n")
        f.write("#       KS Dickson 1990, KF Perez & Fraga 1987, "
                "BT Lee 2010 (432.6*S/35)\n")
        f.write("# Note: Lueker 2000 K1/K2 valid for S=19-43. "
                "Low-S cases are extrapolations.\n")
        f.write("# Columns: T(degC), S(PSU), DIC(mol/kg), TA(mol/kg), "
                "pH_total, fCO2(atm)\n")

        for r in results:
            f.write(f"{r['T']},{r['S']},{r['DIC_molkg']:.6f},"
                    f"{r['TA_molkg']:.6f},{r['pH_total']:.10f},"
                    f"{r['pCO2_atm']:.10e}\n")

    print()
    print(f"Wrote {len(results)} test cases to {output_file}")

    return results


def verify_format(csv_file='reference_cases.csv'):
    """Verify the CSV file can be parsed correctly."""

    print(f"\nVerifying {csv_file}...")

    count = 0
    with open(csv_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split(',')
            if len(parts) != 6:
                print(f"ERROR: Expected 6 columns, got {len(parts)}: {line}")
                continue

            try:
                float(parts[0])
                float(parts[1])
                float(parts[2])
                float(parts[3])
                float(parts[4])
                float(parts[5])
                count += 1
            except ValueError as e:
                print(f"ERROR: Failed to parse line: {line}")
                print(f"       {e}")

    print(f"Successfully verified {count} test cases")
    return count > 0


if __name__ == '__main__':
    results = generate_reference_cases()
    verify_format()
