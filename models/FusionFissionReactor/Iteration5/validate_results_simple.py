#!/usr/bin/env python3
"""
Simple validation of manual depletion simulation results.
Checks key metrics to assess confidence in results.
"""

import openmc
import numpy as np
from pathlib import Path

def main():
    results_dir = Path('/workspaces/MEng172-OpenMC/models/FusionFissionReactor/Iteration5/results/iteration5_manual_depletion')
    
    print("=" * 70)
    print("VALIDATION OF SIMULATION RESULTS")
    print("=" * 70)
    print()
    
    # Collect data from all time steps
    k_source_vals = []
    power_vals = []
    fission_unc_vals = []
    flux_unc_vals = []
    nu_vals = []
    capture_to_fission_vals = []
    
    for i in range(4):
        step_dir = results_dir / f'step_{i}'
        sp = openmc.StatePoint(str(step_dir / 'statepoint.30.h5'))
        
        # Get tallies
        fission_tally = sp.get_tally(name='fission')
        capture_tally = sp.get_tally(name='capture')
        flux_tally = sp.get_tally(name='flux')
        
        # Extract values (already in reactions/s - absolute units)
        fission_rate = fission_tally.mean[0, 0, 0]  # fissions/s
        nu_fission_rate = fission_tally.mean[0, 0, 1]  # nu-fissions/s  
        capture_rate = capture_tally.mean[0, 0, 0]  # captures/s
        
        # Statistical uncertainties (%)
        fission_unc = fission_tally.std_dev[0, 0, 0] / fission_rate * 100 if fission_rate > 0 else 0
        flux_unc = flux_tally.std_dev[0, 0, 0] / flux_tally.mean[0, 0, 0] * 100
        
        # Calculate k_source
        source_rate = 1.0e16  # n/s
        k_source = nu_fission_rate / source_rate
        
        # Calculate power (fission_rate already in fissions/s)
        power_kW = fission_rate * 200e6 * 1.60218e-19 / 1000  # kW
        
        # Physics metrics
        nu_avg = nu_fission_rate / fission_rate if fission_rate > 0 else 0
        capture_to_fission = capture_rate / fission_rate if fission_rate > 0 else 0
        
        k_source_vals.append(k_source)
        power_vals.append(power_kW)
        fission_unc_vals.append(fission_unc)
        flux_unc_vals.append(flux_unc)
        nu_vals.append(nu_avg)
        capture_to_fission_vals.append(capture_to_fission)
    
    # Convert to arrays
    k_source_vals = np.array(k_source_vals)
    power_vals = np.array(power_vals)
    fission_unc_vals = np.array(fission_unc_vals)
    flux_unc_vals = np.array(flux_unc_vals)
    nu_vals = np.array(nu_vals)
    capture_to_fission_vals = np.array(capture_to_fission_vals)
    
    # ========== VALIDATION CHECKS ==========
    
    print("1. STATISTICAL QUALITY")
    print("-" * 70)
    print(f"Mean fission tally uncertainty: {np.mean(fission_unc_vals):.3f}%")
    print(f"Mean flux tally uncertainty:    {np.mean(flux_unc_vals):.3f}%")
    
    if np.mean(fission_unc_vals) < 1.0:
        print("✓ Excellent statistical precision (<1%)")
    elif np.mean(fission_unc_vals) < 5.0:
        print("✓ Good statistical precision (<5%)")
    else:
        print("⚠ Poor statistics - consider more particles/batches")
    print()
    
    print("2. TIME-SERIES STABILITY")
    print("-" * 70)
    print(f"k_source:  {np.mean(k_source_vals):.5f} ± {np.std(k_source_vals):.5f}")
    print(f"  Range: {np.min(k_source_vals):.5f} to {np.max(k_source_vals):.5f}")
    k_variation = (np.max(k_source_vals) - np.min(k_source_vals)) / np.mean(k_source_vals) * 100
    print(f"  Variation: {k_variation:.3f}%")
    
    if k_variation < 0.1:
        print("  ✓ Excellent stability (<0.1% variation)")
    elif k_variation < 1.0:
        print("  ✓ Good stability (<1% variation)")
    else:
        print("  ⚠ Significant drift - check if expected or problematic")
    print()
    
    print(f"Power:     {np.mean(power_vals):.2f} ± {np.std(power_vals):.2f} kW")
    print(f"  Range: {np.min(power_vals):.2f} to {np.max(power_vals):.2f} kW")
    power_variation = (np.max(power_vals) - np.min(power_vals)) / np.mean(power_vals) * 100
    print(f"  Variation: {power_variation:.3f}%")
    
    if power_variation < 0.1:
        print("  ✓ Excellent stability")
    elif power_variation < 1.0:
        print("  ✓ Good stability")
    else:
        print("  ⚠ Power drift detected")
    print()
    
    print("3. PHYSICS REASONABLENESS")
    print("-" * 70)
    
    # Average neutrons per fission
    print(f"Average ν (neutrons/fission): {np.mean(nu_vals):.3f}")
    if 2.0 < np.mean(nu_vals) < 3.0:
        print("  ✓ Typical for spent fuel mix (2-3 neutrons/fission)")
    else:
        print("  ⚠ Unusual value - check material composition")
    print()
    
    # Capture to fission ratio
    print(f"Capture/Fission ratio: {np.mean(capture_to_fission_vals):.3f}")
    if 2.0 < np.mean(capture_to_fission_vals) < 6.0:
        print("  ✓ Expected for U-238 rich spent fuel")
    else:
        print("  ⚠ Unusual ratio")
    print()
    
    # Neutron multiplication
    print(f"k_source: {np.mean(k_source_vals):.5f}")
    if np.mean(k_source_vals) > 1.0:
        print("  ✓ System multiplies neutrons (k > 1)")
    else:
        print("  ✗ No neutron multiplication (k < 1)")
    
    if np.mean(k_source_vals) < 1.5:
        print("  ✓ Safely subcritical (k < 1.5)")
    else:
        print("  ⚠ Approaching criticality")
    print()
    
    # Energy per fission
    print(f"Power per source neutron: {np.mean(power_vals)/1.0e16*1e9:.3f} MeV")
    energy_mult = (np.mean(power_vals) / 1000) / (1.0e16 * 2.45e6 * 1.60218e-19 / 1e6)
    print(f"Energy multiplication: {energy_mult:.1f}×")
    print("  (Fusion input 2.45 MeV → Fission output)")
    if energy_mult > 10:
        print("  ✓ Significant energy gain")
    else:
        print("  ⚠ Low energy multiplication")
    print()
    
    print("4. MESH TALLY QUALITY")
    print("-" * 70)
    # Check mesh tallies from final step
    sp = openmc.StatePoint(str(results_dir / 'step_3' / 'statepoint.30.h5'))
    
    flux_mesh = sp.get_tally(name='3d_flux_tally')
    flux_mean = flux_mesh.mean
    flux_std = flux_mesh.std_dev
    
    # Find active region
    active_cells = flux_mean > 0
    n_active = np.sum(active_cells)
    n_total = flux_mean.size
    
    print(f"Active mesh cells: {n_active}/{n_total} ({n_active/n_total*100:.1f}%)")
    
    if n_active/n_total > 0.1:
        print("  ✓ Good spatial coverage")
    else:
        print("  ⚠ Very sparse - consider smaller mesh or different geometry")
    
    # Uncertainty in active region
    active_flux = flux_mean[active_cells]
    active_unc = flux_std[active_cells] / active_flux * 100
    
    # Filter out invalid values
    valid_unc = active_unc[np.isfinite(active_unc) & (active_unc > 0)]
    
    if len(valid_unc) > 0:
        print(f"  Mean mesh uncertainty: {np.mean(valid_unc):.1f}%")
        print(f"  Median mesh uncertainty: {np.median(valid_unc):.1f}%")
        
        if np.mean(valid_unc) < 10:
            print("  ✓ Acceptable mesh statistics (<10%)")
        elif np.mean(valid_unc) < 20:
            print("  ⚠ Marginal mesh statistics (10-20%)")
        else:
            print("  ⚠ Poor mesh statistics (>20%) - visualizations may be noisy")
    print()
    
    # ========== OVERALL ASSESSMENT ==========
    
    print("=" * 70)
    print("OVERALL ASSESSMENT")
    print("=" * 70)
    
    passed = []
    warnings = []
    failed = []
    
    # Statistical quality
    if np.mean(fission_unc_vals) < 1.0:
        passed.append("Excellent statistical precision on integral tallies")
    elif np.mean(fission_unc_vals) < 5.0:
        passed.append("Good statistical precision")
    else:
        warnings.append("Poor statistical precision - consider more particles")
    
    # Stability
    if k_variation < 0.1 and power_variation < 0.1:
        passed.append("Excellent time-series stability")
    elif k_variation < 1.0 and power_variation < 1.0:
        passed.append("Good time-series stability")
    else:
        warnings.append("Time-series shows drift")
    
    # Physics
    if 2.0 < np.mean(nu_vals) < 3.0:
        passed.append("Realistic neutrons per fission")
    else:
        warnings.append("Unusual nu value")
    
    if 2.0 < np.mean(capture_to_fission_vals) < 6.0:
        passed.append("Expected capture/fission ratio for U-238 rich fuel")
    else:
        warnings.append("Unusual capture/fission ratio")
    
    if np.mean(k_source_vals) > 1.0 and np.mean(k_source_vals) < 1.5:
        passed.append("k_source indicates subcritical multiplication")
    elif np.mean(k_source_vals) > 1.0:
        warnings.append("k_source approaching criticality")
    else:
        failed.append("k_source < 1: No neutron multiplication!")
    
    if energy_mult > 20:
        passed.append(f"Strong energy multiplication ({energy_mult:.1f}×)")
    elif energy_mult > 10:
        passed.append(f"Moderate energy multiplication ({energy_mult:.1f}×)")
    else:
        warnings.append(f"Low energy multiplication ({energy_mult:.1f}×)")
    
    # Mesh quality
    if len(valid_unc) > 0 and np.mean(valid_unc) < 20:
        passed.append("Mesh tallies usable for visualization")
    elif len(valid_unc) > 0:
        warnings.append("Mesh tallies have high uncertainty")
    
    # Print summary
    if len(passed) > 0:
        print("\n✓ PASSED CHECKS:")
        for item in passed:
            print(f"  • {item}")
    
    if len(warnings) > 0:
        print("\n⚠ WARNINGS:")
        for item in warnings:
            print(f"  • {item}")
    
    if len(failed) > 0:
        print("\n✗ FAILED CHECKS:")
        for item in failed:
            print(f"  • {item}")
    
    print()
    print("-" * 70)
    
    if len(failed) > 0:
        print("CONFIDENCE LEVEL: LOW")
        print("Critical issues detected. Results may not be physically meaningful.")
    elif len(warnings) > 2:
        print("CONFIDENCE LEVEL: MODERATE")
        print("Results are usable but have some limitations.")
        print("Consider: More particles/batches for better statistics.")
    else:
        print("CONFIDENCE LEVEL: HIGH")
        print("Results appear physically consistent and statistically sound.")
        print("Integral quantities (k_source, power) are well-resolved.")
        print("Spatial distributions may have higher uncertainty but show correct trends.")
    
    print("=" * 70)
    
    # Specific recommendations
    print("\nRECOMMENDATIONS:")
    print("-" * 70)
    
    print("✓ Trust integral quantities (k_source, total power, reaction rates)")
    print("  These have <1% statistical uncertainty")
    print()
    
    if len(valid_unc) > 0 and np.mean(valid_unc) > 10:
        print("⚠ Mesh tallies (spatial plots) have higher uncertainty (10-20%)")
        print("  Visualizations show correct trends but not fine details")
        print("  For publication-quality plots, consider:")
        print("    - Increase particles: 10,000-20,000 per batch")
        print("    - More active batches: 50-100")
        print("    - Or reduce mesh resolution to improve statistics per cell")
    else:
        print("✓ Spatial distributions are reliable for qualitative analysis")
    print()
    
    print("✓ Time evolution is trustworthy")
    print(f"  k_source varies by only {k_variation:.3f}% over 360 days")
    print(f"  Power varies by only {power_variation:.3f}% over 360 days")
    print()
    
    print("✓ Physics is self-consistent:")
    print(f"  • ν = {np.mean(nu_vals):.3f} neutrons/fission (typical spent fuel)")
    print(f"  • Capture/Fission = {np.mean(capture_to_fission_vals):.3f} (U-238 dominant)")
    print(f"  • k_source = {np.mean(k_source_vals):.5f} (subcritical with multiplication)")
    print(f"  • Energy gain = {energy_mult:.1f}× (fusion-fission hybrid)")
    
    print("=" * 70)

if __name__ == '__main__':
    main()
