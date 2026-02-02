#!/usr/bin/env python3
"""
Validation script for manual depletion results.
Checks statistical quality, physics consistency, and convergence.
"""

import openmc
import numpy as np
import os
from pathlib import Path

def check_statistical_uncertainties(step_dir):
    """Extract statistical uncertainties from statepoint."""
    sp_path = step_dir / 'statepoint.30.h5'
    sp = openmc.StatePoint(str(sp_path))
    
    results = {}
    
    # Get tallies and their uncertainties
    for tally in sp.tallies.values():
        if tally.name:
            mean = tally.mean.flatten()
            std_dev = tally.std_dev.flatten()
            # Calculate relative uncertainty (%)
            rel_unc = np.where(mean > 0, (std_dev / mean) * 100, 0)
            
            # Get statistics for non-zero values
            nonzero = rel_unc[rel_unc > 0]
            if len(nonzero) > 0:
                results[tally.name] = {
                    'mean_uncertainty': np.mean(nonzero),
                    'max_uncertainty': np.max(nonzero),
                    'median_uncertainty': np.median(nonzero)
                }
    
    return results, sp

def check_energy_balance(sp):
    """Verify energy balance from heating tally."""
    heating_tally = sp.get_tally(name='3d_heating_tally')
    heating_rate = heating_tally.mean.sum()  # Sum all mesh cells, eV/source
    
    # Convert to power
    source_rate = 1.0e16  # neutrons/s
    heating_power = heating_rate * source_rate * 1.60218e-19 / 1000  # kW
    
    # Get fission power independently
    fission_tally = sp.get_tally(name='fission')
    fission_rate = fission_tally.get_slice(scores=['fission']).mean[0][0][0] * source_rate  # fissions/s
    fission_power = fission_rate * 200e6 * 1.60218e-19 / 1000  # kW
    
    return {
        'heating_power_kW': heating_power,
        'fission_power_kW': fission_power,
        'energy_discrepancy_%': abs(heating_power - fission_power) / fission_power * 100
    }

def check_neutron_balance(sp):
    """Verify neutron balance makes physical sense."""
    fission_tally = sp.get_tally(name='fission')
    capture_tally = sp.get_tally(name='capture')
    
    source_rate = 1.0e16
    
    # Get fission and nu-fission from the same tally (has both scores)
    fission_rate = fission_tally.get_slice(scores=['fission']).mean[0][0][0] * source_rate
    nu_fission_rate = fission_tally.get_slice(scores=['nu-fission']).mean[0][0][0] * source_rate
    capture_rate = capture_tally.mean[0][0][0] * source_rate
    absorption_rate = fission_rate + capture_rate
    
    # Average neutrons per fission
    nu_avg = nu_fission_rate / fission_rate if fission_rate > 0 else 0
    
    # k_source
    k_source = nu_fission_rate / (source_rate + absorption_rate)
    
    # Capture rate (absorption - fission)
    capture_rate = absorption_rate - fission_rate
    
    return {
        'nu_avg': nu_avg,
        'k_source': k_source,
        'fission_rate': fission_rate,
        'capture_rate': capture_rate,
        'absorption_rate': absorption_rate,
        'capture_to_fission_ratio': capture_rate / fission_rate if fission_rate > 0 else 0
    }

def check_mesh_convergence(step_dir):
    """Check if mesh tally values are well-resolved."""
    sp_path = step_dir / 'statepoint.30.h5'
    sp = openmc.StatePoint(str(sp_path))
    
    # Check 3D flux tally
    flux_tally = sp.get_tally(name='3d_flux_tally')
    flux_mean = flux_tally.mean
    flux_std = flux_tally.std_dev
    
    # Find fuel region (non-zero flux)
    nonzero_mask = flux_mean > 0
    if np.any(nonzero_mask):
        flux_values = flux_mean[nonzero_mask]
        flux_unc = flux_std[nonzero_mask]
        rel_unc = (flux_unc / flux_values) * 100
        
        return {
            'active_mesh_cells': np.sum(nonzero_mask),
            'total_mesh_cells': flux_mean.size,
            'fraction_active': np.sum(nonzero_mask) / flux_mean.size * 100,
            'mean_flux_uncertainty_%': np.mean(rel_unc),
            'max_flux_uncertainty_%': np.max(rel_unc)
        }
    
    return None

def check_time_consistency(base_dir):
    """Check consistency across time steps."""
    results_dir = Path(base_dir) / 'results' / 'iteration5_manual_depletion'
    
    k_source_values = []
    power_values = []
    uncertainties = []
    
    for i in range(4):
        step_dir = results_dir / f'step_{i}'
        if step_dir.exists():
            sp_path = step_dir / 'statepoint.30.h5'
            sp = openmc.StatePoint(str(sp_path))
            
            # Get k_source
            fission = sp.get_tally(name='fission')
            capture = sp.get_tally(name='capture')
            
            source_rate = 1.0e16
            nu_fission_rate = fission.get_slice(scores=['nu-fission']).mean[0][0][0] * source_rate
            fission_rate = fission.get_slice(scores=['fission']).mean[0][0][0] * source_rate
            capture_rate = capture.mean[0][0][0] * source_rate
            absorption_rate = fission_rate + capture_rate
            k_source = nu_fission_rate / (source_rate + absorption_rate)
            
            k_source_values.append(k_source)
            
            # Get power
            fission = sp.get_tally(name='fission')
            fission_rate = fission.get_slice(scores=['fission']).mean[0][0][0] * source_rate
            power = fission_rate * 200e6 * 1.60218e-19 / 1000  # kW
            power_values.append(power)
            
            # Get representative uncertainty
            flux = sp.get_tally(name='flux')
            flux_unc = flux.std_dev[0][0][0] / flux.mean[0][0][0] * 100
            uncertainties.append(flux_unc)
    
    k_source_values = np.array(k_source_values)
    power_values = np.array(power_values)
    uncertainties = np.array(uncertainties)
    
    return {
        'k_source_mean': np.mean(k_source_values),
        'k_source_std': np.std(k_source_values),
        'k_source_variation_%': (np.max(k_source_values) - np.min(k_source_values)) / np.mean(k_source_values) * 100,
        'power_mean_kW': np.mean(power_values),
        'power_std_kW': np.std(power_values),
        'power_variation_%': (np.max(power_values) - np.min(power_values)) / np.mean(power_values) * 100,
        'mean_statistical_uncertainty_%': np.mean(uncertainties)
    }

def validate_physics_reasonableness(neutron_balance):
    """Check if results are physically reasonable."""
    checks = []
    
    # Nu should be between 2.2-2.9 for typical spent fuel
    if 2.0 < neutron_balance['nu_avg'] < 3.0:
        checks.append(("✓", f"Neutrons per fission: {neutron_balance['nu_avg']:.3f} (reasonable for spent fuel)"))
    else:
        checks.append(("✗", f"Neutrons per fission: {neutron_balance['nu_avg']:.3f} (unusual!)"))
    
    # Capture-to-fission ratio should be 2-5 for spent fuel with U-238
    ratio = neutron_balance['capture_to_fission_ratio']
    if 2.0 < ratio < 6.0:
        checks.append(("✓", f"Capture/fission ratio: {ratio:.2f} (typical for U-238 rich fuel)"))
    else:
        checks.append(("✗", f"Capture/fission ratio: {ratio:.2f} (unusual!)"))
    
    # k_source should be > 1 for multiplication
    if neutron_balance['k_source'] > 1.0:
        checks.append(("✓", f"k_source: {neutron_balance['k_source']:.5f} (system multiplies neutrons)"))
    else:
        checks.append(("✗", f"k_source: {neutron_balance['k_source']:.5f} (no multiplication!)"))
    
    # k_source should be < 1.5 for subcritical system
    if neutron_balance['k_source'] < 1.5:
        checks.append(("✓", f"k_source < 1.5 (safely subcritical)"))
    else:
        checks.append(("⚠", f"k_source > 1.5 (approaching criticality)"))
    
    return checks

def main():
    base_dir = Path('/workspaces/MEng172-OpenMC/models/FusionFissionReactor/Iteration5')
    results_dir = base_dir / 'results' / 'iteration5_manual_depletion'
    
    print("=" * 70)
    print("VALIDATION OF SIMULATION RESULTS")
    print("=" * 70)
    print()
    
    # Check time consistency across all steps
    print("1. TIME-SERIES CONSISTENCY")
    print("-" * 70)
    consistency = check_time_consistency(base_dir)
    print(f"k_source:  {consistency['k_source_mean']:.5f} ± {consistency['k_source_std']:.5f}")
    print(f"  Variation: {consistency['k_source_variation_%']:.3f}%")
    if consistency['k_source_variation_%'] < 1.0:
        print("  ✓ k_source is stable over time")
    else:
        print("  ⚠ k_source shows significant variation")
    print()
    print(f"Power:     {consistency['power_mean_kW']:.2f} ± {consistency['power_std_kW']:.2f} kW")
    print(f"  Variation: {consistency['power_variation_%']:.3f}%")
    if consistency['power_variation_%'] < 1.0:
        print("  ✓ Power is stable over time")
    else:
        print("  ⚠ Power shows significant variation")
    print()
    print(f"Mean statistical uncertainty: {consistency['mean_statistical_uncertainty_%']:.3f}%")
    if consistency['mean_statistical_uncertainty_%'] < 1.0:
        print("  ✓ Excellent statistical precision")
    elif consistency['mean_statistical_uncertainty_%'] < 5.0:
        print("  ✓ Good statistical precision")
    else:
        print("  ⚠ Consider increasing particle count")
    print()
    
    # Analyze final step in detail
    final_step = results_dir / 'step_3'
    print("2. STATISTICAL UNCERTAINTY (Final Step)")
    print("-" * 70)
    uncertainties, sp = check_statistical_uncertainties(final_step)
    for tally_name, unc_data in uncertainties.items():
        print(f"{tally_name}:")
        print(f"  Mean uncertainty: {unc_data['mean_uncertainty']:.3f}%")
        print(f"  Max uncertainty:  {unc_data['max_uncertainty']:.3f}%")
        if unc_data['mean_uncertainty'] < 1.0:
            print("  ✓ Excellent precision")
        elif unc_data['mean_uncertainty'] < 5.0:
            print("  ✓ Good precision")
        else:
            print("  ⚠ Consider more particles")
        print()
    
    # Energy balance
    print("3. ENERGY BALANCE")
    print("-" * 70)
    energy = check_energy_balance(sp)
    print(f"Heating tally power:  {energy['heating_power_kW']:.2f} kW")
    print(f"Fission power (200 MeV/fission): {energy['fission_power_kW']:.2f} kW")
    print(f"Discrepancy: {energy['energy_discrepancy_%']:.3f}%")
    if energy['energy_discrepancy_%'] < 5.0:
        print("✓ Energy balance is consistent")
    else:
        print("⚠ Energy discrepancy may indicate missing physics")
    print()
    
    # Neutron balance
    print("4. NEUTRON BALANCE")
    print("-" * 70)
    neutron = check_neutron_balance(sp)
    print(f"Average ν (neutrons/fission): {neutron['nu_avg']:.3f}")
    print(f"k_source: {neutron['k_source']:.5f}")
    print(f"Fission rate:    {neutron['fission_rate']:.4e} /s")
    print(f"Capture rate:    {neutron['capture_rate']:.4e} /s")
    print(f"Absorption rate: {neutron['absorption_rate']:.4e} /s")
    print(f"Capture/Fission: {neutron['capture_to_fission_ratio']:.3f}")
    print()
    
    # Physics checks
    print("5. PHYSICS REASONABLENESS")
    print("-" * 70)
    physics_checks = validate_physics_reasonableness(neutron)
    for status, message in physics_checks:
        print(f"{status} {message}")
    print()
    
    # Mesh convergence
    print("6. SPATIAL MESH RESOLUTION")
    print("-" * 70)
    mesh_data = check_mesh_convergence(final_step)
    if mesh_data:
        print(f"Active mesh cells: {mesh_data['active_mesh_cells']}/{mesh_data['total_mesh_cells']} ({mesh_data['fraction_active']:.1f}%)")
        print(f"Mean flux uncertainty: {mesh_data['mean_flux_uncertainty_%']:.3f}%")
        print(f"Max flux uncertainty:  {mesh_data['max_flux_uncertainty_%']:.3f}%")
        
        if mesh_data['fraction_active'] > 5:
            print("✓ Sufficient spatial resolution")
        else:
            print("⚠ Very few active mesh cells")
        
        if mesh_data['mean_flux_uncertainty_%'] < 10:
            print("✓ Mesh tallies have good statistics")
        else:
            print("⚠ Mesh tallies need more particles")
    print()
    
    # Overall assessment
    print("=" * 70)
    print("OVERALL ASSESSMENT")
    print("=" * 70)
    
    issues = []
    if consistency['k_source_variation_%'] > 1.0:
        issues.append("k_source shows >1% variation")
    if consistency['mean_statistical_uncertainty_%'] > 5.0:
        issues.append("Statistical uncertainty >5%")
    if energy['energy_discrepancy_%'] > 5.0:
        issues.append("Energy balance discrepancy >5%")
    if not (2.0 < neutron['nu_avg'] < 3.0):
        issues.append("Unusual nu value")
    if mesh_data and mesh_data['mean_flux_uncertainty_%'] > 10:
        issues.append("Mesh tally uncertainty >10%")
    
    if len(issues) == 0:
        print("✓ All validation checks passed!")
        print("✓ Results appear to be high-quality and physically consistent")
        print("✓ Statistical uncertainties are acceptable")
        print("✓ Energy and neutron balances are consistent")
        print()
        print("CONFIDENCE LEVEL: HIGH")
    else:
        print("⚠ Some concerns identified:")
        for issue in issues:
            print(f"  - {issue}")
        print()
        print("CONFIDENCE LEVEL: MODERATE")
        print("Consider: More particles, more batches, or longer active cycles")
    
    print("=" * 70)

if __name__ == '__main__':
    main()
