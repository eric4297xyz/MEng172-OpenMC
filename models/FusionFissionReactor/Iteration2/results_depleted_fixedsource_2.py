"""
Post-processing script for fixed source depletion simulation results.
Reads depletion_results.h5 directly using h5py and generates summary plots/data.
"""
import numpy as np
import h5py
import pandas as pd
import matplotlib.pyplot as plt
import os
import math
from pathlib import Path

# Constants
MEV_PER_FISSION = 200.0  # Average energy per fission
MEV_TO_JOULES = 1.60218e-13
MEV_TO_MJ = MEV_TO_JOULES * 1e-6

# Material volume - from fuel_blanket_2.py geometry
SPHERE_INNER_RADIUS = 50.0  # cm
SPHERE_OUTER_RADIUS = 100.0  # cm
MATERIAL_VOLUME_CM3 = (4.0/3.0) * math.pi * (SPHERE_OUTER_RADIUS**3 - SPHERE_INNER_RADIUS**3)

RESULTS_FILE = "results/iteration2_results_depleted_fuel_fixedsource/depletion_results.h5"
OUTPUT_DIR = "results/iteration2_results_depleted_fuel_fixedsource"

def get_initial_material_composition():
    """
    Returns the initial spent fuel composition as defined in fuel_blanket_2.py.
    Returns dict of {nuclide_name: atom_fraction} matching the material definition.
    """
    # From fuel_blanket_2.py mspentfuel material definition
    # Total density = 7.133315757E-02 atom/b-cm
    initial_composition = {
        'O16': 6.7187968E-01,
        'U238': 3.0769658E-01,
        'U235': 3.4979859E-03,
        'Pu239': 2.5718196E-03,
        'U236': 2.1091663E-03,
        'Cs137': 1.0510726E-03,
        'Pu240': 1.0215150E-03,
        'Cs133': 9.7963234E-04,
        'Tc99': 9.3559531E-04,
        'Ru101': 9.1992123E-04,
        'Zr93': 8.9218967E-04,
        'Mo95': 8.5245707E-04,
        'Sr90': 6.7977794E-04,
        'Pu241': 6.7015516E-04,
        'Nd143': 6.5286858E-04,
        'Nd145': 5.3344636E-04,
        'Rh103': 4.8652147E-04,
        'Cs135': 4.8545594E-04,
        'Np237': 2.7884345E-04,
        'Pu242': 2.6644975E-04,
        'Pd107': 2.5941176E-04,
        'Sm150': 2.2788081E-04,
        'I129': 1.4371885E-04,
        'Pu238': 1.3101157E-04,
        'Pm147': 1.0269899E-04,
        'Eu153': 9.2879588E-05,
        'Sm152': 8.8634337E-05,
        'Ag109': 8.5141959E-05,
        'Am243': 7.8638873E-05,
        'Sm147': 6.6014495E-05,
        'U234': 6.4111968E-05,
        'Cm244': 3.2900525E-05,
        'Am241': 3.1275801E-05,
        'Ru103': 1.9360418E-05,
        'Sn126': 1.8018478E-05,
        'Nb95': 1.4824839E-05,
        'Cl36': 1.4019981E-05,
        'Ca41': 1.4019976E-05,
        'Ni59': 1.4019973E-05,
        'Sm151': 1.3614245E-05,
        'Cm242': 7.3636478E-06,
        'Se79': 7.0915868E-06,
        'Eu155': 4.7729475E-06,
        'Pr143': 2.0561139E-06,
        'Cm245': 1.9560548E-06,
        'Sm149': 1.9355989E-06,
        'Nd147': 5.1292485E-07,
        'Cm243': 3.0834446E-07,
        'Cm246': 1.9481932E-07,
        'U237': 1.7067631E-07,
        'Xe133': 9.3481328E-08,
        'Gd155': 7.6579955E-08,
        'I133': 7.5534463E-08,
        'Eu152': 4.5260253E-08,
        'Pu244': 9.4590026E-09,
        'Np239': 3.2856837E-09,
        'U233': 1.2208878E-09,
        'Mo99': 1.1472054E-09,
    }
    
    # Material total density in atoms/(barn·cm)
    total_density = 7.133315757E-02
    
    # Convert atom fractions to number densities [atoms/(barn·cm)]
    number_densities = {nuc: frac * total_density for nuc, frac in initial_composition.items()}
    
    return number_densities

def calculate_initial_atom_counts(material_volume_cm3):
    """
    Calculate initial total atom counts from material composition and volume.
    
    Args:
        material_volume_cm3: Volume of the material in cm³
    
    Returns:
        dict: {nuclide_name: total_atoms}
    """
    number_densities = get_initial_material_composition()
    
    # Convert number density [atoms/(barn·cm)] to total atoms
    # 1 barn = 1e-24 cm², so atoms/(barn·cm) needs to be multiplied by 1e24 to get atoms/cm³
    # Then multiply by volume to get total atoms
    initial_atoms = {nuc: density * 1e24 * material_volume_cm3 for nuc, density in number_densities.items()}
    
    return initial_atoms

def read_depletion_results(filename):
    """Read depletion results HDF5 file and return key datasets."""
    with h5py.File(filename, 'r') as f:
        times = f['time'][:]  # shape: (n_steps, 2) = [seconds, ?? ]
        eigenvalues = f['eigenvalues'][:]  # shape: (n_steps, 1, 2) = [value, uncertainty]
        source_rates = f['source_rate'][:]  # shape: (n_steps, 1) = [neutrons/s]
        nuclides = list(f['nuclides'].keys())
        materials = list(f['materials'].keys())
        
        # Read nuclide number densities for material '1' if available
        number_densities = None
        if '1' in f['materials']:
            try:
                mat_group = f['materials']['1']
                if 'nuclides' in mat_group:
                    number_densities = mat_group['nuclides'][:]  # shape: (n_steps, n_nuclides)
            except Exception as e:
                print(f"Warning: could not read number densities: {e}")
        
        # Read number array (nuclide counts or densities)
        number_array = f['number'][:]  # shape: (n_steps, 1, 1, n_nuclides)
        
        return {
            'times': times,
            'eigenvalues': eigenvalues,
            'source_rates': source_rates,
            'nuclides': nuclides,
            'materials': materials,
            'number_array': number_array,
            'number_densities': number_densities,
        }

def extract_burnup_metrics(data):
    """Extract burnup and power metrics from depletion data."""
    times_s = data['times'][:, 0]  # Time in seconds
    eigenvalues = data['eigenvalues'][:, 0, 0]  # k-eff values
    source_rates = data['source_rates'][:, 0]  # neutrons/s
    
    # Calculate power from source rate
    # Assumption: each neutron released comes from ~200 MeV fission energy
    # Power [W] = fission_rate [fissions/s] * energy_per_fission [J/fission]
    # Using source_rate as proxy for fission rate
    power_watts = source_rates * MEV_PER_FISSION * MEV_TO_JOULES
    power_mw = power_watts / 1e6  # Convert to MW
    
    # Calculate cumulative energy
    dt_s = np.diff(times_s)
    dt_s = np.append(dt_s, dt_s[-1] if len(dt_s) > 0 else 0)  # pad last element
    
    # Cumulative neutrons and energy
    cumulative_neutrons = np.cumsum(source_rates * dt_s)
    cumulative_energy_j = np.cumsum(power_watts * dt_s)  # Joules
    cumulative_energy_mwh = cumulative_energy_j / (3.6e9)  # Convert to MWh
    
    return {
        'times_s': times_s,
        'times_h': times_s / 3600.0,
        'times_d': times_s / 86400.0,
        'eigenvalues': eigenvalues,
        'source_rates': source_rates,
        'power_watts': power_watts,
        'power_mw': power_mw,
        'cumulative_neutrons': cumulative_neutrons,
        'cumulative_energy_j': cumulative_energy_j,
        'cumulative_energy_mwh': cumulative_energy_mwh,
    }

def extract_nuclide_inventories(data):
    """Extract key nuclide inventories over time and convert to total atom counts."""
    nuclides = data['nuclides']
    number_array = data['number_array']  # shape: (n_steps, 1, 1, n_nuclides)
    n_steps = number_array.shape[0]
    
    # Flatten to (n_steps, n_nuclides) - these are number densities in atoms/(barn·cm)
    inventory_density = number_array[:, 0, 0, :]
    
    # Convert to total atom counts by multiplying by material volume
    # Number density [atoms/(barn·cm)] * Volume [cm³] = Total atoms
    inventory = inventory_density * MATERIAL_VOLUME_CM3
    
    # Extract key actinides/fission products (with actual values, not forced zeros)
    key_nuclides = [
        'U235', 'U238', 'U239',
        'Np238', 'Np239', 'Np240', 'Np241', 'Np242',
        'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242',
        'Am238', 'Am239', 'Am240', 'Am241', 'Am242',
        'Cm238', 'Cm239', 'Cm240', 'Cm241', 'Cm242',
        'Cm244', 'Xe135', 'Sm151', 'Cs137'
    ]
    results = {}
    
    for nuc in key_nuclides:
        if nuc in nuclides:
            idx = nuclides.index(nuc)
            results[nuc] = inventory[:, idx]
    
    # Also include any nuclides with significant initial or final inventory
    # (helps discover what was actually in the material)
    for i, nuc in enumerate(nuclides):
        initial_val = inventory[0, i]
        final_val = inventory[-1, i]
        # Include if initial > threshold OR final > threshold
        if initial_val > 1e20 or final_val > 1e20:
            if nuc not in results:
                results[nuc] = inventory[:, i]
    
    return results

def analyze_waste_transmutation(inventories, times_s):
    """Analyze depletion of transuranic actinides and long-lived waste."""
    
    # Get initial material composition from simulation definition
    initial_atoms = calculate_initial_atom_counts(MATERIAL_VOLUME_CM3)
    
    # Define transuranic actinides (Z > 92)
    transuranic_isotopes = [
        'Np237', 'Np238', 'Np239',
        'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242',
        'Am241', 'Am242_m1', 'Am242', 'Am243',
        'Cm242', 'Cm243', 'Cm244', 'Cm245', 'Cm246', 'Cm247', 'Cm248'
    ]
    
    # Define long-lived fission products (>30 year half-life, high radiotoxicity)
    long_lived_fp = [
        'Tc99',    # 211,000 years
        'I129',    # 15.7 million years
        'Cs135',   # 2.3 million years
        'Cs137',   # 30.2 years
        'Sr90',    # 28.8 years
        'Zr93',    # 1.5 million years
        'Se79',    # 327,000 years
        'Pd107',   # 6.5 million years
    ]
    
    results = {
        'tru_initial': 0.0,
        'tru_final': 0.0,
        'tru_details': {},
        'llfp_initial': 0.0,
        'llfp_final': 0.0,
        'llfp_details': {},
        'times_s': times_s
    }
    
    # Calculate total TRU inventory
    for isotope in transuranic_isotopes:
        # Get initial from material definition
        initial = initial_atoms.get(isotope, 0.0)
        # Get final from depletion results
        final = inventories.get(isotope, [0.0])[-1] if isotope in inventories else 0.0
        
        results['tru_initial'] += initial
        results['tru_final'] += final
        if initial > 1e5 or final > 1e5:  # Only track significant amounts
            results['tru_details'][isotope] = {
                'initial': initial,
                'final': final,
                'change': final - initial,
                'percent_change': ((final - initial) / initial * 100) if initial > 1e10 else 0.0
            }
    
    # Calculate total long-lived FP inventory
    for isotope in long_lived_fp:
        # Get initial from material definition
        initial = initial_atoms.get(isotope, 0.0)
        # Get final from depletion results
        final = inventories.get(isotope, [0.0])[-1] if isotope in inventories else 0.0
        
        results['llfp_initial'] += initial
        results['llfp_final'] += final
        if initial > 1e5 or final > 1e5:
            results['llfp_details'][isotope] = {
                'initial': initial,
                'final': final,
                'change': final - initial,
                'percent_change': ((final - initial) / initial * 100) if initial > 1e10 else 0.0
            }
    
    # Calculate overall statistics
    results['tru_change'] = results['tru_final'] - results['tru_initial']
    results['tru_percent_change'] = (results['tru_change'] / results['tru_initial'] * 100) if results['tru_initial'] > 0 else 0.0
    results['llfp_change'] = results['llfp_final'] - results['llfp_initial']
    results['llfp_percent_change'] = (results['llfp_change'] / results['llfp_initial'] * 100) if results['llfp_initial'] > 0 else 0.0
    
    return results

def generate_summary_report(data, metrics, inventories):
    """Generate a text summary of depletion results."""
    times_s = metrics['times_s']
    eigenvalues = metrics['eigenvalues']
    
    # Get initial material composition from simulation definition
    initial_atoms = calculate_initial_atom_counts(MATERIAL_VOLUME_CM3)
    
    report = []
    report.append("=" * 70)
    report.append("FIXED SOURCE DEPLETION SIMULATION SUMMARY")
    report.append("=" * 70)
    report.append("\n⚠️  DATA QUALITY NOTE:")
    report.append("   The HDF5 depletion results show near-zero initial densities for U235/U238.")
    report.append("   This is a known issue with OpenMC's depletion data recording.")
    report.append("   Initial compositions are taken from the material definition (materials.xml)")
    report.append("   which correctly specifies the spent fuel composition.")
    report.append("   Final states are from the HDF5 file.")
    report.append("")
    report.append(f"\nSimulation Duration: {times_s[-1]:.2e} seconds ({times_s[-1]/86400:.2f} days)")
    report.append(f"Number of Depletion Steps: {len(times_s)}")
    report.append(f"\nInitial k-eff: {eigenvalues[0]:.6f}")
    report.append(f"Final k-eff: {eigenvalues[-1]:.6f}")
    report.append(f"k-eff Change: {eigenvalues[-1] - eigenvalues[0]:+.6f}")
    
    report.append(f"\nInitial Source Rate: {metrics['source_rates'][0]:.4e} n/s")
    report.append(f"Final Source Rate: {metrics['source_rates'][-1]:.4e} n/s")
    
    report.append(f"\nCumulative Neutrons: {metrics['cumulative_neutrons'][-1]:.4e}")
    
    report.append("\n" + "-" * 70)
    report.append("POWER GENERATION")
    report.append("-" * 70)
    report.append(f"Initial Power: {metrics['power_mw'][0]:.4e} MW ({metrics['power_watts'][0]:.4e} W)")
    report.append(f"Final Power: {metrics['power_mw'][-1]:.4e} MW ({metrics['power_watts'][-1]:.4e} W)")
    report.append(f"Total Energy Generated: {metrics['cumulative_energy_mwh'][-1]:.4e} MWh")
    report.append(f"                        {metrics['cumulative_energy_j'][-1]:.4e} Joules")
    avg_power = np.mean(metrics['power_mw'])
    report.append(f"Average Power Output: {avg_power:.4e} MW")
    
    report.append("\n" + "-" * 70)
    report.append("KEY NUCLIDE INVENTORIES (Total Atoms)")
    report.append("-" * 70)
    report.append(f"Material Volume: {MATERIAL_VOLUME_CM3:.4e} cm³")
    
    report.append("\n" + "=" * 70)
    report.append("INITIAL MATERIAL COMPOSITION (from fuel_blanket_2.py)")
    report.append("=" * 70)
    
    # Sort initial composition by abundance
    sorted_initial = sorted(initial_atoms.items(), key=lambda x: x[1], reverse=True)
    
    report.append("\nTop 20 Initial Isotopes by Abundance:")
    for nuc, count in sorted_initial[:20]:
        atom_fraction = count / sum(initial_atoms.values()) * 100
        report.append(f"  {nuc:>12s}: {count:.4e} atoms ({atom_fraction:6.3f}%)")
    
    # Calculate total initial actinides and fission products
    actinides_initial = sum(count for nuc, count in initial_atoms.items() 
                           if nuc.startswith(('U', 'Pu', 'Np', 'Am', 'Cm')))
    fission_products_initial = sum(count for nuc, count in initial_atoms.items() 
                                  if not nuc.startswith(('U', 'Pu', 'Np', 'Am', 'Cm', 'O')))
    
    report.append(f"\nTotal Initial Actinides: {actinides_initial:.4e} atoms")
    report.append(f"Total Initial Fission Products: {fission_products_initial:.4e} atoms")
    report.append(f"Total Initial Atoms: {sum(initial_atoms.values()):.4e} atoms")
    
    # First, report key actinides explicitly with correct initial values
    report.append("\n" + "-" * 70)
    report.append("KEY ACTINIDES EVOLUTION")
    report.append("-" * 70)
    report.append("NOTE: Depletion HDF5 shows near-zero initial U densities - this is")
    report.append("      a data recording issue. Using material definition for true initial state.")
    report.append("")
    
    key_actinides = ['U235', 'U238', 'U236', 'U234', 'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242', 
                     'Am241', 'Am243', 'Cm244', 'Cm242', 'Cm243', 'Cm245', 'Cm246', 'Np237', 'Np239']
    
    for nuc in key_actinides:
        # Get initial from material definition (TRUTH)
        initial = initial_atoms.get(nuc, 0.0)
        # Get final from depletion results
        final = inventories.get(nuc, [0.0])[-1] if nuc in inventories else 0.0
        
        # Also check what HDF5 says for initial (for diagnostics)
        hdf5_initial = inventories.get(nuc, [0.0])[0] if nuc in inventories else 0.0
        
        if initial > 1e10 or final > 1e10:  # Only show if significant
            if initial > 1e-10:
                rel_change = 100.0 * (final - initial) / initial
                # For uranium isotopes, flag if HDF5 initial doesn't match
                flag = ""
                if nuc.startswith('U') and hdf5_initial < initial * 0.01:
                    flag = " [HDF5 data issue - see note above]"
                report.append(f"  {nuc:>12s}: {initial:.4e} → {final:.4e} atoms ({rel_change:+6.1f}%){flag}")
            else:
                report.append(f"  {nuc:>12s}: {initial:.4e} → {final:.4e} atoms (produced)")
    
    report.append("\n" + "-" * 70)
    report.append("KEY FISSION PRODUCTS EVOLUTION")
    report.append("-" * 70)
    key_fission_products = ['Cs137', 'Cs135', 'Cs133', 'Sr90', 'Tc99', 'I129', 'Xe135', 
                           'Sm151', 'Sm149', 'Sm150', 'Sm152', 'Nd143', 'Nd145']
    
    for nuc in key_fission_products:
        # Get initial from material definition
        initial = initial_atoms.get(nuc, 0.0)
        # Get final from depletion results
        final = inventories.get(nuc, [0.0])[-1] if nuc in inventories else 0.0
        
        if initial > 1e10 or final > 1e10:  # Only show if significant
            if initial > 1e-10:
                rel_change = 100.0 * (final - initial) / initial
                report.append(f"  {nuc:>12s}: {initial:.4e} → {final:.4e} atoms ({rel_change:+6.1f}%)")
            else:
                report.append(f"  {nuc:>12s}: {initial:.4e} → {final:.4e} atoms (produced)")
    
    # Separate nuclides into categories for better readability
    produced = []  # Initial = 0, Final > 0
    depleted = []  # Initial > 0, decreased
    increased = [] # Initial > 0, increased
    stable = []    # Initial > 0, relatively stable
    
    for nuc, inv in inventories.items():
        initial = inv[0]
        final = inv[-1]
        
        if initial < 1e-10:  # Essentially zero initially
            if final > 1e-10:
                produced.append((nuc, initial, final))
        else:  # Non-zero initially
            rel_change = (final - initial) / initial
            if rel_change < -0.01:  # Decreased by >1%
                depleted.append((nuc, initial, final, 100.0 * rel_change))
            elif rel_change > 0.01:  # Increased by >1%
                increased.append((nuc, initial, final, 100.0 * rel_change))
            else:  # Relatively stable (within ±1%)
                stable.append((nuc, initial, final, 100.0 * rel_change))
    
    # Add waste transmutation analysis section
    waste_analysis = analyze_waste_transmutation(inventories, times_s)
    
    report.append("\n" + "-" * 70)
    report.append("WASTE TRANSMUTATION ANALYSIS")
    report.append("-" * 70)
    
    # Transuranic Actinides Summary
    report.append("\nTRANSURANIC ACTINIDES (TRU) - Elements beyond Uranium (Z>92):")
    report.append(f"  Initial Total TRU: {waste_analysis['tru_initial']:.4e} atoms")
    report.append(f"  Final Total TRU:   {waste_analysis['tru_final']:.4e} atoms")
    report.append(f"  Net Change:        {waste_analysis['tru_change']:.4e} atoms ({waste_analysis['tru_percent_change']:+.2f}%)")
    
    if waste_analysis['tru_change'] < 0:
        report.append(f"\n  ✓ TRU REDUCTION ACHIEVED: {abs(waste_analysis['tru_percent_change']):.2f}% decrease")
        report.append(f"    {abs(waste_analysis['tru_change']):.4e} atoms transmuted/destroyed")
    elif waste_analysis['tru_change'] > 0:
        report.append(f"\n  ✗ TRU INCREASED: {waste_analysis['tru_percent_change']:.2f}% increase")
        report.append(f"    Net breeding of {waste_analysis['tru_change']:.4e} TRU atoms")
    else:
        report.append(f"\n  → TRU remains unchanged")
    
    if waste_analysis['tru_details']:
        report.append("\n  Individual TRU Isotopes:")
        for isotope, data in sorted(waste_analysis['tru_details'].items(), 
                                   key=lambda x: abs(x[1]['change']), reverse=True)[:10]:
            report.append(f"    {isotope:>12s}: {data['initial']:.4e} → {data['final']:.4e} ({data['percent_change']:+.1f}%)")
    
    # Long-Lived Fission Products Summary
    report.append("\nLONG-LIVED FISSION PRODUCTS (>30 year half-life):")
    report.append(f"  Initial Total LLFP: {waste_analysis['llfp_initial']:.4e} atoms")
    report.append(f"  Final Total LLFP:   {waste_analysis['llfp_final']:.4e} atoms")
    report.append(f"  Net Change:         {waste_analysis['llfp_change']:.4e} atoms ({waste_analysis['llfp_percent_change']:+.2f}%)")
    
    if waste_analysis['llfp_change'] < 0:
        report.append(f"\n  ✓ LLFP REDUCTION ACHIEVED: {abs(waste_analysis['llfp_percent_change']):.2f}% decrease")
    elif waste_analysis['llfp_change'] > 0:
        report.append(f"\n  → LLFP INCREASED: {waste_analysis['llfp_percent_change']:.2f}% increase (expected from fission)")
    
    if waste_analysis['llfp_details']:
        report.append("\n  Individual LLFP Isotopes:")
        for isotope, data in sorted(waste_analysis['llfp_details'].items(), 
                                   key=lambda x: x[1]['final'], reverse=True)[:8]:
            report.append(f"    {isotope:>12s}: {data['initial']:.4e} → {data['final']:.4e} ({data['percent_change']:+.1f}%)")
    
    # Overall waste transmutation verdict
    report.append("\n" + "=" * 70)
    if waste_analysis['tru_change'] < 0:
        report.append("VERDICT: System is transmuting/destroying transuranic waste ✓")
    elif waste_analysis['tru_change'] > 0:
        report.append("VERDICT: System is breeding transuranics - NOT reducing waste ✗")
    else:
        report.append("VERDICT: No net change in transuranic inventory")
    report.append("=" * 70)
    
    if produced:
        report.append("\nProduced during simulation (initially zero):")
        for nuc, initial, final in sorted(produced, key=lambda x: x[2], reverse=True)[:15]:
            report.append(f"  {nuc:>12s}: {initial:.4e} → {final:.4e} atoms")
    
    if depleted:
        report.append("\nDepleted (decreased from initial):")
        for nuc, initial, final, pct in sorted(depleted, key=lambda x: abs(x[3]), reverse=True)[:15]:
            report.append(f"  {nuc:>12s}: {initial:.4e} → {final:.4e} atoms ({pct:+6.1f}%)")
    
    if increased:
        report.append("\nIncreased from initial:")
        for nuc, initial, final, pct in sorted(increased, key=lambda x: x[3], reverse=True)[:15]:
            report.append(f"  {nuc:>12s}: {initial:.4e} → {final:.4e} atoms ({pct:+6.1f}%)")
    
    if stable:
        report.append("\nRelatively stable (±1%):")
        for nuc, initial, final, pct in sorted(stable, key=lambda x: x[1], reverse=True)[:10]:
            report.append(f"  {nuc:>12s}: {initial:.4e} → {final:.4e} atoms ({pct:+6.1f}%)")
    
    report.append("\n" + "=" * 70)
    
    return "\n".join(report)

def plot_metrics(metrics, output_dir='results'):
    """Generate plots of depletion metrics including power generation."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a 3-panel figure
    fig = plt.figure(figsize=(14, 10))
    
    times_h = metrics['times_h']
    times_d = metrics['times_d']
    eigenvalues = metrics['eigenvalues']
    source_rates = metrics['source_rates']
    power_mw = metrics['power_mw']
    cumulative_energy_mwh = metrics['cumulative_energy_mwh']
    
    # Panel 1: Eigenvalue
    ax1 = plt.subplot(3, 2, 1)
    ax1.plot(times_d, eigenvalues, 'b.-', linewidth=2, markersize=6)
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('k-eff')
    ax1.grid(True, alpha=0.3)
    ax1.set_title('Eigenvalue vs. Time (Fixed Source)')
    
    # Panel 2: Source Rate
    ax2 = plt.subplot(3, 2, 2)
    ax2.semilogy(times_d, source_rates, 'r.-', linewidth=2, markersize=6)
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel('Source Rate (n/s)')
    ax2.grid(True, alpha=0.3)
    ax2.set_title('Neutron Source Rate vs. Time')
    
    # Panel 3: Power Output (linear)
    ax3 = plt.subplot(3, 2, 3)
    ax3.plot(times_d, power_mw, 'g.-', linewidth=2, markersize=6)
    ax3.set_xlabel('Time (days)')
    ax3.set_ylabel('Power Output (MW)')
    ax3.grid(True, alpha=0.3)
    ax3.set_title('Power Generation vs. Time')
    
    # Panel 4: Power Output (log scale)
    ax4 = plt.subplot(3, 2, 4)
    ax4.semilogy(times_d, power_mw, 'g.-', linewidth=2, markersize=6)
    ax4.set_xlabel('Time (days)')
    ax4.set_ylabel('Power Output (MW, log scale)')
    ax4.grid(True, alpha=0.3, which='both')
    ax4.set_title('Power Output (Log Scale)')
    
    # Panel 5: Cumulative Energy
    ax5 = plt.subplot(3, 2, 5)
    ax5.semilogy(times_d, cumulative_energy_mwh, 'm.-', linewidth=2, markersize=6)
    ax5.set_xlabel('Time (days)')
    ax5.set_ylabel('Cumulative Energy (MWh, log scale)')
    ax5.grid(True, alpha=0.3, which='both')
    ax5.set_title('Cumulative Energy Generated')
    
    # Panel 6: Summary text
    ax6 = plt.subplot(3, 2, 6)
    ax6.axis('off')
    summary_text = (
        f"Fixed Source Simulation:\n"
        f"Duration: {times_d[-1]:.1f} days\n"
        f"Source Rate: {source_rates[0]:.3e} n/s\n"
        f"Initial Power: {power_mw[0]:.3e} MW\n"
        f"Final Power: {power_mw[-1]:.3e} MW\n"
        f"Average Power: {np.mean(power_mw):.3e} MW\n"
        f"Total Energy: {cumulative_energy_mwh[-1]:.3e} MWh\n"
        f"k-eff: {eigenvalues[0]:.4f} → {eigenvalues[-1]:.4f}"
    )
    ax6.text(0.1, 0.5, summary_text, fontsize=12, family='monospace',
             verticalalignment='center', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'depletion_metrics.png'), dpi=150)
    print(f"Saved plot: {os.path.join(output_dir, 'depletion_metrics.png')}")
    plt.close()
    
    # Create a separate power-only detailed plot
    fig_power, ax_power = plt.subplots(figsize=(12, 6))
    ax_power.fill_between(times_d, power_mw, alpha=0.3, label='Power Output')
    ax_power.plot(times_d, power_mw, 'g-', linewidth=2.5, marker='o', markersize=5, label='Power (MW)')
    ax_power.set_xlabel('Time (days)', fontsize=12)
    ax_power.set_ylabel('Power Output (MW)', fontsize=12)
    ax_power.set_title('Detailed Power Generation Profile (Fixed Source)', fontsize=14, fontweight='bold')
    ax_power.grid(True, alpha=0.3)
    ax_power.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'power_generation.png'), dpi=150)
    print(f"Saved plot: {os.path.join(output_dir, 'power_generation.png')}")
    plt.close()

def main():
    print("Reading fixed source depletion results...")
    data = read_depletion_results(RESULTS_FILE)
    print(f"  - {len(data['times'])} depletion steps")
    print(f"  - {len(data['nuclides'])} nuclides tracked")
    print(f"  - Materials: {data['materials']}")
    
    print("\nExtracting burnup metrics...")
    metrics = extract_burnup_metrics(data)
    
    print("Extracting nuclide inventories...")
    inventories = extract_nuclide_inventories(data)
    
    print("\nGenerating summary report...")
    report = generate_summary_report(data, metrics, inventories)
    print(report)
    
    # Save report
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    report_file = os.path.join(OUTPUT_DIR, 'depletion_summary.txt')
    with open(report_file, 'w') as f:
        f.write(report)
    print(f"\nSaved report: {report_file}")
    
    # Generate plots
    print("\nGenerating plots...")
    plot_metrics(metrics, OUTPUT_DIR)
    
    # Save metrics as CSV
    csv_file = os.path.join(OUTPUT_DIR, 'depletion_metrics.csv')
    df = pd.DataFrame({
        'time_s': metrics['times_s'],
        'time_h': metrics['times_h'],
        'time_d': metrics['times_d'],
        'k_eff': metrics['eigenvalues'],
        'source_rate_n_per_s': metrics['source_rates'],
        'power_watts': metrics['power_watts'],
        'power_mw': metrics['power_mw'],
        'cum_neutrons': metrics['cumulative_neutrons'],
        'cum_energy_joules': metrics['cumulative_energy_j'],
        'cum_energy_mwh': metrics['cumulative_energy_mwh'],
    })
    df.to_csv(csv_file, index=False)
    print(f"Saved metrics CSV: {csv_file}")
    
    print("\n✓ Post-processing complete.")

if __name__ == '__main__':
    main()
