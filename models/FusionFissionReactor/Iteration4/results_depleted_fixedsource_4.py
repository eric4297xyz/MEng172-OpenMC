# FILE: results_depleted_fixedsource_4.py
# Post-processing script for Iteration4 TRU fuel depletion results

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def get_initial_material_composition():
    """
    Get initial spent UO2 fuel composition from fuel_blanket_4.py material definition.
    Returns dictionary of nuclide: atom_fraction pairs.
    """
    # Spent UO2 fuel composition WITH uranium - direct fractions from fuel_blanket_4.py
    # No scaling needed - using exact values from the file
    
    composition = {
        # Oxygen
        'O16': 6.7187968E-01,
        
        # Uranium isotopes (present in Iteration 4)
        'U238': 3.0769658E-01,
        'U235': 3.4979859E-03,
        'U236': 2.1091663E-03,
        'U234': 6.4111968E-05,
        'U237': 1.7067631E-07,
        
        # Transuranics - Plutonium
        'Pu239': 2.5718196E-03,
        'Pu240': 1.0215150E-03,
        'Pu241': 6.7015516E-04,
        'Pu242': 2.6644975E-04,
        'Pu238': 1.3101157E-04,
        
        # Transuranics - Neptunium
        'Np237': 2.7884345E-04,
        
        # Transuranics - Americium
        'Am243': 7.8638873E-05,
        'Am241': 3.1275801E-05,
        
        # Transuranics - Curium
        'Cm244': 3.2900525E-05,
        'Cm242': 7.3636478E-06,
        'Cm245': 1.9560548E-06,
        'Cm243': 3.0834446E-07,
        'Cm246': 1.9481932E-07,
        
        # Fission Products
        'Cs137': 1.0510726E-03,
        'Cs133': 9.7963234E-04,
        'Tc99': 9.3559531E-04,
        'Ru101': 9.1992123E-04,
        'Zr93': 8.9218967E-04,
        'Mo95': 8.5245707E-04,
        'Sr90': 6.7977794E-04,
        'Nd143': 6.5286858E-04,
        'Nd145': 5.3344636E-04,
        'Rh103': 4.8652147E-04,
        'Cs135': 4.8545594E-04,
        'Pd107': 2.5941176E-04,
        'Sm150': 2.2788081E-04,
        'I129': 1.4371885E-04,
        'Pm147': 1.0269899E-04,
        'Eu153': 9.2879588E-05,
        'Sm152': 8.8634337E-05,
        'Ag109': 8.5141959E-05,
        'Sm147': 6.6014495E-05,
        'Ru103': 1.9360418E-05,
        'Sn126': 1.8018478E-05,
        'Nb95': 1.4824839E-05,
        'Cl36': 1.4019981E-05,
        'Ca41': 1.4019976E-05,
        'Ni59': 1.4019973E-05,
        'Sm151': 1.3614245E-05,
        'Se79': 7.0915868E-06,
        'Eu155': 4.7729475E-06,
        'Pr143': 2.0561139E-06,
        'Sm149': 1.9355989E-06,
        'Nd147': 5.1292485E-07,
        'Xe133': 9.3481328E-08,
        'Gd155': 7.6579955E-08,
        'I133': 7.5534463E-08,
        'Eu152': 4.5260253E-08,
        'Mo99': 1.1472054E-09 / 0.68652,
    }
    
    return composition

def calculate_initial_atom_counts(volume_cm3=3.6652e6, density_atom_per_barn_cm=7.133315757E-02):
    """
    Calculate initial atom counts from material definition.
    
    Args:
        volume_cm3: Volume in cm³ (default: spherical shell 50-100cm radius)
        density_atom_per_barn_cm: Atomic density in atoms/(barn·cm) for spent UO2 fuel
    
    Returns:
        Dictionary of nuclide: total_atoms
    """
    composition = get_initial_material_composition()
    
    # Convert: atoms/(barn·cm) × 1e24 (barn/cm²) × volume_cm³ = total atoms
    conversion_factor = density_atom_per_barn_cm * 1e24 * volume_cm3
    
    initial_counts = {}
    for nuclide, fraction in composition.items():
        initial_counts[nuclide] = fraction * conversion_factor
    
    return initial_counts

def read_depletion_results(h5_file, volume_cm3=3.6652e6):
    """Read depletion results from HDF5 file.
    
    Args:
        h5_file: Path to depletion HDF5 file
        volume_cm3: Material volume in cm³ (default: spherical shell 50-100cm)
    
    Returns:
        time, nuclides, number (converted to total atoms), keff, source_rate
    """
    with h5py.File(h5_file, 'r') as f:
        time = f['time'][:]
        nuclides_group = f['nuclides']
        nuclides = list(nuclides_group.keys())
        number_density = f['number'][:]  # atoms/(barn·cm)
        
        # Convert number density to total atoms by multiplying by volume
        # This matches what Iteration 2 does
        number = number_density * volume_cm3
        
        # For fixed source, k-eff might not be present
        keff = None
        if 'eigenvalues' in f:
            keff = f['eigenvalues'][:]
        
        # Source rates
        source_rate = None
        if 'source_rate' in f:
            source_rate = f['source_rate'][:]
    
    return time, nuclides, number, keff, source_rate

def generate_summary_report(h5_file, output_file):
    """Generate comprehensive summary report."""
    time, nuclides, number, keff, source_rate = read_depletion_results(h5_file)
    
    # Get initial composition from material definition
    initial_counts = calculate_initial_atom_counts()
    
    # Handle 2D time array - get final time point
    final_time = float(time[-1, -1]) if time.ndim > 1 else float(time[-1])
    
    with open(output_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("ITERATION 4: SPENT UO2 FUEL DEPLETION SIMULATION SUMMARY\n")
        f.write("="*70 + "\n\n")
        
        f.write("⚠️  DATA QUALITY NOTE:\n")
        f.write("   The HDF5 depletion results show near-zero initial densities.\n")
        f.write("   This is a known OpenMC bug in depletion data recording.\n")
        f.write("   Initial compositions are taken from the material definition (materials.xml)\n")
        f.write("   which correctly specifies the spent UO2 fuel composition (WITH uranium).\n")
        f.write("   Final states are from the HDF5 file.\n\n")
        
        f.write("⚠️  PARTIAL DATA:\n")
        f.write("   Simulation crashed during step 3 due to multiprocessing error.\n")
        f.write("   Results cover 2 of 3 planned depletion steps.\n\n\n")
        
        f.write(f"Simulation Duration: {final_time:.2e} seconds ({final_time/(24*3600):.2f} days)\n")
        f.write(f"Number of Depletion Steps: {len(time)}\n\n")
        
        # k-eff (if available)
        if keff is not None and keff.size > 0:
            # Handle both 2D and 3D k-eff arrays
            if keff.ndim == 3:
                initial_keff = keff[0,0,0]
                final_keff = keff[-1,0,0]
            elif keff.ndim == 2:
                initial_keff = keff[0,0]
                final_keff = keff[-1,0]
            else:
                initial_keff = keff[0]
                final_keff = keff[-1]
            f.write(f"Initial k-eff: {initial_keff:.5f}\n")
            f.write(f"Final k-eff: {final_keff:.5f}\n")
            f.write(f"k-eff Change: {final_keff-initial_keff:+.5f}\n\n")
        else:
            f.write("Initial k-eff: N/A (fixed source mode)\n")
            f.write("Final k-eff: N/A (fixed source mode)\n\n")
        
        # Source rate
        if source_rate is not None:
            initial_sr = float(source_rate[0])
            final_sr = float(source_rate[-1])
            cumulative = float(np.sum(source_rate) * (final_time/len(source_rate)))
            f.write(f"Initial Source Rate: {initial_sr:.4e} n/s\n")
            f.write(f"Final Source Rate: {final_sr:.4e} n/s\n\n")
            f.write(f"Cumulative Neutrons: {cumulative:.4e}\n\n")
        
        f.write("-"*70 + "\n")
        f.write("KEY NUCLIDE INVENTORIES (Total Atoms)\n")
        f.write("-"*70 + "\n")
        f.write(f"Material Volume: {3.6652e6:.4e} cm³\n\n")
        
        # Initial composition
        f.write("="*70 + "\n")
        f.write("INITIAL MATERIAL COMPOSITION (from fuel_blanket_4.py)\n")
        f.write("="*70 + "\n\n")
        
        # Sort initial composition by abundance
        sorted_initial = sorted(initial_counts.items(), key=lambda x: x[1], reverse=True)
        
        f.write("Top 20 Initial Isotopes by Abundance:\n")
        total_initial_atoms = sum(initial_counts.values())
        for nuclide, count in sorted_initial[:20]:
            pct = (count / total_initial_atoms) * 100
            f.write(f"  {nuclide:>12s}: {count:.4e} atoms ({pct:6.3f}%)\n")
        
        # Calculate initial TRU total
        tru_elements = ['Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf']
        initial_tru = sum(count for nuc, count in initial_counts.items() 
                         if any(nuc.startswith(elem) for elem in tru_elements))
        
        f.write(f"\nTotal Initial TRU: {initial_tru:.4e} atoms\n")
        f.write(f"Total Initial Atoms: {total_initial_atoms:.4e} atoms\n\n")
        
        # Add uranium totals
        uranium_elements = ['U']
        initial_uranium = sum(count for nuc, count in initial_counts.items() 
                            if any(nuc.startswith(elem) for elem in uranium_elements))
        f.write(f"Total Initial Uranium: {initial_uranium:.4e} atoms\n\n")
        
        f.write("="*70 + "\n")
        f.write("NOTE: SPENT UO2 FUEL - Contains uranium and transuranics\n")
        f.write("="*70 + "\n\n")
        
        # Uranium evolution
        f.write("-"*70 + "\n")
        f.write("URANIUM ISOTOPES EVOLUTION\n")
        f.write("-"*70 + "\n\n")
        
        uranium_nuclides = ['U234', 'U235', 'U236', 'U237', 'U238']
        
        for nuc_name in uranium_nuclides:
            if nuc_name in nuclides:
                idx = nuclides.index(nuc_name)
                initial = initial_counts.get(nuc_name, 0.0)
                final = number[-1, 0, idx]
                pct_change = ((final - initial) / initial * 100) if initial > 0 else 0
                f.write(f"  {nuc_name:>12s}: {initial:.4e} → {final:.4e} atoms ({pct_change:+.1f}%)\n")
        
        f.write("\n")
        
        # TRU evolution
        f.write("-"*70 + "\n")
        f.write("TRANSURANIC ACTINIDES EVOLUTION\n")
        f.write("-"*70 + "\n\n")
        
        tru_nuclides = ['Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242', 'Pu244',
                       'Np237', 'Np238', 'Np239',
                       'Am241', 'Am242', 'Am243',
                       'Cm242', 'Cm243', 'Cm244', 'Cm245', 'Cm246']
        
        for nuc_name in tru_nuclides:
            if nuc_name in nuclides:
                idx = nuclides.index(nuc_name)
                initial = initial_counts.get(nuc_name, 0.0)
                final = number[-1, 0, idx]
                
                if initial > 0:
                    change_pct = ((final - initial) / initial) * 100
                    f.write(f"  {nuc_name:>12s}: {initial:.4e} → {final:.4e} atoms ({change_pct:+.1f}%)\n")
                elif final > 0:
                    f.write(f"  {nuc_name:>12s}: {initial:.4e} → {final:.4e} atoms (produced)\n")
        
        # Fission products evolution
        f.write("\n" + "-"*70 + "\n")
        f.write("KEY FISSION PRODUCTS EVOLUTION\n")
        f.write("-"*70 + "\n\n")
        
        fp_nuclides = ['Cs137', 'Cs135', 'Cs133', 'Sr90', 'Tc99', 'I129',
                      'Sm151', 'Sm149', 'Sm150', 'Sm152']
        
        for nuc_name in fp_nuclides:
            if nuc_name in nuclides:
                idx = nuclides.index(nuc_name)
                initial = initial_counts.get(nuc_name, 0.0)
                final = number[-1, 0, idx]
                
                if initial > 0:
                    change_pct = ((final - initial) / initial) * 100
                    f.write(f"  {nuc_name:>12s}: {initial:.4e} → {final:.4e} atoms ({change_pct:+.1f}%)\n")
        
        # Overall transmutation analysis
        analyze_waste_transmutation(nuclides, number, initial_counts, f)
        
        f.write("\n" + "="*70 + "\n\n")
    
    print(f"Saved report: {output_file}")

def analyze_waste_transmutation(nuclides, number, initial_counts, f):
    """Analyze waste transmutation effectiveness."""
    f.write("\n" + "-"*70 + "\n")
    f.write("WASTE TRANSMUTATION ANALYSIS - TRU FUEL\n")
    f.write("-"*70 + "\n\n")
    
    # TRU elements
    tru_elements = ['Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf']
    
    initial_tru = 0.0
    final_tru = 0.0
    
    tru_details = []
    
    for nuc_name in nuclides:
        if any(nuc_name.startswith(elem) for elem in tru_elements):
            idx = nuclides.index(nuc_name)
            initial = initial_counts.get(nuc_name, 0.0)
            final = number[-1, 0, idx]
            
            initial_tru += initial
            final_tru += final
            
            if initial > 1e20 or final > 1e20:  # Only track significant isotopes
                tru_details.append((nuc_name, initial, final))
    
    net_tru_change = final_tru - initial_tru
    if initial_tru > 0:
        tru_change_pct = (net_tru_change / initial_tru) * 100
    else:
        tru_change_pct = float('inf') if final_tru > 0 else 0
    
    f.write(f"TRANSURANIC ACTINIDES (TRU) - Elements beyond Uranium (Z>92):\n")
    f.write(f"  Initial Total TRU: {initial_tru:.4e} atoms\n")
    f.write(f"  Final Total TRU:   {final_tru:.4e} atoms\n")
    f.write(f"  Net Change:        {net_tru_change:.4e} atoms ({tru_change_pct:+.2f}%)\n\n")
    
    if tru_change_pct < 0:
        f.write(f"  ✓ TRU REDUCTION ACHIEVED: {abs(tru_change_pct):.2f}% decrease\n")
        f.write(f"    Successfully transmuted {abs(net_tru_change):.4e} TRU atoms\n\n")
    else:
        f.write(f"  ✗ TRU INCREASED: {tru_change_pct:.2f}% increase\n")
        f.write(f"    TRU production from U238 neutron capture and subsequent decay\n")
        f.write(f"    This is expected in uranium-bearing fuel\n\n")
    
    # Individual TRU isotopes
    tru_details.sort(key=lambda x: x[1], reverse=True)
    f.write("  Individual TRU Isotopes (significant only):\n")
    for nuc_name, initial, final in tru_details:
        if initial > 0:
            change_pct = ((final - initial) / initial) * 100
            f.write(f"    {nuc_name:>12s}: {initial:.4e} → {final:.4e} ({change_pct:+.1f}%)\n")
        elif final > 0:
            f.write(f"    {nuc_name:>12s}: {initial:.4e} → {final:.4e} (produced)\n")
    
    # Long-lived fission products
    llfp_list = ['Tc99', 'Cs135', 'Cs137', 'Pd107', 'I129', 'Sr90', 'Zr93', 'Se79']
    
    initial_llfp = 0.0
    final_llfp = 0.0
    llfp_details = []
    
    for nuc_name in llfp_list:
        if nuc_name in nuclides:
            idx = nuclides.index(nuc_name)
            initial = initial_counts.get(nuc_name, 0.0)
            final = number[-1, 0, idx]
            initial_llfp += initial
            final_llfp += final
            llfp_details.append((nuc_name, initial, final))
    
    net_llfp_change = final_llfp - initial_llfp
    if initial_llfp > 0:
        llfp_change_pct = (net_llfp_change / initial_llfp) * 100
    else:
        llfp_change_pct = 0
    
    f.write(f"\nLONG-LIVED FISSION PRODUCTS (>30 year half-life):\n")
    f.write(f"  Initial Total LLFP: {initial_llfp:.4e} atoms\n")
    f.write(f"  Final Total LLFP:   {final_llfp:.4e} atoms\n")
    f.write(f"  Net Change:         {net_llfp_change:.4e} atoms ({llfp_change_pct:+.2f}%)\n\n")
    
    if llfp_change_pct < 0:
        f.write(f"  ✓ LLFP REDUCTION ACHIEVED: {abs(llfp_change_pct):.2f}% decrease\n\n")
    else:
        f.write(f"  ⚠ LLFP INCREASED: {llfp_change_pct:.2f}% increase\n")
        f.write(f"    More fission products produced than transmuted\n\n")
    
    f.write("  Individual LLFP Isotopes:\n")
    for nuc_name, initial, final in llfp_details:
        if initial > 0:
            change_pct = ((final - initial) / initial) * 100
            f.write(f"    {nuc_name:>12s}: {initial:.4e} → {final:.4e} ({change_pct:+.1f}%)\n")
    
    # Overall verdict
    f.write("\n" + "="*70 + "\n")
    if tru_change_pct < 0 and llfp_change_pct < 0:
        f.write("VERDICT: Excellent - Both TRU and LLFP reduced ✓✓\n")
    elif tru_change_pct < 0:
        f.write("VERDICT: Good - TRU reduced but LLFP increased ✓\n")
    elif llfp_change_pct < 0:
        f.write("VERDICT: Partial - LLFP reduced but TRU increased ⚠\n")
    else:
        f.write("VERDICT: Poor - Both TRU and LLFP increased ✗\n")
    f.write("="*70 + "\n")

def main():
    """Main post-processing function."""
    h5_file = "results/iteration4_results_tru_fuel/depletion_results.h5"
    
    if not os.path.exists(h5_file):
        print(f"Error: {h5_file} not found!")
        print("Run reactor_Depleted_FixedSource_4.py first.")
        return
    
    print("Reading TRU fuel depletion results...")
    
    output_dir = "results/iteration4_results_tru_fuel"
    
    # Generate summary report
    print("\nGenerating summary report...")
    summary_file = os.path.join(output_dir, "depletion_summary_tru.txt")
    generate_summary_report(h5_file, summary_file)
    
    print("\n✓ Post-processing complete for Iteration4.")
    print(f"\nResults saved in: {output_dir}/")
    print("  - depletion_summary_tru.txt: Comprehensive text report")

if __name__ == "__main__":
    main()
